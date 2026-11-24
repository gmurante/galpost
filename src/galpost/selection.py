# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Selection of the galaxies to analyse.

Reproduces the three ExtractGalaxies modes (FoF groups, SubFind subhaloes,
explicit list) with the same rules and the same output radii, plus an optional
check of the low-resolution contamination for zoom-in runs.

All positions and radii of :class:`Galaxy` are in ExtractGalaxies units:
physical kpc (divided by h if Htrue).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np

from .catalogs import FofCatalog, SubfindCatalog
from .units import Units

__all__ = ["Galaxy", "select_from_fof", "select_from_subfind", "select_from_list",
           "lowres_fraction", "virial_radius_eg"]


@dataclass
class Galaxy:
    gal_id: int
    com: np.ndarray               # centre (kpc), possibly refined later
    rmin: float
    rmax: float
    rgal: float
    rvir: float
    ralign: float
    log_bins: bool = False
    compute_cm: bool = False
    track_cm: bool = False
    align: bool = True
    pot: bool = False
    central: int = 1              # EG: 1 for FoF; FirstSubOfHalo index for SubFind
    fof_mass: float = 0.0         # Msun
    group: int = -1               # catalogue indices (galpost additions)
    subhalo: int = -1
    lowres_fraction: float = float("nan")
    extra: dict = field(default_factory=dict)


def virial_radius_eg(mass_msun: float, header) -> float:
    """Virial radius (kpc) exactly as computed in ExtractGalaxies from a FoF mass."""
    h = header.hubble_param
    z = header.redshift
    ez2 = header.omega0 * (1.0 + z) ** 3 + header.omega_lambda
    # NB: EG uses the exponent 0.333333 (not 1/3): kept for bit-compatibility
    # (relative effect ~5e-6 on rvir for Milky-Way-sized haloes)
    return (4.296e-6 * mass_msun / h / h / ez2) ** 0.333333


def _wrap(x_kpc: np.ndarray, box_kpc: float, periodic: bool) -> np.ndarray:
    if not periodic:
        return x_kpc
    return np.mod(x_kpc, box_kpc)


def _select_by_threshold(values: np.ndarray, threshold: float, eligible: np.ndarray):
    """EG rule: threshold>0 -> values >= threshold; <0 -> the -threshold largest."""
    if threshold > 0:
        return eligible & (values >= threshold), threshold
    n = int(-threshold)
    v = np.where(eligible, values, -np.inf)
    if n <= 0 or not np.any(eligible):
        return np.zeros_like(eligible), threshold
    cut = np.sort(v)[::-1][min(n, eligible.sum()) - 1]
    return eligible & (values >= cut), cut


# --------------------------------------------------------------------------
def select_from_fof(fof: FofCatalog, units: Units, header, min_mass: float,
                    rmin: float, rmax: float, rgal: float, ralign: float,
                    log_bins: bool, do_not_align: bool) -> List[Galaxy]:
    """ExtractGalaxies ``ReadFofHalos`` mode."""
    mass = fof.mass_type.sum(axis=1) * units.mass_msun      # EG: sum of MassType
    sel, cut = _select_by_threshold(mass, min_mass, np.ones(len(mass), bool))
    idx = np.nonzero(sel)[0]
    box = header.boxsize * units.length_kpc
    gals = []
    for i in idx:
        rv = virial_radius_eg(mass[i], header)
        rg = rgal if rgal > 0 else 0.1 * rv
        gals.append(Galaxy(
            gal_id=0 if len(idx) == 1 else int(i),
            com=_wrap(fof.cm[i].astype(np.float64) * units.length_kpc, box, units.periodic),
            rmin=rmin, rvir=rv, rmax=rmax if rmax > 0 else rv, rgal=rg,
            ralign=ralign if ralign > 0 else rg, log_bins=log_bins,
            compute_cm=False, track_cm=True, align=not do_not_align, pot=False,
            central=1, fof_mass=float(mass[i]), group=int(i)))
    return gals


def select_from_subfind(cat: SubfindCatalog, units: Units, header,
                        min_stellar_mass: float, rmin: float, rmax: float,
                        ralign: float, log_bins: bool, do_not_align: bool,
                        lowres_types: Sequence[int] = (2, 3)) -> List[Galaxy]:
    """ExtractGalaxies ``ReadSubFindGals`` mode.

    Subhaloes containing low-resolution particles are excluded (as in EG).
    """
    smst = cat.subhalo("SMST").astype(np.float64)
    mstar = smst[:, 4] * units.mass_msun
    clean = np.all(smst[:, list(lowres_types)] == 0, axis=1)
    sel, cut = _select_by_threshold(mstar, min_stellar_mass, clean)
    parent = cat.subhalo("GRNR")
    rcrit = cat.group("RCRI").astype(np.float64)
    first = cat.group("FSUB")
    spos = cat.subhalo("SPOS").astype(np.float64)
    box = header.boxsize * units.length_kpc
    gals = []
    for i in np.nonzero(sel)[0]:
        g = int(parent[i])
        rv = rcrit[g] * units.length_kpc
        rg = 0.1 * rv
        gals.append(Galaxy(
            gal_id=int(i), com=_wrap(spos[i] * units.length_kpc, box, units.periodic),
            rmin=rmin, rvir=rv, rmax=rmax if rmax > 0 else rv, rgal=rg,
            ralign=ralign if ralign > 0 else rg, log_bins=log_bins,
            compute_cm=False, track_cm=False, align=not do_not_align, pot=False,
            central=int(first[g]), group=g, subhalo=int(i),
            fof_mass=float(cat.group("MTOT")[g] * units.mass_msun)))
    return gals


def select_from_list(lines, units: Units, header) -> List[Galaxy]:
    """``ListOfGals`` / ``Galfile`` mode: positions and radii are given in kpc."""
    gals = []
    for i, l in enumerate(lines):
        gals.append(Galaxy(
            gal_id=i, com=np.array([l.x, l.y, l.z], dtype=np.float64),
            rmin=l.rmin, rmax=l.rmax, rgal=l.ralign, rvir=l.rmax, ralign=l.ralign,
            log_bins=l.log, compute_cm=l.compute_cm, track_cm=l.track_cm,
            align=l.align, pot=l.pot, central=1))
    return gals


# --------------------------------------------------------------------------
def lowres_fraction(snapshot, galaxies: List[Galaxy], units: Units,
                    lowres_types: Sequence[int] = (2, 3), radius: str = "rvir") -> np.ndarray:
    """Mass fraction of low-resolution particles within `radius` of each galaxy.

    FoF groups never contain low-resolution particles (they are not linked),
    so contamination must be measured on the particle distribution.
    Sets ``Galaxy.lowres_fraction`` and returns the fractions.
    """
    from scipy.spatial import cKDTree

    box = snapshot.header.boxsize * units.length_kpc
    pos_lr, m_lr = [], []
    for t in lowres_types:
        if snapshot.npart[t] > 0:
            pos_lr.append(snapshot.read("POS", t).astype(np.float64) * units.length_kpc)
            m_lr.append(snapshot.read("MASS", t).astype(np.float64) * units.mass_msun)
    out = np.zeros(len(galaxies))
    if not pos_lr:
        for g in galaxies:
            g.lowres_fraction = 0.0
        return out
    pos_lr = np.concatenate(pos_lr)
    m_lr = np.concatenate(m_lr)
    tree = cKDTree(np.mod(pos_lr, box) if units.periodic else pos_lr,
                   boxsize=box if units.periodic else None)
    for k, g in enumerate(galaxies):
        r = getattr(g, radius)
        idx = tree.query_ball_point(np.mod(g.com, box) if units.periodic else g.com, r)
        mlow = m_lr[idx].sum()
        mtot = g.fof_mass if g.fof_mass > 0 else np.nan
        out[k] = mlow / (mtot + mlow) if mtot == mtot else np.nan
        g.lowres_fraction = float(out[k])
        g.extra["lowres_mass"] = float(mlow)
    return out
