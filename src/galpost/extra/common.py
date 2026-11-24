# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Common tools of the galpost extra analyses.

Units: physical kpc, km/s, Msun, Gyr, yr (h included if Htrue), as in the
ExtractGalaxies outputs.
"""

from __future__ import annotations

import math
from typing import Dict, Optional

import numpy as np

__all__ = ["G", "KPC_KM", "YR_S", "GalaxyView", "radius_value", "profile_edges",
           "write_table", "solar_log_ratio", "ATOMIC_WEIGHT", "LOG_EPS_SUN"]

G = 4.30091e-6                 # kpc (km/s)^2 / Msun
KPC_KM = 3.085678e16           # km per kpc
YR_S = 3.15576e7               # s per yr
KMS_KPC_TO_PER_YR = YR_S / KPC_KM   # (km/s)/kpc -> 1/yr

# Asplund et al. (2009) photospheric abundances, log eps (H = 12)
LOG_EPS_SUN = {"He": 10.93, "C": 8.43, "N": 7.83, "O": 8.69, "Ne": 7.93, "Na": 6.24,
               "Mg": 7.60, "Al": 6.45, "Si": 7.51, "S": 7.12, "Ar": 6.40, "Ca": 6.34,
               "Fe": 7.50, "Ni": 6.22}
ATOMIC_WEIGHT = {"H": 1.008, "He": 4.0026, "C": 12.011, "N": 14.007, "O": 15.999,
                 "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085,
                 "S": 32.06, "Ar": 39.948, "Ca": 40.078, "Fe": 55.845, "Ni": 58.693}


def solar_log_ratio(el: str) -> float:
    """log10 of the solar mass ratio X_el / X_H."""
    return LOG_EPS_SUN[el] - 12.0 + math.log10(ATOMIC_WEIGHT[el] / ATOMIC_WEIGHT["H"])


def radius_value(tag, g) -> float:
    s = str(tag)
    try:
        return float(s)
    except ValueError:
        return float(getattr(g, s))


def profile_edges(params, g):
    rmax = radius_value(params["ExtraRmax"], g)
    n = int(params["ExtraBins"])
    if params["ExtraLog"]:
        rmin = max(rmax * 1.0e-3, 0.05)
        return np.concatenate([[0.0], np.logspace(math.log10(rmin), math.log10(rmax), n)])
    return np.linspace(0.0, rmax, n + 1)


def write_table(path: str, columns, data, header_lines=(), fmt="%14.6e"):
    """ASCII table: '# ' header lines, numbered column descriptions, rows."""
    with open(path, "w") as f:
        for h in header_lines:
            f.write("# " + h + "\n")
        for i, c in enumerate(columns):
            f.write("# %2d: %s\n" % (i + 1, c))
        data = np.atleast_2d(np.asarray(data, dtype=np.float64))
        for row in data:
            f.write(" ".join(fmt % v for v in row) + "\n")


class GalaxyView:
    """Particles of a galaxy in its frame, within a sphere of radius `rmax_load`.

    `restrict` applies the LoadParticleList / SubUsePartAssignment restrictions
    (for 'galaxy' quantities); halo-scale analyses use restrict=False.
    """

    def __init__(self, analyser, resolver, igal: int, radius: float, restrict: bool = True):
        d = analyser.d
        g = analyser.gals[igal]
        self.d, self.a, self.g, self.res = d, analyser, g, resolver
        dp = d.rel_pos(g.com)
        sel = np.einsum("ij,ij->i", dp, dp) <= radius * radius
        if restrict:
            if analyser.sf_member is not None:
                sel &= analyser.sf_member == g.subhalo
            if analyser.particle_flag is not None:
                sel &= analyser.particle_flag
        self.idx = np.nonzero(sel)[0]
        self.pos, self.vel = analyser.frame(g, self.idx)
        self.mass = d.mass[self.idx].astype(np.float64)
        self.type = d.ptype()[self.idx]
        self.r = np.sqrt(np.einsum("ij,ij->i", self.pos, self.pos))
        self.R = np.sqrt(self.pos[:, 0] ** 2 + self.pos[:, 1] ** 2)
        self.z = self.pos[:, 2]
        self._c: Dict[str, np.ndarray] = {}

    def of(self, t: int) -> np.ndarray:
        return self.type == t

    # ---- gas ------------------------------------------------------------
    def _gas_index(self):
        m = self.of(0)
        return m, self.idx[m] - self.d.nstart[0]

    def gas_array(self, name: str) -> np.ndarray:
        """Gas quantity over the view (0 for non-gas): mcld, mhot, sfr, temp, mmol,
        rho, u, thot."""
        key = "gas:" + name
        if key not in self._c:
            out = np.zeros(len(self.idx))
            m, gi = self._gas_index()
            d = self.d
            if name in ("mmol", "thot"):
                full = self.res.derived_block(name.upper())
                out[m] = full[self.idx[m]].astype(np.float64)
            else:
                src = {"mcld": d.mcld, "mhot": d.mhot, "sfr": d.sfr, "temp": d.temp,
                       "rho": d.rho, "u": d.u}[name]
                out[m] = src[gi].astype(np.float64)
            self._c[key] = out
        return self._c[key]

    # ---- metals ----------------------------------------------------------
    def metals(self) -> Optional[np.ndarray]:
        """Element masses (Msun) for gas and stars (NaN rows for other types)."""
        if "metals" in self._c:
            return self._c["metals"]
        d = self.d
        if d.metals is None:
            return None
        out = np.full((len(self.idx), d.nmet), np.nan)
        ng = int(d.npart[0])
        mg = self.of(0)
        out[mg] = d.metals[self.idx[mg] - d.nstart[0]]
        ms = self.of(4)
        out[ms] = d.metals[ng + self.idx[ms] - d.nstart[4]]
        self._c["metals"] = out
        return out

    def reference_mass(self) -> np.ndarray:
        """Mass to which element masses refer: gas mass, initial stellar mass."""
        ref = self.mass.copy()
        ms = self.of(4)
        if self.d.im is not None:
            ref[ms] = self.d.im[self.idx[ms] - self.d.nstart[4]].astype(np.float64)
        return ref

    def hydrogen_mass(self) -> Optional[np.ndarray]:
        met = self.metals()
        if met is None:
            return None
        return self.reference_mass() - np.nansum(met, axis=1)

    def element(self, names, el: str) -> Optional[np.ndarray]:
        met = self.metals()
        if met is None or el not in names:
            return None
        return met[:, names.index(el)]

    # ---- stars -------------------------------------------------------------
    def stellar_age(self) -> Optional[np.ndarray]:
        """Age (Gyr) of stars (NaN for other types)."""
        d = self.d
        if d.age is None:
            return None
        out = np.full(len(self.idx), np.nan)
        ms = self.of(4)
        a = d.age[self.idx[ms] - d.nstart[4]].astype(np.float64)
        tnow = float(self.a.cosmological_time(d.header.time))
        out[ms] = tnow - self.a.cosmological_time(a)
        return out


def binned(x, edges, weights, values=None):
    """Sum of weights (and of weights*values) in bins; returns (sumw, sumwv, count)."""
    ib = np.searchsorted(edges, x, side="right") - 1
    ok = (ib >= 0) & (ib < len(edges) - 1)
    n = len(edges) - 1
    sw = np.bincount(ib[ok], weights=weights[ok], minlength=n)
    cnt = np.bincount(ib[ok], minlength=n)
    swv = None
    if values is not None:
        swv = np.bincount(ib[ok], weights=(weights * values)[ok], minlength=n)
    return sw, swv, cnt


def wmean(sw, swv):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(sw > 0, swv / sw, np.nan)


def linear_fit(x, y, w=None):
    ok = np.isfinite(x) & np.isfinite(y)
    if w is not None:
        ok &= np.isfinite(w) & (w > 0)
    if ok.sum() < 2:
        return np.nan, np.nan
    p = np.polyfit(x[ok], y[ok], 1, w=None if w is None else np.sqrt(w[ok]))
    return p[0], p[1]
