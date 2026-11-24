# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Global quantities: central black hole, scaling relations, satellites.
"""

from __future__ import annotations

import math
import os

import numpy as np

from .common import G, GalaxyView, write_table
from .structure import half_mass_radius_2d

YR_S = 3.15576e7


def central_bh(names, params, view: GalaxyView, mbulge: float):
    g = view.g
    d = view.d
    bh = view.of(5) & (view.r <= g.rgal)
    rows = []
    out = dict(n_bh=int(bh.sum()), mbh=np.nan, mdot_bh=np.nan, dbh=np.nan)
    if bh.any():
        bi = view.idx[bh] - d.nstart[5]
        u = d.units
        mbh = (d.raw_block("BHMA", 5)[bi].astype(np.float64) * u.mass_msun
               if d.has("BHMA") else view.mass[bh])
        mdot = (d.raw_block("BHMD", 5)[bi].astype(np.float64)
                * (u.unit_mass_g / 1.989e33) / (u.unit_length_cm / u.unit_velocity_cms) * YR_S
                if d.has("BHMD") else np.full(bh.sum(), np.nan))
        dist = view.r[bh]
        ids = d.ids[view.idx[bh]]
        o = np.argsort(-mbh)
        for k in o:
            rows.append([ids[k], mbh[k], mdot[k], dist[k], view.mass[bh][k]])
        k = o[0]
        out.update(mbh=mbh[k], mdot_bh=mdot[k], dbh=dist[k])
    write_table(names.per_gal("bh", g.gal_id),
                ["particle ID", "BH mass (BHMA, Msun)", "accretion rate (BHMD, Msun/yr)",
                 "distance from the galaxy centre (kpc)", "dynamical mass (MASS, Msun)"],
                rows if rows else np.zeros((0, 5)),
                header_lines=[f"Black holes within rgal = {g.rgal:f} kpc of galaxy {g.gal_id}, "
                              "most massive first",
                              "M_BH / M_bulge = %.5e (bulge from the EG/kinematic decomposition)"
                              % (out["mbh"] / mbulge if mbulge and mbulge > 0 else float("nan"))],
                fmt="%14.6e")
    return out


def basic_quantities(view: GalaxyView, resolver):
    """Global quantities within rgal used by the scaling relations."""
    g = view.g
    inr = view.r <= g.rgal
    st, gas = view.of(4) & inr, view.of(0) & inr
    ms = view.mass[st].sum()
    o = np.argsort(view.r[view.r <= g.rvir])
    m_in = view.mass[view.r <= g.rvir][o]
    r_in = view.r[view.r <= g.rvir][o]
    cm = np.cumsum(m_in)
    with np.errstate(divide="ignore", invalid="ignore"):
        vc = np.sqrt(G * cm / np.maximum(r_in, 1e-3))
    sfr = view.gas_array("sfr")[gas].sum()
    return dict(mstar=ms, mgas=view.mass[gas].sum(), mcold=view.gas_array("mcld")[gas].sum(),
                mh2=view.gas_array("mmol")[gas].sum(), sfr=sfr,
                ssfr=sfr / ms if ms > 0 else np.nan, mvir=cm[-1] if len(cm) else 0.0,
                vmax=np.nanmax(vc[r_in > 1.0]) if np.any(r_in > 1.0) else np.nan,
                r_half_2d=half_mass_radius_2d(view.R[st], view.mass[st]),
                _rcurve=(r_in, vc))


SCALING_COLUMNS = [
    ("gal_id", "galaxy ID"), ("snap", "snapshot number"), ("redshift", "redshift"),
    ("mstar", "stellar mass within rgal (Msun)"), ("mgas", "gas mass within rgal (Msun)"),
    ("mcold", "cold gas mass within rgal (Msun)"), ("mh2", "H2 mass within rgal (Msun)"),
    ("sfr", "SFR within rgal (Msun/yr)"), ("ssfr", "specific SFR (1/yr)"),
    ("mvir", "total mass within rvir (Msun)"), ("rvir", "rvir (kpc)"), ("rgal", "rgal (kpc)"),
    ("r_half_2d", "projected stellar half-mass radius (kpc)"),
    ("r_half_3d", "3D stellar half-mass radius (kpc)"), ("rd", "disc scale length (kpc)"),
    ("v22", "circular velocity at 2.2 R_d (km/s)"), ("vmax", "maximum circular velocity (km/s)"),
    ("bt_eg", "B/T of the EG decomposition"), ("dt_e", "D/T from eps(E) > threshold"),
    ("bt_e", "B/T = 2 x counter-rotating fraction (eps(E))"), ("kappa_rot", "kappa_rot"),
    ("kappa_co", "kappa_co"), ("zstar", "mass-weighted stellar metallicity Z"),
    ("feh_star", "mass-weighted stellar [Fe/H]"), ("age_star", "mass-weighted stellar age (Gyr)"),
    ("oh_re", "12+log(O/H) of cold gas at R_e"), ("grad_gas", "gas O/H gradient (dex/kpc)"),
    ("mbh", "mass of the most massive BH within rgal (Msun)"),
    ("mdot_bh", "its accretion rate (Msun/yr)"), ("a2max", "bar strength A2max"),
    ("m_cgm_hot", "hot CGM mass (Msun)"), ("fbar", "baryon fraction within rvir"),
    ("lowres", "low-resolution mass fraction within rvir"),
]


def scaling_row(res: dict):
    return [res.get(k, np.nan) for k, _ in SCALING_COLUMNS]


def write_scaling(names, rows, tables: bool, gal_ids):
    cols = [d for _, d in SCALING_COLUMNS]
    write_table(names.per_snap("scaling"), cols, rows,
                header_lines=["Global properties of the analysed galaxies (galpost)",
                              "NaN: quantity not computed (enable the corresponding analysis)"])
    if tables:
        for row, gid in zip(rows, gal_ids):
            path = names.table("table-scaling", gid)
            new = not os.path.exists(path)
            with open(path, "a") as f:
                if new:
                    for i, c in enumerate(cols):
                        f.write("# %2d: %s\n" % (i + 1, c))
                f.write(" ".join("%14.6e" % v for v in row) + "\n")


def satellites(names, params, view: GalaxyView, catalog):
    g = view.g
    d = view.d
    u = d.units
    if catalog is None:
        return None
    spos = catalog.subhalo("SPOS").astype(np.float64) * u.length_kpc
    dp = d.pbc_diff(spos - g.com)
    dist = np.sqrt(np.einsum("ij,ij->i", dp, dp))
    own = g.subhalo
    if own < 0:
        own = int(np.argmin(dist)) if dist.min() < 0.1 * g.rgal + 1.0 else -1
    sel = np.nonzero((dist <= g.rvir) & (np.arange(len(dist)) != own))[0]
    smst = catalog.subhalo("SMST").astype(np.float64) * u.mass_msun
    msub = catalog.subhalo("MSUB").astype(np.float64) * u.mass_msun
    vmax = catalog.subhalo("VMAX").astype(np.float64) if "VMAX" in catalog.snap.blocks else np.full(len(dist), np.nan)
    rhms = (catalog.subhalo("RHMS").astype(np.float64) * u.length_kpc
            if "RHMS" in catalog.snap.blocks else np.full(len(dist), np.nan))
    grnr = catalog.subhalo("GRNR")
    first = catalog.group("FSUB")
    o = sel[np.argsort(-smst[sel, 4])]
    rows = [[i, grnr[i], int(first[grnr[i]] == i), dist[i], smst[i, 4], smst[i, 0], smst[i, 1],
             smst[i, 2] + smst[i, 3], msub[i], vmax[i], rhms[i]] for i in o]
    write_table(names.per_gal("satellites", g.gal_id),
                ["subhalo index", "FoF group", "1 if central of its group", "distance (kpc)",
                 "stellar mass (Msun)", "gas mass (Msun)", "DM mass (Msun)",
                 "low-resolution mass (Msun)", "total subhalo mass (Msun)",
                 "Vmax (km/s, SubFind)", "half-mass radius (kpc, SubFind)"],
                rows if rows else np.zeros((0, 11)),
                header_lines=[f"SubFind subhaloes within rvir = {g.rvir:f} kpc of galaxy {g.gal_id} "
                              f"(own subhalo {own} excluded), most massive in stars first",
                              f"catalogue: {catalog.group_dir}"])
    return dict(n_sat=len(rows), m_sat_star=float(smst[sel, 4].sum()))
