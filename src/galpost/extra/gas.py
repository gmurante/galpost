# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Gas analyses: spherical flows and mass loading, resolved Kennicutt-Schmidt
relation, metallicity gradients, dust and circumgalactic medium.
"""

from __future__ import annotations

import math

import numpy as np

from ..units import read_gadget_params
from .common import (ATOMIC_WEIGHT, KMS_KPC_TO_PER_YR, GalaxyView, binned, linear_fit,
                     profile_edges, solar_log_ratio, wmean, write_table)


def _vr(view):
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.einsum("ij,ij->i", view.pos, view.vel) / view.r
    return np.nan_to_num(v)


# --------------------------------------------------------------------------
def spherical_flows(names, params, view: GalaxyView, sfr_gal: float, metal_names):
    g = view.g
    fr = [float(x) for x in str(params["FlowRadii"]).split(",") if x.strip()]
    wfrac = float(params["FlowShellWidth"])
    vr = _vr(view)
    gas = view.of(0)
    mcld, mhot = view.gas_array("mcld"), view.gas_array("mhot")
    met = view.metals()
    mz = np.nansum(met[:, 1:], axis=1) if met is not None else np.zeros(len(vr))
    zfrac = np.where(gas, mz / np.maximum(view.reference_mass(), 1e-300), 0.0) if met is not None else 0.0
    comps = [("gas", gas, view.mass), ("cold gas", gas, mcld), ("hot gas", gas, mhot),
             ("stars", view.of(4), view.mass), ("DM", view.of(1), view.mass),
             ("gas metals", gas, view.mass * zfrac)]
    rows = []
    for f in fr:
        r0 = f * g.rvir
        dr = wfrac * r0
        shell = np.abs(view.r - r0) < 0.5 * dr
        row = [f, r0, dr]
        for _, msk, w in comps:
            s = shell & msk
            flux = w[s] * vr[s] / dr * KMS_KPC_TO_PER_YR          # Msun/yr
            row += [flux[flux > 0].sum(), -flux[flux < 0].sum(), flux.sum()]
        out_gas = row[3]
        row.append(out_gas / sfr_gal if sfr_gal > 0 else np.nan)
        rows.append(row)
    cols = ["r / rvir", "shell radius (kpc)", "shell thickness (kpc)"]
    for lab, _, _ in comps:
        cols += [f"{lab}: outflow rate (Msun/yr)", f"{lab}: inflow rate (Msun/yr)",
                 f"{lab}: net rate (Msun/yr, >0 outwards)"]
    cols.append("mass loading eta = gas outflow / SFR(rgal)")
    write_table(names.per_gal("flows", g.gal_id), cols, rows,
                header_lines=[f"Mass flows through spherical shells, galaxy {g.gal_id}",
                              f"rvir = {g.rvir:f} kpc, SFR within rgal = {sfr_gal:.5e} Msun/yr",
                              "rate = sum(m v_r) / thickness"])


# --------------------------------------------------------------------------
def ks_relation(names, params, view: GalaxyView):
    g = view.g
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    gas = view.of(0) & (np.abs(view.z) < zmax)
    R = view.R[gas]
    area = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2) * 1.0e6    # pc^2
    area_kpc = area / 1.0e6
    mtot, _, cnt = binned(R, edges, view.mass[gas])
    mc, _, _ = binned(R, edges, view.gas_array("mcld")[gas])
    mh2, _, _ = binned(R, edges, view.gas_array("mmol")[gas])
    sfr, _, _ = binned(R, edges, view.gas_array("sfr")[gas])
    with np.errstate(divide="ignore", invalid="ignore"):
        tdep = np.where(sfr > 0, mc / sfr / 1e9, np.nan)
        tdep2 = np.where(sfr > 0, mh2 / sfr / 1e9, np.nan)
    write_table(names.per_gal("ks", g.gal_id),
                ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)", "Sigma_gas total (Msun/pc^2)",
                 "Sigma_cold (Msun/pc^2)", "Sigma_H2 (Msun/pc^2)", "Sigma_HI+ = cold - H2 (Msun/pc^2)",
                 "Sigma_SFR (Msun/yr/kpc^2)", "depletion time of cold gas (Gyr)",
                 "depletion time of H2 (Gyr)", "number of gas particles"],
                np.column_stack([edges[:-1], edges[1:], 0.5 * (edges[1:] + edges[:-1]),
                                 mtot / area, mc / area, mh2 / area, (mc - mh2) / area,
                                 sfr / area_kpc, tdep, tdep2, cnt]),
                header_lines=[f"Resolved Kennicutt-Schmidt relation of galaxy {g.gal_id}, "
                              f"|z| < {zmax:g} kpc",
                              "H2 from the MUPPI pressure law (derived block MMOL)"])
    inr = view.of(0) & (view.r <= g.rgal)
    return dict(mcold=view.gas_array("mcld")[inr].sum(), mh2=view.gas_array("mmol")[inr].sum(),
                mgas=view.mass[inr].sum())


# --------------------------------------------------------------------------
def metallicity_gradients(names, params, view: GalaxyView, metal_names, r_e: float):
    g = view.g
    met = view.metals()
    if met is None or "O" not in metal_names:
        return {}
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    rmid = 0.5 * (edges[1:] + edges[:-1])
    mh = view.hydrogen_mass()
    o = met[:, metal_names.index("O")]
    with np.errstate(divide="ignore", invalid="ignore"):
        oh = 12.0 + np.log10((o / ATOMIC_WEIGHT["O"]) / (mh / ATOMIC_WEIGHT["H"]))
        feh = (np.log10(met[:, metal_names.index("Fe")] / mh) - solar_log_ratio("Fe")
               if "Fe" in metal_names else np.full(len(o), np.nan))
    gas = view.of(0) & (np.abs(view.z) < zmax) & np.isfinite(oh)
    st = view.of(4) & (np.abs(view.z) < zmax) & np.isfinite(feh)
    wc = view.gas_array("mcld")
    sw, s, _ = binned(view.R[gas], edges, wc[gas], oh[gas])
    oh_c = wmean(sw, s)
    sws, s2, _ = binned(view.R[gas], edges, view.gas_array("sfr")[gas], oh[gas])
    oh_sfr = wmean(sws, s2)
    sw3, s3, _ = binned(view.R[st], edges, view.mass[st], feh[st])
    feh_s = wmean(sw3, s3)

    def opt(tag, default):
        x = str(params[tag])
        return default if x == "auto" else float(x)
    r1, r2 = opt("GradientRmin", 0.5 * r_e), opt("GradientRmax", 2.0 * r_e)
    fit = (rmid >= r1) & (rmid <= r2)
    s_g, _ = linear_fit(rmid[fit], oh_c[fit])
    s_s, _ = linear_fit(rmid[fit], feh_s[fit])
    write_table(names.per_gal("gradients", g.gal_id),
                ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)", "12+log(O/H) cold gas (cold-mass weighted)",
                 "12+log(O/H) gas (SFR weighted)", "[Fe/H] stars (mass weighted)"],
                np.column_stack([edges[:-1], edges[1:], rmid, oh_c, oh_sfr, feh_s]),
                header_lines=[f"Metallicity gradients of galaxy {g.gal_id}, |z| < {zmax:g} kpc",
                              f"fit range {r1:.3f}-{r2:.3f} kpc (R_e = {r_e:.4f} kpc)",
                              f"gas slope = {s_g:.6f} dex/kpc = {s_g * r_e:.6f} dex/R_e",
                              f"stellar slope = {s_s:.6f} dex/kpc = {s_s * r_e:.6f} dex/R_e",
                              "Solar abundances: Asplund et al. (2009)"])
    i_e = np.argmin(np.abs(rmid - r_e))
    return dict(oh_re=oh_c[i_e], grad_gas=s_g, grad_star=s_s)


# --------------------------------------------------------------------------
def dust(names, params, view: GalaxyView):
    g = view.g
    d = view.d
    if not (d.has("DSTL") or d.has("DSTS")):
        return {}
    gas = view.of(0)
    gi = view.idx[gas] - d.nstart[0]
    md = {}
    for b in ("DSTL", "DSTS"):
        if d.has(b):
            raw = d.raw_block(b, 0)[gi].astype(np.float64)
            md[b] = (raw.sum(axis=1) if raw.ndim > 1 else raw) * d.units.mass_msun
        else:
            md[b] = np.zeros(gas.sum())
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    R, z = view.R[gas], view.z[gas]
    sl = np.abs(z) < zmax
    met = view.metals()
    mz = np.nansum(met[gas][:, 1:], axis=1) if met is not None else np.zeros(gas.sum())
    mg = view.mass[gas]
    mL, _, _ = binned(R[sl], edges, md["DSTL"][sl])
    mS, _, _ = binned(R[sl], edges, md["DSTS"][sl])
    mG, _, _ = binned(R[sl], edges, mg[sl])
    mZ, _, _ = binned(R[sl], edges, mz[sl])
    area = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2) * 1.0e6
    with np.errstate(divide="ignore", invalid="ignore"):
        dg = (mL + mS) / mG
        dz = (mL + mS) / mZ
    inr = view.r[gas] <= g.rgal
    tot = dict(dust_large=md["DSTL"][inr].sum(), dust_small=md["DSTS"][inr].sum(),
               gas=mg[inr].sum(), metals=mz[inr].sum())
    write_table(names.per_gal("dust", g.gal_id),
                ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)", "Sigma dust large grains (Msun/pc^2)",
                 "Sigma dust small grains (Msun/pc^2)", "dust-to-gas ratio",
                 "dust / gas-phase metals (Zs block)"],
                np.column_stack([edges[:-1], edges[1:], 0.5 * (edges[1:] + edges[:-1]),
                                 mL / area, mS / area, dg, dz]),
                header_lines=[f"Dust of galaxy {g.gal_id} (blocks DSTL/DSTS summed over species), "
                              f"|z| < {zmax:g} kpc",
                              "within rgal: M_dust,large = %.5e  M_dust,small = %.5e  M_gas = %.5e  "
                              "M_metals,gas = %.5e Msun" % (tot["dust_large"], tot["dust_small"],
                                                           tot["gas"], tot["metals"])])
    return tot


# --------------------------------------------------------------------------
def hot_halo(names, params, view: GalaxyView, gadget_paramfile: str):
    g = view.g
    rin = float(params["CGMRmin"]) * g.rvir
    shell = (view.r > rin) & (view.r <= g.rvir)
    gas = view.of(0) & shell
    mg = view.mass[gas]
    mhot = view.gas_array("mhot")[gas]
    mcld = view.gas_array("mcld")[gas]
    temp = view.gas_array("temp")[gas]
    thot = np.where(view.d.mp[view.idx[gas] - view.d.nstart[0]] if gas.any() else np.zeros(0, bool),
                    view.gas_array("thot")[gas], temp)
    met = view.metals()
    zg = (np.nansum(met[gas][:, 1:], axis=1) / mg) if met is not None and gas.any() else np.full(gas.sum(), np.nan)
    hotsel = mhot > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        t_hot = np.average(thot[hotsel], weights=mhot[hotsel]) if hotsel.any() else np.nan
        z_hot = np.average(zg[hotsel], weights=mhot[hotsel]) if hotsel.any() else np.nan
    inv = view.r <= g.rvir
    mb = view.mass[inv & (view.of(0) | view.of(4) | view.of(5))].sum()
    mt = view.mass[inv].sum()
    try:
        ob = float(read_gadget_params(gadget_paramfile).get("OmegaBaryon", "nan"))
    except OSError:
        ob = float("nan")
    fcos = ob / view.d.header.omega0
    vals = [rin, g.rvir, mg.sum(), mhot.sum(), mcld.sum(), t_hot, z_hot,
            view.mass[view.of(4) & shell].sum(), mb, mt, mb / mt if mt else np.nan, fcos]
    write_table(names.per_gal("cgm", g.gal_id),
                ["inner radius (kpc)", "rvir (kpc)", "gas mass (Msun)", "hot phase mass (Msun)",
                 "cold phase mass (Msun)", "hot-mass-weighted temperature (K)",
                 "hot-mass-weighted metallicity Z", "stellar mass in the shell (Msun)",
                 "baryon mass within rvir (Msun)", "total mass within rvir (Msun)",
                 "baryon fraction within rvir", "cosmic baryon fraction Omega_b/Omega_m"],
                [vals],
                header_lines=[f"Circumgalactic medium of galaxy {g.gal_id}: "
                              f"{params['CGMRmin']:g} rvir < r < rvir",
                              "hot phase = hot gas of single-phase particles + hot phase of MUPPI particles"])
    return dict(m_cgm=mg.sum(), m_cgm_hot=mhot.sum(), t_cgm=t_hot, fbar=vals[10])
