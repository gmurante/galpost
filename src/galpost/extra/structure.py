# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Stellar structure and kinematics: binding-energy circularities, kappa_rot,
stellar populations, disc structure, kinematic profiles, Toomre Q, bar
strength and maps.
"""

from __future__ import annotations

import math
import os

import numpy as np

from .common import (G, GalaxyView, binned, linear_fit, profile_edges, radius_value,
                     solar_log_ratio, wmean, write_table)

NB_EPS = 80          # circularity histogram bins in [-2, 2]
NB_EPS_R = 20        # radial bins of the 2D histogram


# --------------------------------------------------------------------------
# spherical potential and jcirc(E)
# --------------------------------------------------------------------------
def spherical_potential(r, m, rgrid):
    """Potential of a spherical mass distribution (particles r, m) on rgrid,
    and the enclosed mass. Mass beyond max(r) is ignored (constant shift)."""
    o = np.argsort(r)
    rs, ms = r[o], m[o]
    cm = np.cumsum(ms)
    # sum_{r_j > R} m_j / r_j
    inv = np.cumsum((ms / np.maximum(rs, 1e-6))[::-1])[::-1]
    k = np.searchsorted(rs, rgrid, side="right")
    menc = np.where(k > 0, cm[np.maximum(k - 1, 0)], 0.0)
    outer = np.where(k < len(rs), inv[np.minimum(k, len(rs) - 1)], 0.0)
    phi = -G * (menc / rgrid + outer)
    return phi, menc


def circularity_energy(view: GalaxyView, rpot: float):
    """eps = jz / jcirc(E) for the stars of the view (Abadi et al. 2003), with a
    spherically averaged potential from all particles within rpot."""
    sel_pot = view.r <= rpot
    rgrid = np.logspace(math.log10(max(1e-3, 1e-4 * rpot)), math.log10(rpot), 2000)
    phi_g, menc_g = spherical_potential(view.r[sel_pot], view.mass[sel_pot], rgrid)
    vc2 = G * menc_g / rgrid
    ec = 0.5 * vc2 + phi_g
    jc = rgrid * np.sqrt(vc2)
    # make E_c(r) monotonic for the inversion
    ec_m = np.maximum.accumulate(ec)
    st = view.of(4)
    rs = view.r[st]
    phi_s = np.interp(rs, rgrid, phi_g)
    v = view.vel[st]
    e = 0.5 * np.einsum("ij,ij->i", v, v) + phi_s
    jcirc = np.interp(e, ec_m, jc, left=jc[0], right=jc[-1])
    p = view.pos[st]
    jz = p[:, 0] * v[:, 1] - p[:, 1] * v[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        eps = jz / jcirc
    return eps, st, e


def kinematic_decomposition(names, params, view: GalaxyView, rgal: float):
    eps, st, energy = circularity_energy(view, view.g.rvir)
    m = view.mass[st]
    r = view.r[st]
    inr = r <= rgal
    thr = float(params["DiscCircularity"])
    me, mm = eps[inr], m[inr]
    mtot = mm.sum()
    m_disc = mm[me > thr].sum()
    m_counter = mm[me < 0].sum()
    m_sph2 = min(2.0 * m_counter, mtot)
    res = dict(mstar=mtot, mdisc_thr=m_disc, dt_thr=m_disc / mtot if mtot else np.nan,
               mspheroid_2counter=m_sph2, bt_2counter=m_sph2 / mtot if mtot else np.nan,
               mthin=mm[me > 0.8].sum(), mthick=mm[(me > 0.5) & (me <= 0.8)].sum(),
               eps_mean=np.average(me, weights=mm) if mtot else np.nan)
    g = view.g
    write_table(names.per_gal("kindec", g.gal_id),
                ["stellar mass within rgal (Msun)",
                 f"mass with eps > {thr:g} (Msun)", f"D/T with eps > {thr:g}",
                 "spheroid mass = 2 x counter-rotating mass (Msun)", "B/T = 2 x counter-rotating fraction",
                 "thin disc mass, eps > 0.8 (Msun)", "thick/warm disc mass, 0.5 < eps <= 0.8 (Msun)",
                 "mass-weighted mean eps"],
                [[res[k] for k in ("mstar", "mdisc_thr", "dt_thr", "mspheroid_2counter",
                                   "bt_2counter", "mthin", "mthick", "eps_mean")]],
                header_lines=[f"Kinematic decomposition of galaxy {g.gal_id}: eps = jz/jcirc(E), "
                              "spherically averaged potential within rvir",
                              f"rgal = {rgal:f} kpc"])
    # histograms (same layout as jhist.*: 1D + 2D in radial bins)
    eb = np.linspace(-2.0, 2.0, NB_EPS + 1)
    rb = np.linspace(0.0, rgal, NB_EPS_R + 1)
    ec = np.clip(eps, -2.0 + 1e-9, 2.0 - 1e-9)
    h1, _ = np.histogram(ec[inr], bins=eb, weights=m[inr])
    h2, _, _ = np.histogram2d(ec[inr], r[inr], bins=[eb, rb], weights=m[inr])
    cols = ["eps = jz/jcirc(E) (bin centre)", "stellar mass in bin (Msun)"] + \
           ["mass in radial bin %d (centre %.3f kpc)" % (k + 1, 0.5 * (rb[k] + rb[k + 1]))
            for k in range(NB_EPS_R)]
    write_table(names.per_gal("jhistE", g.gal_id), cols,
                np.column_stack([0.5 * (eb[1:] + eb[:-1]), h1, h2]))
    return res, eps, st


def kappa_rot(view: GalaxyView, radius: float):
    st = view.of(4) & (view.r <= radius)
    m = view.mass[st]
    p, v = view.pos[st], view.vel[st]
    R = view.R[st]
    jz = p[:, 0] * v[:, 1] - p[:, 1] * v[:, 0]
    ke = np.sum(m * np.einsum("ij,ij->i", v, v))
    with np.errstate(divide="ignore", invalid="ignore"):
        krot = np.where(R > 0, m * (jz / R) ** 2, 0.0)
    kr = krot.sum() / ke if ke > 0 else np.nan
    kco = krot[jz > 0].sum() / ke if ke > 0 else np.nan
    return kr, kco


# --------------------------------------------------------------------------
def stellar_populations(names, params, view: GalaxyView, metal_names):
    g = view.g
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    st = view.of(4) & (np.abs(view.z) < zmax)
    m = view.mass[st]
    R = view.R[st]
    cols, data = ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)", "stellar mass in annulus (Msun)"], []
    sw, _, cnt = binned(R, edges, m)
    data += [edges[:-1], edges[1:], 0.5 * (edges[1:] + edges[:-1]), sw]
    age = view.stellar_age()
    if age is not None:
        a = age[st]
        _, swa, _ = binned(R, edges, m, a)
        data.append(wmean(sw, swa))
        cols.append("mass-weighted age (Gyr)")
        data.append(cnt.astype(float))
        cols.append("number of star particles")
    feh = ofe = mgfe = zt = None
    met = view.metals()
    if met is not None:
        ref = view.reference_mass()[st]
        mh = view.hydrogen_mass()[st]
        mz = np.nansum(met[st][:, 1:], axis=1)          # all metals except He
        zt = mz / ref
        _, swz, _ = binned(R, edges, m, zt)
        data.append(wmean(sw, swz))
        cols.append("mass-weighted total metallicity Z")
        with np.errstate(divide="ignore", invalid="ignore"):
            if "Fe" in metal_names:
                fe = met[st][:, metal_names.index("Fe")]
                feh = np.log10(fe / mh) - solar_log_ratio("Fe")
                _, s, _ = binned(R, edges, m, np.where(np.isfinite(feh), feh, np.nan))
                data.append(wmean(sw, s)); cols.append("mass-weighted [Fe/H]")
                for el in ("O", "Mg"):
                    if el in metal_names:
                        x = met[st][:, metal_names.index(el)]
                        xfe = (np.log10(x / fe) - (solar_log_ratio(el) - solar_log_ratio("Fe")))
                        if el == "O":
                            ofe = xfe
                        else:
                            mgfe = xfe
                        _, s, _ = binned(R, edges, m, xfe)
                        data.append(wmean(sw, s)); cols.append(f"mass-weighted [{el}/Fe]")
    write_table(names.per_gal("stellarpop", g.gal_id), cols, np.column_stack(data),
                header_lines=[f"Stellar populations of galaxy {g.gal_id} in annuli of the disc plane, "
                              f"|z| < {zmax:g} kpc", "Solar abundances: Asplund et al. (2009)"])
    rg = view.r[st] <= g.rgal
    if feh is not None:
        fb = np.linspace(-3.0, 1.0, 81)
        ok = rg & np.isfinite(feh)
        h, _ = np.histogram(feh[ok], bins=fb, weights=m[ok])
        write_table(names.per_gal("mdf", g.gal_id),
                    ["[Fe/H] (bin centre)", "stellar mass (Msun)", "normalised dM/d[Fe/H]"],
                    np.column_stack([0.5 * (fb[1:] + fb[:-1]), h,
                                     h / max(h.sum(), 1e-300) / np.diff(fb)]),
                    header_lines=[f"Metallicity distribution of stars within rgal = {g.rgal:f} kpc"])
        if ofe is not None:
            ob = np.linspace(-0.5, 1.0, 61)
            ok2 = ok & np.isfinite(ofe)
            h2, _, _ = np.histogram2d(feh[ok2], ofe[ok2], bins=[fb, ob], weights=m[ok2])
            with open(names.per_gal("ofe_feh", g.gal_id), "w") as f:
                f.write("# Stellar mass (Msun) in the [O/Fe]-[Fe/H] plane, stars within rgal\n")
                f.write("# rows: [Fe/H] bin centres from %g to %g (%d bins); columns: [O/Fe] "
                        "bin centres from %g to %g (%d bins)\n"
                        % (0.5 * (fb[0] + fb[1]), 0.5 * (fb[-2] + fb[-1]), len(fb) - 1,
                           0.5 * (ob[0] + ob[1]), 0.5 * (ob[-2] + ob[-1]), len(ob) - 1))
                for row in h2:
                    f.write(" ".join("%11.4e" % x for x in row) + "\n")
    return dict(zstar=np.average(zt[rg], weights=m[rg]) if zt is not None and rg.any() else np.nan,
                feh_star=np.average(np.where(np.isfinite(feh[rg]), feh[rg], 0), weights=m[rg] * np.isfinite(feh[rg]))
                if feh is not None and (rg & np.isfinite(feh)).any() else np.nan,
                age_star=np.average(age[st][rg], weights=m[rg]) if age is not None and rg.any() else np.nan)


# --------------------------------------------------------------------------
def half_mass_radius_2d(R, m):
    if not len(R):
        return np.nan
    o = np.argsort(R)
    cm = np.cumsum(m[o])
    return float(np.interp(0.5 * cm[-1], cm, R[o]))


def disc_structure(names, params, view: GalaxyView):
    g = view.g
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    st = view.of(4) & (np.abs(view.z) < zmax)
    m, R, z = view.mass[st], view.R[st], view.z[st]
    ing = view.of(4) & (view.r <= g.rgal)
    r_hm2d = half_mass_radius_2d(view.R[ing], view.mass[ing])
    o = np.argsort(view.r[ing])
    cm = np.cumsum(view.mass[ing][o])
    r_hm3d = float(np.interp(0.5 * cm[-1], cm, view.r[ing][o])) if len(cm) else np.nan
    sw, _, cnt = binned(R, edges, m)
    area = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2) * 1.0e6       # pc^2
    sigma = sw / area
    rmid = 0.5 * (edges[1:] + edges[:-1])

    def opt(tag, default):
        s = str(params[tag])
        return default if s == "auto" else float(s)
    r1, r2 = opt("DiscFitRmin", r_hm2d), opt("DiscFitRmax", 3.0 * r_hm2d)
    ok = (rmid >= r1) & (rmid <= r2) & (sigma > 0)
    slope, icpt = linear_fit(rmid[ok], np.log(sigma[ok]))
    rd = -1.0 / slope if slope < 0 else np.nan
    # vertical structure per annulus
    zmed, zrms, hz = (np.full(len(rmid), np.nan) for _ in range(3))
    ib = np.searchsorted(edges, R, side="right") - 1
    for k in range(len(rmid)):
        s = ib == k
        if s.sum() < 20:
            continue
        az = np.abs(z[s])
        w = m[s]
        o = np.argsort(az)
        c = np.cumsum(w[o])
        zmed[k] = np.interp(0.5 * c[-1], c, az[o])
        zrms[k] = math.sqrt(np.average(z[s] ** 2, weights=w))
        hb = np.linspace(0, min(zmax, 4 * zmed[k] + 1e-3), 13)
        hh, _ = np.histogram(az, bins=hb, weights=w)
        zc = 0.5 * (hb[1:] + hb[:-1])
        good = hh > 0
        sl, _ = linear_fit(zc[good], np.log(hh[good]))
        hz[k] = -1.0 / sl if sl < 0 else np.nan
    write_table(names.per_gal("disc", g.gal_id),
                ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)", "stellar surface density (Msun/pc^2)",
                 "number of star particles", "median |z| (kpc)", "rms z (kpc)",
                 "exponential scale height (kpc)"],
                np.column_stack([edges[:-1], edges[1:], rmid, sigma, cnt, zmed, zrms, hz]),
                header_lines=[f"Stellar disc structure of galaxy {g.gal_id}, |z| < {zmax:g} kpc",
                              f"projected half-mass radius (stars within rgal) = {r_hm2d:.5f} kpc",
                              f"3D half-mass radius (stars within rgal) = {r_hm3d:.5f} kpc",
                              f"exponential scale length R_d = {rd:.5f} kpc (fit {r1:.3f}-{r2:.3f} kpc), "
                              f"central surface density = {math.exp(icpt) if icpt == icpt else float('nan'):.5e} Msun/pc^2"])
    return dict(r_half_2d=r_hm2d, r_half_3d=r_hm3d, rd=rd)


# --------------------------------------------------------------------------
def _cyl_velocities(view, mask):
    p, v = view.pos[mask], view.vel[mask]
    R = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        vR = (p[:, 0] * v[:, 0] + p[:, 1] * v[:, 1]) / R
        vphi = (p[:, 0] * v[:, 1] - p[:, 1] * v[:, 0]) / R
    return R, np.nan_to_num(vR), np.nan_to_num(vphi), v[:, 2]


def kinematic_stats(R, vR, vphi, vz, w, edges):
    sw, _, cnt = binned(R, edges, w)
    out = {"n": cnt}
    for name, x in (("vR", vR), ("vphi", vphi), ("vz", vz)):
        _, s1, _ = binned(R, edges, w, x)
        _, s2, _ = binned(R, edges, w, x * x)
        mu = wmean(sw, s1)
        out["mean_" + name] = mu
        out["sig_" + name] = np.sqrt(np.maximum(wmean(sw, s2) - mu * mu, 0.0))
    s3 = np.sqrt((out["sig_vR"] ** 2 + out["sig_vphi"] ** 2 + out["sig_vz"] ** 2) / 3.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["v_over_sigma"] = out["mean_vphi"] / s3
    out["mass"] = sw
    return out


def kinematic_profiles(names, params, view: GalaxyView, kappa=None, kappa_radius=30.0):
    g = view.g
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    slab = np.abs(view.z) < zmax
    st = view.of(4) & slab
    gas = view.of(0) & slab
    res = {}
    cols = ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)"]
    data = [edges[:-1], edges[1:], 0.5 * (edges[1:] + edges[:-1])]
    for lab, mask, w in (("stars", st, view.mass[st]),
                         ("cold gas", gas, view.gas_array("mcld")[gas])):
        k = kinematic_stats(*_cyl_velocities(view, mask), w, edges)
        res[lab] = k
        for key, desc in (("mean_vphi", "mean v_phi (km/s)"), ("sig_vR", "sigma_R (km/s)"),
                          ("sig_vphi", "sigma_phi (km/s)"), ("sig_vz", "sigma_z (km/s)"),
                          ("v_over_sigma", "V/sigma"), ("n", "number of particles")):
            data.append(k[key])
            cols.append(f"{lab}: {desc}")
    hdr = [f"Kinematic profiles of galaxy {g.gal_id} (cylindrical, galaxy frame), |z| < {zmax:g} kpc",
           "cold gas weighted by the cold phase mass"]
    if kappa is not None:
        hdr.append("kappa_rot = %.5f  kappa_co = %.5f  (stars within %g kpc)"
                   % (kappa[0], kappa[1], kappa_radius))
    write_table(names.per_gal("kinematics", g.gal_id), cols, np.column_stack(data),
                header_lines=hdr)
    return res


# --------------------------------------------------------------------------
def toomre_q(names, params, view: GalaxyView):
    g = view.g
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    rmid = 0.5 * (edges[1:] + edges[:-1])
    area = math.pi * (edges[1:] ** 2 - edges[:-1] ** 2)      # kpc^2
    slab = np.abs(view.z) < zmax
    # rotation curve and epicyclic frequency from the spherical enclosed mass
    o = np.argsort(view.r)
    cm = np.cumsum(view.mass[o])
    rr = np.maximum(rmid, 1e-3)
    menc = np.interp(rr, view.r[o], cm)
    vc = np.sqrt(G * menc / rr)
    omega = vc / rr
    dvc = np.gradient(vc, rr)
    kappa2 = 2.0 * omega * (omega + dvc)                      # (km/s/kpc)^2
    kappa = np.sqrt(np.maximum(kappa2, 0.0))
    st = view.of(4) & slab
    gas = view.of(0) & slab
    ks = kinematic_stats(*_cyl_velocities(view, st), view.mass[st], edges)
    kg = kinematic_stats(*_cyl_velocities(view, gas), view.gas_array("mcld")[gas], edges)
    sig_s = ks["mass"] / area
    sig_g = kg["mass"] / area
    cs = 0.0       # MUPPI cold phase (T_c = 300 K): thermal support negligible
    sg = np.sqrt(kg["sig_vz"] ** 2 + cs ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        q_g = sg * kappa / (math.pi * G * sig_g)
        q_s = ks["sig_vR"] * kappa / (3.36 * G * sig_s)
        q_2 = 1.0 / (1.0 / q_g + 1.0 / q_s)                   # Wang & Silk (1994)
    write_table(names.per_gal("toomre", g.gal_id),
                ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)", "circular velocity (km/s)",
                 "epicyclic frequency kappa (km/s/kpc)", "Sigma_cold gas (Msun/kpc^2)",
                 "Sigma_stars (Msun/kpc^2)", "sigma_z cold gas (km/s)", "sigma_R stars (km/s)",
                 "Q_gas", "Q_stars", "Q_two-component (Wang & Silk 1994)"],
                np.column_stack([edges[:-1], edges[1:], rmid, vc, kappa, sig_g, sig_s, sg,
                                 ks["sig_vR"], q_g, q_s, q_2]),
                header_lines=[f"Toomre stability of galaxy {g.gal_id}, |z| < {zmax:g} kpc",
                              "v_c and kappa from the spherically enclosed mass; gas dispersion = "
                              "sigma_z of the cold phase"])


# --------------------------------------------------------------------------
def bar_strength(names, params, view: GalaxyView):
    g = view.g
    edges = profile_edges(params, g)
    zmax = float(params["ExtraZmax"])
    st = view.of(4) & (np.abs(view.z) < zmax)
    p, m = view.pos[st], view.mass[st]
    R = np.sqrt(p[:, 0] ** 2 + p[:, 1] ** 2)
    phi = np.arctan2(p[:, 1], p[:, 0])
    sw, sc, _ = binned(R, edges, m, np.cos(2 * phi))
    _, ss, _ = binned(R, edges, m, np.sin(2 * phi))
    with np.errstate(divide="ignore", invalid="ignore"):
        a2 = np.sqrt(sc ** 2 + ss ** 2) / sw
    ph = 0.5 * np.degrees(np.arctan2(ss, sc))
    rmid = 0.5 * (edges[1:] + edges[:-1])
    if np.all(~np.isfinite(a2)):
        k = 0; a2max = np.nan; lbar = np.nan
    else:
        k = int(np.nanargmax(a2)); a2max = a2[k]
        j = k
        while j + 1 < len(a2) and np.isfinite(a2[j + 1]) and a2[j + 1] > 0.5 * a2max and \
                abs(((ph[j + 1] - ph[k] + 90) % 180) - 90) < 10.0:
            j += 1
        lbar = edges[j + 1]
    # global A2 within the projected half-mass radius
    rh = half_mass_radius_2d(R[view.r[st] <= g.rgal], m[view.r[st] <= g.rgal])
    w = R <= rh
    a2g = abs(np.sum(m[w] * np.exp(2j * phi[w]))) / max(m[w].sum(), 1e-300)
    write_table(names.per_gal("bar", g.gal_id),
                ["R_min (kpc)", "R_max (kpc)", "R_mid (kpc)", "A2 = |sum m exp(2i phi)| / sum m",
                 "phase of the m=2 mode (deg)"],
                np.column_stack([edges[:-1], edges[1:], rmid, a2, ph]),
                header_lines=[f"m=2 Fourier analysis of the stellar disc of galaxy {g.gal_id}, "
                              f"|z| < {zmax:g} kpc",
                              f"A2max = {a2max:.5f} at R = {rmid[k]:.4f} kpc; bar length (A2 > A2max/2, "
                              f"phase within 10 deg) = {lbar:.4f} kpc",
                              f"A2 within the projected half-mass radius ({rh:.4f} kpc) = {a2g:.5f}"])
    return dict(a2max=a2max, bar_length=lbar, a2_half=a2g)


# --------------------------------------------------------------------------
def maps(names, params, view: GalaxyView):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    g = view.g
    half = radius_value(params["MapSize"], g)
    npx = int(params["MapPixels"])
    comps = [c.strip() for c in str(params["MapComponents"]).split(",") if c.strip()]
    edges = np.linspace(-half, half, npx + 1)
    pix_area = (2 * half / npx * 1000.0) ** 2        # pc^2
    for comp in comps:
        if comp == "stars":
            msk, w = view.of(4), view.mass
        elif comp == "dm":
            msk, w = view.of(1), view.mass
        elif comp == "gas":
            msk, w = view.of(0), view.mass
        elif comp == "cold":
            msk, w = view.of(0), view.gas_array("mcld")
        elif comp == "hot":
            msk, w = view.of(0), view.gas_array("mhot")
        else:
            continue
        p = view.pos[msk]
        for proj, (i, j) in (("faceon", (0, 1)), ("edgeon", (0, 2))):
            h, _, _ = np.histogram2d(p[:, i], p[:, j], bins=[edges, edges], weights=w[msk])
            sig = (h / pix_area).T                     # Msun/pc^2, [row=y, col=x]
            base = f"{names.outdir}map_{comp}_{proj}.{names.flag}{names.snap}.gal{g.gal_id}"
            np.save(base + ".npy", sig)
            with open(base + ".extent", "w") as fx:
                fx.write("# xmin xmax ymin ymax (kpc) of %s.npy; rows = %s, columns = x\n%g %g %g %g\n"
                         % (os.path.basename(base), "y" if proj == "faceon" else "z",
                            -half, half, -half, half))
            fig, ax = plt.subplots(figsize=(5, 5))
            pos = sig[sig > 0]
            vmin = np.percentile(pos, 1) if pos.size else 1e-3
            im = ax.imshow(np.log10(np.maximum(sig, vmin)), origin="lower",
                           extent=[-half, half, -half, half], cmap="inferno")
            ax.set_xlabel("x (kpc)")
            ax.set_ylabel("y (kpc)" if proj == "faceon" else "z (kpc)")
            ax.set_title(f"{comp}, {proj}, galaxy {g.gal_id}")
            fig.colorbar(im, ax=ax, label=r"log $\Sigma$ (M$_\odot$ pc$^{-2}$)", shrink=0.8)
            fig.tight_layout()
            fig.savefig(base + ".png", dpi=110)
            plt.close(fig)
