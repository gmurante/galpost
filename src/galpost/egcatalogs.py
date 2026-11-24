# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
ASCII particle catalogues of ExtractGalaxies:

* WriteGalaxyCatalog / WriteParticleList: per galaxy, files
  ``{catalog|particle_list}_{mp,sp,st,dm,ds,bl}.<flag><snap>.gal<ID>``
  (mp/sp: multi-/single-phase gas, st: stars, dm: type 1, ds: type 2,
  bl: type 3 -- as in the EG code, which differs from its documentation);
* WriteFullCatalog / WriteStarCatalog: ``fullcatalog_<type>.<flag><snap>``.

Same columns and formats as EG.  Particles are those assigned to the galaxy
for the profiles (sphere / cylinder / cube of radius rmax).
"""

from __future__ import annotations

import numpy as np

from .particles import MU_H, MU_C, T_C, E_OVER_M_TO_T

__all__ = ["write_galaxy_catalogs", "write_full_catalog"]

RHO_TO_N = 40.6321263
NGALSMAX = 5


def _h_line(units, header):
    return ("# ALL QUANTITIES ARE IN h=%f\n" % header.hubble_param if units.htrue
            else "# ALL QUANTITIES ARE IN h=1\n")


def _header_gas(f, g, units, header, metal_names, has_metals):
    f.write(_h_line(units, header))
    f.write("# In output all particles within %8.3f kpc from galaxy centre\n" % g.rmax)
    for s in ["#  1- 3: Position (physical kpc, galaxy system)", "#  4- 6: Velocity (km/s, galaxy system)",
              "#  7   : Mass (Msun)", "#  8   : particle ID", "#  9   : polar radius (kpc)",
              "# 10   : gas density (Msun/pc^3)", "# 11   : circularity (j_z/j_circ)",
              "# 12   : SPH temperature (K)", "# 13   : hot phase temperature (K)",
              "# 14   : internal energy (code units)", "# 15   : hot gas mass (Msun)",
              "# 16   : cold gas mass (Msun)", "# 17   : stellar mass (Msun)",
              "# 18   : molecular gas mass (Msun)", "# 19   : ISM pressure (K cm^-3)",
              "# 20   : SPH pressure (K cm^-3)", "# 21   : hot phase thermal energy (10^50 erg)",
              "# 22   : fraction of hot gas", "# 23   : fraction of cold gas",
              "# 24   : filling factor of hot phase (pc^3)", "# 25   : SPH volume of the gas particle (pc^3)",
              "# 26   : hot phase density (Msun/pc^3)", "# 27   : cold phase density (Msun/pc^3)",
              "# 28   : hot phase particle density (cm^-3)", "# 29   : cold phase particle density (cm^-3)",
              "# 30   : star formation rate (Msun/yr)", "# 31   : outflow rate (Msun/yr)",
              "# 32   : inflow rate (Msun/yr)", "# 33   : electron density (cm^-3)",
              "# 34-36: position (physical kpc, box system)", "# 37-39: velocity (km/s, box system)",
              "# 40-42: density gradient (if present in snapshot)", "# 43   : number of steps in multi-phase"]:
        f.write(s + "\n")
    if has_metals:
        for i, m in enumerate(metal_names):
            f.write("# %2d   : mass in %s\n" % (44 + i, (m + "   ")[:3]))
    f.write("#\n")


def _header_dm(f, g, units, header):
    f.write(_h_line(units, header))
    f.write("# In output all particles within %8.3f kpc from galaxy centre\n" % g.rmax)
    for s in ["#  1- 3: Position (physical kpc, galaxy system)", "#  4- 6: Velocity (km/s, galaxy system)",
              "#  7   : Mass (Msun)", "#  8   : particle ID", "#  9   : polar radius (kpc)",
              "# 10-12: position (physical kpc, box system)", "# 13-15: velocity (km/s, box system)", "#"]:
        f.write(s + "\n")


def _header_star(f, g, units, header, metal_names, has_metals):
    f.write(_h_line(units, header))
    f.write("# In output all particles within %8.3f kpc from galaxy centre\n" % g.rmax)
    for s in ["#  1- 3: Position (physical kpc, galaxy system)", "#  4- 6: Velocity (km/s, galaxy system)",
              "#  7   : Mass (Msun)", "#  8   : particle ID", "#  9   : polar radius (kpc)",
              "# 10   : stellar age (code units)", "# 11   : circularity (j_z/j_circ)",
              "# 12   : component (1=bulge, 2=disc, 3=bar, >3=halo)",
              "# 13-15: position (physical kpc, box system)", "# 16-18: velocity (km/s, box system)"]:
        f.write(s + "\n")
    if has_metals:
        for i, m in enumerate(metal_names):
            f.write("# %2d   : mass in %s\n" % (19 + i, (m + "   ")[:3]))
        f.write("# %2d   : initial stellar mass\n" % (19 + len(metal_names)))
    f.write("#\n")


def write_galaxy_catalogs(names, params, analyser, resolver, metal_names):
    d, an = analyser.d, analyser
    units, header = d.units, d.header
    full = bool(params["WriteGalaxyCatalog"])
    kind = "catalog" if full else "particle_list"
    has_metals = d.metals is not None
    gs = d.type_slice(0)
    ng = int(d.npart[0])
    outf, infl = resolver._flows() if full else (None, None)
    gasq = {}
    if full and ng:
        mcld, mhot = d.mcld.astype(np.float64), d.mhot.astype(np.float64)
        ehot, rho = d.ehot.astype(np.float64), d.rho.astype(np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            mgas = mcld + mhot
            fh, fc = mhot / mgas, mcld / mgas
            thot = np.where(mhot > 0, ehot * MU_H / mhot * E_OVER_M_TO_T, 0.0)
            filh = np.where(mhot > 0, 1.0 / (1.0 + fc / fh * MU_H / MU_C * T_C / thot), 0.0)
            vol = mgas / rho
            rho_h = np.where(mhot > 0, mhot / (vol * filh), 0.0)
            rho_c = np.where(1.0 - filh > 0, mcld / (vol * (1.0 - filh)), 0.0)
            hotn, cldn = rho_h / MU_H * RHO_TO_N, rho_c / MU_C * RHO_TO_N
            pres_sph = d.u.astype(np.float64) * rho * RHO_TO_N * units.temp_K
            pres = np.where(d.mp, T_C * cldn, pres_sph)
            mmol = np.where(d.mp, mcld / (1.0 + float(params["Pressure"]) / pres), 0.0)
        gasq = dict(fh=fh, fc=fc, thot=thot, filh=filh, vol=vol, rho_h=rho_h, rho_c=rho_c,
                    hotn=hotn, cldn=cldn, pres_sph=pres_sph, pres=pres, mmol=mmol)

    for igal, g in enumerate(an.gals):
        fmp = open(f"{names.outdir}{kind}_mp.{names.flag}{names.snap}.gal{g.gal_id}", "w") if d.muppi else None
        fsp = open(f"{names.outdir}{kind}_sp.{names.flag}{names.snap}.gal{g.gal_id}", "w")
        for f in (fmp, fsp):
            if f:
                _header_gas(f, g, units, header, metal_names, has_metals)
        sel = np.nonzero(an.belongs[gs] == igal + 1)[0]
        if len(sel):
            pos, vel = an.frame(g, sel + gs.start)
            rp = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
            for k, i in enumerate(sel):
                f = fmp if (d.mp[i] and d.muppi) else fsp
                ia = i + gs.start
                if not full:
                    f.write(" %u %d   %f %f %f\n" % (d.ids[ia], 0, *d.pos[ia]))
                    continue
                grad = (g.rotation @ d.grad[i].astype(np.float64)) if d.grad is not None else np.zeros(3)
                q = gasq
                line = (" %f %f %f %f %f %f   %e %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e "
                        "%e %e %e %e %e %e %e %e %e  %f %f %f %f %f %f  %e %e %e %d  ") % (
                    *pos[k], *vel[k], d.mass[ia], d.ids[ia], rp[k], d.rho[i], an.circularity[i],
                    d.temp[i], q["thot"][i], d.u[i], d.mhot[i], d.mcld[i], d.msf[i], q["mmol"][i],
                    q["pres"][i], q["pres_sph"][i], d.ehot[i], q["fh"][i], q["fc"][i], q["filh"][i],
                    q["vol"][i], q["rho_h"][i], q["rho_c"][i], q["hotn"][i], q["cldn"][i], d.sfr[i],
                    outf[ia], infl[ia], d.ne[i] if d.ne is not None else 0.0,
                    *d.pos[ia], *d.vel[ia], *grad, d.mpsteps[i])
                if has_metals:
                    line += "".join(" %e" % x for x in d.metals[i])
                f.write(line + "\n")
        for f in (fmp, fsp):
            if f:
                f.close()

        # stars
        ss = d.type_slice(4)
        with open(f"{names.outdir}{kind}_st.{names.flag}{names.snap}.gal{g.gal_id}", "w") as f:
            _header_star(f, g, units, header, metal_names, has_metals)
            sel = np.nonzero(an.belongs[ss] == igal + 1)[0]
            if len(sel):
                pos, vel = an.frame(g, sel + ss.start)
                rp = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
            for k, i in enumerate(sel):
                ia = i + ss.start
                if not full:
                    f.write(" %u %d   %f %f %f\n" % (d.ids[ia], 4, *d.pos[ia]))
                    continue
                line = " %f %f %f %f %f %f   %e %u %f %f %f %d   %f %f %f %f %f %f  " % (
                    *pos[k], *vel[k], d.mass[ia], d.ids[ia], rp[k],
                    d.age[i] if d.age is not None else 0.0, an.circularity[ng + i],
                    an.component[ia], *d.pos[ia], *d.vel[ia])
                if has_metals:
                    line += "".join(" %e" % x for x in d.metals[ng + i])
                    line += " %e" % (d.im[i] if d.im is not None else 0.0)
                f.write(line + "\n")

        # dark matter and low-resolution types
        for t, tag in ((1, "dm"), (2, "ds"), (3, "bl")):
            st = d.type_slice(t)
            with open(f"{names.outdir}{kind}_{tag}.{names.flag}{names.snap}.gal{g.gal_id}", "w") as f:
                _header_dm(f, g, units, header)
                sel = np.nonzero(an.belongs[st] == igal + 1)[0]
                if len(sel):
                    pos, vel = an.frame(g, sel + st.start)
                    rp = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
                for k, i in enumerate(sel):
                    ia = i + st.start
                    if full:
                        f.write(" %f %f %f %f %f %f   %e %u %f   %f %f %f %f %f %f\n" % (
                            *pos[k], *vel[k], d.mass[ia], d.ids[ia], rp[k], *d.pos[ia], *d.vel[ia]))
                    else:
                        f.write(" %u %d   %f %f %f\n" % (d.ids[ia], t, *d.pos[ia]))


def write_full_catalog(names, params, analyser):
    d, an = analyser.d, analyser
    if params["WriteFullCatalog"]:
        types = range(6)
    elif params["WriteStarCatalog"]:
        types = [4]
    else:
        return
    ngo = min(len(an.gals), NGALSMAX)
    for t in types:
        st = d.type_slice(t)
        with open(f"{names.outdir}fullcatalog_{t}.{names.flag}{names.snap}", "w") as f:
            f.write(_h_line(d.units, d.header))
            f.write("#  1- 3: Position (physical kpc, box system)\n")
            f.write("#  4- 6: Velocity (km/s, box system)\n")
            f.write("#  7   : Mass (Msun)\n")
            c = 8
            for i in range(ngo):
                f.write("# %2d-%2d: position in the system of galaxy %d (physical kpc)\n" % (c, c + 2, i))
                f.write("# %2d-%2d: velocity in the system of galaxy %d (km/s)\n" % (c + 3, c + 5, i))
                f.write("#\n")
                c += 6
            frames = [an.frame(an.gals[i], np.arange(st.start, st.stop)) for i in range(ngo)]
            for k in range(st.stop - st.start):
                ia = st.start + k
                line = "%f %f %f  %f %f %f  %e " % (*d.pos[ia], *d.vel[ia], d.mass[ia])
                for (pp, vv) in frames:
                    line += "  %f %f %f  %f %f %f  " % (*pp[k], *vv[k])
                f.write(line + "\n")
