# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Output files of ExtractGalaxies, byte-for-byte in the same format:
profiles.*, info.*, jhist.*, sfr.*, rotation_matrix.*, list_of_selected_gals.*,
list_of_blocks.dat, table-*.gal*, (table-)to_find_gals.*.
"""

from __future__ import annotations

import math
import os
from typing import List

import numpy as np

from .galaxy import NBIN_J, NBIN_SFR, PI_EG
from .profiles import KNOWN_BLOCKS

__all__ = ["OutputNames", "write_profiles", "write_info", "write_jhist", "write_sfr",
           "write_rotation_matrix", "write_selected_gals", "write_list_of_blocks",
           "write_tables", "write_to_find_gals"]


def _pad4(s: str) -> str:
    return (s + "    ")[:4] if len(s) < 4 else s


class OutputNames:
    def __init__(self, outdir: str, fileflag: str, snapshot: str):
        self.outdir = (outdir.rstrip("/") + "/") if outdir else ""
        self.flag = (fileflag + ".") if fileflag and not fileflag.endswith(".") else fileflag
        self.snap = os.path.basename(snapshot)
        self.snapnumber = self.snap[-3:]
        if self.outdir:
            os.makedirs(self.outdir, exist_ok=True)

    def per_gal(self, prefix: str, gal_id: int) -> str:
        return f"{self.outdir}{prefix}.{self.flag}{self.snap}.gal{gal_id}"

    def table(self, prefix: str, gal_id: int) -> str:
        return f"{self.outdir}{prefix}.{self.flag}gal{gal_id}"

    def per_snap(self, prefix: str) -> str:
        return f"{self.outdir}{prefix}.{self.flag}{self.snap}"


# --------------------------------------------------------------------------
def write_profiles(names: OutputNames, params, gals, profiles, counters):
    p = params
    nbin = int(p["NumberOfBins"])
    grid = bool(p["GridDensity"])
    for igal, g in enumerate(gals):
        with open(names.per_gal("profiles", g.gal_id), "w") as f:
            f.write("# Number of bins: ")
            if grid:
                f.write(" %d x %d\n" % (nbin, nbin))
                f.write("#  1: Bin size,  2: X,  3: Y ")
            else:
                f.write(" %d\n" % nbin)
                f.write("#  1: Rmed,  2: Rmin,  3: Rmax ")
            for iq, q in enumerate(p.quantities):
                f.write(" %2d: %4s [Type %d]," % (4 + iq * 2, _pad4(q.quantity), q.ptype))
            f.write("\n")
            f.write("# Gal N. %d, ID: %d\n" % (igal, g.gal_id))
            nb = nbin * nbin if grid else nbin
            if not grid:
                rmin, rmax = g.rmin, g.rmax
                if g.log_bins:
                    if rmin < 1.0e-10:
                        rmin = rmax * 0.01
                    rmin, rmax = math.log10(rmin), math.log10(rmax)
                dx = (rmax - rmin) / nbin
            for ib in range(nb):
                if grid:
                    xmed = 2 * g.rmax / float(nbin)
                    x1 = ((ib % nbin) + 0.5) * xmed - g.rmax
                    x2 = ((ib // nbin) + 0.5) * xmed - g.rmax
                else:
                    x1 = rmin + dx * ib
                    x2 = x1 + dx
                    xmed = 0.5 * (x1 + x2)
                    if g.log_bins:
                        x1, x2, xmed = 10.0 ** x1, 10.0 ** x2, 10.0 ** xmed
                line = " %5f %5f %5f   " % (xmed, x1, x2)
                for iq in range(len(p.quantities)):
                    v = float(np.float32(profiles[igal][iq][ib]))
                    line += " %12e %6d  " % (v, counters[igal][iq][ib])
                f.write(line + "\n")


def write_info(names: OutputNames, params, analyser, igal: int):
    p = params
    d = analyser.d
    g = analyser.gals[igal]
    sets = analyser.analysis[igal]
    h = d.header
    L = []
    L.append("Galaxy No. %d \n" % igal)
    L.append("\n")
    L.append(("All quantities are in h=%f\n" % h.hubble_param) if d.units.htrue
             else "All quantities are in h=1\n")
    L.append("Galaxy radius:  %8f  - Virial radius: %8f\n" % (g.rgal, g.rvir))
    if p["SurfaceDensity"] or p["VolumeDensity"]:
        L.append("Profiles are given from rmin=%8f to rmax=%8f " % (g.rmin, g.rmax))
        L.append(("in %d Log10 intervals\n" if g.log_bins else "in %d linear intervals\n")
                 % int(p["NumberOfBins"]))
    else:
        L.append("Profiles are given in %d x %d cells in a square of side 2 x %8f centred "
                 "on the galaxy" % (p["NumberOfBins"], p["NumberOfBins"], g.rmax))
    L.append("\n")
    if int(p["CmType"]) > 0:
        L.append("The following quantities are computed for particles of type %d within %f kpc\n"
                 % (p["CmType"], g.ralign))
    else:
        L.append("The following quantities are computed for cold gas and star particles "
                 "within %f kpc\n" % g.ralign)
    L.append("\n")
    L.append("  Center of mass position: [ %10f  %10f  %10f ] (kpc) \n" % tuple(g.com))
    L.append("  Center of mass velocity: [ %10f  %10f  %10f ] (km/s) \n" % tuple(g.vel))
    L.append("\n")
    L.append("  Moment of inertia in box coodinates (10^10 Msun kpc^2):\n")
    for j in range(3):
        L.append("    [ %10f  %10f  %10f ]  \n" % tuple(g.inertia[j]))
    L.append("\n")
    L.append("  Eigenvalues (10^10 Msun kpc^2)\n")
    L.append("    [ %10f  %10f  %10f ] \n" % tuple(g.eigen))
    L.append("\n")
    L.append("  Specific angular momentum in the box frame:\n")
    L.append("    [ %10g  %10g  %10g ] kpc km/s" % tuple(g.am))
    L.append("  -  modulus: %10g  kpc km/s\n" % math.sqrt(float(np.dot(g.am, g.am))))
    L.append("\n")
    L.append("  Rotation matrix: \n")
    for j in range(3):
        L.append("    [ %6f  %6f  %6f ]  \n" % tuple(g.rotation[j]))
    L.append("\n")
    vm = g.rotation @ g.am
    L.append("  Specific angular momentum in the box frame:\n")      # (sic) galaxy frame
    L.append("    [ %10f  %10f  %10f ] kpc km/s\n" % tuple(vm))
    L.append("  Angle between angular momentum and z-axis of inertia tensor: %8f degree\n"
             % g.angle)
    L.append("\n\n")
    labels = ["  gas particles:             ", "  DM particles:              ",
              "  cold gas particles:        ", "  hot gas particles:         ",
              "  star particles:            ", "  cold gas + star particles: "]
    for turn in range(2):
        if turn and g.rgal == g.rvir:
            break
        if turn:
            L.append("The following quantities refer to r<rvir, with rvir=%f kpc\n" % g.rvir)
        else:
            L.append("The following quantities refer to r<rgal, with rgal=%f kpc\n" % g.rgal)
        L.append("\n")
        for k in range(6):
            s = sets[turn * 6 + k]
            L.append(labels[k])
            L.append(" %10u particles\n" % s.number)
            L.append("     Mass: %10g Msun\n" % s.mass)
            L.append("     eigenvalues of the inertia tensor:               [  %10g  %10g  %10g  ]"
                     " 10^10 Msun kpc^2\n" % tuple(s.eigenvalues))
            L.append("     specific angular momenum in box coordinates:     [  %10g  %10g  %10g  ]"
                     " kpc km/s\n" % tuple(s.specific_spin))
            vm = g.rotation @ s.specific_spin
            L.append("     specific angular momenum in galaxy coordinates:  [  %10g  %10g  %10g  ]"
                     " kpc km/s" % tuple(vm))
            mod = math.sqrt(float(np.dot(vm, vm)))
            L.append("  -  modulus: %10g  kpc km/s\n" % mod)
            with np.errstate(invalid="ignore", divide="ignore"):
                ang = math.acos(min(1.0, abs(vm[2] / mod))) * 180.0 / PI_EG if mod > 0 else float("nan")
            L.append("     angle between angular momentum and galaxy z-axis: %8f degree\n" % ang)
    a = sets
    L.append("\n")
    bar = a[6].mass + a[10].mass
    tot = bar + a[7].mass
    L.append("Total baryonic mass within virial radius: %10g  (baryon fraction = %f)\n"
             % (bar, bar / tot if tot else float("nan")))
    L.append("Fraction of total mass in the (central) galaxy: %10f\n"
             % (a[5].mass / tot if tot else float("nan")))
    L.append("\n")
    for name, kk in (("dark matter", 1), ("stars", 4), ("hot gas", 3)):
        dm = a[6 + kk].mass - a[kk].mass
        L.append("Specific angular momentum of halo %s in box coordinates (%10g Msun):\n"
                 % (name, dm))
        with np.errstate(invalid="ignore", divide="ignore"):
            vm = (a[6 + kk].spin - a[kk].spin) / dm
        L.append("  [  %10g  %10g  %10g  ]  kpc km/s - modulus %10g  kpc km/s\n"
                 % (vm[0], vm[1], vm[2], math.sqrt(float(np.dot(vm, vm)))))
        L.append("\n")
    with np.errstate(invalid="ignore", divide="ignore"):
        bt = g.mbulge / (g.mbulge + g.mdisc)
    L.append("Simple B/T decomposition gives stellar bulge mass: %10g, stellar disc mass: "
             "%10g, B/T=%10f, stellar halo mass: %10g\n" % (g.mbulge, g.mdisc, bt, g.mhalo))
    with open(names.per_gal("info", g.gal_id), "w") as f:
        f.write("".join(L))


def write_jhist(names: OutputNames, g, jh):
    with open(names.per_gal("jhist", g.gal_id), "w") as f:
        f.write("#  1) jz/jcirc - 2) stellar mass in bin - 3...) 2D histogram\n")
        f.write("# radial bin centers:")
        for j in range(1, NBIN_J + 1):
            f.write(" %12g" % ((j - 0.5) / NBIN_J * g.rgal))
        f.write("\n#\n")
        for i in range(NBIN_J):
            f.write("  %10f  %12g" % ((i + 0.5) / NBIN_J * 4.0 - 2.0, jh[i, 0]))
            for j in range(1, NBIN_J + 1):
                f.write(" %12g" % jh[i, j])
            f.write("\n")


def write_sfr(names: OutputNames, g, igal: int, hist, age_of_u: float, htrue: bool):
    with open(names.per_gal("sfr", g.gal_id), "w") as f:
        f.write("# Star formation rate of stars within rgal of galaxy %d\n" % igal)
        f.write("#  1) time (Gyr) -  2) SFR (Msun/yr)\n" if htrue
                else "#  1) time (Gyr/h) -  2) SFR (Msun/yr)\n")
        for i in range(NBIN_SFR):
            f.write("  %12f  %12f\n" % ((i + 0.5) / NBIN_SFR * age_of_u, hist[i]))


def write_rotation_matrix(names: OutputNames, g):
    with open(names.per_gal("rotation_matrix", g.gal_id), "w") as f:
        f.write("#  Rotation matrix: \n")
        for j in range(3):
            f.write(" %6f  %6f  %6f\n" % tuple(g.rotation[j]))


def _hunit(htrue, s):
    return "" if htrue else s


def write_selected_gals(names: OutputNames, gals, header, units):
    ht = units.htrue
    with open(names.per_snap("list_of_selected_gals"), "w") as f:
        f.write("# Complete list of galaxies\n")
        f.write("#  1   : galaxy ID\n")
        f.write("#  2- 4: snapshot number, time, redshift\n")
        f.write("#  5- 7: centre of mass (kpc%s)\n" % _hunit(ht, "/h"))
        f.write("#  8-10: velocity (km/s)\n")
        f.write("# 11-15: Stellar, gas, cold gas masses within rgal, FoF mass (Msun%s), "
                "SFR (Msun/yr)\n" % _hunit(ht, "/h"))
        f.write("# 16-19: Stellar bulge, disc and halo mass within rgal (Msun%s), stellar B/T\n"
                % _hunit(ht, "/h"))
        f.write("# 20-24: rmin, rmax, rgal, rvir, ralign (kpc%s)\n" % _hunit(ht, "/h"))
        f.write("# 25-27: eigenvalues of inertia tensor (10^10 Msun kpc^2 %s)\n" % _hunit(ht, "/h^3"))
        f.write("# 28-30: specific angular momentum (kpc%s km/s)\n" % _hunit(ht, "/h"))
        f.write("# 31   : angle between angular momentum and inertia tensor (deg)\n")
        f.write("# 32-34: centre of mass (code units)\n")
        f.write("# 35-37: velocity (code units)\n")
        f.write("# 38-44: Central, LogFlag, swap, ComputeCM, TrackCM, Align, Pot\n")
        for g in gals:
            f.write((" %12d  %s %f %f   %8f %8f %8f   %8f %8f %8f   %10e %10e %10e %10e %10e   "
                     "%10e %10e %10e %10e   %8f %8f %8f %8f %8f   %10g %10g %10g   "
                     "%10g %10g %10g   %6f   %8f %8f %8f   %8f %8f %8f   %1d %1d %1d %1d %1d %1d %1d\n")
                    % (g.gal_id, names.snapnumber, header.time, header.redshift,
                       *g.com, *g.vel,
                       g.stellar_mass, g.tot_gas_mass, g.cold_gas_mass, g.fof_mass, g.sfr,
                       g.mbulge, g.mdisc, g.mhalo, g.bovert,
                       g.rmin, g.rmax, g.rgal, g.rvir, g.ralign,
                       *g.eigen, *g.am, g.angle,
                       *(g.com / units.length_kpc), *(g.vel / units.vel_kms),
                       g.central, g.log_bins, g.swap, g.compute_cm, g.track_cm, g.align, g.pot))


def write_list_of_blocks(names: OutputNames, derived: List[str]):
    with open(f"{names.outdir}list_of_blocks.dat", "w") as f:
        f.write("known blocks:\n")
        for i, b in enumerate(KNOWN_BLOCKS):
            f.write("  %d %s\n" % (i, _pad4(b)))
        f.write("Derived blocks:\n")
        for i, b in enumerate(derived):
            f.write("  %d %s\n" % (i + len(KNOWN_BLOCKS), _pad4(b)))


def _append(path, header_lines):
    exists = os.path.exists(path)
    f = open(path, "a" if exists else "w")
    if not exists:
        f.write("".join(header_lines))
    return f


def write_tables(names: OutputNames, analyser, header, units, snapnum: int):
    """WriteTables: one line appended per snapshot (table-mass/specspin/spin/inertia)."""
    ht = units.htrue
    tnow = float(analyser.cosmological_time(header.time))
    for igal, g in enumerate(analyser.gals):
        a = analyser.analysis[igal]
        lead = " %3d %8f %8f   " % (snapnum, tnow, header.redshift)
        f = _append(names.table("table-mass", g.gal_id), [
            "# 1) snapshot number - 2) time (Gyr%s) - 3) redshift\n" % _hunit(ht, "/h"),
            "# within rgal:  4) gas mass -  5) DM mass -  6) cold gas mass -  7) hot gas mass -  "
            "8) star mass -  9) cold + star mass (Msun%s)\n" % _hunit(ht, "/h"),
            "# within rvir: 10) gas mass - 11) DM mass - 12) cold gas mass - 13) hot gas mass - "
            "14) star mass - 15) cold + star mass (Msun%s)\n" % _hunit(ht, "/h"),
            "# 16-18: bulge and disc stellar masses within rgal (Msun%s), B/T ratio\n" % _hunit(ht, "/h"),
            "# within rgal: half-mass radius of 19) gas - 20) DM - 21) cold gas - 22) hot gas - "
            "23) stars - 24) cold + stars (kpc%s)\n" % _hunit(ht, "/h"),
            "# within rvir: half-mass radius of 25) gas - 26) DM - 27) cold gas - 28) hot gas - "
            "29) stars - 30) cold + stars (kpc%s)\n" % _hunit(ht, "/h")])
        with np.errstate(invalid="ignore", divide="ignore"):
            bt = g.mbulge / (g.mbulge + g.mdisc)
        f.write(lead + "".join(" %12g" % s.mass for s in a)
                + "  %12g %12g %12g  " % (g.mbulge, g.mdisc, bt)
                + "".join(" %12g" % s.half_mass_radius for s in a) + "\n")
        f.close()
        f = _append(names.table("table-specspin", g.gal_id), [
            "# 1) snapshot number - 2) time (Gyr%s) - 3) redshift\n" % _hunit(ht, "/h"),
            "# specific angular momentum (components and modulus, kpc%s km/s)\n" % _hunit(ht, "/h"),
            "# of (within rgal)  4- 7) gas -  8-11) DM - 12-15) cold gas - 16-19) hot gas - "
            "20-23) stars - 24-27) cold + stars \n",
            "# of (within rvir) 28-31) gas - 32-35) DM - 36-39) cold gas - 40-43) hot gas - "
            "44-47) stars - 48-51) cold + stars \n"])
        f.write(lead + "".join("  %12g %12g %12g  %12g  " % (*s.specific_spin, s.mod_specific_spin)
                               for s in a) + " \n")
        f.close()
        f = _append(names.table("table-spin", g.gal_id), [
            "# 1) snapshot number - 2) time (Gyr%s) - 3) redshift\n" % _hunit(ht, "/h"),
            "# angular momentum (components and modulus, Msun kpc km/s %s))\n" % _hunit(ht, "/h^2"),
            "# of (within rgal)  4- 7) gas -  8-11) DM - 12-15) cold gas - 16-19) hot gas - "
            "20-23) stars - 24-27) cold + stars \n",
            "# of (within rvir) 28-31) gas - 32-35) DM - 36-39) cold gas - 40-43) hot gas - "
            "44-47) stars - 48-51) cold + stars \n"])
        f.write(lead + "".join("  %12g %12g %12g  %12g  " % (*s.spin, s.mod_spin) for s in a) + " \n")
        f.close()
        f = _append(names.table("table-inertia", g.gal_id), [
            "# 1) snapshot number - 2) time (Gyr%s) - 3) redshift\n" % _hunit(ht, "/h"),
            "# within rgal: 6 x (inertia tensor (6x), eigenvalues (3x), eigenvectors (9x))\n",
            "# within rvir: 6 x (inertia tensor (6x), eigenvalues (3x), eigenvectors (9x))\n"])
        line = lead
        for s in a:
            I, ev, V = s.inertia, s.eigenvalues, s.eigenvectors
            line += ("   %12g %12g %12g %12g %12g %12g   %12g %12g %12g   %12g %12g %12g  "
                     "%12g %12g %12g  %12g %12g %12g     "
                     % (I[0, 0], I[1, 1], I[2, 2], I[0, 1], I[0, 2], I[1, 2], *ev,
                        *V[0], *V[1], *V[2]))
        f.write(line + "\n")
        f.close()

    for g in analyser.gals:
        f = _append(names.table("table-info", g.gal_id), [
            "#  1   : galaxy ID\n", "#  2- 4: snapshot number, time, redshift\n",
            "#  5- 7: centre of mass (kpc%s)\n" % _hunit(ht, "/h"), "#  8-10: velocity (km/s)\n",
            "# 11-14: Stellar, gas, cold gas and FoF masses within rgal (Msun%s)\n" % _hunit(ht, "/h"),
            "# 15-19: rmin, rmax, rgal, rvir, ralign (kpc%s)\n" % _hunit(ht, "/h"),
            "# 20-22: eigenvalues of inertia tensor (10^10 Msun kpc^2 %s)\n" % _hunit(ht, "/h^3"),
            "# 23-25: specific angular momentum (kpc%s km/s)\n" % _hunit(ht, "/h"),
            "# 26   : angle between angular momentum and inertia tensor (deg)\n",
            "# 27-33: Central, LogFlag, swap, ComputeCM, TrackCM, Align, Pot\n"])
        f.write((" %12d  %s %f %f  %8f %8f %8f   %8f %8f %8f   %10e %10e %10e %10e  "
                 "%8f %8f %8f %8f %8f  %10g %10g %10g  %10g %10g %10g  %6f  "
                 "%1d %1d %1d %1d %1d %1d %1d\n")
                % (g.gal_id, names.snapnumber, header.time, header.redshift, *g.com, *g.vel,
                   g.stellar_mass, g.tot_gas_mass, g.cold_gas_mass, g.fof_mass,
                   g.rmin, g.rmax, g.rgal, g.rvir, g.ralign, *g.eigen, *g.am, g.angle,
                   g.central, g.log_bins, g.swap, g.compute_cm, g.track_cm, g.align, g.pot))
        f.close()


def write_to_find_gals(names: OutputNames, gals, header, units, tables: bool):
    head = ["# Complete list of galaxies\n", "#  1   : galaxy ID\n",
            "#  2- 4: snapshot number, time, redshift\n", "#  5- 7: centre of mass (code units)\n",
            "#  8-10: velocity (code units)\n", "# 11-19: rotation matrix\n"]
    fmt = (" %12d  %s %f %f   %12g %12g %12g   %12g %12g %12g    %8f %8f %8f  %8f %8f %8f  "
           "%8f %8f %8f\n")

    def line(g):
        return fmt % (g.gal_id, names.snapnumber, header.time, header.redshift,
                      *(g.com / units.length_kpc), *(g.vel / units.vel_kms),
                      *g.rotation[0], *g.rotation[1], *g.rotation[2])
    if tables:
        for g in gals:
            f = _append(names.table("table-to_find_gals", g.gal_id), head[1:])
            f.write(line(g))
            f.close()
    else:
        with open(names.per_snap("to_find_gals"), "w") as f:
            f.write("".join(head))
            for g in gals:
                f.write(line(g))
