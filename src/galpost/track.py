# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Main-progenitor tracking from particle IDs (no halo catalogues needed).

    galpost-track <parameter file> --snapshots snap_134 snap_130 ... [--tracers 1000]

1. At the latest snapshot the galaxies are selected with the rules of the
   parameter file (FoF, SubFind or list) and, for each of them, the N most
   bound dark-matter particles within rgal are taken as tracers (lowest POT;
   closest to the centre if the snapshot has no POT block). DM IDs do not
   change with time.
2. Going backwards in time, in each snapshot the tracers are found by ID and
   the centre of the main progenitor is the result of a shrinking-sphere
   iteration on them (Power et al. 2003), which follows the densest group of
   tracers, i.e. the main branch through mergers. M200c and R200c are measured
   around that centre, and the N most bound DM particles of the progenitor
   become the tracers for the next step, so the chain survives the dispersal
   of the original tracers.
3. The result is written to track.<Fileflag>.gal<ID> in Outdir: one line per
   snapshot with centre, M200c, R200c, the analysis radii and the fraction of
   tracers found within rgal of the new centre (a low value flags a major
   merger or a lost galaxy).

The track file can then be used by galpost with the tag TrackFile: the galaxy
of each snapshot is taken from it, so the main progenitor can be analysed on
any set of snapshots with a single parameter file.
"""

from __future__ import annotations

import math
import os
import sys
import time
from typing import List, Optional

import numpy as np

from .gadget_io import Snapshot
from .params import read_params
from .particles import region_indices, _wrap_diff
from .units import Units

__all__ = ["track_galaxies", "read_track", "track_galaxy_for_snapshot", "main"]

G = 4.30091e-6           # kpc (km/s)^2 / Msun

TRACK_COLUMNS = [
    "snapshot number", "scale factor", "redshift",
    "centre x (kpc)", "centre y (kpc)", "centre z (kpc)",
    "centre x (code units)", "centre y (code units)", "centre z (code units)",
    "M200c (Msun)", "R200c (kpc)", "rmin (kpc)", "rgal (kpc)", "ralign (kpc)", "rmax (kpc)",
    "log bins (1/0)", "fraction of tracers within rgal of the centre", "number of tracers",
]


# --------------------------------------------------------------------------
def rho_crit(header, units) -> float:
    """Critical density in Msun/kpc^3 (EG units: h=1 units unless Htrue)."""
    hub = units.hubble if units.htrue else 1.0
    z = header.redshift
    e2 = header.omega0 * (1 + z) ** 3 + header.omega_lambda
    H = 0.1 * hub * math.sqrt(e2)               # km/s/kpc
    return 3.0 * H * H / (8.0 * math.pi * G)


def shrinking_sphere(pos: np.ndarray, start: np.ndarray, box: float, periodic: bool,
                     shrink: float = 0.9, nmin: int = 50) -> np.ndarray:
    """Centre of the densest concentration of the points (kpc)."""
    rel = pos - start
    if periodic:
        rel = _wrap_diff(rel, box)
    c = rel.mean(axis=0)
    r = np.sqrt(((rel - c) ** 2).sum(axis=1)).max()
    sel = np.ones(len(rel), bool)
    while sel.sum() > nmin:
        r *= shrink
        d2 = ((rel - c) ** 2).sum(axis=1)
        new = d2 <= r * r
        if new.sum() < nmin:
            break
        sel = new
        c = rel[sel].mean(axis=0)
    centre = start + c
    if periodic:
        centre = np.mod(centre, box)
    return centre


def _region(snap, units, centre, radius):
    """Positions (relative, kpc), masses (Msun), types and indices within radius."""
    idx = region_indices(snap, units, centre[None, :], radius)
    box = snap.header.boxsize * units.length_kpc
    pos, mass, typ, ind = [], [], [], []
    for t in range(6):
        if len(idx[t]) == 0:
            continue
        p = snap.read("POS", t, idx[t]).astype(np.float64) * units.length_kpc - centre
        if units.periodic:
            p = _wrap_diff(p, box)
        pos.append(p)
        mass.append(snap.read("MASS", t, idx[t]).astype(np.float64) * units.mass_msun)
        typ.append(np.full(len(idx[t]), t))
        ind.append(idx[t])
    return np.concatenate(pos), np.concatenate(mass), np.concatenate(typ), np.concatenate(ind)


def so_radius(r: np.ndarray, m: np.ndarray, rho: float, delta: float = 200.0):
    """(R_delta, M_delta) for particles at distances r with masses m."""
    o = np.argsort(r)
    rs = np.maximum(r[o], 1e-6)
    cm = np.cumsum(m[o])
    dens = cm / (4.0 / 3.0 * math.pi * rs ** 3)
    above = np.nonzero(dens >= delta * rho)[0]
    if len(above) == 0:
        return float("nan"), float("nan")
    k = above[-1]                       # outermost radius still above the threshold
    return float(rs[k]), float(cm[k])


def _radii(params, mode, r200):
    """Analysis radii with the rules of the selection mode of the parameter file."""
    p = params
    if mode == "sub":
        rmin, rmax, ralign, logb = p["SubRmin"], p["SubRmax"], p["SubRalign"], p["SubRlog"]
        rgal = 0.1 * r200
    elif mode == "fof":
        rmin, rmax, ralign, logb = p["FofRmin"], p["FofRmax"], p["FofRalign"], p["FofRlog"]
        rgal = p["FofRgal"] if (p["FofRgal"] or 0) > 0 else 0.1 * r200
    else:
        rmin, rmax, ralign, logb = 0.5, 0.0, 0.0, True
        rgal = 0.1 * r200
    rmin = float(rmin) if rmin is not None else 0.5
    rmax = float(rmax) if rmax and float(rmax) > 0 else r200
    ralign = float(ralign) if ralign and float(ralign) > 0 else rgal
    return rmin, rgal, ralign, rmax, bool(logb)


def _tracers(snap, pos, typ, ind, r, rgal, r200, n):
    """IDs of the n most bound DM particles within rgal (within R200 if too few)."""
    dm = typ == 1
    for radius in (rgal, r200, np.inf):
        sel = dm & (r <= radius)
        if sel.sum() >= n or radius == np.inf:
            break
    cand = np.nonzero(sel)[0]
    if "POT" in snap.blocks:
        pot = snap.read("POT", 1, ind[cand]).astype(np.float64)
        order = np.argsort(pot)
    else:
        order = np.argsort(r[cand])
    chosen = cand[order[:n]]
    return snap.read("ID", 1, ind[chosen])


# --------------------------------------------------------------------------
def track_galaxies(paramfile: str, snapshots: List[str], ntracers: int = 1000,
                   outdir: Optional[str] = None, log=print) -> List[str]:
    from .pipeline import _select, disable_transparent_hugepages
    from .egoutput import OutputNames
    disable_transparent_hugepages()
    p = read_params(paramfile, check=False)
    snaps = []
    for s in snapshots:
        sn = Snapshot(s, companions=())
        snaps.append((sn.header.time, s, sn))
    snaps.sort(key=lambda x: -x[0])                  # latest first
    t0, ref_path, ref = snaps[0]
    u0 = Units.from_files(p["GadgetParamfile"], ref.header, bool(p["Htrue"]))
    mode = "sub" if p["ReadSubFindGals"] else ("fof" if p["ReadFofHalos"] else "list")
    gals, _ = _select(p, ref, ref_path, u0, log)
    names = OutputNames(outdir if outdir is not None else p["Outdir"], p["Fileflag"], ref_path)
    written = []
    for g in gals:
        log(f"*** tracking galaxy {g.gal_id} from {os.path.basename(ref_path)}")
        rows = []
        centre = np.asarray(g.com, dtype=np.float64)
        r200_prev = g.rvir
        tracer_ids = None
        for k, (a, path, snap) in enumerate(snaps):
            tt = time.time()
            u = Units.from_files(p["GadgetParamfile"], snap.header, bool(p["Htrue"]))
            box = snap.header.boxsize * u.length_kpc
            frac = 1.0
            if tracer_ids is not None:
                # centre from the tracers found in this snapshot
                ids = snap.read("ID", 1)
                found = np.nonzero(np.isin(ids, tracer_ids))[0]
                del ids
                if len(found) == 0:
                    log(f"   {os.path.basename(path)}: no tracer found, stopping")
                    break
                tpos = snap.read("POS", 1, found).astype(np.float64) * u.length_kpc
                # start from the previous centre (comoving code units -> kpc of this snapshot)
                centre = shrinking_sphere(tpos, centre_code * u.length_kpc, box, u.periodic)
            radius = max(2.5 * r200_prev * (a / snaps[max(k - 1, 0)][0]), 30.0)
            pos, mass, typ, ind = _region(snap, u, centre, radius)
            r = np.sqrt((pos ** 2).sum(axis=1))
            r200, m200 = so_radius(r, mass, rho_crit(snap.header, u))
            if not np.isfinite(r200):
                log(f"   {os.path.basename(path)}: no overdensity above 200 rho_c, stopping")
                break
            rmin, rgal, ralign, rmax, logb = _radii(p, mode, r200)
            if tracer_ids is not None:
                tid_here = snap.read("ID", 1, ind[typ == 1])
                inside = np.isin(tracer_ids, tid_here[r[typ == 1] <= rgal])
                frac = float(inside.mean())
            tracer_ids = _tracers(snap, pos, typ, ind, r, rgal, r200, ntracers)
            centre_code = centre / u.length_kpc
            rows.append([int(os.path.basename(path)[-3:]), snap.header.time, snap.header.redshift,
                         *centre, *(centre / u.length_kpc), m200, r200, rmin, rgal, ralign, rmax,
                         int(logb), frac, len(tracer_ids)])
            log(f"   {os.path.basename(path)}  z={snap.header.redshift:7.3f}  M200c={m200:10.4e}  "
                f"R200c={r200:8.2f} kpc  tracers within rgal={frac:5.2f}  ({time.time() - tt:.1f} s)")
            r200_prev = r200
        fn = names.table("track", g.gal_id)
        _write_track(fn, rows, g.gal_id, bool(p["Htrue"]), ntracers, paramfile)
        written.append(fn)
        log(f"   written {fn}")
    return written


def _write_track(path, rows, gal_id, htrue, ntracers, paramfile):
    rows = sorted(rows, key=lambda r: r[0])
    with open(path, "w") as f:
        f.write(f"# Main-progenitor track of galaxy {gal_id} (galpost-track)\n")
        f.write(f"# parameter file: {os.path.abspath(paramfile)}\n")
        f.write(f"# htrue = {int(htrue)}   tracers: {ntracers} most bound DM particles\n")
        for i, c in enumerate(TRACK_COLUMNS):
            f.write("# %2d: %s\n" % (i + 1, c))
        for r in rows:
            f.write(" %4d %12.8f %12.6f  " % tuple(r[:3])
                    + " ".join("%14.6f" % x for x in r[3:9]) + "  "
                    + " ".join("%14.6e" % x for x in r[9:15])
                    + "  %d  %8.5f  %6d\n" % (r[15], r[16], r[17]))


def read_track(path):
    """Returns (rows as array, htrue flag)."""
    htrue = None
    with open(path) as f:
        for line in f:
            if line.startswith("# htrue"):
                htrue = bool(int(line.split("=")[1].split()[0]))
    data = np.atleast_2d(np.loadtxt(path))
    return data, htrue


def track_galaxy_for_snapshot(path, snapnum, units, gal_id_hint=None):
    """Galaxy (selection.Galaxy) of the track file for snapshot `snapnum`."""
    from .selection import Galaxy
    data, htrue = read_track(path)
    if htrue is not None and htrue != units.htrue:
        raise ValueError(f"{path} was written with htrue={int(htrue)}, the parameter file "
                         f"has Htrue={'on' if units.htrue else 'off'}")
    row = data[data[:, 0].astype(int) == int(snapnum)]
    if not len(row):
        raise ValueError(f"snapshot {snapnum} is not in the track file {path}")
    r = row[0]
    name = os.path.basename(path)
    gid = int(name.rsplit("gal", 1)[1]) if "gal" in name and name.rsplit("gal", 1)[1].isdigit() else 0
    com = r[6:9] * units.length_kpc                 # from code units: exact for this snapshot
    return Galaxy(gal_id=gid, com=com, rmin=r[11], rmax=r[14], rgal=r[12], rvir=r[10],
                  ralign=r[13], log_bins=bool(r[15]), compute_cm=False, track_cm=True,
                  align=True, pot=False, central=1, fof_mass=r[9])


# --------------------------------------------------------------------------
def main(argv=None):
    import argparse
    import glob
    ap = argparse.ArgumentParser(prog="galpost-track",
                                 description="Track the main progenitor of the selected galaxies "
                                             "through a set of snapshots, using particle IDs")
    ap.add_argument("paramfile", help="galpost parameter file (selection, units, Outdir)")
    ap.add_argument("--snapshots", nargs="+", required=True,
                    help="snapshot files (any order; shell patterns allowed); the latest one "
                         "is the reference")
    ap.add_argument("--tracers", type=int, default=1000, help="number of tracer particles")
    ap.add_argument("--outdir", default=None, help="output directory (default: Outdir)")
    a = ap.parse_args(argv)
    snaps = []
    for s in a.snapshots:
        g = sorted(glob.glob(s))
        snaps += [x for x in (g or [s]) if not x.endswith(".muppi")]
    track_galaxies(a.paramfile, snaps, a.tracers, a.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
