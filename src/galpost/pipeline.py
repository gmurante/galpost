# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
The analysis pipeline: the equivalent of the ExtractGalaxies main program.

    galaxies  <- selection (FoF / SubFind / list)
    particles <- region around the galaxies, EG units and precision
    frames    <- centre tracking and alignment
    analysis  <- masses, inertia, spins, B/T, jz/jcirc, SFR history
    profiles  <- derived blocks + binning
    outputs   <- the same files as ExtractGalaxies
"""

from __future__ import annotations

import os
import sys
import time
import warnings
from typing import List, Optional

import numpy as np

from . import egoutput as out
from .egcatalogs import write_full_catalog, write_galaxy_catalogs
from .zoom import write_zoom_snapshots, zoom_radius
from .extra import extra_radius, run_extras
from .catalogs import FofCatalog, SubfindCatalog, find_group_dir
from .galaxy import GalaxyAnalyser
from .gadget_io import Snapshot
from .params import GalaxyLine, Params, read_params, _parse_galaxy
from .particles import ParticleData, region_indices
from .membership import particle_list_flags, subfind_membership
from .units import read_gadget_params
from .profiles import QuantityResolver, compute_profiles, derived_names, read_metal_names
from .selection import (Galaxy, lowres_fraction, select_from_fof, select_from_list,
                        select_from_subfind)
from .units import Units

__all__ = ["run", "load_galaxy_file"]


def load_galaxy_file(path: str) -> List[GalaxyLine]:
    """Galfile: same lines as ListOfGals (without ListOfGals/End)."""
    out_ = []
    with open(path) as f:
        for i, line in enumerate(f, 1):
            s = line.strip()
            if not s or s[0] in "#%" or s.split()[0] in ("ListOfGals", "End"):
                continue
            g = _parse_galaxy(s, f"{path}:{i}")
            if g:
                out_.append(g)
    return out_


def _snapshot_path(p: Params, snap_override: Optional[str]) -> str:
    snap = snap_override or p["Snapshot"]
    return os.path.join(p["Snapdir"], snap) if p["Snapdir"] else snap


def _select(p: Params, snap: Snapshot, snap_path: str, units: Units, log):
    h = snap.header
    cat = None
    num = int(os.path.basename(snap_path)[-3:])
    gdir_base = p["GroupDir"] or (os.path.dirname(snap_path) or ".")
    if p["TrackFile"]:
        from .track import track_galaxy_for_snapshot
        g = track_galaxy_for_snapshot(p["TrackFile"], num, units)
        log(f"   galaxy from track file {p['TrackFile']} (snapshot {num})")
        gals = [g]
    elif p["ReadFofHalos"]:
        fof = FofCatalog(find_group_dir(gdir_base, num), num)
        log(f"   FoF catalogue {fof.group_dir}: {len(fof)} groups, format {fof.fmt}")
        gals = select_from_fof(fof, units, h, p["FofMinMass"], p["FofRmin"], p["FofRmax"],
                               p["FofRgal"], p["FofRalign"], p["FofRlog"], p["FofDoNotAlign"])
    elif p["ReadSubFindGals"]:
        cat = SubfindCatalog(find_group_dir(gdir_base, num), num)
        log(f"   SubFind catalogue {cat.group_dir}: {cat.ngroups} groups, "
            f"{cat.nsubhalos} subhaloes")
        lowres = [int(x) for x in str(p["LowResTypes"]).split(",") if x.strip()]
        gals = select_from_subfind(cat, units, h, p["SubMinStellarMass"], p["SubRmin"],
                                   p["SubRmax"], p["SubRalign"], p["SubRlog"],
                                   p["SubDoNotAlign"], lowres)
    else:
        lines = p.galaxies
        if not lines:
            gf = p["Galfile"] or "list_of_gals.dat"
            lines = load_galaxy_file(gf)
        gals = select_from_list(lines, units, h)
    if not gals:
        raise RuntimeError("no galaxies selected")

    if float(p["MaxLowResFraction"]) >= 0:
        lowres = [int(x) for x in str(p["LowResTypes"]).split(",") if x.strip()]
        lowres_fraction(snap, gals, units, lowres)
        keep = [g for g in gals if not (g.lowres_fraction > float(p["MaxLowResFraction"]))]
        for g in gals:
            if g not in keep:
                log(f"   galaxy {g.gal_id} discarded: low-resolution fraction "
                    f"{g.lowres_fraction:.3f}")
        gals = keep
        if not gals:
            raise RuntimeError("all selected galaxies are contaminated")
    return gals, cat


def _region_radius(p: Params, gals: List[Galaxy]) -> float:
    """Radius of the sphere to load: must contain the profile volume
    (sphere, cylinder of radius and half-height up to rmax, or cube of side
    2 rmax), the analysis spheres and the alignment sphere."""
    geo = np.sqrt(3.0) if p["GridDensity"] else (np.sqrt(2.0) if p["SurfaceDensity"] else 1.0)
    r = 0.0
    for g in gals:
        rr = max(geo * g.rmax, g.rvir, g.rgal, g.ralign)
        if p["AquilaAlignment"]:
            rr = max(rr, 3.0 * g.rgal)
        if p["WriteGalaxySnapshot"]:
            rr = max(rr, geo * zoom_radius(p, g))
        rr = max(rr, extra_radius(p, g))
        r = max(r, rr)
    return 1.2 * r + 30.0


def disable_transparent_hugepages() -> bool:
    """Disable transparent huge pages for this process (Linux only).

    On nodes with THP='always' and memory full of page cache, every large
    allocation triggers memory compaction in the kernel: allocating a few GB
    can take minutes of system time. Disabling THP for the process (no
    privileges needed) makes allocations ~6x faster. Set GALPOST_KEEP_THP=1
    to skip.
    """
    if os.environ.get("GALPOST_KEEP_THP") or not sys.platform.startswith("linux"):
        return False
    try:
        import ctypes
        libc = ctypes.CDLL(None, use_errno=True)
        return libc.prctl(41, 1, 0, 0, 0) == 0          # PR_SET_THP_DISABLE
    except Exception:
        return False


def run(paramfile: str, snap: Optional[str] = None, log=print) -> dict:
    disable_transparent_hugepages()
    t0 = time.time()
    p = read_params(paramfile)
    verbose = bool(p["Verbose"])
    vlog = log if verbose else (lambda *a, **k: None)

    snap_path = _snapshot_path(p, snap)
    table = p["BlockTable"] or None
    snapshot = Snapshot(snap_path, block_table=table)
    log(f"*** galpost: {snapshot}")
    units = Units.from_files(p["GadgetParamfile"], snapshot.header, bool(p["Htrue"]))
    vlog(units.summary())

    gals, catalog = _select(p, snapshot, snap_path, units, vlog)
    for i, g in enumerate(gals):
        vlog("      %d) ID: %5d, CoM=[%f %f %f]  rmin=%f  ralign=%f  rmax=%f LOG:%d  COMP_CM:%d  "
             "TRACK_CM:%d  ALIGN:%d  POT:%d  FoFMass=%g, Central=%d"
             % (i, g.gal_id, *g.com, g.rmin, g.ralign, g.rmax, g.log_bins, g.compute_cm,
                g.track_cm, g.align, g.pot, g.fof_mass, g.central))

    # particles of the region (whole box for COMP_CM)
    full = (any(g.compute_cm for g in gals) or bool(p["WriteFullCatalog"])
            or bool(p["WriteStarCatalog"]))
    radius = _region_radius(p, gals)
    init_com = np.array([g.com for g in gals])
    for attempt in range(3):
        if full:
            index = None
        else:
            index = region_indices(snapshot, units, init_com, radius)
        data = ParticleData(snapshot, units, index)
        vlog(f"   loaded {data.ntot} particles ({'whole box' if full else f'within {radius:.1f} kpc'})")
        work = [Galaxy(**{k: (v.copy() if isinstance(v, np.ndarray) else
                              (dict(v) if isinstance(v, dict) else v))
                          for k, v in g.__dict__.items() if k in Galaxy.__dataclass_fields__})
                for g in gals]
        an = GalaxyAnalyser(data, work, p, log)
        if p["LoadParticleList"]:
            gen = int(float(read_gadget_params(p["GadgetParamfile"]).get("Generations", 4)))
            an.particle_flag = particle_list_flags(data, p["LoadParticleList"], gen, vlog)
        if p["SubUsePartAssignment"]:
            if catalog is None:
                raise ValueError("SubUsePartAssignment requires ReadSubFindGals")
            an.sf_member = subfind_membership(data, catalog, work)
        an.compute_frames()
        if full:
            break
        shift = max(np.linalg.norm(an.d.pbc_diff(w.com - c)) for w, c in zip(work, init_com))
        if shift + _region_radius(p, work) - 30.0 <= radius:
            break
        radius = 2.0 * radius
        vlog(f"   centre moved by {shift:.2f} kpc: reloading a larger region ({radius:.1f} kpc)")
    gals = work

    if p["AnalyzeGalaxies"]:
        an.analyse()
    else:
        an.masses_only()
    an.final_assignment()

    names = out.OutputNames(p["Outdir"], p["Fileflag"], snap_path)
    metal_names = (read_metal_names(p["MetalsFile"], data.nmet) if data.nmet else [])
    out.write_list_of_blocks(names, derived_names(metal_names))

    if p["AnalyzeGalaxies"]:
        for igal, g in enumerate(gals):
            if g.nparticles == 0:
                continue
            out.write_jhist(names, g, an.jhist[igal])
            out.write_info(names, p, an, igal)
            out.write_rotation_matrix(names, g)
            if data.age is not None:
                out.write_sfr(names, g, igal, an.sfr_history[igal], an.age_of_universe,
                              units.htrue)
        if p["WriteTables"]:
            out.write_tables(names, an, snapshot.header, units, int(names.snapnumber))
    if not p["WriteTables"]:
        out.write_selected_gals(names, gals, snapshot.header, units)
    if p["ToFindGalaxy"]:
        out.write_to_find_gals(names, gals, snapshot.header, units, bool(p["WriteTables"]))

    if p["WriteFullCatalog"] or p["WriteStarCatalog"]:
        write_full_catalog(names, p, an)

    resolver = QuantityResolver(data, an, p, metal_names)
    if p["WriteGalaxyCatalog"] or p["WriteParticleList"]:
        write_galaxy_catalogs(names, p, an, resolver, metal_names)
    profiles = counters = None
    if p.quantities:
        profiles, counters, volf = compute_profiles(data, an, p, resolver)
        out.write_profiles(names, p, gals, profiles, counters)

    run_extras(names, p, an, resolver, metal_names, catalog, snap_path, vlog)

    if p["WriteGalaxySnapshot"]:
        for fn in write_zoom_snapshots(names, p, an):
            vlog(f"   written {fn}")

    log(f"*** galpost done in {time.time() - t0:.1f} s")
    return dict(params=p, snapshot=snapshot, units=units, galaxies=gals, data=data,
                analyser=an, profiles=profiles, counters=counters, names=names)
