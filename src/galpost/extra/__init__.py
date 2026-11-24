# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
galpost extra analyses (beyond ExtractGalaxies). Each one is enabled by its
own tag in the parameter file (see galpost.params.TAGS, section 'extra').
"""

from __future__ import annotations

import math
import os
import warnings

import numpy as np

from . import gas as _gas
from . import globalq as _glob
from . import structure as _st
from .common import GalaxyView, radius_value

__all__ = ["run_extras", "extra_enabled", "extra_radius"]

_TAGS = ["KinematicDecomposition", "KappaRot", "StellarPopulations", "DiscStructure",
         "KinematicProfiles", "ToomreQ", "BarStrength", "SphericalFlows", "KSRelation",
         "MetallicityGradients", "DustAnalysis", "HotHaloAnalysis", "CentralBH",
         "ScalingRelations", "Satellites", "Maps"]


def extra_enabled(params) -> bool:
    return any(params[t] for t in _TAGS)


def extra_radius(params, g) -> float:
    """Radius that the loaded region must cover for the extra analyses."""
    if not extra_enabled(params):
        return 0.0
    r = g.rvir
    r = max(r, math.sqrt(2.0) * radius_value(params["ExtraRmax"], g) + float(params["ExtraZmax"]))
    if params["KappaRot"]:
        r = max(r, float(params["KappaRotRadius"]))
    if params["Maps"]:
        r = max(r, math.sqrt(3.0) * radius_value(params["MapSize"], g))
    return r


def _load_subfind(params, snap_path, catalog):
    if catalog is not None:
        return catalog
    from ..catalogs import SubfindCatalog, find_group_dir
    try:
        num = int(os.path.basename(snap_path)[-3:])
        base = params["GroupDir"] or (os.path.dirname(snap_path) or ".")
        return SubfindCatalog(find_group_dir(base, num), num)
    except Exception as e:                      # no SubFind catalogue
        warnings.warn(f"Satellites: no SubFind catalogue available ({e})")
        return None


def run_extras(names, params, analyser, resolver, metal_names, catalog, snap_path, log=print):
    p = params
    if not extra_enabled(p):
        return []
    d = analyser.d
    rows, ids = [], []
    snapnum = int(names.snapnumber) if names.snapnumber.isdigit() else -1
    sub = _load_subfind(p, snap_path, catalog) if p["Satellites"] else None
    if p["ScalingRelations"] and any(not (g.lowres_fraction == g.lowres_fraction)
                                     for g in analyser.gals):
        from ..selection import lowres_fraction
        lr = [int(x) for x in str(p["LowResTypes"]).split(",") if x.strip()]
        lowres_fraction(d.snap, analyser.gals, d.units, lr)
    for igal, g in enumerate(analyser.gals):
        if g.nparticles == 0:
            continue
        log(f"   extra analyses for galaxy {g.gal_id}")
        halo = GalaxyView(analyser, resolver, igal, extra_radius(p, g), restrict=False)
        gal = GalaxyView(analyser, resolver, igal, extra_radius(p, g), restrict=True)
        res = dict(gal_id=g.gal_id, snap=snapnum, redshift=d.header.redshift, rvir=g.rvir,
                   rgal=g.rgal, bt_eg=g.bovert if g.mbulge or g.mdisc else np.nan,
                   lowres=g.lowres_fraction)
        basic = _glob.basic_quantities(gal, resolver)
        rcurve = basic.pop("_rcurve")
        res.update(basic)
        mbulge = g.mbulge if g.mbulge else np.nan
        if p["KinematicDecomposition"]:
            kd, _, _ = _st.kinematic_decomposition(names, p, halo, g.rgal)
            res.update(dt_e=kd["dt_thr"], bt_e=kd["bt_2counter"])
            if not (mbulge == mbulge):
                mbulge = kd["mspheroid_2counter"]
        kappa = None
        if p["KappaRot"]:
            kappa = _st.kappa_rot(gal, float(p["KappaRotRadius"]))
            res.update(kappa_rot=kappa[0], kappa_co=kappa[1])
        if p["StellarPopulations"]:
            res.update(_st.stellar_populations(names, p, gal, metal_names))
        if p["DiscStructure"] or p["MetallicityGradients"] or p["ScalingRelations"]:
            ds = _st.disc_structure(names, p, gal) if p["DiscStructure"] else \
                dict(r_half_2d=res["r_half_2d"], rd=np.nan, r_half_3d=np.nan)
            res.update(ds)
            if res.get("rd", np.nan) == res.get("rd", np.nan):
                r, vc = rcurve
                res["v22"] = float(np.interp(2.2 * res["rd"], r, vc)) if len(r) else np.nan
        if p["KinematicProfiles"]:
            _st.kinematic_profiles(names, p, gal, kappa, float(p["KappaRotRadius"]))
        if p["ToomreQ"]:
            _st.toomre_q(names, p, gal)
        if p["BarStrength"]:
            res.update(_st.bar_strength(names, p, gal))
        if p["SphericalFlows"]:
            _gas.spherical_flows(names, p, halo, res["sfr"], metal_names)
        if p["KSRelation"]:
            _gas.ks_relation(names, p, gal)
        if p["MetallicityGradients"]:
            res.update(_gas.metallicity_gradients(names, p, gal, metal_names, res["r_half_2d"]))
        if p["DustAnalysis"]:
            _gas.dust(names, p, gal)
        if p["HotHaloAnalysis"]:
            res.update(_gas.hot_halo(names, p, halo, p["GadgetParamfile"]))
        if p["CentralBH"]:
            res.update(_glob.central_bh(names, p, gal, mbulge))
        if p["Satellites"]:
            _glob.satellites(names, p, halo, sub)
        if p["Maps"]:
            _st.maps(names, p, gal)
        rows.append(_glob.scaling_row(res))
        ids.append(g.gal_id)
    if p["ScalingRelations"] and rows:
        _glob.write_scaling(names, rows, bool(p["WriteTables"]), ids)
    return rows
