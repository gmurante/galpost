# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Galaxy snapshots ("zoom" files, WriteGalaxySnapshot).

Improvement over ExtractGalaxies: the zoom contains ALL blocks of the
original snapshot (companion .muppi files included), in GADGET format 2 with
an INFO block, so any snapshot reader can open it.

As in EG:
* particles: those within GalaxySnapshotRadius (number in kpc or
  rmin/rmax/rvir/rgal/ralign) with the geometry of the profiles (sphere,
  cylinder or cube);
* header: copy of the original one, NPart of the galaxy, all masses in the
  MASS block (mass table set to 0), num_files = 1;
* AlignGalaxySnapshot: POS centred on the galaxy and rotated in its frame,
  VEL relative to the galaxy and rotated, GRAD rotated (galpost also treats
  VBLK as a velocity and ACCE as a vector);
* without GalaxySnapshotInCodeUnits, POS (kpc), VEL (km/s), MASS, MHOT,
  MCLD, MSF (Msun), RHO (Msun/pc^3), EHOT (1e50 erg) and GRAD are in physical
  units (h included if Htrue), all other blocks in code units.

Difference from EG: MHOT/MCLD/MSF/EHOT/SFR are the original snapshot values,
not those modified by the phase re-assignment of the analysis.
"""

from __future__ import annotations

import struct
import warnings
from typing import Dict, List

import numpy as np

from .gadget_io import BlockInfo, NTYPES

__all__ = ["write_zoom_snapshots", "zoom_radius"]

_VELOCITY_LIKE = {"VEL", "VBLK"}
_VECTORS = {"GRAD", "ACCE"}
_MASS_LIKE = {"MASS", "MHOT", "MCLD", "MSF"}

_INFO_NAME = {("f", 1): "FLOAT", ("f", 0): "FLOATN", ("d", 1): "DOUBLE", ("d", 0): "DOUBLEN",
              ("i", 1): "LONG", ("i", 0): "LONGN", ("l", 1): "LLONG", ("l", 0): "LLONGN"}


def zoom_radius(params, g) -> float:
    tag = str(params["GalaxySnapshotRadius"])
    try:
        return float(tag)
    except ValueError:
        pass
    if tag in ("rmin", "rmax", "rvir", "rgal", "ralign"):
        return float(getattr(g, tag))
    warnings.warn("invalid GalaxySnapshotRadius, assuming the virial radius")
    return float(g.rvir)


def _info_type(bi: BlockInfo) -> str:
    k = bi.dtype.kind
    s = bi.dtype.itemsize
    key = ("f" if (k == "f" and s == 4) else "d" if k == "f" else
           "i" if s == 4 else "l", 1 if bi.ndim == 1 else 0)
    return _INFO_NAME[key]


def _record(f, payload: bytes, e: str):
    n = len(payload) % (1 << 32)
    f.write(struct.pack(e + "I", n))
    f.write(payload)
    f.write(struct.pack(e + "I", n))


def _block(f, name: str, payload: bytes, e: str):
    _record(f, name.ljust(4)[:4].encode() + struct.pack(e + "I", (len(payload) + 8) % (1 << 32)), e)
    _record(f, payload, e)


def _header_bytes(raw: bytes, endian: str, npart) -> bytes:
    h = bytearray(raw[:256].ljust(256, b"\0"))
    e = endian
    h[0:24] = struct.pack(e + "6I", *[int(n) % (1 << 32) for n in npart])
    h[24:72] = struct.pack(e + "6d", *([0.0] * 6))
    h[96:120] = struct.pack(e + "6I", *[int(n) % (1 << 32) for n in npart])
    h[120:128] = h[120:124] + struct.pack(e + "i", 1)            # num_files = 1
    h[168:192] = struct.pack(e + "6I", *[int(n) >> 32 for n in npart])
    return bytes(h)


def write_zoom_snapshots(names, params, analyser) -> List[str]:
    d, an = analyser.d, analyser
    snap = d.snap
    u = d.units
    p = params
    geom = int(bool(p["SurfaceDensity"])) + 2 * int(bool(p["GridDensity"]))
    aligned = bool(p["AlignGalaxySnapshot"])
    code = bool(p["GalaxySnapshotInCodeUnits"])
    h = d.header.hubble_param

    an.belongs[:] = 0
    for igal, g in enumerate(an.gals):
        if g.nparticles == 0:
            continue
        an.assign(igal, zoom_radius(p, g), geom)

    main = snap.files[0][0]
    endian = main.endian
    info = snap.info
    blocks = [b for b in snap.blocks if b in info]
    skipped = [b for g in snap.files[0] for b in g.truncated]
    written = []
    for igal, g in enumerate(an.gals):
        if g.nparticles == 0:
            continue
        loc = {t: np.nonzero(an.belongs[d.type_slice(t)] == igal + 1)[0] for t in range(NTYPES)}
        orig = {t: d.index[t][loc[t]] for t in range(NTYPES)}
        npart = [len(loc[t]) for t in range(NTYPES)]
        fname = (f"{names.outdir}zoom.{names.flag}{names.snap}" if len(an.gals) == 1
                 else f"{names.outdir}zoom.{names.flag}gal{g.gal_id}.{names.snap}")
        with open(fname, "wb") as f:
            hb = _header_bytes(main.header.raw, endian, npart)
            _block(f, "HEAD", hb, endian)
            # INFO block (MASS now present for every type)
            entries = []
            for b in blocks:
                bi = info[b]
                pres = bi.present.copy()
                if b == "MASS":
                    pres[:] = True
                entries.append(b.ljust(4)[:4].encode() + _info_type(bi).ljust(8).encode()
                               + struct.pack(endian + "i6i", bi.ndim, *[int(x) for x in pres]))
            _block(f, "INFO", b"".join(entries), endian)
            for b in blocks:
                bi = info[b]
                types = [t for t in range(NTYPES) if (bi.present[t] or b == "MASS") and npart[t]]
                if not types:
                    continue
                parts = []
                for t in types:
                    gi = loc[t] + d.type_slice(t).start
                    parts.append(_values(b, bi, t, orig[t], gi, g, d, an, aligned, code, h))
                payload = b"".join(x.astype(bi.dtype.newbyteorder(endian), copy=False).tobytes()
                                   for x in parts)
                _block(f, b, payload, endian)
        written.append(fname)
        if skipped:
            warnings.warn(f"{fname}: incomplete blocks not written: {', '.join(skipped)}")
    return written


def _values(b, bi, t, orig_idx, loaded_idx, g, d, an, aligned, code, h):
    """Values of block b for the selected particles of type t."""
    u = d.units
    if b == "POS":
        if aligned:
            pos, _ = an.frame(g, loaded_idx)
            return (pos / u.length_kpc if code else pos).astype(np.float32)
        if code:
            return d.snap.read("POS", t, orig_idx)
        return d.pos[loaded_idx]
    if b in _VELOCITY_LIKE:
        raw = d.snap.read(b, t, orig_idx)
        v = raw.astype(np.float64) * u.vel_kms
        if aligned:
            v = (v - g.vel) @ g.rotation.T
        elif code:
            return raw
        return (v / u.vel_kms if code else v).astype(np.float32)
    if b == "MASS":
        raw = d.snap.read("MASS", t, orig_idx)
        return raw if code else (raw.astype(np.float64) * u.mass_msun).astype(np.float32)
    raw = d.snap.read(b, t, orig_idx)
    if b in _VECTORS:
        x = raw.astype(np.float64)
        if b == "GRAD" and not code:
            x = x * u.dens_msun_pc3 / u.length_kpc
        if aligned:
            x = x @ g.rotation.T
        return x.astype(np.float32)
    if code:
        return raw
    if b in _MASS_LIKE:
        return (raw.astype(np.float64) * u.mass_msun).astype(np.float32)
    if b == "RHO":
        return (raw.astype(np.float64) * u.dens_msun_pc3).astype(np.float32)
    if b == "EHOT":
        return (raw.astype(np.float64) / h).astype(np.float32) if u.htrue else raw
    return raw
