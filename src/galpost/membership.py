# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Particle membership restrictions:

* ``LoadParticleList``: analyses restricted to a list of particle IDs
  (file lines ``ID type``, types in ascending order), as in ExtractGalaxies.
  IDs of gas particles are compared with the most significant bits masked
  (they encode the star-formation generation, GADGET ``Generations``), and
  stars of the list not found among stars are searched among gas particles
  (the list may come from a later snapshot).
* ``SubUsePartAssignment``: particles restricted to those bound to the
  galaxy's subhalo according to SubFind.

Fix with respect to EG: EG compared the subhalo index of each particle with
the *position of the galaxy in the list*, which is correct only when the two
coincide; here the galaxy's own subhalo index is used.
"""

from __future__ import annotations

import warnings
from typing import List, Optional

import numpy as np

__all__ = ["read_particle_list", "particle_list_flags", "subfind_membership", "id_mask_bits"]


def id_mask_bits(generations: int = 4) -> int:
    bits = 0
    while generations > (1 << bits):
        bits += 1
    return bits


def _masked(ids: np.ndarray, bits: int) -> np.ndarray:
    nb = ids.dtype.itemsize * 8
    if bits == 0:
        return ids
    return ids & ids.dtype.type((1 << (nb - bits)) - 1)


def read_particle_list(path: str):
    ids, types = [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#%":
                continue
            w = s.split()
            ids.append(int(w[0]))
            types.append(int(w[1]) if len(w) > 1 else 0)
    ids = np.array(ids, dtype=np.uint64)
    types = np.array(types, dtype=np.int64)
    if np.any(np.diff(types) < 0):
        raise ValueError(f"{path}: particle types are not in ascending order")
    if np.any((types < 0) | (types > 5)):
        raise ValueError(f"{path}: inconsistent particle types")
    return ids, types


def particle_list_flags(data, path: str, generations: int = 4, log=print) -> np.ndarray:
    """Boolean flag over the loaded particles (EG ParticleFlag)."""
    lids, ltypes = read_particle_list(path)
    bits = id_mask_bits(generations)
    flag = np.zeros(data.ntot, dtype=bool)
    ids = data.ids
    gas = data.type_slice(0)
    gas_masked = _masked(ids[gas].astype(np.uint64), bits)
    for t in range(6):
        want = lids[ltypes == t]
        if not len(want):
            continue
        s = data.type_slice(t)
        if t == 0:
            hit = np.isin(gas_masked, _masked(want, bits))
            flag[s] |= hit
            found = np.isin(_masked(want, bits), gas_masked)
        else:
            tid = ids[s].astype(np.uint64)
            flag[s] |= np.isin(tid, want)
            found = np.isin(want, tid)
        lost = want[~found]
        if len(lost) and t == 4:
            log(f"   running again on {len(lost)} unfound stars (searched among gas)")
            flag[gas] |= np.isin(gas_masked, _masked(lost, bits))
            lost = lost[~np.isin(_masked(lost, bits), gas_masked)]
        if len(lost):
            warnings.warn(f"{len(lost)} type-{t} particles of the list not found in the "
                          "loaded region (they are outside the analysed volume)")
    log(f"   particle list {path}: {int(flag.sum())} loaded particles flagged")
    return flag


def subfind_membership(data, catalog, galaxies) -> np.ndarray:
    """Subhalo index (or -1) of each loaded particle, for the subhaloes of the
    selected galaxies."""
    subs = sorted({g.subhalo for g in galaxies if g.subhalo >= 0})
    out = np.full(data.ntot, -1, dtype=np.int64)
    if not subs:
        return out
    ids = [catalog.subhalo_ids(s).astype(np.uint64) for s in subs]
    lab = np.concatenate([np.full(len(x), s, dtype=np.int64) for s, x in zip(subs, ids)])
    ids = np.concatenate(ids)
    o = np.argsort(ids, kind="stable")
    ids, lab = ids[o], lab[o]
    pid = data.ids.astype(np.uint64)
    pos = np.searchsorted(ids, pid)
    ok = pos < len(ids)
    ok[ok] = ids[pos[ok]] == pid[ok]
    out[ok] = lab[pos[ok]]
    return out
