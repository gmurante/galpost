# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Particle data in ExtractGalaxies units, for the region around the galaxies.

Reproduces ``read_snapshot()`` of ExtractGalaxies, including its numerical
precision: arrays are float32 obtained as float32(float64(code value) * factor),
temperatures and MUPPI phase masses follow the same rules (single-phase gas
re-assigned to the cold or hot phase, MUPPI particles with a non-pressurised
hot phase demoted to single phase, SFR of single-phase particles set to zero).

Only particles within a sphere around the galaxies are loaded (the original
snapshot order is preserved, so that sums are done in the same order as in
EG); ``index[t]`` gives the position of each loaded particle in the snapshot.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from .gadget_io import Snapshot
from .units import CONST, Units

__all__ = ["ParticleData", "region_indices", "E_OVER_M_TO_T", "MU_H", "MU_C", "T_C"]

X_H = CONST.HYDROGEN_MASSFRAC
MU_H = 4.0 / (5.0 * X_H + 3.0)
MU_C = 4.0 / (3.0 * X_H + 1.0)
T_C = 300.0
E_OVER_M_TO_T = (CONST.GAMMA_MINUS1 * 1.0e50 * CONST.PROTONMASS / CONST.BOLTZMANN
                 / CONST.SOLAR_MASS)

f32 = np.float32


def _wrap_diff(d: np.ndarray, box: float) -> np.ndarray:
    """PBC_DISTANCE of EG (one wrap)."""
    half = 0.5 * box
    return np.where(d > half, d - box, np.where(d < -half, d + box, d))


def region_indices(snap: Snapshot, units: Units, centers_kpc: np.ndarray,
                   radius_kpc: float, types: Sequence[int] = range(6)) -> Dict[int, np.ndarray]:
    """Indices (per type) of particles within radius of any of the centres."""
    box = snap.header.boxsize * units.length_kpc
    out = {}
    centers_kpc = np.atleast_2d(centers_kpc)
    # selection in float32 code units, in chunks (no huge temporaries); the
    # radius gets a small safety margin, the exact geometry is applied later
    L = units.length_kpc
    boxc = np.float32(snap.header.boxsize)
    rc = np.float32(radius_kpc / L * (1.0 + 1e-5) + 1e-6)
    cc = (centers_kpc / L).astype(np.float32)
    need_wrap = [units.periodic and bool(np.any(c - rc < 0) or np.any(c + rc > boxc)) for c in cc]
    chunk = 2_000_000
    for t in types:
        if snap.npart[t] == 0:
            out[t] = np.zeros(0, dtype=np.int64)
            continue
        parts = []
        for group in snap.files:
            mm = group[0].memmap("POS", t)
            if mm is not None:
                parts.append(mm)
        ntype = sum(len(x) for x in parts)
        keep = np.zeros(ntype, bool)
        off = 0
        chunks = []
        for mm in parts:
            for s0 in range(0, len(mm), chunk):
                chunks.append((off + s0, mm[s0:s0 + chunk]))
            off += len(mm)
        for s0, pm in chunks:
            p = np.asarray(pm, dtype=np.float32)
            k = np.zeros(len(p), bool)
            for c, wrap in zip(cc, need_wrap):
                d = p - c
                if wrap:
                    d = _wrap_diff(d, boxc)
                k |= np.einsum("ij,ij->i", d, d) <= rc * rc
            keep[s0:s0 + len(p)] = k
        out[t] = np.nonzero(keep)[0]
        del parts, chunks
    return out


class ParticleData:
    """
    Arrays over the loaded particles, ordered by type as in the snapshot.

    Attributes
    ----------
    npart, nstart : per-type counts and offsets in the loaded arrays
    pos, vel, mass, ids : all particles (kpc, km/s, Msun)
    gas arrays (length npart[0]): rho, u, ne, sfr, temp, mhot, mcld, msf, ehot,
        mp (MUPPI multi-phase flag), grad (or None)
    metals : (npart[0]+npart[4], nmet) masses in Msun (gas first), or None
    im, age : initial mass and formation time of stars (or None)
    """

    def __init__(self, snap: Snapshot, units: Units,
                 index: Optional[Dict[int, np.ndarray]] = None):
        self.snap = snap
        self.units = units
        h = snap.header
        self.header = h
        self.box = h.boxsize * units.length_kpc
        self.half_box = h.boxsize / 2.0 * units.length_kpc
        if index is None:
            index = {t: None for t in range(6)}
        self.index = {t: (np.arange(int(snap.npart[t])) if index.get(t) is None
                          else np.asarray(index[t], dtype=np.int64)) for t in range(6)}
        self.npart = np.array([len(self.index[t]) for t in range(6)], dtype=np.int64)
        self.nstart = np.concatenate([[0], np.cumsum(self.npart)[:-1]])
        self.ntot = int(self.npart.sum())
        self._load()

    # ------------------------------------------------------------------
    def _read(self, name, t):
        return self.snap.read(name, t, self.index[t])

    def _all(self, name, types=range(6), fill=None):
        parts = []
        for t in types:
            if self.npart[t] == 0:
                continue
            parts.append(self._read(name, t))
        return np.concatenate(parts) if parts else np.zeros(0, dtype="f4")

    def has(self, name) -> bool:
        return name in self.snap.blocks

    def _load(self):
        u = self.units
        L, V, M = u.length_kpc, u.vel_kms, u.mass_msun

        pos = self._all("POS").astype(np.float64) * L
        if u.periodic:
            pos = np.where(pos > self.box, pos - self.box, np.where(pos < 0.0, pos + self.box, pos))
        self.pos = pos.astype(f32)
        self.vel = (self._all("VEL").astype(np.float64) * V).astype(f32)
        mass = []
        for t in range(6):
            if self.npart[t] == 0:
                continue
            if self.header.mass[t] > 0:
                m = np.full(self.npart[t], f32(self.header.mass[t]), dtype=f32)
            else:
                m = self._read("MASS", t)
            mass.append(m)
        self.mass = (np.concatenate(mass).astype(np.float64) * M).astype(f32)
        self.ids = self._all("ID")

        ng = int(self.npart[0])
        if ng:
            self.rho = (self._read("RHO", 0).astype(np.float64) * u.dens_msun_pc3).astype(f32)
            self.u = self._read("U", 0)
            self.ne = self._read("NE", 0) if self.has("NE") else None
            self.sfr = self._read("SFR", 0).copy() if self.has("SFR") else np.zeros(ng, f32)
        else:
            self.rho = self.u = self.sfr = np.zeros(0, f32)
            self.ne = None
        uu = self.u.astype(np.float64)
        if self.ne is not None:
            temp = uu * 4.0 / (3 * X_H + 1 + 4 * X_H * self.ne.astype(np.float64)) * u.temp_K
        else:
            temp = uu * (1.0 / (4.0 - X_H)) * u.temp_K
        self.temp = temp.astype(f32)

        self.muppi = self.has("EHOT")
        if self.muppi and ng:
            self.mhot = (self._read("MHOT", 0).astype(np.float64) * M).astype(f32)
            self.mcld = (self._read("MCLD", 0).astype(np.float64) * M).astype(f32)
            self.msf = (self._read("MSF", 0).astype(np.float64) * M).astype(f32)
            ehot = self._read("EHOT", 0)
            if u.htrue:
                ehot = (ehot.astype(np.float64) / self.header.hubble_param).astype(f32)
            self.ehot = ehot.copy()
            mfst = self._read("MFST", 0)
            self.mp = mfst > 0
            # EG stores MFST in an int array: the raw bits (OpenGadget writes an
            # integer number of multi-phase steps in a block declared FLOAT)
            self.mpsteps = (mfst.view(np.int32) if mfst.dtype.itemsize == 4
                            else mfst.astype(np.int32)).copy()
            self.grad = None
            if self.has("GRAD"):
                self.grad = (self._read("GRAD", 0).astype(np.float64)
                             * u.dens_msun_pc3 / L).astype(f32)
        else:
            z = np.zeros(ng, f32)
            self.mhot, self.mcld, self.msf, self.ehot = z.copy(), z.copy(), z.copy(), z.copy()
            self.mp = np.zeros(ng, bool)
            self.mpsteps = np.zeros(ng, np.int32)
            self.grad = None

        # metals (gas + stars) and initial stellar masses
        self.metals = None
        self.im = None
        self.nmet = 0
        if self.has("Zs"):
            zs = []
            for t in (0, 4):
                if self.npart[t]:
                    zs.append(self._read("Zs", t))
            if zs:
                self.metals = (np.concatenate(zs).astype(np.float64) * M).astype(f32)
                self.nmet = self.metals.shape[1]
        if self.npart[4] and self.has("iM"):
            self.im = (self._read("iM", 4).astype(np.float64) * M).astype(f32)
        self.age = self._read("AGE", 4) if (self.npart[4] and self.has("AGE")) else None

        self._phase_fix()

    def _phase_fix(self):
        """EG: demote MP particles with T_hot < 2000 K; assign single-phase gas."""
        ng = int(self.npart[0])
        if ng == 0:
            return
        mass_g = self.mass[:ng]
        if self.muppi:
            mp = self.mp
            with np.errstate(divide="ignore", invalid="ignore"):
                thot = (self.ehot.astype(np.float64) * MU_H / self.mhot.astype(np.float64)
                        * E_OVER_M_TO_T)
            demote = mp & (thot < 2000.0)
            self._set_cold(demote, mass_g)
            self.sfr[demote] = 0.0
            self.mpsteps[demote] = 0
            self.mp = mp & ~demote
            sp = ~mp
            self.sfr[sp] = 0.0
        else:
            sp = np.ones(ng, bool)
        cold = sp & (self.temp < CONST.T_THRESHOLD)
        hot = sp & ~(self.temp < CONST.T_THRESHOLD)
        self._set_cold(cold, mass_g)
        self.mhot[hot] = mass_g[hot]
        self.ehot[hot] = (self.temp[hot].astype(np.float64) / MU_H
                          * self.mhot[hot].astype(np.float64) / E_OVER_M_TO_T).astype(f32)
        self.mcld[hot] = 0.0
        self.msf[hot] = 0.0

    def _set_cold(self, sel, mass_g):
        self.mcld[sel] = mass_g[sel]
        self.mhot[sel] = 0.0
        self.msf[sel] = 0.0
        self.ehot[sel] = 0.0

    # ------------------------------------------------------------------
    def type_slice(self, t: int) -> slice:
        return slice(int(self.nstart[t]), int(self.nstart[t] + self.npart[t]))

    def ptype(self) -> np.ndarray:
        return np.repeat(np.arange(6), self.npart)

    def cold_gas(self) -> np.ndarray:
        """EG cold gas flag (multi-phase or T < 1e5 K), length npart[0]."""
        return self.mp | (self.temp < CONST.T_THRESHOLD)

    def hot_gas(self) -> np.ndarray:
        """EG hot gas flag (multi-phase or T > 1e5 K), length npart[0]."""
        return self.mp | (self.temp > CONST.T_THRESHOLD)

    def raw_block(self, name: str, t: int) -> np.ndarray:
        """Any snapshot block for type t, in code units (float32)."""
        return self._read(name, t).astype(f32)

    def pbc_diff(self, d: np.ndarray) -> np.ndarray:
        if not self.units.periodic:
            return d
        return _wrap_diff(d, self.box)

    def _no_wrap_needed(self, center) -> bool:
        """True if |pos - center| < box/2 for every loaded particle, so that the
        periodic correction is the identity (compact region, centre inside)."""
        if not self.units.periodic:
            return True
        if not hasattr(self, "_bounds"):
            if self.ntot == 0:
                self._bounds = (np.zeros(3), np.zeros(3))
            else:
                self._bounds = (self.pos.min(axis=0).astype(np.float64),
                                self.pos.max(axis=0).astype(np.float64))
        lo, hi = self._bounds
        c = np.asarray(center, dtype=np.float64)
        return bool(np.all(np.maximum(hi - c, c - lo) < 0.5 * self.box * 0.999))

    def rel_pos(self, center, idx=None) -> np.ndarray:
        """PBC-corrected positions relative to a centre (float64).

        For the full set the last result is cached (the centres do not change
        after the tracking), avoiding repeated large allocations."""
        if idx is None:
            key = tuple(np.asarray(center, dtype=np.float64))
            cache = getattr(self, "_relpos_cache", None)
            if cache is not None and cache[0] == key:
                return cache[1]
            res = self._rel_pos(center, None)
            res.setflags(write=False)
            self._relpos_cache = (key, res)
            return res
        return self._rel_pos(center, idx)

    def _rel_pos(self, center, idx=None) -> np.ndarray:
        p = self.pos if idx is None else self.pos[idx]
        d = p.astype(np.float64) - center
        if self._no_wrap_needed(center):
            return d
        return self.pbc_diff(d)
