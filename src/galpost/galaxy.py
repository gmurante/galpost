# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Galaxy frames and standard analysis: the ``define_galaxies()`` part of
ExtractGalaxies, reproduced step by step (same particle sets, same order of
operations, same quirks, noted in comments).

The particle-to-galaxy assignment ``belongs`` (EG: BelongsToGal) is a global
array: 0 = no galaxy, i+1 = galaxy i. As in EG it accumulates across
galaxies unless explicitly reset.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .particles import ParticleData
from .selection import Galaxy
from .units import CONST

__all__ = ["GalaxyAnalyser", "SetProperties", "NBIN_SFR", "NBIN_J", "NBIN_MP"]

NBIN_SFR = 500
NBIN_J = 50
NBIN_MP = 1000
ONE_MINUS_FREST = 0.8
SQRT_G = 2.07e-3          # sqrt(G) in kpc (km/s)^2 / Msun units, as in EG
PI_EG = 3.141592          # value used by EG in angle conversions

# Scannapieco / Murante decomposition constants (EG)
_Es, _C1, _C2, _C3 = 0.0, 0.1, 0.5, 2.0
_Es1, _C4, _C5, _C6 = 0.6, 0.1, 0.1, 2.0

SET_NAMES = ["gas", "DM", "cold gas", "hot gas", "star", "cold gas + star"]


@dataclass
class SetProperties:
    number: int = 0
    mass: float = 0.0
    inertia: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    eigenvalues: np.ndarray = field(default_factory=lambda: np.zeros(3))
    eigenvectors: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    specific_spin: np.ndarray = field(default_factory=lambda: np.zeros(3))
    mod_specific_spin: float = 0.0
    spin: np.ndarray = field(default_factory=lambda: np.zeros(3))
    mod_spin: float = 0.0
    half_mass_radius: float = 0.0


def _eigh_sorted(m: np.ndarray):
    """Eigen-decomposition with ascending eigenvalues (GSL_EIGEN_SORT_VAL_ASC).

    NB: eigenvector signs are implementation dependent (GSL vs LAPACK), so the
    inertia-tensor frame may differ from EG by axis flips; the angular-momentum
    frame (AlignWithAngularMomentum) is unaffected.
    """
    w, v = np.linalg.eigh(m)
    o = np.argsort(w)
    return w[o], v[:, o]


class GalaxyAnalyser:
    def __init__(self, data: ParticleData, galaxies: List[Galaxy], params, log=print):
        self.d = data
        self.gals = galaxies
        self.p = params
        self.log = log if params["Verbose"] else (lambda *a, **k: None)
        self.belongs = np.zeros(data.ntot, dtype=np.int32)
        self.cut = self._parse_cut(str(params["Cut"]))
        self.circularity = np.zeros(int(data.npart[0] + data.npart[4]))
        self.component = np.zeros(data.ntot, dtype=np.int32)
        self.analysis: List[List[SetProperties]] = []
        self.jhist = []
        self.sfr_history = []
        self.cmtype = int(params["CmType"])
        self.comoving = data.units.comoving
        self.particle_flag: Optional[np.ndarray] = None   # LoadParticleList
        self.sf_member: Optional[np.ndarray] = None       # SubUsePartAssignment
        for g in galaxies:
            g.vel = np.zeros(3)
            g.rotation = np.eye(3)
            g.inertia = np.zeros((3, 3))
            g.eigen = np.zeros(3)
            g.am = np.zeros(3)
            g.angle = 0.0
            g.swap = 0
            g.stellar_mass = g.cold_gas_mass = g.tot_gas_mass = 0.0
            g.mbulge = g.mdisc = g.mhalo = g.bovert = g.sfr = 0.0
            g.nparticles_type = np.zeros(6, dtype=np.int64)
            g.nparticles = 0

    # ------------------------------------------------------------------
    @staticmethod
    def _parse_cut(tag: str) -> float:
        try:
            return float(tag)
        except ValueError:
            return {"rmin": -1.0, "rmax": -2.0, "rvir": -3.0, "rgal": -4.0,
                    "ralign": -5.0}.get(tag, 0.0)

    def rvert(self, g: Galaxy, r: float) -> float:
        c = self.cut
        if c == 0.0:
            return r
        return {-1.0: g.rmin, -2.0: g.rmax, -3.0: g.rvir, -4.0: g.rgal,
                -5.0: g.ralign}.get(c, r if r <= c else c)

    def cut_value(self, g: Galaxy) -> float:
        """'rvert' used for mass flows (derived OUTF/INFL)."""
        c = self.cut
        return {-1.0: g.rmin, -2.0: g.rmax, -3.0: g.rvir, -4.0: g.rgal,
                -5.0: g.ralign}.get(c, c)

    def cosmological_time(self, a) -> np.ndarray:
        h = self.d.header
        hub = h.hubble_param if self.d.units.htrue else 1.0
        if not self.comoving:
            return a
        return (9.776 / hub / 1.5 / math.sqrt(h.omega_lambda)
                * np.arcsinh(math.sqrt(h.omega_lambda / h.omega0) * np.power(a, 1.5)))

    # ---- geometry --------------------------------------------------------
    def frame(self, g: Galaxy, idx=None):
        """Positions and velocities in the galaxy frame (EG position()).
        The full-set result is cached per galaxy state (centre, velocity, rotation)."""
        d = self.d
        if idx is None:
            key = (id(g), tuple(g.com), tuple(np.ravel(g.vel)), tuple(np.ravel(g.rotation)))
            c = getattr(self, "_frame_cache", None)
            if c is not None and c[0] == key:
                return c[1]
        dp = d.rel_pos(g.com, idx)
        v = (d.vel if idx is None else d.vel[idx]).astype(np.float64) - g.vel
        res = (dp @ g.rotation.T, v @ g.rotation.T)
        if idx is None:
            self._frame_cache = (key, res)
        return res

    def assign(self, igal: int, r: float, geometry: int, reset: bool = False):
        """EG assign_BelongsToGal: 0 sphere, 1 cylinder (height from Cut), 2 cube."""
        if reset:
            self.belongs[:] = 0
        g = self.gals[igal]
        rv = self.rvert(g, r)
        if geometry == 0:
            dp = self.d.rel_pos(g.com)
            cond = np.einsum("ij,ij->i", dp, dp) <= r * r
        else:
            pos, _ = self.frame(g)
            if geometry == 1:
                cond = (pos[:, 0] ** 2 + pos[:, 1] ** 2 <= r * r) & (np.abs(pos[:, 2]) < rv)
            else:
                cond = (np.abs(pos[:, 0]) < r) & (np.abs(pos[:, 1]) < r) & (np.abs(pos[:, 2]) < rv)
        if self.sf_member is not None:
            cond &= self.sf_member == g.subhalo
        if self.particle_flag is not None:
            cond &= self.particle_flag
        self.belongs[cond] = igal + 1
        self.log(f"   BelongsToGal: {int(cond.sum())} particles within r={r:f} with "
                 f"geometry {geometry} selected for galaxy {igal}")

    # ---- particle sets -----------------------------------------------------
    def use_cmtype(self):
        """usePart for the galaxy frame and myMass (always the total mass here)."""
        d = self.d
        use = np.zeros(d.ntot, bool)
        ct = self.cmtype
        if 0 <= ct < 6:
            use[d.type_slice(ct)] = True
        else:
            use[d.type_slice(0)] = d.cold_gas()
            if ct != 6:
                use[d.type_slice(4)] = True
        return use, d.mass

    def use_set(self, kind: int):
        """The six EG analysis sets and the corresponding myMass."""
        d = self.d
        use = np.zeros(d.ntot, bool)
        mass = d.mass
        gs = d.type_slice(0)
        if kind == 0:
            use[gs] = True
        elif kind == 1:
            use[d.type_slice(1)] = True
        elif kind == 2:
            use[gs] = d.cold_gas()
            mass = d.mass.copy()
            mass[gs] = d.mcld
        elif kind == 3:
            use[gs] = d.hot_gas()
            mass = d.mass.copy()
            mass[gs] = d.mhot
        elif kind == 4:
            use[d.type_slice(4)] = True
        elif kind == 5:
            use[gs] = d.cold_gas()
            use[d.type_slice(4)] = True
        return use, mass

    # ---- elementary computations --------------------------------------------
    def inertia(self, igal: int, use, mymass):
        """EG compute_inertia: mass, inertia (1e10 Msun kpc^2), specific AM (box frame)."""
        g = self.gals[igal]
        sel = use & (self.belongs == igal + 1)
        n = int(sel.sum())
        m = mymass[sel].astype(np.float64)
        dp = self.d.rel_pos(g.com, sel)
        dv = self.d.vel[sel].astype(np.float64) - g.vel
        mass = float(np.sum(m)) if n else 0.0
        r2 = np.einsum("ij,ij->i", dp, dp)
        I = np.zeros((3, 3))
        w = m * 1.0e-10
        for i in range(3):
            I[i, i] = np.sum(w * (r2 - dp[:, i] * dp[:, i]))
            for j in range(i + 1, 3):
                I[i, j] = I[j, i] = -np.sum(w * dp[:, i] * dp[:, j])
        am = np.array([np.sum(m * (dp[:, 1] * dv[:, 2] - dp[:, 2] * dv[:, 1])),
                       np.sum(m * (dp[:, 2] * dv[:, 0] - dp[:, 0] * dv[:, 2])),
                       np.sum(m * (dp[:, 0] * dv[:, 1] - dp[:, 1] * dv[:, 0]))])
        with np.errstate(invalid="ignore", divide="ignore"):
            am = am / mass
        return n, mass, I, am

    def half_mass_radius(self, igal: int, use, mymass, number=0, mass=0.0) -> float:
        """EG compute_half_mass_radius (including its off-by-one: the radius of
        the particle *after* the one that crosses half of the mass)."""
        g = self.gals[igal]
        sel = np.nonzero(use & (self.belongs == igal + 1))[0]
        if number == 0:
            number = len(sel)
            mass = float(np.sum(mymass[sel].astype(np.float64)))
        if number == 0:
            return 0.0
        dp = self.d.rel_pos(g.com, sel)
        d2 = np.einsum("ij,ij->i", dp, dp)
        o = np.argsort(d2, kind="stable")
        cm = np.cumsum(mymass[sel][o].astype(np.float64))
        k = int(np.searchsorted(cm, mass / 2.0, side="left"))   # first with cum >= half
        j = min(k + 1, number - 1)
        return math.sqrt(d2[o][j])

    # ---- part one: centres and frames ----------------------------------------
    def count_particles(self):
        d = self.d
        types = d.ptype()
        for g in self.gals:
            dp = d.rel_pos(g.com)
            inside = np.einsum("ij,ij->i", dp, dp) <= g.rmax * g.rmax
            if self.particle_flag is not None:
                inside &= self.particle_flag
            g.nparticles_type = np.bincount(types[inside], minlength=6)
            g.nparticles = int(inside.sum())

    def compute_frames(self):
        d = self.d
        p = self.p
        use, mymass = self.use_cmtype()

        if self.gals and self.gals[0].compute_cm:
            if len(self.gals) > 1:
                raise ValueError(f"COMP_CM required but {len(self.gals)} galaxies are given")
            if d.units.periodic:
                raise ValueError("COMP_CM cannot be used with periodic boundaries")
            m = mymass[use].astype(np.float64)
            g = self.gals[0]
            g.com = (d.pos[use].astype(np.float64) * m[:, None]).sum(0) / m.sum()
            g.vel = (d.vel[use].astype(np.float64) * m[:, None]).sum(0) / m.sum()

        self.count_particles()

        for igal, g in enumerate(self.gals):
            self.log(f"   GALAXY {igal}, ID {g.gal_id}")
            if g.nparticles == 0:
                self.log(f"WARNING: galaxy {igal} has no particles, skipping...")
                continue

            if g.pot:
                if not d.has("POT"):
                    print("WARNING: Pot search impossible because POT block is not present")
                else:
                    pot = np.concatenate([d.raw_block("POT", t) for t in range(6) if d.npart[t]])
                    dp = d.rel_pos(g.com)
                    r = np.sqrt(np.einsum("ij,ij->i", dp, dp))
                    cand = use & (r < g.rmax)
                    i = np.nonzero(cand)[0][np.argmin(pot[cand])]
                    g.com = d.pos[i].astype(np.float64)

            if g.track_cm:
                cm = g.com.astype(np.float64).copy()
                r2 = g.ralign * g.ralign
                diff, counter = 1.0e6, 0
                # EG: already assigned particles excluded; candidates kept in
                # snapshot order so that sums are done in the same order
                cand = np.nonzero(use & (self.belongs == 0))[0]
                mc = mymass[cand].astype(np.float64)
                while diff > 1.0e-3:
                    dp = d.rel_pos(cm, cand)
                    sel = np.einsum("ij,ij->i", dp, dp) <= r2
                    m = mc[sel]
                    shift = (dp[sel] * m[:, None]).sum(0) / m.sum()
                    diff = math.sqrt(float(np.dot(shift, shift)))
                    counter += 1
                    cm = cm + shift
                    if d.units.periodic:
                        cm = np.where(cm > d.box, cm - d.box, np.where(cm < 0.0, cm + d.box, cm))
                self.log(f"   CoM recomputed with {counter} iterations:")
                self.log("   Old: [%8f,%8f,%8f]" % tuple(g.com))
                self.log("   New: [%8f,%8f,%8f]" % tuple(cm))
                g.com = cm
                g.track_iterations = counter

            if p["AquilaAlignment"]:
                self.assign(igal, g.rgal, 0)
                g.ralign = 3.0 * self.half_mass_radius(igal, use, mymass)

            self.assign(igal, g.ralign, 0)
            sel = use & (self.belongs == igal + 1)
            m = mymass[sel].astype(np.float64)
            g.vel = (d.vel[sel].astype(np.float64) * m[:, None]).sum(0) / m.sum()

            n, mass, I, am = self.inertia(igal, use, mymass)
            g.inertia, g.am = I, am
            w, v = _eigh_sorted(I)
            g.eigen = w
            rot = np.zeros((3, 3))
            if p["AlignWithAngularMomentum"] or p["AquilaAlignment"]:
                z = am / math.sqrt(float(np.dot(am, am)))
                rot[2] = z
                nrm = math.sqrt(z[2] * z[2] + z[0] * z[0])
                rot[0] = [z[2] / nrm, 0.0, -z[0] / nrm]
                rot[1] = [z[1] * rot[0, 2] - z[2] * rot[0, 1],
                          z[2] * rot[0, 0] - z[0] * rot[0, 2],
                          z[0] * rot[0, 1] - z[1] * rot[0, 0]]
                g.angle = math.acos(abs(float(np.dot(z, v[:, 2])))) * 180.0 / PI_EG
                vm = rot @ am
            else:
                rot = v.T.copy()
                vm = rot @ am
                g.angle = math.acos(abs(vm[2] / math.sqrt(float(np.dot(vm, vm))))) * 180.0 / PI_EG
                g.swap = int(vm[2] < 0)
                if g.swap:
                    rot *= -1
                vp = ((rot[1, 0] * rot[2, 1] - rot[2, 0] * rot[1, 1]) * rot[0, 2]
                      + (rot[2, 0] * rot[0, 1] - rot[0, 0] * rot[2, 1]) * rot[1, 2]
                      + (rot[0, 0] * rot[1, 1] - rot[1, 0] * rot[0, 1]) * rot[2, 2])
                if vp < 0:
                    rot[[0, 1]] = rot[[1, 0]]
                g.swap += 2 * int(vp < 0)
            if p["RotateGalaxyFrame"] == 1:
                rot = rot[[2, 0, 1]]
            elif p["RotateGalaxyFrame"] == -1:
                rot = rot[[1, 2, 0]]
            g.rotation = rot if g.align else np.eye(3)
            g.frame_number = n
            self.log(f"   Computation performed on {n} particles")
            self.log("   Angular momentum in the box frame:   [  %f %f %f  ]" % tuple(am))
            self.log("   Angular momentum in the galaxy frame:[  %f %f %f  ]" % tuple(vm))
            self.log(f"   Inclination of angular momentum with respect to z-axis of "
                     f"inertia tensor: {g.angle:f} degrees")

    # ---- part two: analysis --------------------------------------------------
    def analyse(self):
        d = self.d
        age_of_u = float(self.cosmological_time(1.0))
        self.age_of_universe = age_of_u
        for igal, g in enumerate(self.gals):
            sets = [SetProperties() for _ in range(12)]
            self.analysis.append(sets)
            self.jhist.append(np.zeros((NBIN_J, NBIN_J + 1)))
            self.sfr_history.append(np.zeros(NBIN_SFR))
            if g.nparticles == 0:
                continue
            self.log(f"   Computing masses and inertia tensors for galaxy {igal}, ID {g.gal_id}")
            for turn in range(2):
                if turn and g.rgal == g.rvir:
                    break
                r = g.rgal if turn == 0 else g.rvir
                self.assign(igal, r, 0, reset=True)
                for kind in range(6):
                    use, mymass = self.use_set(kind)
                    sp = sets[turn * 6 + kind]
                    sp.number, sp.mass, sp.inertia, sp.specific_spin = self.inertia(igal, use, mymass)
                    sp.half_mass_radius = self.half_mass_radius(igal, use, mymass, sp.number, sp.mass)
                    sp.spin = sp.specific_spin * sp.mass
                    sp.mod_spin = math.sqrt(float(np.dot(sp.spin, sp.spin)))
                    sp.mod_specific_spin = math.sqrt(float(np.dot(sp.specific_spin, sp.specific_spin)))
                    sp.eigenvalues, sp.eigenvectors = _eigh_sorted(sp.inertia)
                    if turn == 0 and kind == 4 and d.age is not None:
                        self.sfr_history[igal] = self._sfr_history(igal, mymass, age_of_u)
                if turn == 0:
                    mb, md, mh = self.bulge_disc(igal)
                    g.mbulge, g.mdisc, g.mhalo = mb, md, mh
                    with np.errstate(invalid="ignore", divide="ignore"):
                        g.bovert = mb / (mb + md)
                    gs = d.type_slice(0)
                    g.sfr = float(np.sum(d.sfr[self.belongs[gs] == igal + 1].astype(np.float64)))
            g.stellar_mass = sets[4].mass
            g.tot_gas_mass = sets[0].mass
            g.cold_gas_mass = sets[2].mass

    def masses_only(self):
        """EG branch without AnalyzeGalaxies (note: cold gas tested with Mhot>0)."""
        d = self.d
        gs = d.type_slice(0)
        ss = d.type_slice(4)
        for igal, g in enumerate(self.gals):
            if g.nparticles == 0:
                continue
            self.assign(igal, g.rgal, 0)
            b = self.belongs
            ing = b[gs] == igal + 1
            mg = d.mass[gs].astype(np.float64)
            g.tot_gas_mass = float(np.sum(mg[ing]))
            coldsel = ing & ((d.mhot > 0) | (d.temp < CONST.T_THRESHOLD))
            g.cold_gas_mass = float(np.sum(mg[coldsel]))
            g.stellar_mass = float(np.sum(d.mass[ss][b[ss] == igal + 1].astype(np.float64)))

    def final_assignment(self):
        """Particles of each galaxy for profiles and derived blocks."""
        p = self.p
        geom = int(bool(p["SurfaceDensity"])) + 2 * int(bool(p["GridDensity"]))
        self.belongs[:] = 0
        for igal, g in enumerate(self.gals):
            self.assign(igal, g.rmax, geom)

    # ---- details -------------------------------------------------------------
    def _sfr_history(self, igal, mymass, age_of_u):
        d = self.d
        ss = d.type_slice(4)
        sel = self.belongs[ss] == igal + 1
        t = self.cosmological_time(d.age[sel].astype(np.float64))
        b = (t / age_of_u * NBIN_SFR).astype(np.int64)
        b = np.minimum(b, NBIN_SFR - 1)
        if d.nmet > 1 and d.im is not None:
            w = d.im[sel].astype(np.float64)
            factor = 1.0e9 * age_of_u / NBIN_SFR
        else:
            w = mymass[ss][sel].astype(np.float64)
            factor = 1.0e9 * age_of_u / NBIN_SFR * ONE_MINUS_FREST
        h = np.bincount(b, weights=w, minlength=NBIN_SFR)[:NBIN_SFR]
        return h / factor

    def bulge_disc(self, igal):
        """EG compute_stellar_B_over_T: circularities, jz/jcirc histograms and
        bulge/disc/halo masses within rgal (decomposition methods -1, 0, 1, 2)."""
        d = self.d
        g = self.gals[igal]
        mymass = d.mass
        inside = np.nonzero(self.belongs == igal + 1)[0]
        dp = d.rel_pos(g.com, inside)
        r_all = np.sqrt(np.einsum("ij,ij->i", dp, dp))
        mp = np.zeros(NBIN_MP + 1)
        ib = (r_all / g.rgal * NBIN_MP).astype(np.int64)
        np.add.at(mp, np.minimum(ib, NBIN_MP), mymass[inside].astype(np.float64))
        mp = np.cumsum(mp)
        dr = g.rgal / NBIN_MP

        def jcirc_of(r):
            rb = (r / g.rgal * NBIN_MP).astype(np.int64)
            out = np.empty_like(r)
            big = rb >= NBIN_MP - 1
            out[big] = mp[NBIN_MP - 1]
            s = ~big
            rbs = rb[s]
            out[s] = (mp[rbs + 1] * (r[s] - rbs * dr) + mp[rbs] * ((rbs + 1) * dr - r[s])) / dr
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.sqrt(out / r) * r * SQRT_G

        jh = self.jhist[igal]
        jh[:] = 0.0
        # gas circularities
        gs = d.type_slice(0)
        gi = np.nonzero(self.belongs[gs] == igal + 1)[0]
        if len(gi):
            pos, vel = self.frame(g, gi + gs.start)
            r = np.sqrt(np.einsum("ij,ij->i", pos, pos))
            jz = pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0]
            with np.errstate(divide="ignore", invalid="ignore"):
                self.circularity[gi] = jz / jcirc_of(r)

        mb = md = mh = 0.0
        ss = d.type_slice(4)
        si = np.nonzero(self.belongs[ss] == igal + 1)[0]
        method = int(self.p["DecompositionMethod"])
        if len(si):
            pos, vel = self.frame(g, si + ss.start)
            r = np.sqrt(np.einsum("ij,ij->i", pos, pos))
            rc = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
            jx = pos[:, 1] * vel[:, 2] - pos[:, 2] * vel[:, 1]
            jy = pos[:, 2] * vel[:, 0] - pos[:, 0] * vel[:, 2]
            jz = pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0]
            with np.errstate(divide="ignore", invalid="ignore"):
                eps = jz / jcirc_of(r)
            self.circularity[int(d.npart[0]) + si] = eps
            m = mymass[ss][si].astype(np.float64)
            jb = ((2.0 + eps) / 4.0 * NBIN_J)
            jb = np.where(np.isfinite(jb), jb, 0).astype(np.int64)
            jb = np.clip(jb, 0, NBIN_J - 1)
            inr = r <= g.rgal
            np.add.at(jh[:, 0], jb[inr], m[inr])
            rb = (r / g.rgal * NBIN_J).astype(np.int64) + 1
            ok = rb <= NBIN_J
            np.add.at(jh, (jb[ok], rb[ok]), m[ok])

            comp = np.zeros(len(si), dtype=np.int32)
            w = inr
            if method == -1:
                bul = w & (eps < 0) & (r < 0.2 * g.rgal)
                hal = w & (eps < 0) & ~(r < 0.2 * g.rgal)
                dis = w & ~(eps < 0)
                comp[bul], comp[hal], comp[dis] = 1, 4, 2
                mb, mh, md = m[bul].sum(), m[hal].sum(), m[dis].sum()
            elif method == 0:
                bul = w & (eps < 0)
                dis = w & ~(eps < 0)
                comp[bul], comp[dis] = 1, 2
                mb, md = m[bul].sum(), m[dis].sum()
            elif method in (1, 2):
                env = np.exp(-(rc - _C2) / _C3)
                disc = (eps < _Es - _C1 - env) | (eps > _Es + _C1 + env)
                c = np.where(disc, 2, 1)
                if method == 2:
                    env2 = np.exp(-(rc - _C5) / _C6)
                    c = np.where((eps < _Es + _C1 + env) & (eps > _Es1 - _C4 - env2), 3, c)
                    c = np.where((c == 1) & (eps > _Es1 + _C4 + env2), 4, c)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        jn = jz / np.sqrt(jx * jx + jy * jy + jz * jz)
                    c = np.where(((c == 2) | (c == 3)) & (jn < 0.7), 5, c)
                    c = np.where((c == 1) & (r > 10.0), 6, c)
                comp = np.where(w, c, 0).astype(np.int32)
                mb, md = m[w & (comp == 1)].sum(), m[w & (comp == 2)].sum()
            self.component[ss.start + si[w]] = comp[w]
            self.component[ss.start + si[~w]] = 0

        if method in (0, -1):
            if mb + mh < md:
                md -= mb + mh
                mb += mb
                mh += mh
            else:
                mb = mb + md + mh
                md = 0.0
                mh = 0.0
        return float(mb), float(md), float(mh)
