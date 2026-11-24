# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Derived blocks and profiles (ExtractGalaxies ``build_derived_blocks``,
``build_pointers`` and ``ComputeProfiles``).

Name resolution of a profile quantity / weight, as in EG:
  1. blocks known to libreadgad (``KNOWN_BLOCKS``): the preloaded ones are in
     EG units (MASS, RHO, MHOT, ...), the others are read in code units;
  2. derived blocks (``derived_names``);
  3. galpost extension: any other block present in the snapshot (code units).
Values are float32, as in EG.

Fix with respect to EG: OUTF/INFL of particles with zero flow are 0 (EG kept
the value of the previous particle).
"""

from __future__ import annotations

import math
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np

from .particles import ParticleData, MU_H, MU_C, T_C, E_OVER_M_TO_T
from .units import CONST

__all__ = ["KNOWN_BLOCKS", "derived_names", "QuantityResolver", "compute_profiles",
           "read_metal_names"]

f32 = np.float32
KM_PER_KPC = 3.085678e16
RHO_TO_N = 40.6321263
SMALL = 1.0e-10
PI_EG_PROF = CONST.PI

KNOWN_BLOCKS = ["HEAD", "POS", "VEL", "ID", "MASS", "U", "TTTT", "RHO", "NE", "NH", "HSML",
                "SFR", "AGE", "Z", "Zs", "iM", "ZAGE", "ZALV", "CLDX", "TSTP", "POT", "ACCE",
                "ENDT", "CLDX", "IDU", "HOTT", "MHOT", "MCLD", "EHOT", "MSF", "MFST", "NMF",
                "EOUT", "EREC", "EOLD", "TDYN", "SFRo", "CLCK", "Egy0", "GRAD", "BHMA", "BHMD",
                "BHPC", "ACRB"]
_PRELOADED = {"POS", "VEL", "MASS", "RHO", "U", "NE", "SFR", "MHOT", "MCLD", "EHOT", "MSF", "AGE"}

_DERIVED_BASE = ["TEMP", "THOT", "HOTN", "CLDN", "FILH", "PRES", "MMOL", "OUTF", "INFL",
                 "VX", "VY", "VZ", "SVX", "SVY", "SVZ", "MVX", "MVY", "MVZ",
                 "VR", "VT", "VP", "SVR", "SVT", "SVP", "MVR", "MVT", "MVP",
                 "VRp", "VTp", "VZp", "SVRp", "SVTp", "SVZp", "MVRp", "MVTp", "MVZp",
                 "SV", "MV", "SVp", "MVp", "HZ", "SHZ", "MHZ", "ENTR"]
_GAS_ONLY = {"TEMP", "THOT", "HOTN", "CLDN", "FILH", "PRES", "MMOL", "ENTR"}


def read_metal_names(path: str, nmet: int) -> List[str]:
    names = []
    with open(path) as f:
        for line in f:
            w = line.split()
            if w:
                names.append(w[0][:3])
            if len(names) == nmet:
                break
    if len(names) < nmet:
        raise ValueError(f"{path}: {len(names)} elements listed, {nmet} in the snapshot")
    return names


def derived_names(metal_names: List[str]) -> List[str]:
    out = list(_DERIVED_BASE)
    if metal_names:
        out += ["EMet", "ZMet"]
        out += ["E" + m for m in metal_names]
        out += [("Y" if i == 0 else "Z") + m for i, m in enumerate(metal_names)]
    return out


class QuantityResolver:
    """Builds the per-particle arrays needed by the profiles (lazily)."""

    def __init__(self, data: ParticleData, analyser, params, metal_names: List[str]):
        self.d = data
        self.a = analyser
        self.p = params
        self.metal_names = metal_names
        self.derived = derived_names(metal_names)
        self._cache: Dict[str, np.ndarray] = {}
        self.surface = bool(params["SurfaceDensity"] or params["GridDensity"])

    # ---- derived blocks over the galaxy particles -------------------------
    def _frame_quantities(self):
        """Velocities/positions in the galaxy frames for all assigned particles."""
        if "_frame" in self._cache:
            return self._cache["_frame"]
        d, a = self.d, self.a
        n = d.ntot
        pos = np.zeros((n, 3))
        vel = np.zeros((n, 3))
        for igal, g in enumerate(a.gals):
            sel = np.nonzero(a.belongs == igal + 1)[0]
            if len(sel):
                pos[sel], vel[sel] = a.frame(g, sel)
        self._cache["_frame"] = (pos, vel)
        return pos, vel

    def _velocity_block(self, name):
        pos, vel = self._frame_quantities()
        x, y, z = pos.T
        vx, vy, vz = vel.T
        with np.errstate(divide="ignore", invalid="ignore"):
            radius = np.sqrt(x * x + y * y + z * z)
            rp = np.sqrt(x * x + y * y)
            vr = (x * vx + y * vy + z * vz) / radius
            vthe = (-x * z * vx - y * z * vy + (x * x + y * y) * vz) / radius / rp
            vphi = (-y * vx + x * vy) / rp
            vrp = (x * vx + y * vy) / rp
        table = {
            "VX": vx, "VY": vy, "VZ": vz, "SVX": vx * vx, "SVY": vy * vy, "SVZ": vz * vz,
            "MVX": np.abs(vx), "MVY": np.abs(vy), "MVZ": np.abs(vz),
            "VR": vr, "VT": vthe, "VP": vphi, "SVR": vr * vr, "SVT": vthe * vthe,
            "SVP": vphi * vphi, "MVR": np.abs(vr), "MVT": np.abs(vthe), "MVP": np.abs(vphi),
            "VRp": vrp, "VTp": vphi, "VZp": vz, "SVRp": vrp * vrp, "SVTp": vphi * vphi,
            "SVZp": vz * vz, "MVRp": np.abs(vrp), "MVTp": np.abs(vphi), "MVZp": np.abs(vz),
            "SV": vx * vx + vy * vy + vz * vz, "MV": np.sqrt(vx * vx + vy * vy + vz * vz),
            "SVp": vx * vx + vy * vy, "MVp": np.sqrt(vx * vx + vy * vy),
            "HZ": z, "SHZ": z * z, "MHZ": np.abs(z),
        }
        out = table[name]
        out = np.where(self.a.belongs > 0, out, 0.0)
        return out.astype(f32)

    def _flows(self):
        if "_flows" in self._cache:
            return self._cache["_flows"]
        d, a = self.d, self.a
        pos, vel = self._frame_quantities()
        outf = np.zeros(d.ntot)
        infl = np.zeros(d.ntot)
        for igal, g in enumerate(a.gals):
            sel = np.nonzero(a.belongs == igal + 1)[0]
            rv = a.cut_value(g)
            if rv <= 0 or not len(sel):
                continue
            x, y, z = pos[sel].T
            if self.surface:
                inside = (z >= rv) & (z < 2.0 * rv)
                v = vel[sel, 2]
            else:
                radius = np.sqrt(x * x + y * y + z * z)
                inside = (radius >= rv) & (radius < 2.0 * rv)
                with np.errstate(divide="ignore", invalid="ignore"):
                    v = (x * vel[sel, 0] + y * vel[sel, 1] + z * vel[sel, 2]) / radius
            flow = np.where(inside, d.mass[sel].astype(np.float64) * v / rv / KM_PER_KPC
                            * CONST.SEC_PER_YEAR, 0.0)
            outf[sel] = np.where(flow > 0, flow, 0.0)
            infl[sel] = np.where(flow < 0, -flow, 0.0)
        self._cache["_flows"] = (outf.astype(f32), infl.astype(f32))
        return self._cache["_flows"]

    def _gas_block(self, name):
        d = self.d
        mp = d.mp
        mcld = d.mcld.astype(np.float64)
        mhot = d.mhot.astype(np.float64)
        ehot = d.ehot.astype(np.float64)
        rho = d.rho.astype(np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            mgas = mcld + mhot
            frac_h = mhot / mgas
            frac_c = mcld / mgas
            thot = np.where(mhot > 0, ehot * MU_H / mhot * E_OVER_M_TO_T, 0.0)
            filh = np.where(mhot > 0, 1.0 / (1.0 + frac_c / frac_h * MU_H / MU_C * T_C / thot), 0.0)
            fill_c = 1.0 - filh
            vol = mgas / rho
            rho_h = np.where(mhot > 0, mhot / (vol * filh), 0.0)
            rho_c = np.where(fill_c > 0, mcld / (vol * fill_c), 0.0)
            hotn = rho_h / MU_H * RHO_TO_N
            cldn = rho_c / MU_C * RHO_TO_N
            pres_sph = d.u.astype(np.float64) * rho * RHO_TO_N * d.units.temp_K
            pres = np.where(mp, T_C * cldn, pres_sph)
            mmol = np.where(mp, mcld / (1.0 + float(self.p["Pressure"]) / pres), 0.0)
            entr = d.temp.astype(np.float64) / np.power(rho, 0.66666666666667)
        table = {"THOT": thot, "HOTN": hotn, "CLDN": cldn, "FILH": filh, "PRES": pres,
                 "MMOL": mmol, "ENTR": entr, "TEMP": d.temp.astype(np.float64)}
        return table[name].astype(f32)

    def _metal_block(self, name):
        d = self.d
        if d.metals is None:
            raise KeyError(name)
        ng, ns = int(d.npart[0]), int(d.npart[4])
        n = d.ntot
        out = np.zeros(n, dtype=f32)
        assigned = self.a.belongs > 0
        idx_all = np.r_[np.arange(ng), d.nstart[4] + np.arange(ns)]
        denom = np.r_[d.mass[:ng].astype(np.float64),
                      d.im.astype(np.float64) if d.im is not None else np.ones(ns)]
        met = d.metals.astype(np.float64)
        if name in ("EMet", "ZMet"):
            val = met[:, 1:].sum(axis=1)
            if name == "ZMet":
                val = val / denom
        else:
            el = name[1:]
            k = self.metal_names.index(el)
            val = met[:, k]
            if name[0] in "YZ":
                val = val / denom
        out[idx_all] = val.astype(f32)
        out[~assigned] = 0.0
        return out

    def derived_block(self, name) -> np.ndarray:
        """Full-length (all loaded particles) float32 array; gas-only blocks are
        returned with zeros for the other types."""
        if name in self._cache:
            return self._cache[name]
        d = self.d
        if name in _GAS_ONLY:
            arr = np.zeros(d.ntot, dtype=f32)
            arr[d.type_slice(0)] = self._gas_block(name)
        elif name in ("OUTF", "INFL"):
            arr = self._flows()[0 if name == "OUTF" else 1]
        elif name in _DERIVED_BASE:
            arr = self._velocity_block(name)
        else:
            arr = self._metal_block(name)
        self._cache[name] = arr
        return arr

    # ---- resolution ----------------------------------------------------------
    def resolve(self, label: str, ptype: int, is_weight: bool
                ) -> Tuple[Optional[np.ndarray], Optional[float]]:
        """Array (over particles of ptype) and volume factor override (or None)."""
        d = self.d
        s = d.type_slice(ptype)
        if label in KNOWN_BLOCKS:
            if label in _PRELOADED:
                volf = 1.0 if (label == "SFR" and not is_weight) else None
                if label == "MASS":
                    return d.mass[s], None
                if label in ("POS", "VEL"):
                    return None, None       # EG leaves these unlinked
                if label == "AGE":
                    return (d.age if ptype == 4 else None), None
                gas = {"RHO": d.rho, "U": d.u, "NE": d.ne, "SFR": d.sfr, "MHOT": d.mhot,
                       "MCLD": d.mcld, "EHOT": d.ehot, "MSF": d.msf}[label]
                if ptype != 0:
                    warnings.warn(f"{label} is defined only for gas (type 0)")
                    return None, None
                return gas, volf
            if label in d.snap.blocks and d.snap.info.get(label) is not None \
                    and d.snap.info[label].present[ptype]:
                a = d.raw_block(label, ptype)
                return (a if a.ndim == 1 else None), None
            warnings.warn(f"block {label} not available for type {ptype}")
            return None, None
        if label in self.derived:
            if label in _GAS_ONLY and ptype != 0:
                warnings.warn(f"derived block {label} is defined only for gas (type 0)")
                return None, None
            return self.derived_block(label)[s], None
        if label in d.snap.blocks and d.snap.info.get(label) is not None \
                and d.snap.info[label].present[ptype]:
            a = d.raw_block(label, ptype)
            if a.ndim == 1:
                warnings.warn(f"{label}: snapshot block used in code units (galpost extension)")
                return a, None
        warnings.warn(f"quantity not recognized: {label}")
        return None, None


# --------------------------------------------------------------------------
def compute_profiles(data: ParticleData, analyser, params, resolver: QuantityResolver):
    """Returns (profiles[gal][q] float64 arrays, counters[gal][q] int arrays, bins)."""
    p = params
    nbin = int(p["NumberOfBins"])
    gals = analyser.gals
    quants = p.quantities
    grid = bool(p["GridDensity"])
    surface = 1 if p["SurfaceDensity"] else 0
    ngal, nq = len(gals), len(quants)
    nb = nbin * nbin if grid else nbin

    # geometry of the bins (EG mutates rmin/rmax to log10 when LogFlag)
    bins = []
    for g in gals:
        if grid:
            dx = 2 * g.rmax / nbin
            bins.append(dict(dx=dx, vbin=np.full(nb, dx * dx)))
            continue
        rmin, rmax = g.rmin, g.rmax
        if g.log_bins:
            if rmin < SMALL:
                rmin = rmax * 0.01
                print(f"WARNING: found very small rmin with Log spacing, assuming: "
                      f"0.01 * rmax = {rmin:f}")
            rmin, rmax = math.log10(rmin), math.log10(rmax)
        dx = (rmax - rmin) / nbin
        f = 1.0 / dx
        a = -rmin * f
        edges = rmin + dx * np.arange(nbin + 1)
        if g.log_bins:
            edges = np.power(10.0, edges)
        xr, xo = edges[1:], edges[:-1]
        if surface:
            vbin = PI_EG_PROF * (xr * xr - xo * xo)
            vol = PI_EG_PROF * xr * xr
        else:
            vbin = 4.0 * PI_EG_PROF * (xr * xr * xr - xo * xo * xo) / 3.0
            vol = 4.0 * PI_EG_PROF * xr * xr * xr / 3.0
        bins.append(dict(rmin_b=rmin, rmax_b=rmax, dx=dx, f=f, a=a, vbin=vbin, vol=vol))

    flags = [q.flag & (1 if grid else 3) for q in quants]
    profiles = [[np.zeros(nb) for _ in range(nq)] for _ in range(ngal)]
    weights = [[np.zeros(nb) for _ in range(nq)] for _ in range(ngal)]
    counters = [[np.zeros(nb, dtype=np.int64) for _ in range(nq)] for _ in range(ngal)]
    volfact = []
    arrays = []
    for q in quants:
        vf = 1.0e6 if (p["SurfaceDensity"] or p["GridDensity"]) else 1.0e9
        qa, over = resolver.resolve(q.quantity, q.ptype, False)
        if over is not None:
            vf = over
        wa = None
        if q.weight != "NULL":
            wa, _ = resolver.resolve(q.weight, q.ptype, True)
        volfact.append(vf)
        arrays.append((qa, wa))

    belongs = analyser.belongs
    for t in range(6):
        qs = [iq for iq, q in enumerate(quants) if q.ptype == t]
        if not qs:
            continue
        s = data.type_slice(t)
        bt = belongs[s]
        for igal, g in enumerate(gals):
            loc = np.nonzero(bt == igal + 1)[0]
            if not len(loc):
                continue
            pos, _ = analyser.frame(g, loc + s.start)
            if not grid:
                dist = np.sqrt((pos[:, :3 - surface] ** 2).sum(axis=1))
                b = bins[igal]
                if g.log_bins:
                    dist = np.log10(np.maximum(dist, SMALL))
                ok = ~(dist < b["rmin_b"])
                ibin = (dist * b["f"] + b["a"] + 1.0e-6).astype(np.int64)
            else:
                rmax = g.rmax
                ok = ~((np.abs(pos[:, 0]) >= rmax) | (np.abs(pos[:, 1]) >= rmax)
                       | (np.abs(pos[:, 2]) >= rmax))
                dx = bins[igal]["dx"]
                ibin = ((pos[:, 0] + rmax) / dx).astype(np.int64) \
                    + ((pos[:, 1] + rmax) / dx).astype(np.int64) * nbin
            ok &= (ibin >= 0) & (ibin < nb)
            for iq in qs:
                qa, wa = arrays[iq]
                if qa is None or quants[iq].quantity == "VonS":
                    continue
                qv = qa[loc[ok]].astype(np.float64)
                wv = wa[loc[ok]].astype(np.float64) if wa is not None else np.ones(len(qv))
                ib = ibin[ok]
                profiles[igal][iq] += np.bincount(ib, weights=qv * wv, minlength=nb)
                weights[igal][iq] += np.bincount(ib, weights=wv, minlength=nb)
                counters[igal][iq] += np.bincount(ib, minlength=nb)

    for igal in range(ngal):
        b = bins[igal]
        for iq in range(nq):
            prof, w, cnt = profiles[igal][iq], weights[igal][iq], counters[igal][iq]
            if flags[iq] & 2:
                prof[:] = np.cumsum(prof)
                w[:] = np.cumsum(w)
            pos_w = w > 0
            prof[pos_w] /= w[pos_w]
            if flags[iq] & 1:
                if flags[iq] & 2:
                    prof *= np.cumsum(cnt).astype(np.float64) / b["vol"] / volfact[iq]
                else:
                    prof *= cnt.astype(np.float64) / b["vbin"] / volfact[iq]
    return profiles, counters, volfact
