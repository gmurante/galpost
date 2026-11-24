# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Halo catalogues: FoF group tables and SubFind catalogues.

* :class:`FofCatalog` reads ``groups_NNN/group_tab_NNN.*`` written either by the
  standalone P-FoF code ("old" format) or by Gadget3/OpenGadget on the fly or
  in post-processing ("G3" format). Format and endianness are detected as in
  ExtractGalaxies (the number of files is stored in the header).
* :class:`SubfindCatalog` reads ``groups_NNN/sub_NNN.*`` written by OpenGadget
  SubFind with WRITE_SUB_IN_SNAP_FORMAT (GADGET format 2), through the generic
  snapshot reader, so every block present in the files is available.

All quantities are returned in code units; conversions are done elsewhere.
"""

from __future__ import annotations

import glob
import os
import struct
from typing import Dict, Optional, Sequence

import numpy as np

from .gadget_io import Snapshot

__all__ = ["FofCatalog", "SubfindCatalog", "find_group_dir"]

# per-group optional fields appended at the end of a G3 group_tab file,
# in the order they are written by fof.c (SFR, then BLACK_HOLES)
_G3_EXTRAS = [("sfr", "f4"), ("bh_mass", "f4"), ("bh_mdot", "f4")]
# P-FoF writes one extra float per group (GroupSfr, identically zero)
_OLD_EXTRAS = [("sfr", "f4")]


def find_group_dir(base_dir: str, num: int) -> str:
    """Return ``base_dir/groups_NNN`` (also accepts 4-digit numbering)."""
    for fmt in ("groups_%03d", "groups_%04d"):
        d = os.path.join(base_dir, fmt % num)
        if os.path.isdir(d):
            return d
    raise FileNotFoundError(f"no groups_{num:03d} directory in {base_dir}")


def _tag(group_dir: str, num: int) -> str:
    """Snapshot number as written in the file names (from 'groups_NNN', or from
    `num` if the directory has another name)."""
    t = os.path.basename(group_dir.rstrip("/")).split("_", 1)[-1]
    if t.isdigit():
        return t
    for fmt in ("%03d", "%04d"):
        if glob.glob(os.path.join(group_dir, ("*_" + fmt + ".*") % num)):
            return fmt % num
    return "%03d" % num


def _numbered_files(pattern: str):
    files = [f for f in glob.glob(pattern + ".*") if f.rsplit(".", 1)[1].isdigit()]
    files.sort(key=lambda s: int(s.rsplit(".", 1)[1]))
    return files


# --------------------------------------------------------------------------
# FoF
# --------------------------------------------------------------------------
class FofCatalog:
    """
    FoF group table.

    Attributes (arrays over all groups, in file order)
    ----------
    length, offset : int        particles per group, offset in the ID list
    len_type       : int[N,6]   particles per type
    mass_type      : float[N,6] mass per type
    mass           : float      total mass (G3: as written; old: sum of mass_type)
    cm             : float[N,3] centre of mass
    vel            : float[N,3] or None (old format has no velocity)
    extra          : dict       optional per-group fields (sfr, bh_mass, bh_mdot)
    fmt            : 'old' or 'g3';  endian: '<' or '>'
    """

    def __init__(self, group_dir: str, num: int):
        self.group_dir = group_dir
        self.num = num
        tag = _tag(group_dir, num)
        self.tab_files = _numbered_files(os.path.join(group_dir, f"group_tab_{tag}"))
        self.ids_files = _numbered_files(os.path.join(group_dir, f"group_ids_{tag}"))
        if not self.tab_files:
            raise FileNotFoundError(f"no group_tab files in {group_dir}")
        self.fmt, self.endian = self._guess_format()
        self._read()

    @classmethod
    def from_dir(cls, base_dir: str, num: int) -> "FofCatalog":
        return cls(find_group_dir(base_dir, num), num)

    # ---- format detection (same logic as ExtractGalaxies) ----------------
    def _guess_format(self):
        n = len(self.tab_files)
        raw = open(self.tab_files[0], "rb").read(24)
        for e in ("<", ">"):
            if len(raw) >= 16 and struct.unpack(e + "4i", raw[:16])[3] == n:
                return "old", e
        for e in ("<", ">"):
            if len(raw) >= 24 and struct.unpack(e + "6i", raw[:24])[5] == n:
                return "g3", e
        raise IOError(f"{self.group_dir}: unrecognised FoF table format")

    def _header(self, raw: bytes):
        e = self.endian
        if self.fmt == "old":
            ng, nids, totng, ntask = struct.unpack(e + "4i", raw[:16])
            return ng, nids, totng, ntask, 16
        ng, totng, nids = struct.unpack(e + "3i", raw[:12])
        ntask, = struct.unpack(e + "i", raw[20:24])
        return ng, nids, totng, ntask, 24

    # ---- reading ---------------------------------------------------------
    def _read(self):
        e = self.endian
        cols: Dict[str, list] = {}
        known = (_OLD_EXTRAS if self.fmt == "old" else _G3_EXTRAS)
        self.ngroups_per_file, self.nids_per_file = [], []
        for path in self.tab_files:
            raw = open(path, "rb").read()
            ng, nids, totng, ntask, off = self._header(raw)
            self.ngroups_per_file.append(ng)
            self.nids_per_file.append(nids)
            self.total_groups = totng
            if ng == 0:
                continue
            if self.fmt == "old":
                fields = [("length", "i4", 1), ("offset", "i4", 1), ("len_type", "i4", 6),
                          ("mass_type", "f8", 6), ("cm", "f4", 3)]
            else:
                fields = [("length", "i4", 1), ("offset", "i4", 1), ("mass", "f4", 1),
                          ("cm", "f4", 3), ("vel", "f4", 3), ("len_type", "i4", 6),
                          ("mass_type", "f4", 6)]
            for name, code, nd in fields:
                a = np.frombuffer(raw, dtype=e + code, count=ng * nd, offset=off)
                cols.setdefault(name, []).append(a.reshape(ng, nd) if nd > 1 else a)
                off += ng * nd * np.dtype(code).itemsize
            # trailing optional per-group fields
            rest = len(raw) - off
            if rest % (4 * ng):
                raise IOError(f"{path}: {rest} trailing bytes not a multiple of 4*Ngroups")
            nextra = rest // (4 * ng)
            for k in range(nextra):
                name = known[k][0] if k < len(known) else f"extra{k}"
                a = np.frombuffer(raw, dtype=e + "f4", count=ng, offset=off)
                cols.setdefault("x:" + name, []).append(a)
                off += 4 * ng

        def cat(name, nd=1, dt="f4"):
            if name in cols:
                return np.concatenate(cols[name]).astype(np.dtype(cols[name][0].dtype).newbyteorder("="))
            return np.zeros((0, nd) if nd > 1 else 0, dtype=dt)

        self.length = cat("length", 1, "i4")
        self.offset = cat("offset", 1, "i4")
        self.len_type = cat("len_type", 6, "i4")
        self.mass_type = cat("mass_type", 6).astype(np.float64)
        self.cm = cat("cm", 3)
        self.vel = cat("vel", 3) if self.fmt == "g3" else None
        self.mass = cat("mass") if self.fmt == "g3" else self.mass_type.sum(axis=1)
        self.extra = {k[2:]: cat(k) for k in cols if k.startswith("x:")}

    def __len__(self):
        return len(self.length)

    def __repr__(self):
        return (f"<FofCatalog {self.group_dir}: {len(self)} groups, format={self.fmt}, "
                f"files={len(self.tab_files)}>")

    def contamination(self, lowres_types: Sequence[int] = (2, 3)) -> np.ndarray:
        """Mass fraction in low-resolution particle types (zoom simulations)."""
        tot = self.mass_type.sum(axis=1)
        low = self.mass_type[:, list(lowres_types)].sum(axis=1)
        return np.divide(low, tot, out=np.zeros_like(tot), where=tot > 0)

    def group_ids(self, igroup: int) -> np.ndarray:
        """Particle IDs of a group (needs group_ids files)."""
        if not self.ids_files:
            raise FileNotFoundError(f"no group_ids files in {self.group_dir}")
        if not hasattr(self, "_ids"):
            e, parts = self.endian, []
            for path in self.ids_files:
                raw = open(path, "rb").read()
                if self.fmt == "old":
                    nids = struct.unpack(e + "i", raw[4:8])[0]
                    hdr = 16
                else:
                    nids = struct.unpack(e + "i", raw[8:12])[0]
                    hdr = 28
                if nids:
                    isz = (len(raw) - hdr) // nids
                    parts.append(np.frombuffer(raw, dtype=e + ("u8" if isz == 8 else "u4"),
                                               count=nids, offset=hdr))
            self._ids = np.concatenate(parts) if parts else np.zeros(0, "u4")
        o, n = int(self.offset[igroup]), int(self.length[igroup])
        return self._ids[o:o + n]


# --------------------------------------------------------------------------
# SubFind
# --------------------------------------------------------------------------
class SubfindCatalog:
    """
    SubFind catalogue in snapshot format (``sub_NNN.*``).

    Groups are "particle type" 0, subhaloes type 1, the ID list type 2.
    Any block can be read with :meth:`group` / :meth:`subhalo`; frequently used
    ones are cached.  The FoF table of the same directory, if present, is
    attached as ``self.fof`` (it has the same group order) and provides the
    mass per particle type, used for the contamination of zoom runs.

    Notable blocks: groups GLEN GOFF MTOT GPOS M200/R200 (mean), MCRI/RCRI
    (critical), MVIR/RVIR, M500/R500, NSUB, FSUB (first = central subhalo);
    subhaloes SLEN SOFF MSUB SPOS SVEL SCM SPIN VMAX RMAX RHMS MBID GRNR SMST
    (mass per type).
    """

    def __init__(self, group_dir: str, num: int, block_table=None):
        self.group_dir = group_dir
        self.num = num
        tag = _tag(group_dir, num)
        base = os.path.join(group_dir, f"sub_{tag}")
        self.snap = Snapshot(base, companions=(), block_table=block_table)
        self.header = self.snap.header
        self.ngroups = int(self.snap.npart[0])
        self.nsubhalos = int(self.snap.npart[1])
        self._cache: Dict[tuple, np.ndarray] = {}
        try:
            self.fof: Optional[FofCatalog] = FofCatalog(group_dir, num)
            if len(self.fof) != self.ngroups:
                self.fof = None
        except (FileNotFoundError, IOError):
            self.fof = None

    @classmethod
    def from_dir(cls, base_dir: str, num: int, **kw) -> "SubfindCatalog":
        return cls(find_group_dir(base_dir, num), num, **kw)

    def __repr__(self):
        return (f"<SubfindCatalog {self.group_dir}: {self.ngroups} groups, "
                f"{self.nsubhalos} subhaloes>")

    @property
    def blocks(self):
        return {t: self.snap.blocks_for(t) for t in (0, 1, 2)}

    def _get(self, name: str, ptype: int) -> np.ndarray:
        key = (name, ptype)
        if key not in self._cache:
            self._cache[key] = self.snap.read(name, ptype)
        return self._cache[key]

    def group(self, name: str) -> np.ndarray:
        return self._get(name, 0)

    def subhalo(self, name: str) -> np.ndarray:
        return self._get(name, 1)

    def central_subhalo(self, igroup: int) -> int:
        return int(self.group("FSUB")[igroup])

    def subhalos_of_group(self, igroup: int) -> np.ndarray:
        f, n = int(self.group("FSUB")[igroup]), int(self.group("NSUB")[igroup])
        return np.arange(f, f + n)

    def _ids(self) -> np.ndarray:
        return self._get("PID", 2)

    def subhalo_ids(self, isub: int) -> np.ndarray:
        o, n = int(self.subhalo("SOFF")[isub]), int(self.subhalo("SLEN")[isub])
        return self._ids()[o:o + n]

    def group_ids(self, igroup: int) -> np.ndarray:
        o, n = int(self.group("GOFF")[igroup]), int(self.group("GLEN")[igroup])
        return self._ids()[o:o + n]

    def contamination(self, lowres_types: Sequence[int] = (2, 3)) -> np.ndarray:
        """Low-resolution mass fraction of each group (needs the FoF table)."""
        if self.fof is None:
            raise RuntimeError("no matching group_tab files: contamination unknown")
        return self.fof.contamination(lowres_types)
