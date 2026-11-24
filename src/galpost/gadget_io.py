# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Reader for GADGET / OpenGadget binary snapshots (format 2, and format 1 with a
block table), including SubFind catalogues written with WRITE_SUB_IN_SNAP_FORMAT.

How the layout of each block (data type, components, particle types) is found,
in order of priority:

1. a user **block table** (text file, see :func:`read_block_table`);
2. the **INFO block** written by OpenGadget (WRITE_INFO_BLOCK), if present;
3. **inference** from the block size and the header particle numbers, helped by
   a built-in table of common blocks. Since the true layout is always among the
   candidates, an inference is accepted only when it is unique; otherwise the
   block is marked ambiguous and reading it raises an error asking for a table.

Other features: multi-file snapshots (snap_NNN.0, .1, ...), companion files
(snap_NNN.muppi), endianness detection, blocks larger than 4 GB, detection of
truncated files, reading only the requested particle types.

Only numpy is required.
"""

from __future__ import annotations

import glob
import os
import struct
import warnings
from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Optional, Sequence

import numpy as np

__all__ = ["Header", "BlockInfo", "GadgetFile", "Snapshot", "read_header",
           "read_block_table", "write_block_table", "AmbiguousBlockError"]

NTYPES = 6
_TWO32 = 1 << 32

# INFO-block / block-table type strings -> numpy dtype code
_TYPE_CODES = {
    "FLOAT": "f4", "FLOATN": "f4", "DOUBLE": "f8", "DOUBLEN": "f8",
    "LONG": "i4", "LONGN": "i4", "LLONG": "i8", "LLONGN": "i8",
    "INT": "i4", "INTN": "i4",
    # numpy-style codes are accepted in block tables
    "f4": "f4", "f8": "f8", "i4": "i4", "i8": "i8", "u4": "u4", "u8": "u8",
}

# integer blocks that GADGET writes as unsigned although INFO says LONG/LLONG
_UNSIGNED_BLOCKS = {"ID", "PID", "MBID"}

# ---------------------------------------------------------------------------
# Built-in hints used only by the size inference (layer 3).
# name -> (kind, ndim or None, particle types that usually carry the block)
# kind: 'f' float, 'u' unsigned int, 'i' signed int; ndim None = any (e.g. metals)
# ---------------------------------------------------------------------------
_ALL = (0, 1, 2, 3, 4, 5)
_GAS = (0,)
_KNOWN = {
    "POS": ("f", 3, _ALL), "VEL": ("f", 3, _ALL), "ID": ("u", 1, _ALL),
    "MASS": ("f", 1, None),          # special: types with zero header mass
    "POT": ("f", 1, _ALL), "ACCE": ("f", 3, _ALL), "TSTP": ("f", 1, _ALL),
    "U": ("f", 1, _GAS), "RHO": ("f", 1, _GAS), "NE": ("f", 1, _GAS),
    "NH": ("f", 1, _GAS), "HSML": ("f", 1, _GAS), "SFR": ("f", 1, _GAS),
    "ENDT": ("f", 1, _GAS), "TEMP": ("f", 1, _GAS),
    "AGE": ("f", 1, (4, 5)), "iM": ("f", 1, (4,)), "ACRS": ("f", 1, (4,)),
    "Z": ("f", 1, (0, 4)), "Zs": ("f", None, (0, 4)),
    "BHMA": ("f", 1, (5,)), "BHMD": ("f", 1, (5,)), "BHPC": ("i", 1, (5,)),
    "ACRB": ("f", 1, (5,)),
    # MUPPI / OpenGadget gas blocks
    "MHOT": ("f", 1, _GAS), "MCLD": ("f", 1, _GAS), "MMOL": ("f", 1, _GAS),
    "EHOT": ("f", 1, _GAS), "MSF": ("f", 1, _GAS), "MFST": ("f", 1, _GAS),
    "NMF": ("f", 1, _GAS), "EOUT": ("f", 1, _GAS), "EREC": ("f", 1, _GAS),
    "CLCK": ("f", 1, _GAS), "Egy0": ("f", 1, _GAS), "TDYN": ("f", 1, _GAS),
    "TSMP": ("f", 1, _GAS), "EKRC": ("f", 3, _GAS), "GRAD": ("f", 3, _GAS),
    "VBLK": ("f", 3, _GAS), "ABVC": ("f", 1, _GAS), "ACVC": ("f", 1, _GAS),
    "VRMS": ("f", 1, _GAS), "DETI": ("f", 1, _GAS),
    "DSTL": ("f", None, _GAS), "DSTS": ("f", None, _GAS),
    # SubFind catalogue (WRITE_SUB_IN_SNAP_FORMAT): integer fields; the particle
    # types are left to the size inference (None = no preference)
    "GLEN": ("i", 1, None), "GOFF": ("i", 1, None), "NCON": ("i", 1, None),
    "NSUB": ("i", 1, None), "FSUB": ("i", 1, None), "SLEN": ("i", 1, None),
    "SOFF": ("i", 1, None), "SSUB": ("i", 1, None), "GRNR": ("i", 1, None),
    "MBID": ("u", 1, None), "PID": ("u", 1, None),
}


class AmbiguousBlockError(KeyError):
    """Layout of a block cannot be determined without a block table."""


# --------------------------------------------------------------------------
# Header
# --------------------------------------------------------------------------
@dataclass
class Header:
    npart: np.ndarray            # particles of each type in this file
    mass: np.ndarray             # mass table (0 -> masses in MASS block)
    time: float                  # scale factor (comoving) or time
    redshift: float
    flag_sfr: int
    flag_feedback: int
    npart_total: np.ndarray      # total particles of each type (all files, 64 bit)
    flag_cooling: int
    num_files: int
    boxsize: float
    omega0: float
    omega_lambda: float
    hubble_param: float
    flag_age: int
    flag_metals: int
    raw: bytes = field(repr=False, default=b"")

    @classmethod
    def from_bytes(cls, d: bytes, endian: str) -> "Header":
        e = endian
        npart = np.array(struct.unpack(e + "6I", d[0:24]), dtype=np.int64)
        mass = np.array(struct.unpack(e + "6d", d[24:72]))
        time, z = struct.unpack(e + "2d", d[72:88])
        fsfr, ffb = struct.unpack(e + "2i", d[88:96])
        ntot = np.array(struct.unpack(e + "6I", d[96:120]), dtype=np.int64)
        fcool, nfiles = struct.unpack(e + "2i", d[120:128])
        box, om, ol, h = struct.unpack(e + "4d", d[128:160])
        fage, fmet = struct.unpack(e + "2i", d[160:168])
        ntot_hw = np.array(struct.unpack(e + "6I", d[168:192]), dtype=np.int64)
        ntot = ntot + (ntot_hw << 32)
        return cls(npart, mass, time, z, fsfr, ffb, ntot, fcool, nfiles,
                   box, om, ol, h, fage, fmet, raw=d)


# --------------------------------------------------------------------------
# Block descriptions
# --------------------------------------------------------------------------
@dataclass
class BlockInfo:
    name: str
    dtype: np.dtype
    ndim: int
    present: np.ndarray          # bool[6], which particle types carry the block
    source: str = "info"         # 'table' | 'info' | 'inferred'

    @property
    def itemsize(self) -> int:
        return self.dtype.itemsize * self.ndim

    def with_endian(self, endian: str) -> "BlockInfo":
        return BlockInfo(self.name, self.dtype.newbyteorder(endian), self.ndim,
                         self.present.copy(), self.source)


def _make_dtype(name: str, code: str, endian: str) -> np.dtype:
    if name in _UNSIGNED_BLOCKS and code[0] == "i":
        code = "u" + code[1:]                 # GADGET IDs are unsigned (MyIDType)
    return np.dtype(endian + code)


def _parse_info(data: bytes, endian: str) -> Dict[str, BlockInfo]:
    out = {}
    for i in range(len(data) // 40):
        c = data[i * 40:(i + 1) * 40]
        lab = c[0:4].decode("ascii", "replace").strip()
        typ = c[4:12].decode("ascii", "replace").replace("\x00", "").strip()
        ndim = struct.unpack(endian + "i", c[12:16])[0]
        pres = np.array(struct.unpack(endian + "6i", c[16:40]), dtype=bool)
        if typ not in _TYPE_CODES:
            raise ValueError(f"unknown INFO type '{typ}' for block {lab}")
        out[lab] = BlockInfo(lab, _make_dtype(lab, _TYPE_CODES[typ], endian),
                             ndim, pres, "info")
    return out


def read_block_table(path: str) -> Dict[str, BlockInfo]:
    """Read a block table: one line per block, in file order (needed for format 1)::

        # NAME  TYPE    NDIM  p0 p1 p2 p3 p4 p5
        POS     FLOATN  3     1  1  1  1  1  1
        ID      LONG    1     1  1  1  1  1  1
        Zs      FLOATN  21    1  0  0  0  1  0

    TYPE is an INFO type (FLOAT, FLOATN, DOUBLE, LONG, LLONG, ...) or a numpy
    code (f4, f8, i4, u4, i8, u8). Lines starting with '#' or '%' are comments.
    The returned dict keeps the file order. Endianness is fixed when used.
    """
    out: Dict[str, BlockInfo] = {}
    with open(path) as f:
        for ln, line in enumerate(f, 1):
            line = line.split("#")[0].strip()
            if not line or line.startswith("%"):
                continue
            w = line.split()
            if len(w) != 9 or w[1] not in _TYPE_CODES:
                raise ValueError(f"{path}:{ln}: expected 'NAME TYPE NDIM p0..p5'")
            pres = np.array([int(x) for x in w[3:9]], dtype=bool)
            out[w[0]] = BlockInfo(w[0], _make_dtype(w[0], _TYPE_CODES[w[1]], "="),
                                  int(w[2]), pres, "table")
    return out


def write_block_table(info: Dict[str, BlockInfo], path: str):
    """Write block descriptions (e.g. ``Snapshot.info``) as a block table."""
    with open(path, "w") as f:
        f.write("# NAME  TYPE  NDIM  p0 p1 p2 p3 p4 p5   (source)\n")
        for b in info.values():
            code = b.dtype.str[1:]
            f.write(f"{b.name:<5} {code:<4} {b.ndim:4d}   "
                    + " ".join(str(int(p)) for p in b.present)
                    + f"   # {b.source}\n")


# --------------------------------------------------------------------------
# A single binary file
# --------------------------------------------------------------------------
class GadgetFile:
    """One physical file: locates blocks and their layout without reading data."""

    def __init__(self, path: str, info: Optional[Dict[str, BlockInfo]] = None,
                 header: Optional[Header] = None,
                 table: Optional[Dict[str, BlockInfo]] = None,
                 double_floats: bool = False):
        self.path = path
        self.size = os.path.getsize(path)
        self.header: Optional[Header] = header
        self.double_floats = double_floats
        self.format, self.endian = self._detect_format()
        self.table = {k: v.with_endian(self.endian) for k, v in (table or {}).items()}
        self.info: Dict[str, BlockInfo] = dict(info) if info else {}
        self.blocks: Dict[str, tuple] = {}     # name -> (data offset, nbytes)
        self.truncated: set = set()            # blocks cut short at end of file
        self.ambiguous: Dict[str, list] = {}   # name -> candidate layouts
        if self.format == 2:
            self._scan_format2()
        else:
            self._scan_format1()

    # ---- format / endianness --------------------------------------------
    def _detect_format(self):
        with open(self.path, "rb") as f:
            raw = f.read(4)
        for e in ("<", ">"):
            m = struct.unpack(e + "i", raw)[0]
            if m == 8:
                return 2, e
            if m == 256:
                return 1, e
        raise IOError(f"{self.path}: not a GADGET file (first record is neither "
                      "an 8-byte label nor a 256-byte header)")

    # ---- layout helpers ---------------------------------------------------
    def _describe(self, name: str) -> Optional[BlockInfo]:
        """Table has priority over INFO."""
        return self.table.get(name) or self.info.get(name)

    def _nelem_types(self, types) -> int:
        return int(sum(self.header.npart[t] for t in types))

    def _nelem(self, bi: BlockInfo) -> int:
        h = self.header
        return int(sum(h.npart[t] for t in range(NTYPES)
                       if bi.present[t] and not (bi.name == "MASS" and h.mass[t] != 0)))

    def _record_ok(self, f, data_off: int, nbytes: int, marker: int) -> bool:
        """True if a record of nbytes ends with the right marker and is followed
        by EOF or by another label record."""
        end = data_off + nbytes
        if end + 4 > self.size:
            return False
        f.seek(end)
        if struct.unpack(self.endian + "I", f.read(4))[0] != marker:
            return False
        if end + 4 == self.size:
            return True
        nxt = f.read(4)
        return len(nxt) == 4 and struct.unpack(self.endian + "i", nxt)[0] in (8, 256)

    def _infer(self, name: str, nbytes: int) -> Optional[BlockInfo]:
        h = self.header
        nonzero = [t for t in range(NTYPES) if h.npart[t] > 0]
        hint = _KNOWN.get(name)
        if name == "MASS":
            subsets = [tuple(t for t in nonzero if h.mass[t] == 0)]
        else:
            subsets = [c for r in range(1, len(nonzero) + 1)
                       for c in combinations(nonzero, r)]
        cands = []
        for s in subsets:
            n = self._nelem_types(s)
            if n == 0 or nbytes % n:
                continue
            isz = nbytes // n
            if isz % 4 == 0 and 4 <= isz <= 4 * 256:
                cands.append((s, isz))
        if hint is not None and hint[2] is not None:
            pref = tuple(t for t in hint[2] if t in nonzero)
            match = [c for c in cands if c[0] == pref]
            if match:
                cands = match
        if not cands:
            self.ambiguous[name] = []
            return None
        if len(cands) > 1:
            self.ambiguous[name] = cands
            return None
        types, isz = cands[0]
        kind, ndim = (hint[0], hint[1]) if hint else ("f", None)
        if ndim is not None:
            if isz % ndim:
                self.ambiguous[name] = cands
                return None
            esz = isz // ndim
        else:
            esz = 8 if (kind == "f" and self.double_floats and isz % 8 == 0) else 4
            if kind in "ui" and isz == 8:
                esz = 8
            ndim = isz // esz
        code = f"{kind}{esz}"
        pres = np.zeros(NTYPES, dtype=bool)
        pres[list(types)] = True
        return BlockInfo(name, np.dtype(self.endian + code), ndim, pres, "inferred")

    def _register(self, name, data_off, nbytes):
        self.blocks[name] = (data_off, nbytes)

    # ---- format 2 -----------------------------------------------------------
    def _scan_format2(self):
        e = self.endian
        with open(self.path, "rb") as f:
            pos = 0
            while pos < self.size:
                f.seek(pos)
                lab_rec = f.read(16)
                if len(lab_rec) < 16:
                    break
                if struct.unpack(e + "i", lab_rec[0:4])[0] != 8:
                    raise IOError(f"{self.path}: corrupted label record at byte {pos}")
                name = lab_rec[4:8].decode("ascii", "replace").strip()
                marker, = struct.unpack(e + "I", f.read(4))
                data_off = pos + 20

                if name == "HEAD" and self.header is None:
                    self.header = Header.from_bytes(f.read(marker), e)
                if name == "INFO":
                    f.seek(data_off)
                    self.info.update(_parse_info(f.read(marker), e))

                bi = self._describe(name) if name not in ("HEAD", "INFO") else None
                if bi is not None:
                    nbytes = self._nelem(bi) * bi.itemsize
                    if nbytes % _TWO32 != marker:
                        raise IOError(f"{self.path}: block {name} size {marker} does "
                                      f"not match its {bi.source} description ({nbytes})")
                else:
                    # size from the markers; for >4 GB blocks try marker + k*2^32
                    nbytes = None
                    for k in range(16):
                        trial = marker + k * _TWO32
                        if data_off + trial + 4 > self.size:
                            break
                        if self._record_ok(f, data_off, trial, marker):
                            nbytes = trial
                            break
                    if nbytes is None:
                        nbytes = marker            # -> reported as truncated below
                    if name not in ("HEAD", "INFO") and self.header is not None:
                        inf = self._infer(name, nbytes)
                        if inf is not None:
                            self.info[name] = inf

                if data_off + nbytes + 4 > self.size:
                    missing = data_off + nbytes + 4 - self.size
                    warnings.warn(f"{self.path} is truncated: block {name} misses "
                                  f"{missing} bytes and will not be readable")
                    self.truncated.add(name)
                    break
                self._register(name, data_off, nbytes)
                pos = data_off + nbytes + 4

        for name, c in self.ambiguous.items():
            warnings.warn(f"{self.path}: layout of block {name} is ambiguous "
                          f"({len(c)} candidates); provide a block table to read it")

    # ---- format 1 -----------------------------------------------------------
    def _scan_format1(self):
        if not self.table:
            raise IOError(f"{self.path}: format-1 file (no block labels): a block "
                          "table listing the blocks in file order is required")
        e = self.endian
        names = list(self.table)
        with open(self.path, "rb") as f:
            m, = struct.unpack(e + "I", f.read(4))
            if self.header is None:
                self.header = Header.from_bytes(f.read(m), e)
            pos = 4 + m + 4
            i = 0
            while pos < self.size:
                if i >= len(names):
                    warnings.warn(f"{self.path}: more blocks in file than in the table")
                    break
                bi = self.table[names[i]]
                nbytes = self._nelem(bi) * bi.itemsize
                f.seek(pos)
                marker, = struct.unpack(e + "I", f.read(4))
                if nbytes % _TWO32 != marker:
                    raise IOError(f"{self.path}: block #{i} ({bi.name}) has {marker} "
                                  f"bytes, table implies {nbytes}")
                if pos + 4 + nbytes + 4 > self.size:
                    warnings.warn(f"{self.path} is truncated in block {bi.name}")
                    self.truncated.add(bi.name)
                    break
                self._register(bi.name, pos + 4, nbytes)
                pos += 4 + nbytes + 4
                i += 1

    # ---- data access ----------------------------------------------------------
    def layout(self, name: str) -> BlockInfo:
        bi = self._describe(name)
        if bi is None:
            if name in self.ambiguous:
                cands = ", ".join(f"types {s} x {isz} B" for s, isz in self.ambiguous[name])
                raise AmbiguousBlockError(
                    f"block {name} in {self.path}: layout not described (no INFO entry) "
                    f"and not uniquely inferable [{cands or 'no candidate'}]; "
                    "provide a block table (galpost.gadget_io.read_block_table)")
            raise KeyError(f"block {name} has no layout description")
        return bi

    def type_offsets(self, name: str) -> Dict[int, tuple]:
        """{ptype: (first element index inside block, n elements)}"""
        bi = self.layout(name)
        h = self.header
        out, start = {}, 0
        for t in range(NTYPES):
            if not bi.present[t] or (name == "MASS" and h.mass[t] != 0):
                continue
            out[t] = (start, int(h.npart[t]))
            start += int(h.npart[t])
        return out

    def read(self, name: str, ptype: int, index: Optional[np.ndarray] = None
             ) -> Optional[np.ndarray]:
        """Read block `name` for particle type `ptype` (None if absent).

        `index` (sorted int array, positions within this file's particles of
        that type) reads only a subset, through a memory map.
        """
        if name not in self.blocks:
            return None
        bi = self.layout(name)
        offs = self.type_offsets(name)
        if ptype not in offs:
            return None
        first, n = offs[ptype]
        data_off, _ = self.blocks[name]
        shape = (n, bi.ndim) if bi.ndim > 1 else (n,)
        if index is None:
            arr = np.fromfile(self.path, dtype=bi.dtype, count=n * bi.ndim,
                              offset=data_off + first * bi.itemsize).reshape(shape)
        else:
            if n == 0 or len(index) == 0:
                return np.zeros((0,) + shape[1:], dtype=bi.dtype.newbyteorder("="))
            arr = self._read_subset(data_off + first * bi.itemsize, bi, n, index)
        return arr.astype(bi.dtype.newbyteorder("="), copy=False)

    def _read_subset(self, offset: int, bi: BlockInfo, n: int, index: np.ndarray) -> np.ndarray:
        """Sorted subset read through a memory map: data come straight from the
        page cache and only the output array is allocated."""
        shape = (n, bi.ndim) if bi.ndim > 1 else (n,)
        mm = np.memmap(self.path, dtype=bi.dtype, mode="r", offset=offset, shape=shape)
        out = np.array(mm[index])
        del mm
        return out

    def memmap(self, name: str, ptype: int):
        """Read-only memory map of block `name` for type `ptype` (or None)."""
        if name not in self.blocks:
            return None
        bi = self.layout(name)
        offs = self.type_offsets(name)
        if ptype not in offs or offs[ptype][1] == 0:
            return None
        first, n = offs[ptype]
        data_off, _ = self.blocks[name]
        shape = (n, bi.ndim) if bi.ndim > 1 else (n,)
        return np.memmap(self.path, dtype=bi.dtype, mode="r",
                         offset=data_off + first * bi.itemsize, shape=shape)

    def read_raw(self, name: str, ptype: int) -> Optional[bytes]:
        """Raw bytes of block `name` for type `ptype` (as stored on disk)."""
        if name not in self.blocks:
            return None
        bi = self.layout(name)
        offs = self.type_offsets(name)
        if ptype not in offs:
            return None
        first, n = offs[ptype]
        data_off, _ = self.blocks[name]
        with open(self.path, "rb") as f:
            f.seek(data_off + first * bi.itemsize)
            return f.read(n * bi.itemsize)


# --------------------------------------------------------------------------
# Snapshot = all files (multi-file + companions)
# --------------------------------------------------------------------------
def _expand_files(base: str) -> List[str]:
    if os.path.isfile(base):
        return [base]
    files = [f for f in glob.glob(base + ".*") if f.rsplit(".", 1)[1].isdigit()]
    files.sort(key=lambda s: int(s.rsplit(".", 1)[1]))
    if not files:
        raise FileNotFoundError(f"no file '{base}' nor '{base}.N'")
    return files


class Snapshot:
    """
    Lazy access to a (possibly multi-file) GADGET snapshot or SubFind catalogue.

    Parameters
    ----------
    base : str
        e.g. 'snap_134' (single file) or 'groups_134/sub_134' (-> .0, .1, ...).
    companions : sequence of str
        suffixes of companion files holding additional blocks with the same
        particle ordering. Default ('.muppi',); missing companions are ignored.
    block_table : str or dict, optional
        block table file (see :func:`read_block_table`) or already-parsed dict.
        It overrides the INFO block; required for format-1 files.
    double_floats : bool
        for *inferred* float blocks of unknown dimension, assume double precision.

    Example
    -------
    >>> s = Snapshot('snap_134')
    >>> pos = s.read('POS', 4)       # stellar positions
    >>> s.blocks_for(0)               # blocks available for gas
    """

    def __init__(self, base: str, companions: Sequence[str] = (".muppi",),
                 block_table=None, double_floats: bool = False):
        self.base = base
        table = read_block_table(block_table) if isinstance(block_table, str) else block_table
        self.files: List[List[GadgetFile]] = []   # per file number: [main, companions...]
        first: Optional[GadgetFile] = None
        for path in _expand_files(base):
            main = GadgetFile(path, table=table, double_floats=double_floats,
                              info=first.info if first else None)
            first = first or main
            group = [main]
            stem, num = (path, "") if path == base else path.rsplit(".", 1)
            for suf in companions:
                cpath = stem + suf + ("." + num if num else "")
                if not os.path.isfile(cpath) and os.path.islink(path):
                    # snapshot given through a symlink: look next to the real file
                    real = os.path.realpath(path)
                    rstem = real if not num else real.rsplit(".", 1)[0]
                    cpath = rstem + suf + ("." + num if num else "")
                if os.path.isfile(cpath):
                    group.append(GadgetFile(cpath, info=main.info, header=main.header,
                                            table=table, double_floats=double_floats))
            self.files.append(group)
        self.header: Header = self.files[0][0].header

    # ---- metadata ------------------------------------------------------
    @property
    def info(self) -> Dict[str, BlockInfo]:
        """Resolved layout of every block (table > INFO > inferred)."""
        out: Dict[str, BlockInfo] = {}
        for gf in self.files[0]:
            for b in gf.blocks:
                if b in ("HEAD", "INFO") or b in out:
                    continue
                d = gf._describe(b)
                if d is not None:
                    out[b] = d
        return out

    @property
    def blocks(self) -> List[str]:
        names = []
        for g in self.files[0]:
            names += [b for b in g.blocks if b not in ("HEAD", "INFO") and b not in names]
        return names

    def blocks_for(self, ptype: int) -> List[str]:
        info = self.info
        return [b for b in self.blocks if b in info and info[b].present[ptype]]

    @property
    def ambiguous(self) -> List[str]:
        return sorted({b for g in self.files for gf in g for b in gf.ambiguous})

    @property
    def npart(self) -> np.ndarray:
        return self.header.npart_total

    def __repr__(self):
        h = self.header
        return (f"<Snapshot {self.base}: z={h.redshift:.4f} a={h.time:.4f} "
                f"npart={list(h.npart_total)} files={len(self.files)}>")

    # ---- data ----------------------------------------------------------
    def _locate(self, name: str):
        for group in self.files:
            gf = next((g for g in group if name in g.blocks), None)
            yield gf

    def read(self, name: str, ptype: int, index: Optional[np.ndarray] = None) -> np.ndarray:
        """Concatenate block `name` of type `ptype` over all files.

        `index`: optional sorted array of particle positions (over all files,
        within type `ptype`) to read only a subset.
        MASS for types with a fixed mass in the header is synthesised.
        """
        if name == "MASS" and self.header.mass[ptype] != 0:
            n = int(self.npart[ptype]) if index is None else len(index)
            return np.full(n, self.header.mass[ptype], dtype="f4")
        parts, found = [], False
        start = 0
        for group, gf in zip(self.files, self._locate(name)):
            nfile = int(group[0].header.npart[ptype])
            sub = None
            if index is not None:
                lo, hi = np.searchsorted(index, [start, start + nfile])
                sub = index[lo:hi] - start
            start += nfile
            if gf is None:
                continue
            found = True
            if sub is not None and len(sub) == 0:
                continue
            a = gf.read(name, ptype, sub)
            if a is not None:
                parts.append(a)
        if not found:
            if any(name in gf.truncated for g in self.files for gf in g):
                raise IOError(f"block '{name}' is incomplete (truncated file)")
            raise KeyError(f"block '{name}' not found in {self.base} "
                           f"(available: {', '.join(self.blocks)})")
        if not parts:
            bi = self.info.get(name)
            if index is not None and bi is not None and bi.present[ptype]:
                shape = (0, bi.ndim) if bi.ndim > 1 else (0,)
                return np.zeros(shape, dtype=bi.dtype.newbyteorder("="))
            raise KeyError(f"block '{name}' not present for particle type {ptype}")
        return parts[0] if len(parts) == 1 else np.concatenate(parts)


def read_header(path: str) -> Header:
    return GadgetFile(_expand_files(path)[0]).header
