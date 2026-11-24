# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Self-contained tests: write a tiny synthetic format-2 snapshot and read it back."""
import struct
import numpy as np
import pytest

from galpost.gadget_io import Snapshot, AmbiguousBlockError


def _rec(f, payload):
    f.write(struct.pack("<I", len(payload)) + payload + struct.pack("<I", len(payload)))


def _block(f, name, payload):
    _rec(f, name.ljust(4).encode() + struct.pack("<i", len(payload) + 8))
    _rec(f, payload)


def _info(entries):
    out = b""
    for name, typ, ndim, pres in entries:
        out += name.ljust(4).encode() + typ.ljust(8).encode() + struct.pack("<i6i", ndim, *pres)
    return out


@pytest.fixture
def snap(tmp_path):
    npart = [3, 2, 0, 0, 4, 1]
    mass = [0, 0.5, 0, 0, 0, 0]
    head = struct.pack("<6I6d2d2i6I2i4d2i6I", *npart, *mass, 1.0, 0.0, 1, 1, *npart,
                       1, 1, 100.0, 0.25, 0.75, 0.73, 1, 1, 0, 0, 0, 0, 0, 0)
    head = head.ljust(256, b"\0")
    n = sum(npart)
    rng = np.random.default_rng(1)
    data = {
        "POS": rng.random((n, 3), dtype=np.float32),
        "ID": np.arange(n, dtype=np.uint32) + np.uint32(2**31),
        "MASS": rng.random(3 + 4 + 1, dtype=np.float32),       # no DM (fixed mass)
        "Zs": rng.random((3 + 4, 2), dtype=np.float32),        # gas + stars
    }
    info = _info([("POS", "FLOATN", 3, [1] * 6), ("ID", "LONG", 1, [1] * 6),
                  ("MASS", "FLOAT", 1, [1, 0, 1, 1, 1, 1]),
                  ("Zs", "FLOATN", 2, [1, 0, 0, 0, 1, 0]),
                  ("EOUT", "FLOAT", 1, [1, 0, 0, 0, 0, 0])])
    p = tmp_path / "snap_001"
    with open(p, "wb") as f:
        _block(f, "HEAD", head)
        _block(f, "INFO", info)
        for k, v in data.items():
            _block(f, k, v.tobytes())
    eout = np.array([7, 8, 9], dtype=np.float32)
    with open(str(p) + ".muppi", "wb") as f:
        _block(f, "EOUT", eout.tobytes())
    return str(p), data, eout


def test_read(snap):
    path, data, eout = snap
    s = Snapshot(path)
    assert list(s.npart) == [3, 2, 0, 0, 4, 1]
    np.testing.assert_array_equal(s.read("POS", 4), data["POS"][5:9])
    np.testing.assert_array_equal(s.read("ID", 5), data["ID"][9:])
    assert s.read("ID", 5).dtype.kind == "u"
    np.testing.assert_array_equal(s.read("MASS", 1), [0.5, 0.5])          # from header
    np.testing.assert_array_equal(s.read("MASS", 4), data["MASS"][3:7])
    np.testing.assert_array_equal(s.read("Zs", 4), data["Zs"][3:])
    np.testing.assert_array_equal(s.read("EOUT", 0), eout)                # companion file
    assert "Zs" in s.blocks_for(4) and "Zs" not in s.blocks_for(1)
    with pytest.raises(KeyError):
        s.read("NOPE", 0)


def test_truncated_companion(snap):
    path, data, eout = snap
    cpath = path + ".muppi"
    raw = open(cpath, "rb").read()
    open(cpath, "wb").write(raw[:-6])          # cut the last block short
    with pytest.warns(UserWarning):
        s = Snapshot(path)
    np.testing.assert_array_equal(s.read("POS", 4), data["POS"][5:9])   # main file fine
    with pytest.raises(IOError):
        s.read("EOUT", 0)


# ---------------------------------------------------------------------------
# files without INFO block: inference, ambiguity, block table, format 1
# ---------------------------------------------------------------------------
def _noinfo_data():
    npart = [3, 2, 0, 0, 4, 1]
    mass = [0, 0.5, 0, 0, 0, 0]
    head = struct.pack("<6I6d2d2i6I2i4d2i6I", *npart, *mass, 1.0, 0.0, 1, 1, *npart,
                       1, 1, 100.0, 0.25, 0.75, 0.73, 1, 1, 0, 0, 0, 0, 0, 0).ljust(256, b"\0")
    rng = np.random.default_rng(2)
    data = {
        "POS": rng.random((10, 3), dtype=np.float32),
        "ID": np.arange(10, dtype=np.uint32),
        "MASS": rng.random(8, dtype=np.float32),
        "Zs": rng.random((7, 2), dtype=np.float32),
        "XYZW": rng.random(3, dtype=np.float32),      # unknown gas block: ambiguous
    }
    return head, data


@pytest.fixture
def noinfo(tmp_path):
    head, data = _noinfo_data()
    p = tmp_path / "snap_002"
    with open(p, "wb") as f:
        _block(f, "HEAD", head)
        for k, v in data.items():
            _block(f, k, v.tobytes())
    return str(p), data


def test_inference_without_info(noinfo):
    path, data = noinfo
    with pytest.warns(UserWarning):                    # XYZW is ambiguous
        s = Snapshot(path)
    np.testing.assert_array_equal(s.read("POS", 5), data["POS"][9:])
    np.testing.assert_array_equal(s.read("ID", 4), data["ID"][5:9])
    np.testing.assert_array_equal(s.read("MASS", 5), data["MASS"][7:])
    np.testing.assert_array_equal(s.read("Zs", 4), data["Zs"][3:])
    assert s.info["Zs"].source == "inferred" and s.info["Zs"].ndim == 2
    assert s.ambiguous == ["XYZW"]
    with pytest.raises(AmbiguousBlockError):
        s.read("XYZW", 0)


def test_block_table_resolves(noinfo, tmp_path):
    path, data = noinfo
    t = tmp_path / "blocks.txt"
    t.write_text("# NAME TYPE NDIM p0..p5\nXYZW FLOAT 1 1 0 0 0 0 0\n")
    s = Snapshot(path, block_table=str(t))
    np.testing.assert_array_equal(s.read("XYZW", 0), data["XYZW"])
    assert s.info["XYZW"].source == "table"


def test_format1_with_table(tmp_path):
    head, data = _noinfo_data()
    p = tmp_path / "snap_003"
    with open(p, "wb") as f:
        _rec(f, head)
        for v in data.values():
            _rec(f, v.tobytes())
    t = tmp_path / "blocks.txt"
    t.write_text("POS FLOATN 3 1 1 1 1 1 1\nID LONG 1 1 1 1 1 1 1\n"
                 "MASS FLOAT 1 1 0 1 1 1 1\nZs FLOATN 2 1 0 0 0 1 0\n"
                 "XYZW FLOAT 1 1 0 0 0 0 0\n")
    with pytest.raises(IOError):
        Snapshot(str(p))                                # no table -> refused
    s = Snapshot(str(p), block_table=str(t))
    np.testing.assert_array_equal(s.read("Zs", 0), data["Zs"][:3])
    np.testing.assert_array_equal(s.read("XYZW", 0), data["XYZW"])
