# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Synthetic tests for FoF tables (P-FoF 'old' and G3 formats) and paramfiles."""
import struct
import numpy as np
import pytest

from galpost.catalogs import FofCatalog
from galpost.params import read_params


def _write_old(d, ng_per_file):
    nfiles = len(ng_per_file)
    rng = np.random.default_rng(3)
    allg = []
    for k, ng in enumerate(ng_per_file):
        L = rng.integers(32, 1000, ng).astype("<i4"); O = np.zeros(ng, "<i4")
        LT = rng.integers(0, 100, (ng, 6)).astype("<i4"); MT = rng.random((ng, 6)).astype("<f8")
        CM = rng.random((ng, 3)).astype("<f4"); SFR = np.zeros(ng, "<f4")
        with open(d / f"group_tab_005.{k}", "wb") as f:
            f.write(struct.pack("<4i", ng, int(L.sum()), sum(ng_per_file), nfiles))
            for a in (L, O, LT, MT, CM, SFR):
                f.write(a.tobytes())
        allg.append((L, MT, CM))
    return allg


def _write_g3(d, ng_per_file):
    nfiles = len(ng_per_file)
    rng = np.random.default_rng(4)
    allg = []
    for k, ng in enumerate(ng_per_file):
        L = rng.integers(32, 1000, ng).astype("<i4"); O = np.zeros(ng, "<i4")
        M = rng.random(ng).astype("<f4"); CM = rng.random((ng, 3)).astype("<f4")
        V = rng.random((ng, 3)).astype("<f4"); LT = rng.integers(0, 9, (ng, 6)).astype("<i4")
        MT = rng.random((ng, 6)).astype("<f4")
        ex = [rng.random(ng).astype("<f4") for _ in range(3)]       # sfr, bh_mass, bh_mdot
        with open(d / f"group_tab_005.{k}", "wb") as f:
            f.write(struct.pack("<3iqi", ng, sum(ng_per_file), int(L.sum()), 10**6, nfiles))
            for a in (L, O, M, CM, V, LT, MT, *ex):
                f.write(a.tobytes())
        allg.append((L, M, CM, ex[1]))
    return allg


def test_fof_old(tmp_path):
    d = tmp_path / "groups_005"; d.mkdir()
    ref = _write_old(d, [0, 3, 2])
    f = FofCatalog.from_dir(str(tmp_path), 5)
    assert f.fmt == "old" and len(f) == 5
    np.testing.assert_array_equal(f.length, np.concatenate([r[0] for r in ref]))
    np.testing.assert_allclose(f.mass, np.concatenate([r[1] for r in ref]).sum(1))
    assert f.vel is None and "sfr" in f.extra


def test_fof_g3(tmp_path):
    d = tmp_path / "groups_005"; d.mkdir()
    ref = _write_g3(d, [2, 0, 4])
    f = FofCatalog.from_dir(str(tmp_path), 5)
    assert f.fmt == "g3" and len(f) == 6
    np.testing.assert_array_equal(f.cm, np.concatenate([r[2] for r in ref]))
    np.testing.assert_array_equal(f.extra["bh_mass"], np.concatenate([r[3] for r in ref]))


def test_params(tmp_path):
    p = tmp_path / "EGparams"
    p.write_text("""Snapshot snap_010
NumberOfBins 80
VolumeDensity
RotateGalaxyFrame 0
CmType -1
Cut 0
Pressure 20000    ; comment
DecompositionMethod 0
gadgetparamfile ../paramfile.par
#AquilaAlignment
ReadFofHalos
FofMinMass -1
FofRmin 0.5
FofRmax rvir
FofRgal 0
FofRalign 8.0
Quantities
# Quantity Weight Type Density Cumulate
MASS NULL 4 DENS CUMU    # 4
TEMP MCLD 0 AVE  DIFF
End
""")
    par = read_params(str(p))
    assert par["GadgetParamfile"] == "../paramfile.par" and par["FofRmax"] == 0.0
    assert not par["AquilaAlignment"] and par["NumberOfBins"] == 80
    assert [(q.quantity, q.ptype, q.flag) for q in par.quantities] == [("MASS", 4, 3), ("TEMP", 0, 0)]
    assert "VolumeDensity" in par.to_text()


def test_params_missing(tmp_path):
    p = tmp_path / "bad"
    p.write_text("Snapshot s\nVolumeDensity\n")
    with pytest.raises(ValueError):
        read_params(str(p))
