# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Synthetic tests of the extra analyses."""
import math
import numpy as np

from galpost.extra.common import G, solar_log_ratio
from galpost.extra.structure import circularity_energy, spherical_potential


class _FakeView:
    """Point mass at the centre + test stars on circular orbits in the x-y plane."""
    def __init__(self, n=200, mpoint=1e11, sign=1.0):
        rng = np.random.default_rng(5)
        r = rng.uniform(1.0, 20.0, n)
        phi = rng.uniform(0, 2 * np.pi, n)
        vc = np.sqrt(G * mpoint / r)
        pos = np.column_stack([r * np.cos(phi), r * np.sin(phi), np.zeros(n)])
        vel = sign * np.column_stack([-vc * np.sin(phi), vc * np.cos(phi), np.zeros(n)])
        self.pos = np.vstack([[0, 0, 0], pos])
        self.vel = np.vstack([[0, 0, 0], vel])
        self.mass = np.r_[mpoint, np.full(n, 1e-6 * mpoint / n)]
        self.type = np.r_[1, np.full(n, 4)]
        self.r = np.sqrt((self.pos ** 2).sum(1))

        class G_:
            rvir = 50.0
        self.g = G_()

    def of(self, t):
        return self.type == t


def test_spherical_potential_point_mass():
    r = np.array([0.0])
    m = np.array([1e10])
    rg = np.array([1.0, 10.0])
    phi, menc = spherical_potential(r, m, rg)
    np.testing.assert_allclose(phi, -G * 1e10 / rg, rtol=1e-12)


def test_circularity_energy_circular_orbits():
    eps, st, _ = circularity_energy(_FakeView(sign=+1.0), 50.0)
    np.testing.assert_allclose(eps, 1.0, atol=5e-3)
    eps, st, _ = circularity_energy(_FakeView(sign=-1.0), 50.0)
    np.testing.assert_allclose(eps, -1.0, atol=5e-3)


def test_solar_ratio():
    # [Fe/H] of the Sun is 0 by construction; mass ratio Fe/H ~ 1.3e-3 / 0.74
    assert abs(10 ** solar_log_ratio("Fe") - 1.29e-3 / 0.738) / (1.29e-3 / 0.738) < 0.05


def test_shrinking_sphere_finds_dense_clump():
    from galpost.track import shrinking_sphere
    rng = np.random.default_rng(7)
    clump = rng.normal([10.0, 20.0, 30.0], 0.5, (300, 3))        # main clump
    other = rng.normal([60.0, 20.0, 30.0], 3.0, (150, 3))        # smaller, diffuse one
    c = shrinking_sphere(np.vstack([clump, other]), np.array([30.0, 20.0, 30.0]), 1e5, False)
    assert np.linalg.norm(c - [10.0, 20.0, 30.0]) < 0.3


def test_so_radius_uniform_sphere():
    from galpost.track import so_radius
    rng = np.random.default_rng(8)
    n, R, M = 200000, 100.0, 1e12
    r = R * rng.random(n) ** (1.0 / 3.0)                          # uniform sphere
    rho = M / (4.0 / 3.0 * math.pi * R ** 3)
    # threshold = half the mean density -> R_delta ~ R (density drops to 0 outside)
    rd, md = so_radius(r, np.full(n, M / n), rho / 400.0, 200.0)
    assert abs(rd - R) / R < 0.01 and abs(md - M) / M < 0.01
