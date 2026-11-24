# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Units and physical constants.

Constants and conversions are those of ExtractGalaxies / libreadgad, so that
results can be compared with the old outputs to numerical precision.

Output units (as in ExtractGalaxies):
  lengths   physical kpc     (code length / (1+z) [/h if Htrue])
  masses    Msun             (code mass [/h if Htrue])
  velocities km/s            (code velocity / sqrt(1+z), i.e. peculiar)
  densities Msun pc^-3
  times     Gyr
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict

__all__ = ["CONST", "read_gadget_params", "Units"]


class CONST:
    GRAVITY = 6.672e-8
    BOLTZMANN = 1.3806e-16
    PROTONMASS = 1.6726e-24
    GAMMA = 5.0 / 3
    GAMMA_MINUS1 = GAMMA - 1
    HYDROGEN_MASSFRAC = 0.76
    SOLAR_MASS = 1.989e33
    PI = 3.141592654
    CM_PER_KM = 1e5
    CM_PER_KPC = 3.085678e21
    SEC_PER_YEAR = 3.155e7
    T_THRESHOLD = 1.0e5          # cold / hot gas separation (K)


def read_gadget_params(path: str) -> Dict[str, str]:
    """Read a GADGET parameter file into {tag: value-string} (first token)."""
    out: Dict[str, str] = {}
    with open(path) as f:
        for line in f:
            w = line.split()
            if len(w) < 2 or w[0][0] in "%#;":
                continue
            out.setdefault(w[0], w[1])
    return out


@dataclass
class Units:
    unit_length_cm: float
    unit_mass_g: float
    unit_velocity_cms: float
    hubble: float                 # h of the snapshot header
    redshift: float
    htrue: bool
    periodic: bool = True
    comoving: bool = True

    @classmethod
    def from_files(cls, gadget_params: str, header, htrue: bool) -> "Units":
        p = read_gadget_params(gadget_params)
        return cls(float(p["UnitLength_in_cm"]), float(p["UnitMass_in_g"]),
                   float(p["UnitVelocity_in_cm_per_s"]), header.hubble_param,
                   header.redshift, htrue,
                   bool(int(p.get("PeriodicBoundariesOn", "1"))),
                   bool(int(p.get("ComovingIntegrationOn", "1"))))

    @property
    def h(self) -> float:
        """h used in conversions (1 unless Htrue)."""
        return self.hubble if self.htrue else 1.0

    @property
    def mass_msun(self) -> float:
        return self.unit_mass_g / CONST.SOLAR_MASS / self.h

    @property
    def length_kpc(self) -> float:
        return self.unit_length_cm / CONST.CM_PER_KPC / (1.0 + self.redshift) / self.h

    @property
    def vel_kms(self) -> float:
        return self.unit_velocity_cms / CONST.CM_PER_KM / math.sqrt(1.0 + self.redshift)

    @property
    def dens_msun_pc3(self) -> float:
        return self.mass_msun / (1000.0 * self.length_kpc) ** 3

    @property
    def temp_K(self) -> float:
        """Multiply the specific internal energy (code units) to get T for mu=1."""
        return (CONST.GAMMA_MINUS1 * CONST.PROTONMASS / CONST.BOLTZMANN
                * self.unit_velocity_cms ** 2)

    def summary(self) -> str:
        return ("Unit transformation constants:\n"
                f"Mass_from_codeu_to_msun: {self.mass_msun:g}\n"
                f"Length_from_codeu_to_kpc: {self.length_kpc:g}\n"
                f"Vel_from_codeu_to_kms: {self.vel_kms:g}\n"
                f"Dens_from_codeu_to_msun_per_pc3: {self.dens_msun_pc3:g}\n"
                f"Temp_from_codeu_to_K: {self.temp_K:g}\n")
