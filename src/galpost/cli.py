# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Command line: ``galpost <paramfile> [-snap <snapshot>]`` (same syntax as
ExtractGalaxies.x)."""

from __future__ import annotations

import argparse
import sys

from .pipeline import run


def main(argv=None):
    ap = argparse.ArgumentParser(prog="galpost",
                                 description="Galaxy post-processing of OpenGadget/MUPPI "
                                             "snapshots (ExtractGalaxies replacement)")
    ap.add_argument("paramfile", help="parameter file (ExtractGalaxies syntax)")
    ap.add_argument("-snap", "--snap", default=None,
                    help="snapshot file (overrides the Snapshot tag)")
    a = ap.parse_args(argv)
    run(a.paramfile, a.snap)
    return 0


if __name__ == "__main__":
    sys.exit(main())
