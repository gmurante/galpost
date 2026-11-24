# galpost

**Galaxy post-processing for OpenGadget / MUPPI simulations, in pure Python.**

galpost analyses the galaxies of a GADGET/OpenGadget snapshot: it finds them in
FoF or SubFind catalogues, defines their centre and reference frame, computes
radial profiles of any quantity and a detailed set of global properties, and
writes galaxy-centred snapshots.

It is a re-implementation of the C code *ExtractGalaxies* (P. Monaco and
collaborators): it reads the same parameter files and writes the same output
files, column by column, so existing plotting scripts keep working. On top of
that it adds a set of modern analyses (kinematic decomposition from the binding
energy, stellar populations, disc structure, Toomre stability, bars, gas flows,
Kennicutt–Schmidt relation, metallicity gradients, dust, circumgalactic medium,
black holes, satellites, maps, scaling relations), each one switched on by its
own parameter.

* only `numpy`, `scipy` and `matplotlib` are required;
* runs on any Unix system, WSL2 and macOS;
* reads GADGET format 2 (with or without INFO block) and format 1, multi-file
  snapshots, MUPPI companion files (`.muppi`), files larger than 4 GB;
* validated against ExtractGalaxies outputs (see
  [docs/compatibility.md](docs/compatibility.md)).

## Install

```bash
git clone https://github.com/gmurante/galpost.git
cd galpost
./install.sh                 # creates a 'galpost' environment and installs the package
```

or by hand, with any Python ≥ 3.9:

```bash
python3 -m venv ~/.venvs/galpost && . ~/.venvs/galpost/bin/activate
pip install -e .
```

See [docs/installation.md](docs/installation.md) for conda/micromamba, clusters
and WSL2.

## Quick start

```bash
galpost EGparams.surf -snap snap_134
```

`EGparams.surf` is an ExtractGalaxies parameter file, for instance

```
NumberOfBins        80
SurfaceDensity
RotateGalaxyFrame   0
CmType              -1
Cut                 0
Pressure            20000
DecompositionMethod 0
Htrue
AnalyzeGalaxies
AlignWithAngularMomentum
GadgetParamfile     ../paramfile.par
MetalsFile          ../metals.dat
Outdir              EG-official
ReadFofHalos
FofMinMass          -1
FofRmin             0.5
FofRmax             rvir
FofRgal             0
FofRalign           8.0
FofRlog
Quantities
MASS NULL 4 DENS DIFF
SFR  NULL 0 DENS DIFF
End
```

Add, for example, `KinematicDecomposition`, `StellarPopulations` and
`ScalingRelations` to get the new analyses. See
[docs/quickstart.md](docs/quickstart.md).

## Documentation

| | |
|---|---|
| [Installation](docs/installation.md) | environments, clusters, WSL2 |
| [Quick start](docs/quickstart.md) | a complete run, from SubFind to plots |
| [The parameter file](docs/parameter_file.md) | syntax, quantities, galaxy lists |
| [Parameter reference](docs/parameters.md) | every tag (generated from the code) |
| [Output files](docs/outputs.md) | every file and column |
| [Extra analyses](docs/extra_analyses.md) | methods and definitions of the new analyses |
| [Halo finding with OpenGadget](docs/subfind_opengadget.md) | FoF + SubFind in post-processing |
| [Snapshot reader](docs/snapshot_reader.md) | block tables, INFO, formats |
| [Progenitor tracking](docs/tracking.md) | follow a galaxy back in time with `galpost-track` |
| [Graphical interface](docs/gui.md) | configure, run and plot with `galpost-gui` |
| [Python API](docs/python_api.md) | using galpost from Python |
| [Compatibility with ExtractGalaxies](docs/compatibility.md) | validation, fixed bugs, differences |
| [Performance](docs/performance.md) | memory, timing, cluster tips |

## Requirements

Python ≥ 3.9, numpy ≥ 1.20, scipy ≥ 1.7, matplotlib ≥ 3.4. Running FoF/SubFind
requires a working OpenGadget installation, which is not part of this package.

## Authors

Giuseppe Murante (INAF – Osservatorio Astrofisico di Torino) and
Pierluigi Monaco (Università degli Studi di Trieste).

galpost re-implements and extends ExtractGalaxies, written by P. Monaco with
contributions by G. Murante, L. Tornatore and collaborators.

## License

Copyright (C) 2026 INAF – Istituto Nazionale di Astrofisica and
Università degli Studi di Trieste.

galpost is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See [LICENSE](LICENSE) for details.

If you use galpost in a publication, please cite it (see
[CITATION.cff](CITATION.cff)).
