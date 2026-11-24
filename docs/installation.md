# Installation

galpost is a pure-Python package. It needs Python ≥ 3.9 and three libraries:
`numpy`, `scipy` and `matplotlib`. Nothing has to be compiled.

## With the installer

```bash
git clone https://github.com/gmurante/galpost.git
cd galpost
./install.sh
```

The script creates an environment called `galpost` and installs the package in
it, using the first tool it finds among `micromamba`, `mamba`, `conda`; if none
is available it creates a Python virtual environment in `~/.venvs/galpost`.
Run `./install.sh --help` for the options (environment name, location, Python
version, venv instead of conda).

Afterwards:

```bash
micromamba activate galpost        # or: conda activate galpost / . ~/.venvs/galpost/bin/activate
galpost --help
```

## By hand

### Virtual environment (any system with Python ≥ 3.9)

```bash
python3 -m venv ~/.venvs/galpost
. ~/.venvs/galpost/bin/activate
pip install -e /path/to/galpost          # -e: changes to the code are picked up immediately
```

### conda / mamba / micromamba

```bash
micromamba create -f environment.yml     # environment 'galpost'
micromamba activate galpost
pip install -e .
```

### Without installing

For a quick test, the code can be used directly from the source tree with the
system Python, provided numpy/scipy/matplotlib are available:

```bash
export PYTHONPATH=/path/to/galpost/src
python3 -m galpost.cli EGparams -snap snap_134
```

## Graphical interface

`galpost-gui` uses Tkinter, which is part of Python. conda/micromamba
environments include it; with a system Python on Linux it may require the
distribution package (`python3-tkinter` on RHEL/Fedora, `python3-tk` on
Debian/Ubuntu). From a cluster, connect with `ssh -X` (or `-Y`); under WSL2 the
window opens directly with WSLg.

## Running the tests

```bash
pip install -e .[test]
pytest
```

The tests build small synthetic snapshots and catalogues and do not need any
simulation data.

## Notes for clusters

* Install the environment in your home or project space; no root privileges
  are needed. `micromamba` is a single static executable that can be placed in
  `~/.local/bin`.
* Run the analysis on compute nodes for large snapshots (see
  [performance](performance.md)); a single node is enough, galpost is serial.
* On login nodes with *transparent huge pages* set to `always` large memory
  allocations can be very slow; galpost disables them for its own process
  automatically (see [performance](performance.md)).

## Notes for WSL2 and macOS

galpost works unchanged. Snapshots on Windows drives (`/mnt/c/...`) are read
much more slowly than files inside the WSL2 file system: copy large snapshots
to your Linux home first. Plots are written as PNG files and do not require a
display.

## OpenGadget

Running FoF and SubFind on the snapshots requires OpenGadget, which must be
installed separately. See [halo finding with OpenGadget](subfind_opengadget.md).
