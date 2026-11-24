#!/usr/bin/env bash
# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later
# Installer of galpost: creates an environment and installs the package in it.
#
#   ./install.sh                      conda-like tool if available, otherwise a venv
#   ./install.sh --venv               force a Python virtual environment
#   ./install.sh --name NAME          environment name (default: galpost)
#   ./install.sh --prefix DIR         venv directory (default: ~/.venvs/<name>)
#   ./install.sh --python 3.12        Python version for conda environments
#   ./install.sh --no-test            do not run the test suite at the end
set -euo pipefail

NAME=galpost
PREFIX=""
PYVER=3.12
USE_VENV=0
RUN_TESTS=1
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [ $# -gt 0 ]; do
    case "$1" in
        --venv) USE_VENV=1 ;;
        --name) NAME="$2"; shift ;;
        --prefix) PREFIX="$2"; shift ;;
        --python) PYVER="$2"; shift ;;
        --no-test) RUN_TESTS=0 ;;
        -h|--help) sed -n '2,9p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
    shift
done

TOOL=""
if [ $USE_VENV -eq 0 ]; then
    for t in micromamba mamba conda; do
        if command -v $t >/dev/null 2>&1; then TOOL=$t; break; fi
    done
fi

if [ -n "$TOOL" ]; then
    echo ">>> creating environment '$NAME' with $TOOL (Python $PYVER)"
    $TOOL create -y -n "$NAME" -c conda-forge "python=$PYVER" numpy scipy matplotlib pytest pip
    RUN="$TOOL run -n $NAME"
    $RUN python -m pip install -e "$HERE"
    ACTIVATE="$TOOL activate $NAME"
else
    PY=$(command -v python3 || true)
    [ -n "$PY" ] || { echo "python3 not found" >&2; exit 1; }
    "$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)' \
        || { echo "Python >= 3.9 is required" >&2; exit 1; }
    PREFIX=${PREFIX:-$HOME/.venvs/$NAME}
    echo ">>> creating virtual environment in $PREFIX"
    "$PY" -m venv "$PREFIX"
    "$PREFIX/bin/python" -m pip install --upgrade pip
    "$PREFIX/bin/python" -m pip install -e "$HERE[test]"
    RUN="$PREFIX/bin/"
    ACTIVATE=". $PREFIX/bin/activate"
fi

if [ $RUN_TESTS -eq 1 ]; then
    echo ">>> running the tests"
    if [ -n "$TOOL" ]; then (cd "$HERE" && $RUN pytest -q); else (cd "$HERE" && "${RUN}pytest" -q); fi
fi

echo
echo "galpost installed. Activate the environment with:"
echo "    $ACTIVATE"
echo "and run:"
echo "    galpost <parameter file> -snap <snapshot>"
