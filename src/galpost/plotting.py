# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Reading and plotting galpost / ExtractGalaxies output files.

Independent of any GUI: used by the galpost GUI and by the command

    galpost-plot <file> [<file> ...] [--x COL] [--y COL ...] [--logx] [--logy]
                 [--out figure.pdf|png]

Every output file is turned into a :class:`Table` (column names + data); a
registry of file kinds gives sensible default plots, and any pair of columns
can be plotted. Several files of the same kind are overlaid.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

__all__ = ["Table", "read_table", "file_kind", "KINDS", "plot_files", "plot_map",
           "list_output_files", "main"]


# --------------------------------------------------------------------------
@dataclass
class Table:
    path: str
    kind: str
    columns: List[str]
    data: np.ndarray
    header: List[str] = field(default_factory=list)
    label: str = ""

    def col(self, name_or_index) -> np.ndarray:
        if isinstance(name_or_index, (int, np.integer)):
            return self.data[:, int(name_or_index)]
        return self.data[:, self.columns.index(name_or_index)]


@dataclass
class Kind:
    name: str
    description: str
    x: Optional[int] = 0              # default x column (index)
    y: Sequence[int] = (1,)           # default y columns
    logx: bool = False
    logy: bool = False
    plot: str = "line"                # line | step | image | text | map | scatter


# file prefix -> kind (longest prefix wins)
KINDS: Dict[str, Kind] = {
    "profiles": Kind("profiles", "Profiles (ExtractGalaxies)", 0, (3,), True, True),
    "sfr": Kind("sfr", "Star formation history", 0, (1,), False, False, "step"),
    "jhist": Kind("jhist", "Circularity jz/jcirc(r) histogram", 0, (1,), False, False, "step"),
    "jhistE": Kind("jhistE", "Circularity jz/jcirc(E) histogram", 0, (1,), False, False, "step"),
    "info": Kind("info", "Galaxy information", None, (), plot="text"),
    "rotation_matrix": Kind("rotation_matrix", "Rotation matrix", None, (), plot="text"),
    "list_of_selected_gals": Kind("list_of_selected_gals", "Selected galaxies", None, (), plot="text"),
    "list_of_blocks": Kind("list_of_blocks", "Available quantities", None, (), plot="text"),
    "kindec": Kind("kindec", "Kinematic decomposition", None, (), plot="text"),
    "kinematics": Kind("kinematics", "Kinematic profiles", 2, (3, 4, 5, 6)),
    "stellarpop": Kind("stellarpop", "Stellar populations", 2, (4,)),
    "mdf": Kind("mdf", "Metallicity distribution", 0, (2,), plot="step"),
    "ofe_feh": Kind("ofe_feh", "[O/Fe]-[Fe/H] plane", None, (), plot="image"),
    "disc": Kind("disc", "Disc structure", 2, (3,), False, True),
    "toomre": Kind("toomre", "Toomre Q", 2, (9, 10, 11), False, True),
    "bar": Kind("bar", "Bar strength A2(R)", 2, (3,)),
    "flows": Kind("flows", "Spherical flows", 0, (3, 4), False, True),
    "ks": Kind("ks", "Kennicutt-Schmidt", 4, (7,), True, True, "scatter"),
    "gradients": Kind("gradients", "Metallicity gradients", 2, (3, 5)),
    "dust": Kind("dust", "Dust", 2, (5,), False, True),
    "cgm": Kind("cgm", "Circumgalactic medium", None, (), plot="text"),
    "bh": Kind("bh", "Black holes", None, (), plot="text"),
    "satellites": Kind("satellites", "Satellites", 3, (4,), False, True, "scatter"),
    "scaling": Kind("scaling", "Scaling relations", 3, (7,), True, True, "scatter"),
    "table-scaling": Kind("table-scaling", "Scaling relations vs time", 2, (3,), False, True),
    "table-mass": Kind("table-mass", "Masses vs time", 2, (7,), False, True),
    "table-info": Kind("table-info", "Galaxy vs time", None, (), plot="text"),
    "table-spin": Kind("table-spin", "Spins vs time", None, (), plot="text"),
    "table-specspin": Kind("table-specspin", "Specific spins vs time", None, (), plot="text"),
    "table-inertia": Kind("table-inertia", "Inertia vs time", None, (), plot="text"),
    "track": Kind("track", "Main-progenitor track", 2, (9,), False, True),
    "map_": Kind("map", "Surface density map", None, (), plot="map"),
    "zoom": Kind("zoom", "Galaxy snapshot", None, (), plot="none"),
    "catalog_": Kind("catalog", "Particle catalogue", 8, (6,), True, True, "scatter"),
    "particle_list_": Kind("particle_list", "Particle list", None, (), plot="text"),
    "fullcatalog_": Kind("fullcatalog", "Full catalogue", None, (), plot="none"),
}


def file_kind(path: str) -> Optional[Kind]:
    base = os.path.basename(path)
    best = None
    for pref, k in KINDS.items():
        if base.startswith(pref) and (len(base) == len(pref) or not base[len(pref)].isalpha()
                                      or pref.endswith("_")):
            if best is None or len(pref) > len(best[0]):
                best = (pref, k)
    return best[1] if best else None


def short_label(path: str) -> str:
    """'profiles.vol.snap_134.gal0' -> 'vol snap_134 gal0' (directory if ambiguous)."""
    base = os.path.basename(path)
    parts = base.split(".")[1:]
    parts = [p for p in parts if p not in ("npy", "png")]
    lab = " ".join(parts) if parts else base
    return lab


# --------------------------------------------------------------------------
_NUMCOL = re.compile(r"#\s*(\d+)\s*:\s*(.*)$")


def _numeric_rows(lines):
    rows = []
    for l in lines:
        s = l.strip()
        if not s or s.startswith("#"):
            continue
        try:
            rows.append([float(x) for x in s.split()])
        except ValueError:
            return None
    if not rows:
        return np.zeros((0, 0))
    n = max(len(r) for r in rows)
    out = np.full((len(rows), n), np.nan)
    for i, r in enumerate(rows):
        out[i, :len(r)] = r
    return out


def read_table(path: str) -> Table:
    """Parse any galpost/ExtractGalaxies text output into a Table."""
    kind = file_kind(path)
    kname = kind.name if kind else "table"
    with open(path, errors="replace") as f:
        lines = f.read().splitlines()
    header = [l for l in lines if l.startswith("#")]
    data = _numeric_rows(lines)
    if data is None:                      # free text (info.*)
        return Table(path, kname, [], np.zeros((0, 0)), lines, short_label(path))
    ncol = data.shape[1] if data.size else 0
    cols: List[str] = []
    if kname == "profiles":
        h = header[1] if len(header) > 1 else ""
        names = ["Rmed (kpc)", "Rmin (kpc)", "Rmax (kpc)"]
        for m in re.finditer(r"(\d+):\s*([A-Za-z0-9]+)\s*\[Type (\d)\]", h):
            q = f"{m.group(2)} [Type {m.group(3)}]"
            names += [q, f"N particles ({q})"]
        if "Bin size" in h:
            names[:3] = ["cell size (kpc)", "x (kpc)", "y (kpc)"]
        cols = names
    elif kname == "sfr":
        cols = ["time (Gyr)", "SFR (Msun/yr)"]
    elif kname in ("jhist", "jhistE"):
        rb = []
        for l in header:
            if "radial bin centers" in l:
                rb = l.split(":", 1)[1].split()
        cols = ["circularity", "mass (all radii)"] + [f"mass at R={x}" for x in rb]
        if kname == "jhistE":
            cols = []                       # numbered header, parsed below
    if not cols:
        for l in header:
            m = _NUMCOL.match(l)
            if m:
                cols.append(m.group(2).strip())
    if kname == "list_of_selected_gals" or (cols and len(cols) != ncol):
        cols = []
    while len(cols) < ncol:
        cols.append(f"column {len(cols) + 1}")
    return Table(path, kname, cols[:ncol] if ncol else cols, data, header, short_label(path))


def list_output_files(directory: str) -> Dict[str, List[str]]:
    """Output files of a directory grouped by kind description."""
    groups: Dict[str, List[str]] = {}
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return groups
    for n in names:
        p = os.path.join(directory, n)
        if not os.path.isfile(p) or n.endswith((".png", ".extent")) or n.startswith("."):
            continue
        k = file_kind(n)
        desc = k.description if k else "Other files"
        if k and k.plot == "none":
            continue
        groups.setdefault(desc, []).append(p)
    return groups


# --------------------------------------------------------------------------
def plot_map(ax, path: str, log: bool = True):
    """Show a map (.npy array in Msun/pc^2, or a PNG image)."""
    import matplotlib.image as mpimg
    if path.endswith(".png"):
        ax.imshow(mpimg.imread(path))
        ax.axis("off")
        return None
    arr = np.load(path)
    pos = arr[arr > 0]
    vmin = np.percentile(pos, 1) if pos.size else 1e-3
    extent = None
    ext_file = path[:-4] + ".extent" if path.endswith(".npy") else None
    if ext_file and os.path.exists(ext_file):
        try:
            extent = [float(x) for x in open(ext_file).read().split("\n")[1].split()[:4]]
        except (ValueError, IndexError):
            extent = None
    im = ax.imshow(np.log10(np.maximum(arr, vmin)) if log else arr, origin="lower",
                   cmap="inferno", extent=extent)
    if extent:
        ax.set_xlabel("x (kpc)")
        ax.set_ylabel("z (kpc)" if "edgeon" in path else "y (kpc)")
    ax.set_title(os.path.basename(path).replace(".npy", ""), fontsize=9)
    return im


def plot_files(fig, paths: Sequence[str], x=None, y=None, logx=None, logy=None,
               labels: Optional[Sequence[str]] = None, title: Optional[str] = None,
               style: Optional[str] = None):
    """Plot one or more output files of the same kind on a matplotlib figure.

    x, y: column indices (default: those of the file kind). Returns the axes,
    or the text to show for text-only files.
    """
    fig.clear()
    tables = [read_table(p) for p in paths]
    t0 = tables[0]
    kind = file_kind(t0.path) or Kind("table", "Table", 0, (1,))
    if kind.plot == "map" or t0.path.endswith(".npy"):
        n = len(paths)
        ncol = min(n, 3)
        nrow = (n + ncol - 1) // ncol
        for i, p in enumerate(paths):
            ax = fig.add_subplot(nrow, ncol, i + 1)
            im = plot_map(ax, p)
            if im is not None:
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=r"log $\Sigma$ (M$_\odot$ pc$^{-2}$)")
        fig.tight_layout()
        return fig.axes
    if kind.plot == "text" and x is None and not y:
        return "\n\n".join(open(p, errors="replace").read() for p in paths)
    if kind.plot == "image" and x is None and not y:
        ax = fig.add_subplot(111)
        arr = t0.data
        ext = None
        nums = re.findall(r"from ([-\d.e]+) to ([-\d.e]+)", " ".join(t0.header))
        if len(nums) >= 2:
            ext = [float(nums[1][0]), float(nums[1][1]), float(nums[0][0]), float(nums[0][1])]
        im = ax.imshow(np.log10(np.maximum(arr, arr[arr > 0].min() if (arr > 0).any() else 1)),
                       origin="lower", aspect="auto", extent=ext, cmap="viridis")
        ax.set_xlabel("[O/Fe]"); ax.set_ylabel("[Fe/H]")
        fig.colorbar(im, ax=ax, label="log stellar mass")
        ax.set_title(title or t0.label, fontsize=9)
        return [ax]
    if not t0.data.size:
        return "\n".join(open(paths[0], errors="replace").read().splitlines())

    xi = kind.x if x is None else x
    yi = list(kind.y) if not y else list(y)
    lx = kind.logx if logx is None else logx
    ly = kind.logy if logy is None else logy
    st = style or kind.plot
    if st == "text":
        st = "line"
    ax = fig.add_subplot(111)
    multi_y = len(yi) > 1
    for it, t in enumerate(tables):
        xx = t.col(xi) if xi is not None else np.arange(len(t.data))
        for j in yi:
            if j >= t.data.shape[1]:
                continue
            yy = t.col(j)
            lab = []
            if len(tables) > 1:
                lab.append(labels[it] if labels else t.label)
            if multi_y:
                lab.append(t.columns[j])
            lab = " – ".join(lab) if lab else None
            if st == "step":
                ax.step(xx, yy, where="mid", label=lab)
            elif st == "scatter":
                ax.plot(xx, yy, "o", ms=4, label=lab)
            else:
                ax.plot(xx, yy, "-", marker="." if len(xx) < 40 else None, label=lab)
    if lx:
        ax.set_xscale("log")
    if ly:
        ax.set_yscale("log")
    ax.set_xlabel(t0.columns[xi] if xi is not None else "row")
    if not multi_y:
        ax.set_ylabel(t0.columns[yi[0]] if yi and yi[0] < len(t0.columns) else "")
    if (len(tables) > 1 or multi_y) and ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)
    ax.set_title(title or (t0.label if len(tables) == 1 else kind.description), fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return [ax]


# --------------------------------------------------------------------------
def main(argv=None):
    import argparse
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(prog="galpost-plot",
                                 description="Plot galpost/ExtractGalaxies output files")
    ap.add_argument("files", nargs="+")
    ap.add_argument("--x", type=int, default=None, help="x column (1-based)")
    ap.add_argument("--y", type=int, nargs="*", default=None, help="y columns (1-based)")
    ap.add_argument("--logx", action="store_true", default=None)
    ap.add_argument("--logy", action="store_true", default=None)
    ap.add_argument("--linx", action="store_true", help="force linear x")
    ap.add_argument("--liny", action="store_true", help="force linear y")
    ap.add_argument("--columns", action="store_true", help="list the columns and exit")
    ap.add_argument("--out", default=None, help="output file (.png, .pdf, ...)")
    ap.add_argument("--dpi", type=int, default=150)
    a = ap.parse_args(argv)
    if a.columns:
        t = read_table(a.files[0])
        for i, c in enumerate(t.columns):
            print(f"{i + 1:3d}  {c}")
        return 0
    fig = plt.figure(figsize=(7, 5))
    res = plot_files(fig, a.files,
                     x=None if a.x is None else a.x - 1,
                     y=None if not a.y else [k - 1 for k in a.y],
                     logx=False if a.linx else a.logx, logy=False if a.liny else a.logy)
    if isinstance(res, str):
        print(res)
        return 0
    out = a.out or (os.path.basename(a.files[0]) + ".png")
    fig.savefig(out, dpi=a.dpi)
    print(f"written {out}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
