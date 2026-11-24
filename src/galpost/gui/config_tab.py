# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""'Configure' tab: edit a parameter file.

The widgets are generated from galpost.params.TAGS, so new parameters appear
automatically.
"""

from __future__ import annotations

import os
import tkinter as tk
import warnings
from tkinter import filedialog, messagebox, ttk

from ..params import (GalaxyLine, LOGICAL, Params, Quantity, TAGS, _SECTION_TITLES,
                      read_params)
from ..profiles import KNOWN_BLOCKS, derived_names
from .widgets import ScrollFrame, Tooltip

RADIUS_NAMES = ["rmin", "rmax", "rvir", "rgal", "ralign"]
FILE_TAGS = {"GadgetParamfile", "MetalsFile", "Galfile", "BlockTable", "LoadParticleList",
             "Snapshot"}
DIR_TAGS = {"Outdir", "Snapdir", "GroupDir"}
COMBO_TAGS = {"GalaxySnapshotRadius": RADIUS_NAMES, "ExtraRmax": RADIUS_NAMES,
              "MapSize": RADIUS_NAMES, "Cut": ["0"] + RADIUS_NAMES, "FofRmax": ["rvir", "0"],
              "DiscFitRmin": ["auto"], "DiscFitRmax": ["auto"], "GradientRmin": ["auto"],
              "GradientRmax": ["auto"]}
DEFAULT_METALS = ["He", "C", "Ca", "O", "N", "Ne", "Mg", "S", "Si", "Fe", "Na", "Al", "Ar",
                  "Ni"]


class ConfigTab(ttk.Frame):
    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.path = None
        self.vars = {}
        self.quantities = []
        self.galaxies = []
        self._build()
        self.new()

    # ---- layout ------------------------------------------------------------
    def _build(self):
        bar = ttk.Frame(self)
        bar.pack(fill="x", padx=6, pady=4)
        for txt, cmd in (("New", self.new), ("Open…", self.open), ("Save", self.save),
                         ("Save as…", self.save_as), ("Check", self.check),
                         ("Show text", self.show_text)):
            ttk.Button(bar, text=txt, command=cmd).pack(side="left", padx=2)
        self.path_label = ttk.Label(bar, text="(new parameter file)", foreground="#555")
        self.path_label.pack(side="left", padx=12)

        self.nb = ttk.Notebook(self)
        self.nb.pack(fill="both", expand=True, padx=6, pady=4)
        for sec, title in _SECTION_TITLES.items():
            sf = ScrollFrame(self.nb)
            self.nb.add(sf, text=title)
            self._section(sf.inner, sec)
        self._quantities_page()
        self._galaxies_page()

    def _section(self, parent, sec):
        parent.columnconfigure(1, weight=1)
        row = 0
        for t in TAGS.values():
            if t.section != sec:
                continue
            fg = "#b00000" if t.required else "black"
            if t.kind == LOGICAL:
                v = tk.BooleanVar(value=False)
                w = ttk.Checkbutton(parent, variable=v, text=t.name)
                w.grid(row=row, column=0, columnspan=2, sticky="w", padx=6, pady=2)
            else:
                v = tk.StringVar()
                lab = tk.Label(parent, text=t.name, fg=fg, anchor="w")
                lab.grid(row=row, column=0, sticky="w", padx=6, pady=2)
                fr = ttk.Frame(parent)
                fr.grid(row=row, column=1, sticky="ew", padx=4)
                fr.columnconfigure(0, weight=1)
                if t.name in COMBO_TAGS:
                    w = ttk.Combobox(fr, textvariable=v, values=COMBO_TAGS[t.name])
                else:
                    w = ttk.Entry(fr, textvariable=v)
                w.grid(row=0, column=0, sticky="ew")
                if t.name in FILE_TAGS or t.name in DIR_TAGS:
                    ttk.Button(fr, text="…", width=3,
                               command=lambda n=t.name: self._browse(n)).grid(row=0, column=1)
                Tooltip(lab, t.help)
            Tooltip(w, t.help + ("" if t.default in (None, "", False) else f"\n(default: {t.default})"))
            dflt = "" if t.default in (None, "", False) or t.kind == LOGICAL else f"  [default {t.default}]"
            req = "  REQUIRED" if t.required else ""
            ttk.Label(parent, text=t.help + dflt + req, foreground="#666", wraplength=520,
                      justify="left").grid(row=row, column=2, sticky="w", padx=8)
            self.vars[t.name] = v
            row += 1

    def _browse(self, name):
        start = self.base_dir()
        if name in DIR_TAGS:
            p = filedialog.askdirectory(initialdir=start, parent=self)
        else:
            p = filedialog.askopenfilename(initialdir=start, parent=self)
        if p:
            self.vars[name].set(self._relative(p))

    def _relative(self, p):
        try:
            rel = os.path.relpath(p, self.base_dir())
            return rel if not rel.startswith("../../..") else p
        except ValueError:
            return p

    def base_dir(self):
        return os.path.dirname(os.path.abspath(self.path)) if self.path else os.getcwd()

    # ---- quantities ----------------------------------------------------------
    def _quantities_page(self):
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="Quantities (profiles)")
        cols = ("quantity", "weight", "type", "dens", "cumu")
        self.qtree = ttk.Treeview(f, columns=cols, show="headings", height=18)
        for c, w, txt in zip(cols, (110, 110, 60, 90, 90),
                             ("Quantity", "Weight", "Type", "DENS/AVE", "DIFF/CUMU")):
            self.qtree.heading(c, text=txt)
            self.qtree.column(c, width=w, anchor="center")
        self.qtree.grid(row=0, column=0, sticky="nsew", padx=6, pady=6)
        self.qtree.bind("<Double-1>", lambda e: self._edit_quantity())
        f.columnconfigure(0, weight=1)
        f.rowconfigure(0, weight=1)
        b = ttk.Frame(f)
        b.grid(row=0, column=1, sticky="n", pady=6)
        for txt, cmd in (("Add…", self._add_quantity), ("Edit…", self._edit_quantity),
                         ("Remove", self._remove_quantity), ("Up", lambda: self._move_q(-1)),
                         ("Down", lambda: self._move_q(1))):
            ttk.Button(b, text=txt, command=cmd, width=10).pack(pady=2)
        ttk.Label(f, foreground="#666", wraplength=700, justify="left",
                  text="Each line is one profile. DENS: surface/volume density of the quantity "
                       "(Msun/pc^2 or /pc^3 for masses); AVE: weighted average. CUMU: cumulative "
                       "within the outer edge of the bin. Double-click to edit."
                  ).grid(row=1, column=0, columnspan=2, sticky="w", padx=6, pady=4)

    def _quantity_names(self):
        metals = DEFAULT_METALS
        mf = self.vars["MetalsFile"].get().strip()
        if mf:
            p = mf if os.path.isabs(mf) else os.path.join(self.base_dir(), mf)
            try:
                with open(p) as fh:
                    names = [l.split()[0][:3] for l in fh if l.split()]
                metals = [n for n in names if n != "H"] or metals
            except OSError:
                pass
        known = [b for b in KNOWN_BLOCKS if b not in ("HEAD", "POS", "VEL", "ID", "Zs")]
        return sorted(set(known + derived_names(metals)), key=str.lower)

    def _refresh_quantities(self):
        self.qtree.delete(*self.qtree.get_children())
        for q in self.quantities:
            self.qtree.insert("", "end", values=(q.quantity, q.weight, q.ptype,
                                                 "DENS" if q.dens else "AVE",
                                                 "CUMU" if q.cumulative else "DIFF"))

    def _quantity_dialog(self, q=None):
        d = tk.Toplevel(self)
        d.title("Profile quantity")
        d.transient(self.winfo_toplevel())
        names = self._quantity_names()
        vq = tk.StringVar(value=q.quantity if q else "MASS")
        vw = tk.StringVar(value=q.weight if q else "NULL")
        vt = tk.IntVar(value=q.ptype if q else 0)
        vd = tk.StringVar(value=("DENS" if q.dens else "AVE") if q else "DENS")
        vc = tk.StringVar(value=("CUMU" if q.cumulative else "DIFF") if q else "DIFF")
        ttk.Label(d, text="Quantity").grid(row=0, column=0, sticky="w", padx=6, pady=3)
        ttk.Combobox(d, textvariable=vq, values=names, width=12).grid(row=0, column=1, padx=6)
        ttk.Label(d, text="Weight").grid(row=1, column=0, sticky="w", padx=6, pady=3)
        ttk.Combobox(d, textvariable=vw, values=["NULL"] + names, width=12).grid(row=1, column=1, padx=6)
        ttk.Label(d, text="Particle type").grid(row=2, column=0, sticky="w", padx=6, pady=3)
        ttk.Spinbox(d, from_=0, to=5, textvariable=vt, width=5).grid(row=2, column=1, sticky="w", padx=6)
        fr = ttk.Frame(d)
        fr.grid(row=3, column=0, columnspan=2, pady=3)
        for txt in ("DENS", "AVE"):
            ttk.Radiobutton(fr, text=txt, value=txt, variable=vd).pack(side="left", padx=4)
        fr2 = ttk.Frame(d)
        fr2.grid(row=4, column=0, columnspan=2, pady=3)
        for txt in ("DIFF", "CUMU"):
            ttk.Radiobutton(fr2, text=txt, value=txt, variable=vc).pack(side="left", padx=4)
        ttk.Label(d, text="types: 0 gas, 1 DM, 2-3 low-res DM, 4 stars, 5 BH",
                  foreground="#666").grid(row=5, column=0, columnspan=2, padx=6)
        res = {}

        def ok():
            res["q"] = Quantity(vq.get().strip(), vw.get().strip() or "NULL", int(vt.get()),
                                vd.get() == "DENS", vc.get() == "CUMU")
            d.destroy()
        b = ttk.Frame(d)
        b.grid(row=6, column=0, columnspan=2, pady=6)
        ttk.Button(b, text="OK", command=ok).pack(side="left", padx=4)
        ttk.Button(b, text="Cancel", command=d.destroy).pack(side="left", padx=4)
        d.grab_set()
        self.wait_window(d)
        return res.get("q")

    def _sel_index(self, tree):
        s = tree.selection()
        return tree.index(s[0]) if s else None

    def _add_quantity(self):
        q = self._quantity_dialog()
        if q:
            self.quantities.append(q)
            self._refresh_quantities()

    def _edit_quantity(self):
        i = self._sel_index(self.qtree)
        if i is None:
            return
        q = self._quantity_dialog(self.quantities[i])
        if q:
            self.quantities[i] = q
            self._refresh_quantities()

    def _remove_quantity(self):
        i = self._sel_index(self.qtree)
        if i is not None:
            del self.quantities[i]
            self._refresh_quantities()

    def _move_q(self, step):
        i = self._sel_index(self.qtree)
        if i is None or not 0 <= i + step < len(self.quantities):
            return
        q = self.quantities
        q[i], q[i + step] = q[i + step], q[i]
        self._refresh_quantities()
        self.qtree.selection_set(self.qtree.get_children()[i + step])

    # ---- galaxy list ----------------------------------------------------------
    def _galaxies_page(self):
        f = ttk.Frame(self.nb)
        self.nb.add(f, text="Galaxy list")
        cols = ("x", "y", "z", "rmin", "ralign", "rmax", "flags")
        self.gtree = ttk.Treeview(f, columns=cols, show="headings", height=12)
        for c in cols:
            self.gtree.heading(c, text=c)
            self.gtree.column(c, width=110 if c != "flags" else 220, anchor="center")
        self.gtree.grid(row=0, column=0, sticky="nsew", padx=6, pady=6)
        self.gtree.bind("<Double-1>", lambda e: self._edit_galaxy())
        f.columnconfigure(0, weight=1)
        b = ttk.Frame(f)
        b.grid(row=0, column=1, sticky="n", pady=6)
        for txt, cmd in (("Add…", self._add_galaxy), ("Edit…", self._edit_galaxy),
                         ("Remove", self._remove_galaxy)):
            ttk.Button(b, text=txt, command=cmd, width=10).pack(pady=2)
        ttk.Label(f, foreground="#666", wraplength=700, justify="left",
                  text="Explicit galaxies (ListOfGals block), used only when neither ReadFofHalos "
                       "nor ReadSubFindGals is set. Positions and radii in kpc."
                  ).grid(row=1, column=0, columnspan=2, sticky="w", padx=6)

    def _refresh_galaxies(self):
        self.gtree.delete(*self.gtree.get_children())
        for g in self.galaxies:
            flags = " ".join(n for n, on in (("LOG", g.log), ("COMP_CM", g.compute_cm),
                                             ("POT", g.pot), ("TRACK_CM", g.track_cm),
                                             ("ALIGN", g.align)) if on)
            self.gtree.insert("", "end", values=(g.x, g.y, g.z, g.rmin, g.ralign, g.rmax, flags))

    def _galaxy_dialog(self, g=None):
        d = tk.Toplevel(self)
        d.title("Galaxy")
        d.transient(self.winfo_toplevel())
        fields = ["x", "y", "z", "rmin", "ralign", "rmax"]
        vs = {k: tk.StringVar(value=str(getattr(g, k)) if g else "") for k in fields}
        for i, k in enumerate(fields):
            ttk.Label(d, text=k + " (kpc)").grid(row=i, column=0, sticky="w", padx=6, pady=2)
            ttk.Entry(d, textvariable=vs[k], width=14).grid(row=i, column=1, padx=6)
        flags = {"log": "LOG", "compute_cm": "COMP_CM", "pot": "POT", "track_cm": "TRACK_CM",
                 "align": "ALIGN"}
        fv = {k: tk.BooleanVar(value=getattr(g, k) if g else k in ("log", "track_cm", "align"))
              for k in flags}
        fr = ttk.Frame(d)
        fr.grid(row=6, column=0, columnspan=2, pady=4)
        for k, lab in flags.items():
            ttk.Checkbutton(fr, text=lab, variable=fv[k]).pack(side="left", padx=2)
        res = {}

        def ok():
            try:
                vals = [float(vs[k].get()) for k in fields]
            except ValueError:
                messagebox.showerror("Galaxy", "all positions and radii must be numbers", parent=d)
                return
            gl = GalaxyLine(*vals)
            for k in flags:
                setattr(gl, k, fv[k].get())
            res["g"] = gl
            d.destroy()
        b = ttk.Frame(d)
        b.grid(row=7, column=0, columnspan=2, pady=6)
        ttk.Button(b, text="OK", command=ok).pack(side="left", padx=4)
        ttk.Button(b, text="Cancel", command=d.destroy).pack(side="left", padx=4)
        d.grab_set()
        self.wait_window(d)
        return res.get("g")

    def _add_galaxy(self):
        g = self._galaxy_dialog()
        if g:
            self.galaxies.append(g)
            self._refresh_galaxies()

    def _edit_galaxy(self):
        i = self._sel_index(self.gtree)
        if i is None:
            return
        g = self._galaxy_dialog(self.galaxies[i])
        if g:
            self.galaxies[i] = g
            self._refresh_galaxies()

    def _remove_galaxy(self):
        i = self._sel_index(self.gtree)
        if i is not None:
            del self.galaxies[i]
            self._refresh_galaxies()

    # ---- params <-> widgets ----------------------------------------------------
    def to_params(self) -> Params:
        p = Params(path=self.path or "")
        for name, v in self.vars.items():
            t = TAGS[name]
            if t.kind == LOGICAL:
                if v.get():
                    p.values[name] = True
            else:
                s = v.get().strip()
                if s:
                    if t.kind == "int":
                        p.values[name] = int(float(s))
                    elif t.kind == "double":
                        try:
                            p.values[name] = float(s)
                        except ValueError:
                            p.values[name] = 0.0 if s.lower() == "rvir" else s
                    else:
                        p.values[name] = s
        p.quantities = list(self.quantities)
        p.galaxies = list(self.galaxies)
        return p

    def from_params(self, p: Params):
        for name, v in self.vars.items():
            if TAGS[name].kind == LOGICAL:
                v.set(bool(p.values.get(name, False)))
            else:
                val = p.values.get(name)
                v.set("" if val is None else (f"{val:g}" if isinstance(val, float) else str(val)))
        self.quantities = list(p.quantities)
        self.galaxies = list(p.galaxies)
        self._refresh_quantities()
        self._refresh_galaxies()

    # ---- file actions ------------------------------------------------------------
    def _set_path(self, path):
        self.path = path
        self.path_label.config(text=path or "(new parameter file)")
        self.app.on_params_path(path)

    def new(self):
        p = Params()
        for k, v in (("Htrue", True), ("CmType", -1), ("RotateGalaxyFrame", 0), ("Cut", "0"),
                     ("Pressure", 20000.0), ("DecompositionMethod", 0), ("NumberOfBins", 50),
                     ("SurfaceDensity", True), ("AnalyzeGalaxies", True),
                     ("AlignWithAngularMomentum", True), ("Snapshot", "snap_000"),
                     ("Outdir", "galpost-out")):
            p.values[k] = v
        p.quantities = [Quantity("MASS", "NULL", 4, True, False),
                        Quantity("MCLD", "NULL", 0, True, False),
                        Quantity("SFR", "NULL", 0, True, False)]
        self.from_params(p)
        self._set_path(None)

    def open(self, path=None):
        path = path or filedialog.askopenfilename(
            parent=self, title="Open parameter file",
            initialdir=self.app.prefs.get("last_param_dir", os.getcwd()))
        if not path:
            return
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                p = read_params(path, check=False)
            except Exception as e:
                messagebox.showerror("Open", f"Cannot read {path}:\n{e}", parent=self)
                return
        self.from_params(p)
        self._set_path(os.path.abspath(path))
        self.app.prefs["last_param_dir"] = os.path.dirname(os.path.abspath(path))
        self.app.status(f"opened {path}" + (f" ({len(w)} warnings: " + "; ".join(
            str(x.message) for x in w[:3]) + ")" if w else ""))

    def save(self):
        if not self.path:
            return self.save_as()
        try:
            text = self.to_params().to_text()
        except ValueError as e:
            messagebox.showerror("Save", str(e), parent=self)
            return False
        with open(self.path, "w") as f:
            f.write("# galpost parameter file (written by galpost-gui)\n" + text)
        self.app.status(f"saved {self.path}")
        return True

    def save_as(self):
        path = filedialog.asksaveasfilename(
            parent=self, title="Save parameter file",
            initialdir=self.app.prefs.get("last_param_dir", os.getcwd()),
            initialfile=os.path.basename(self.path) if self.path else "galpost_params.txt")
        if not path:
            return False
        self._set_path(os.path.abspath(path))
        self.app.prefs["last_param_dir"] = os.path.dirname(os.path.abspath(path))
        return self.save()

    def check(self, quiet=False) -> bool:
        try:
            p = self.to_params()
            p.check()
        except ValueError as e:
            messagebox.showerror("Check", str(e), parent=self)
            return False
        msgs = []
        if not p.quantities:
            msgs.append("no profile quantities: profiles.* will not be written")
        if not (p["ReadFofHalos"] or p["ReadSubFindGals"] or p.galaxies or p["Galfile"]):
            msgs.append("no galaxy selection: ListOfGals, Galfile, ReadFofHalos or ReadSubFindGals")
        if not quiet:
            messagebox.showinfo("Check", "Parameters are consistent." +
                                ("\n\nNotes:\n- " + "\n- ".join(msgs) if msgs else ""), parent=self)
        return True

    def show_text(self):
        try:
            text = self.to_params().to_text()
        except ValueError as e:
            text = f"error: {e}"
        d = tk.Toplevel(self)
        d.title("Parameter file")
        t = tk.Text(d, width=90, height=40, font="TkFixedFont")
        t.insert("1.0", text)
        t.config(state="disabled")
        t.pack(fill="both", expand=True)
