# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""'Plots' tab: browse an output directory, plot files (overlaying several of
the same kind), save the figure as PNG or PDF."""

from __future__ import annotations

import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from ..plotting import file_kind, list_output_files, plot_files, read_table  # noqa: E402
from .widgets import Tooltip  # noqa: E402


class PlotTab(ttk.Frame):
    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.table = None
        self._build()

    # ---- layout ----------------------------------------------------------------
    def _build(self):
        pw = ttk.PanedWindow(self, orient="horizontal")
        pw.pack(fill="both", expand=True, padx=4, pady=4)

        left = ttk.Frame(pw)
        pw.add(left, weight=1)
        top = ttk.Frame(left)
        top.pack(fill="x")
        self.dir_var = tk.StringVar(value=self.app.prefs.get("plot_dir", ""))
        e = ttk.Entry(top, textvariable=self.dir_var)
        e.pack(side="left", fill="x", expand=True)
        e.bind("<Return>", lambda ev: self.refresh())
        ttk.Button(top, text="…", width=3, command=self._browse).pack(side="left")
        ttk.Button(top, text="⟳", width=3, command=self.refresh).pack(side="left")
        self.tree = ttk.Treeview(left, show="tree", selectmode="extended")
        ys = ttk.Scrollbar(left, command=self.tree.yview)
        self.tree.configure(yscrollcommand=ys.set)
        self.tree.pack(side="left", fill="both", expand=True, pady=4)
        ys.pack(side="left", fill="y", pady=4)
        self.tree.bind("<<TreeviewSelect>>", lambda e: self._on_select())
        self.tree.bind("<Double-1>", lambda e: self.plot())
        Tooltip(self.tree, "Select one file, or several files of the same kind to overlay them "
                           "(Ctrl/Shift-click). Double-click to plot.")

        mid = ttk.Frame(pw)
        pw.add(mid, weight=1)
        ttk.Label(mid, text="x axis").pack(anchor="w")
        self.x_var = tk.StringVar()
        self.x_box = ttk.Combobox(mid, textvariable=self.x_var, state="readonly", width=40)
        self.x_box.pack(fill="x")
        ttk.Label(mid, text="y axis (one or more)").pack(anchor="w", pady=(6, 0))
        yf = ttk.Frame(mid)
        yf.pack(fill="both", expand=True)
        self.y_list = tk.Listbox(yf, selectmode="extended", exportselection=False, height=12)
        ysb = ttk.Scrollbar(yf, command=self.y_list.yview)
        self.y_list.configure(yscrollcommand=ysb.set)
        self.y_list.pack(side="left", fill="both", expand=True)
        ysb.pack(side="left", fill="y")
        of = ttk.Frame(mid)
        of.pack(fill="x", pady=4)
        self.logx, self.logy = tk.BooleanVar(), tk.BooleanVar()
        ttk.Checkbutton(of, text="log x", variable=self.logx).pack(side="left")
        ttk.Checkbutton(of, text="log y", variable=self.logy).pack(side="left", padx=6)
        ttk.Label(of, text="style").pack(side="left")
        self.style = tk.StringVar(value="auto")
        ttk.Combobox(of, textvariable=self.style, values=["auto", "line", "step", "scatter"],
                     width=8, state="readonly").pack(side="left", padx=3)
        tf = ttk.Frame(mid)
        tf.pack(fill="x")
        ttk.Label(tf, text="title").pack(side="left")
        self.title_var = tk.StringVar()
        ttk.Entry(tf, textvariable=self.title_var).pack(side="left", fill="x", expand=True, padx=3)
        bf = ttk.Frame(mid)
        bf.pack(fill="x", pady=6)
        ttk.Button(bf, text="Plot", command=self.plot).pack(side="left", padx=2)
        ttk.Button(bf, text="Defaults", command=self._defaults).pack(side="left", padx=2)
        sf = ttk.LabelFrame(mid, text="Save figure")
        sf.pack(fill="x", pady=4)
        ttk.Button(sf, text="PNG…", command=lambda: self.save("png")).pack(side="left", padx=2, pady=4)
        ttk.Button(sf, text="PDF…", command=lambda: self.save("pdf")).pack(side="left", padx=2)
        ttk.Label(sf, text="dpi").pack(side="left", padx=(8, 2))
        self.dpi = tk.IntVar(value=150)
        ttk.Spinbox(sf, from_=50, to=600, increment=50, textvariable=self.dpi, width=5).pack(side="left")
        self.info = ttk.Label(mid, text="", foreground="#555", wraplength=300, justify="left")
        self.info.pack(anchor="w", pady=4)

        right = ttk.Frame(pw)
        pw.add(right, weight=4)
        self.fig = Figure(figsize=(7, 5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.toolbar = NavigationToolbar2Tk(self.canvas, right, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")
        self.canvas_w = self.canvas.get_tk_widget()
        self.canvas_w.pack(fill="both", expand=True)
        self.text = tk.Text(right, font="TkFixedFont", wrap="none")
        self.showing_text = False

    # ---- directory -------------------------------------------------------------
    def _browse(self):
        d = filedialog.askdirectory(parent=self, initialdir=self.dir_var.get() or os.getcwd())
        if d:
            self.set_dir(d)

    def set_dir(self, d):
        self.dir_var.set(d)
        self.app.prefs["plot_dir"] = d
        self.refresh()

    def refresh(self):
        self.tree.delete(*self.tree.get_children())
        d = self.dir_var.get()
        groups = list_output_files(d)
        if not groups:
            self.info.config(text=f"no galpost output files in {d}")
            return
        for desc in sorted(groups):
            node = self.tree.insert("", "end", text=f"{desc} ({len(groups[desc])})", open=False)
            for p in groups[desc]:
                self.tree.insert(node, "end", text=os.path.basename(p), values=(p,), iid=p)
        self.info.config(text=f"{sum(len(v) for v in groups.values())} files")

    def _selected_files(self):
        return [i for i in self.tree.selection() if os.path.isfile(i)]

    # ---- columns ------------------------------------------------------------------
    def _on_select(self):
        files = self._selected_files()
        if not files:
            return
        kinds = {getattr(file_kind(f), "name", None) for f in files}
        if len(kinds) > 1:
            self.info.config(text="selected files of different kinds: only the kind of the "
                                  "first one is used")
        try:
            self.table = read_table(files[0]) if not files[0].endswith(".npy") else None
        except Exception as e:
            self.table = None
            self.info.config(text=f"cannot read {os.path.basename(files[0])}: {e}")
        cols = self.table.columns if self.table is not None else []
        self.x_box.config(values=[f"{i + 1}: {c}" for i, c in enumerate(cols)])
        self.y_list.delete(0, "end")
        for i, c in enumerate(cols):
            self.y_list.insert("end", f"{i + 1}: {c}")
        self._defaults()

    def _defaults(self):
        files = self._selected_files()
        if not files:
            return
        k = file_kind(files[0])
        self.y_list.selection_clear(0, "end")
        if k is None or self.table is None or not self.table.columns:
            self.x_var.set("")
            return
        if k.x is not None and k.x < len(self.table.columns):
            self.x_box.current(k.x)
        for j in k.y:
            if j < len(self.table.columns):
                self.y_list.selection_set(j)
        self.logx.set(k.logx)
        self.logy.set(k.logy)
        self.style.set("auto")

    # ---- plot ---------------------------------------------------------------------
    def _show_text(self, text):
        if not self.showing_text:
            self.canvas_w.pack_forget()
            self.toolbar.pack_forget()
            self.text.pack(fill="both", expand=True)
            self.showing_text = True
        self.text.config(state="normal")
        self.text.delete("1.0", "end")
        self.text.insert("1.0", text)
        self.text.config(state="disabled")

    def _show_canvas(self):
        if self.showing_text:
            self.text.pack_forget()
            self.toolbar.pack(side="bottom", fill="x")
            self.canvas_w.pack(fill="both", expand=True)
            self.showing_text = False

    def plot(self):
        files = self._selected_files()
        if not files:
            return
        k = file_kind(files[0])
        x = self.x_box.current() if self.x_var.get() else None
        x = None if x is not None and x < 0 else x
        y = list(self.y_list.curselection()) or None
        is_text = k is not None and k.plot == "text" and x is None and not y
        try:
            res = plot_files(self.fig, files, x=x, y=y, logx=self.logx.get(), logy=self.logy.get(),
                             title=self.title_var.get() or None,
                             style=None if self.style.get() == "auto" else self.style.get())
        except Exception as e:
            messagebox.showerror("Plot", f"{type(e).__name__}: {e}", parent=self)
            return
        if isinstance(res, str) or is_text:
            self._show_text(res if isinstance(res, str) else "")
            return
        self._show_canvas()
        self.canvas.draw()

    def save(self, fmt):
        if self.showing_text:
            messagebox.showinfo("Save", "Nothing to save: the current view is a text file.",
                                parent=self)
            return
        files = self._selected_files()
        base = os.path.basename(files[0]) if files else "figure"
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension="." + fmt, initialfile=f"{base}.{fmt}",
            initialdir=self.app.prefs.get("save_dir", self.dir_var.get() or os.getcwd()),
            filetypes=[(fmt.upper(), "*." + fmt), ("all files", "*")])
        if not path:
            return
        self.fig.savefig(path, dpi=self.dpi.get(), bbox_inches="tight")
        self.app.prefs["save_dir"] = os.path.dirname(path)
        self.app.status(f"saved {path}")
