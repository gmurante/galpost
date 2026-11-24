# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""Main window of the galpost GUI."""

from __future__ import annotations

import os
import sys
import tkinter as tk
from tkinter import messagebox, ttk

from .. import __version__
from .config_tab import ConfigTab
from .plot_tab import PlotTab
from .run_tab import RunTab
from .widgets import Prefs


class App:
    def __init__(self, root: tk.Tk, paramfile=None, plotdir=None):
        self.root = root
        self.prefs = Prefs()
        root.title(f"galpost {__version__}")
        root.geometry(self.prefs.get("geometry", "1400x880"))
        root.protocol("WM_DELETE_WINDOW", self.quit)
        try:
            ttk.Style().theme_use("clam")
        except tk.TclError:
            pass
        self._menu()
        self.nb = ttk.Notebook(root)
        self.nb.pack(fill="both", expand=True)
        self.statusbar = ttk.Label(root, text="", anchor="w", relief="sunken")
        self.statusbar.pack(fill="x", side="bottom")
        self.run_tab = None
        self.config_tab = ConfigTab(self.nb, self)
        self.run_tab = RunTab(self.nb, self)
        self.plot_tab = PlotTab(self.nb, self)
        self.nb.add(self.config_tab, text="  1 · Configure  ")
        self.nb.add(self.run_tab, text="  2 · Run  ")
        self.nb.add(self.plot_tab, text="  3 · Plots  ")
        if paramfile:
            self.config_tab.open(paramfile)
        if plotdir:
            self.show_plots(plotdir)
        elif self.prefs.get("plot_dir"):
            self.plot_tab.refresh()
        self.status("ready")

    def _menu(self):
        m = tk.Menu(self.root)
        f = tk.Menu(m, tearoff=False)
        f.add_command(label="New parameter file", command=lambda: self.config_tab.new())
        f.add_command(label="Open parameter file…", command=lambda: self.config_tab.open())
        f.add_command(label="Save", command=lambda: self.config_tab.save())
        f.add_command(label="Save as…", command=lambda: self.config_tab.save_as())
        f.add_separator()
        f.add_command(label="Quit", command=self.quit)
        m.add_cascade(label="File", menu=f)
        h = tk.Menu(m, tearoff=False)
        h.add_command(label="Documentation", command=self._docs)
        h.add_command(label="About", command=self._about)
        m.add_cascade(label="Help", menu=h)
        self.root.config(menu=m)

    def status(self, msg):
        self.statusbar.config(text=" " + msg)

    def on_params_path(self, path):
        if self.run_tab is not None and path:
            self.run_tab.set_param_file(path)

    def show_plots(self, directory):
        self.plot_tab.set_dir(directory)
        self.nb.select(self.plot_tab)

    def _docs(self):
        root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))))
        doc = os.path.join(root, "docs")
        messagebox.showinfo("Documentation",
                            f"The documentation is in\n{doc}\n\n(index.md, parameters.md, "
                            "outputs.md, ...)" if os.path.isdir(doc) else
                            "See the docs/ directory of the galpost repository.")

    def _about(self):
        messagebox.showinfo("About galpost",
                            f"galpost {__version__}\nGalaxy post-processing for OpenGadget/MUPPI "
                            "simulations\n\nGiuseppe Murante (INAF-OATo), Pierluigi Monaco "
                            "(Università di Trieste)\nCopyright (C) 2026 INAF — GPL-3.0-or-later")

    def quit(self):
        self.prefs["geometry"] = self.root.geometry()
        self.prefs.save()
        self.root.destroy()


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog="galpost-gui", description="galpost graphical interface")
    ap.add_argument("paramfile", nargs="?", help="parameter file to open")
    ap.add_argument("--plots", default=None, help="output directory to show in the Plots tab")
    a = ap.parse_args(argv)
    try:
        root = tk.Tk()
    except tk.TclError as e:
        print(f"galpost-gui: cannot open a window ({e}).\n"
              "Is a display available? (ssh -X / -Y, WSLg, or a local X server)", file=sys.stderr)
        return 1
    App(root, a.paramfile, a.plots)
    root.mainloop()
    return 0
