# galpost - galaxy post-processing for OpenGadget/MUPPI simulations
# Copyright (C) 2026 INAF - Istituto Nazionale di Astrofisica and
#                    Università degli Studi di Trieste
# Authors: Giuseppe Murante (INAF - Osservatorio Astrofisico di Torino),
#          Pierluigi Monaco (Università degli Studi di Trieste)
# SPDX-License-Identifier: GPL-3.0-or-later

"""'Run' tab: run the analysis locally (in the background, with a live log)
or through a SLURM job array (one task per snapshot)."""

from __future__ import annotations

import glob
import os
import queue
import shlex
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from .widgets import Tooltip

SLURM_FIELDS = [("account", "Account", ""), ("partition", "Partition", ""),
                ("time", "Time limit", "02:00:00"), ("mem", "Memory", "32G"),
                ("cpus", "CPUs per task", "4"), ("max_parallel", "Max tasks at once", "")]


class RunTab(ttk.Frame):
    def __init__(self, master, app):
        super().__init__(master)
        self.app = app
        self.proc = None
        self.stop_flag = False
        self.q = queue.Queue()
        self.snaps = []
        self._build()
        self.after(150, self._poll)

    # ---- layout ------------------------------------------------------------
    def _build(self):
        top = ttk.Frame(self)
        top.pack(fill="x", padx=6, pady=4)
        ttk.Label(top, text="Parameter file:").grid(row=0, column=0, sticky="w")
        self.param_var = tk.StringVar()
        ttk.Entry(top, textvariable=self.param_var).grid(row=0, column=1, sticky="ew", padx=4)
        ttk.Button(top, text="…", width=3, command=self._browse_param).grid(row=0, column=2)
        ttk.Button(top, text="Use the Configure tab (save)", command=self._use_config).grid(
            row=0, column=3, padx=4)
        ttk.Label(top, text="Working directory:").grid(row=1, column=0, sticky="w")
        self.cwd_var = tk.StringVar()
        e = ttk.Entry(top, textvariable=self.cwd_var)
        e.grid(row=1, column=1, sticky="ew", padx=4)
        Tooltip(e, "Relative paths of the parameter file (GadgetParamfile, Outdir, ...) are "
                   "relative to this directory. Default: the directory of the parameter file.")
        ttk.Button(top, text="…", width=3, command=self._browse_cwd).grid(row=1, column=2)
        top.columnconfigure(1, weight=1)

        mid = ttk.Frame(self)
        mid.pack(fill="both", expand=True, padx=6, pady=4)
        # snapshots
        lf = ttk.LabelFrame(mid, text="Snapshots")
        lf.pack(side="left", fill="y", padx=(0, 6))
        self.snap_list = tk.Listbox(lf, width=46, height=18, selectmode="extended")
        self.snap_list.pack(fill="both", expand=True, padx=4, pady=4)
        b = ttk.Frame(lf)
        b.pack(fill="x", padx=4, pady=2)
        ttk.Button(b, text="Add files…", command=self._add_files).pack(side="left", padx=1)
        ttk.Button(b, text="Add range…", command=self._add_range).pack(side="left", padx=1)
        ttk.Button(b, text="Remove", command=self._remove).pack(side="left", padx=1)
        ttk.Button(b, text="Clear", command=self._clear).pack(side="left", padx=1)

        # how to run
        rf = ttk.LabelFrame(mid, text="Run")
        rf.pack(side="left", fill="both", expand=True)
        self.mode = tk.StringVar(value=self.app.prefs.get("run_mode", "local"))
        mf = ttk.Frame(rf)
        mf.pack(fill="x", padx=4, pady=4)
        ttk.Radiobutton(mf, text="Here, in the background", value="local",
                        variable=self.mode).pack(side="left", padx=4)
        ttk.Radiobutton(mf, text="SLURM job array", value="slurm",
                        variable=self.mode).pack(side="left", padx=4)
        sf = ttk.Frame(rf)
        sf.pack(fill="x", padx=4)
        self.slurm = {}
        saved = self.app.prefs.get("slurm", {})
        for i, (k, lab, dflt) in enumerate(SLURM_FIELDS):
            ttk.Label(sf, text=lab).grid(row=i // 3, column=2 * (i % 3), sticky="w", padx=3, pady=2)
            v = tk.StringVar(value=saved.get(k, dflt))
            ttk.Entry(sf, textvariable=v, width=14).grid(row=i // 3, column=2 * (i % 3) + 1, padx=3)
            self.slurm[k] = v
        ttk.Label(rf, text="Commands executed before galpost in the job (modules, environment):"
                  ).pack(anchor="w", padx=6)
        self.setup_text = tk.Text(rf, height=3, font="TkFixedFont")
        self.setup_text.insert("1.0", saved.get("setup", self._default_setup()))
        self.setup_text.pack(fill="x", padx=6)
        bf = ttk.Frame(rf)
        bf.pack(fill="x", padx=4, pady=6)
        self.run_btn = ttk.Button(bf, text="▶ Run", command=self.run)
        self.run_btn.pack(side="left", padx=2)
        self.stop_btn = ttk.Button(bf, text="■ Stop", command=self.stop, state="disabled")
        self.stop_btn.pack(side="left", padx=2)
        ttk.Button(bf, text="Preview SLURM script", command=self._preview).pack(side="left", padx=2)
        ttk.Button(bf, text="Save SLURM script…", command=self._save_script).pack(side="left", padx=2)
        ttk.Button(bf, text="Show results", command=self._show_results).pack(side="right", padx=2)
        self.progress = ttk.Label(rf, text="", foreground="#555")
        self.progress.pack(anchor="w", padx=6)

        tf = ttk.LabelFrame(rf, text="Progenitor tracking (galpost-track)")
        tf.pack(fill="x", padx=4, pady=6)
        ttk.Label(tf, text="Tracers").pack(side="left", padx=(6, 2))
        self.ntracers = tk.IntVar(value=int(self.app.prefs.get("tracers", 1000)))
        sp = ttk.Spinbox(tf, from_=100, to=100000, increment=100, textvariable=self.ntracers, width=7)
        sp.pack(side="left")
        Tooltip(sp, "Number of most-bound dark-matter particles used as tracers.")
        self.track_btn = ttk.Button(tf, text="▶ Track", command=self.track)
        self.track_btn.pack(side="left", padx=6, pady=4)
        Tooltip(self.track_btn, "Follow the galaxies selected at the latest snapshot of the list "
                                "back through all the listed snapshots (here or as one SLURM "
                                "job, according to the run mode). Writes track.gal<ID> in Outdir.")
        b2 = ttk.Button(tf, text="Use track file in Configure", command=self._use_track)
        b2.pack(side="left", padx=2)
        Tooltip(b2, "Set TrackFile in the Configure tab: galpost will then analyse the main "
                    "progenitor in every snapshot of the track.")
        self.last_track = None

        lg = ttk.LabelFrame(self, text="Log")
        lg.pack(fill="both", expand=True, padx=6, pady=4)
        self.log = tk.Text(lg, height=14, font="TkFixedFont", wrap="none")
        sb = ttk.Scrollbar(lg, command=self.log.yview)
        self.log.configure(yscrollcommand=sb.set)
        self.log.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

    def _default_setup(self):
        env = os.environ.get("VIRTUAL_ENV") or os.environ.get("CONDA_PREFIX")
        if os.environ.get("VIRTUAL_ENV"):
            return f"source {env}/bin/activate"
        if os.environ.get("CONDA_DEFAULT_ENV"):
            return f"# conda activate {os.environ['CONDA_DEFAULT_ENV']}"
        return f"export PYTHONPATH={os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))}:$PYTHONPATH"

    # ---- inputs ----------------------------------------------------------------
    def set_param_file(self, path):
        if path:
            self.param_var.set(path)
            if not self.cwd_var.get():
                self.cwd_var.set(os.path.dirname(path))

    def _use_config(self):
        cfg = self.app.config_tab
        if cfg.save():
            self.param_var.set(cfg.path)
            self.cwd_var.set(os.path.dirname(cfg.path))

    def _browse_param(self):
        p = filedialog.askopenfilename(parent=self, initialdir=self._cwd())
        if p:
            self.param_var.set(p)
            self.cwd_var.set(os.path.dirname(p))

    def _browse_cwd(self):
        p = filedialog.askdirectory(parent=self, initialdir=self._cwd())
        if p:
            self.cwd_var.set(p)

    def _cwd(self):
        return self.cwd_var.get() or os.path.dirname(self.param_var.get()) or os.getcwd()

    def _refresh(self):
        self.snap_list.delete(0, "end")
        for s in self.snaps:
            self.snap_list.insert("end", s)

    def _add_files(self):
        ps = filedialog.askopenfilenames(parent=self, initialdir=self.app.prefs.get(
            "last_snap_dir", self._cwd()), title="Snapshots")
        for p in ps:
            p = p[:-2] if p.endswith(".0") else p      # multi-file snapshot: base name
            if p not in self.snaps and not p.endswith(".muppi"):
                self.snaps.append(p)
        if ps:
            self.app.prefs["last_snap_dir"] = os.path.dirname(ps[0])
        self._refresh()

    def _add_range(self):
        d = tk.Toplevel(self)
        d.title("Add a range of snapshots")
        d.transient(self.winfo_toplevel())
        vb = tk.StringVar(value=os.path.join(self.app.prefs.get("last_snap_dir", self._cwd()),
                                             "snap_"))
        v1, v2, v3 = tk.StringVar(value="0"), tk.StringVar(value="0"), tk.StringVar(value="1")
        for i, (lab, v) in enumerate((("Base name (path + prefix)", vb), ("First", v1),
                                      ("Last", v2), ("Step", v3))):
            ttk.Label(d, text=lab).grid(row=i, column=0, sticky="w", padx=6, pady=2)
            ttk.Entry(d, textvariable=v, width=50 if i == 0 else 8).grid(row=i, column=1,
                                                                          sticky="w", padx=6)

        def ok():
            missing = 0
            for n in range(int(v1.get()), int(v2.get()) + 1, max(1, int(v3.get()))):
                p = f"{vb.get()}{n:03d}"
                if os.path.exists(p) or os.path.exists(p + ".0"):
                    if p not in self.snaps:
                        self.snaps.append(p)
                else:
                    missing += 1
            self._refresh()
            if missing:
                self.app.status(f"{missing} snapshots of the range do not exist and were skipped")
            d.destroy()
        ttk.Button(d, text="Add", command=ok).grid(row=4, column=0, columnspan=2, pady=6)

    def _remove(self):
        for i in reversed(self.snap_list.curselection()):
            del self.snaps[i]
        self._refresh()

    def _clear(self):
        self.snaps = []
        self._refresh()

    # ---- local run ----------------------------------------------------------
    def _check_inputs(self):
        pf = self.param_var.get().strip()
        if not pf or not os.path.exists(pf):
            messagebox.showerror("Run", "Choose (or save) a parameter file first.", parent=self)
            return None
        if not self.snaps:
            messagebox.showerror("Run", "Add at least one snapshot.", parent=self)
            return None
        return pf

    def run(self):
        pf = self._check_inputs()
        if not pf:
            return
        self._save_prefs()
        if self.mode.get() == "slurm":
            return self._submit(pf)
        self.stop_flag = False
        self.run_btn.config(state="disabled")
        self.track_btn.config(state="disabled")
        self.stop_btn.config(state="normal")
        self.log.delete("1.0", "end")
        snaps = list(self.snaps)
        threading.Thread(target=self._worker, args=(pf, self._cwd(), snaps), daemon=True).start()

    def _worker(self, pf, cwd, snaps):
        ok = 0
        for i, s in enumerate(snaps):
            if self.stop_flag:
                break
            self.q.put(("progress", f"snapshot {i + 1}/{len(snaps)}: {os.path.basename(s)}"))
            self.q.put(("log", f"\n$ galpost {pf} -snap {s}\n"))
            cmd = [sys.executable, "-u", "-m", "galpost.cli", pf, "-snap", s]
            try:
                self.proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT, text=True, bufsize=1,
                                             env=self._child_env())
                for line in self.proc.stdout:
                    self.q.put(("log", line))
                rc = self.proc.wait()
            except OSError as e:
                self.q.put(("log", f"cannot start galpost: {e}\n"))
                rc = -1
            self.q.put(("log", f"[exit status {rc}]\n"))
            ok += rc == 0
        self.proc = None
        self.q.put(("done", f"finished: {ok}/{len(snaps)} snapshots OK"
                    + (" (stopped)" if self.stop_flag else "")))

    def _child_env(self):
        env = dict(os.environ)
        src = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        env["PYTHONPATH"] = src + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        env["PYTHONWARNINGS"] = env.get("PYTHONWARNINGS", "ignore")
        return env

    # ---- progenitor tracking --------------------------------------------------
    def _outdir(self, pf):
        try:
            from ..params import read_params
            out = read_params(pf, check=False)["Outdir"] or "."
        except Exception:
            out = "."
        return out if os.path.isabs(out) else os.path.join(self._cwd(), out)

    def _track_cmd(self, pf, python="python3"):
        return [python, "-u", "-m", "galpost.track", pf, "--tracers", str(int(self.ntracers.get())),
                "--snapshots"] + list(self.snaps)

    def track(self):
        pf = self._check_inputs()
        if not pf:
            return
        if len(self.snaps) < 2:
            messagebox.showerror("Track", "Add at least two snapshots.", parent=self)
            return
        self.app.prefs["tracers"] = int(self.ntracers.get())
        self._save_prefs()
        if self.mode.get() == "slurm":
            return self._submit_track(pf)
        self.stop_flag = False
        for b in (self.run_btn, self.track_btn):
            b.config(state="disabled")
        self.stop_btn.config(state="normal")
        self.log.delete("1.0", "end")
        threading.Thread(target=self._track_worker, args=(pf, self._cwd()), daemon=True).start()

    def _track_worker(self, pf, cwd):
        cmd = self._track_cmd(pf, sys.executable)
        self.q.put(("progress", f"tracking through {len(self.snaps)} snapshots…"))
        self.q.put(("log", "$ galpost-track " + " ".join(shlex.quote(c) for c in cmd[4:]) + "\n"))
        rc = -1
        try:
            self.proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, text=True, bufsize=1,
                                         env=self._child_env())
            for line in self.proc.stdout:
                self.q.put(("log", line))
                if "written " in line and "track" in line:
                    self.last_track = line.split("written", 1)[1].strip()
            rc = self.proc.wait()
        except OSError as e:
            self.q.put(("log", f"cannot start galpost-track: {e}\n"))
        self.proc = None
        if self.last_track and not os.path.isabs(self.last_track):
            self.last_track = os.path.join(cwd, self.last_track)
        self.q.put(("done", f"tracking finished (exit status {rc})"
                    + (f": {self.last_track}" if rc == 0 and self.last_track else "")))

    def track_script(self, pf):
        s = {k: v.get().strip() for k, v in self.slurm.items()}
        lines = ["#!/bin/bash", "#SBATCH --job-name=galpost-track", "#SBATCH --nodes=1",
                 "#SBATCH --ntasks=1"]
        for key, opt in (("cpus", "cpus-per-task"), ("mem", "mem"), ("time", "time"),
                         ("account", "account"), ("partition", "partition")):
            if s[key]:
                lines.append(f"#SBATCH --{opt}={s[key]}")
        lines += ["#SBATCH --output=galpost_track_%j.log", "", "# written by galpost-gui",
                  self.setup_text.get("1.0", "end").strip(), "",
                  f"cd {shlex.quote(self._cwd())}",
                  " ".join(shlex.quote(c) for c in self._track_cmd(pf)), ""]
        return "\n".join(lines)

    def _submit_track(self, pf):
        path = os.path.join(self._cwd(), "galpost_track_job.sh")
        if not messagebox.askokcancel("Track", f"Write {path} and submit it with sbatch?",
                                      parent=self):
            return
        with open(path, "w") as f:
            f.write(self.track_script(pf))
        os.chmod(path, 0o755)
        try:
            r = subprocess.run(["sbatch", path], cwd=self._cwd(), capture_output=True, text=True)
            self.log.insert("end", (r.stdout or "") + (r.stderr or ""))
            self.app.status(r.stdout.strip() or r.stderr.strip())
        except OSError as e:
            messagebox.showerror("Track", f"sbatch not available: {e}", parent=self)

    def _use_track(self):
        path = self.last_track
        if not path or not os.path.exists(path):
            pf = self.param_var.get().strip()
            start = self._outdir(pf) if pf else self._cwd()
            path = filedialog.askopenfilename(parent=self, initialdir=start, title="Track file",
                                              filetypes=[("track files", "track*"), ("all", "*")])
        if not path:
            return
        cfg = self.app.config_tab
        cfg.vars["TrackFile"].set(cfg._relative(os.path.abspath(path)))
        self.app.status(f"TrackFile = {path} set in the Configure tab (save the parameter file "
                        "to use it)")
        self.app.nb.select(cfg)
        sec = [cfg.nb.tab(i, "text") for i in cfg.nb.tabs()].index("Galaxy selection")
        cfg.nb.select(sec)

    def stop(self):
        self.stop_flag = True
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
        self.progress.config(text="stopping…")

    def _poll(self):
        try:
            while True:
                kind, msg = self.q.get_nowait()
                if kind == "log":
                    self.log.insert("end", msg)
                    self.log.see("end")
                elif kind == "progress":
                    self.progress.config(text=msg)
                elif kind == "done":
                    self.progress.config(text=msg)
                    self.run_btn.config(state="normal")
                    self.track_btn.config(state="normal")
                    self.stop_btn.config(state="disabled")
                    self.app.status(msg)
        except queue.Empty:
            pass
        self.after(150, self._poll)

    # ---- SLURM ------------------------------------------------------------------
    def slurm_script(self, pf):
        s = {k: v.get().strip() for k, v in self.slurm.items()}
        n = len(self.snaps)
        lines = ["#!/bin/bash", "#SBATCH --job-name=galpost", "#SBATCH --nodes=1",
                 "#SBATCH --ntasks=1"]
        if s["cpus"]:
            lines.append(f"#SBATCH --cpus-per-task={s['cpus']}")
        if s["mem"]:
            lines.append(f"#SBATCH --mem={s['mem']}")
        if s["time"]:
            lines.append(f"#SBATCH --time={s['time']}")
        if s["account"]:
            lines.append(f"#SBATCH --account={s['account']}")
        if s["partition"]:
            lines.append(f"#SBATCH --partition={s['partition']}")
        arr = f"0-{n - 1}" + (f"%{s['max_parallel']}" if s["max_parallel"] else "")
        lines += [f"#SBATCH --array={arr}", "#SBATCH --output=galpost_%A_%a.log", "",
                  "# written by galpost-gui", self.setup_text.get("1.0", "end").strip(), "",
                  "SNAPS=("] + [f"  {shlex.quote(x)}" for x in self.snaps] + [")", "",
                  f"cd {shlex.quote(self._cwd())}",
                  f"python3 -m galpost.cli {shlex.quote(pf)} -snap \"${{SNAPS[$SLURM_ARRAY_TASK_ID]}}\"", ""]
        return "\n".join(lines)

    def _preview(self):
        pf = self._check_inputs()
        if not pf:
            return
        d = tk.Toplevel(self)
        d.title("SLURM script")
        t = tk.Text(d, width=100, height=30, font="TkFixedFont")
        t.insert("1.0", self.slurm_script(pf))
        t.pack(fill="both", expand=True)

    def _save_script(self, path=None):
        pf = self._check_inputs()
        if not pf:
            return None
        path = path or filedialog.asksaveasfilename(parent=self, initialdir=self._cwd(),
                                                    initialfile="galpost_job.sh")
        if not path:
            return None
        with open(path, "w") as f:
            f.write(self.slurm_script(pf))
        os.chmod(path, 0o755)
        self.app.status(f"written {path}")
        return path

    def _submit(self, pf):
        path = os.path.join(self._cwd(), "galpost_job.sh")
        if not messagebox.askokcancel(
                "Submit", f"Write {path} and submit it with sbatch\n"
                          f"({len(self.snaps)} array tasks)?", parent=self):
            return
        self._save_script(path)
        try:
            r = subprocess.run(["sbatch", path], cwd=self._cwd(), capture_output=True, text=True)
            self.log.insert("end", (r.stdout or "") + (r.stderr or ""))
            self.app.status(r.stdout.strip() or r.stderr.strip())
        except OSError as e:
            messagebox.showerror("Submit", f"sbatch not available: {e}", parent=self)

    def _save_prefs(self):
        d = {k: v.get() for k, v in self.slurm.items()}
        d["setup"] = self.setup_text.get("1.0", "end").strip()
        self.app.prefs["slurm"] = d
        self.app.prefs["run_mode"] = self.mode.get()
        self.app.prefs.save()

    # ---- results ----------------------------------------------------------------
    def _show_results(self):
        out = None
        try:
            from ..params import read_params
            p = read_params(self.param_var.get(), check=False)
            out = p["Outdir"] or "."
            if not os.path.isabs(out):
                out = os.path.join(self._cwd(), out)
        except Exception:
            out = self._cwd()
        self.app.show_plots(out)
