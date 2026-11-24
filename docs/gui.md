# Graphical interface

```bash
galpost-gui [parameter file] [--plots <output directory>]
```

(or `python -m galpost.gui`). A desktop application (Tkinter + matplotlib,
no browser) with three tabs. It needs a display: a local desktop, WSLg under
WSL2, or X forwarding (`ssh -X` / `ssh -Y`) from a cluster.

## 1 · Configure

Edit a parameter file. The fields are generated from the parameter table of the
code, grouped as in the [parameter reference](parameters.md): logical switches
are check boxes, the other parameters text fields (with suggestions for radii,
and a browse button for files and directories). Required parameters are in
red; hovering on a field shows its description and default.

* **Quantities (profiles)**: the `Quantities` block, one line per profile; Add /
  Edit (double-click) / Remove / Up / Down. Quantity and weight names are chosen
  from the known and derived blocks (element names are read from `MetalsFile`).
* **Galaxy list**: the `ListOfGals` block.
* **Check** verifies required parameters and consistency; **Show text** shows
  the parameter file that will be written; **Save** writes it in the
  ExtractGalaxies syntax (old `EGparams` files can be opened and saved).

## 2 · Run

* **Parameter file** and **working directory** (relative paths in the parameter
  file are relative to it; default: the directory of the parameter file).
  "Use the Configure tab" saves the current configuration and selects it.
* **Snapshots**: add files, or a range (base name, first, last, step); multi-file
  snapshots are added by their base name.
* **Here, in the background**: snapshots are analysed one after the other in a
  separate process, with a live log; **Stop** interrupts the run.
* **SLURM job array**: one array task per snapshot. Account, partition, time,
  memory, CPUs, maximum number of simultaneous tasks and the commands to run
  before galpost (modules, environment activation) are remembered between
  sessions. "Preview" / "Save SLURM script" show or write the script; **Run**
  writes `galpost_job.sh` in the working directory and submits it with `sbatch`.
* **Show results** opens the output directory in the Plots tab.
* **Progenitor tracking**: **Track** runs [`galpost-track`](tracking.md) on the
  listed snapshots (the latest one is the reference) with the chosen number of
  tracers, here in the background or as a single SLURM job; **Use track file in
  Configure** sets `TrackFile` to the track just produced, so that a following
  **Run** analyses the main progenitor in every listed snapshot (save the
  parameter file first).

## 3 · Plots

* Choose an output directory: its files are listed by kind (profiles, star
  formation history, every extra analysis, maps, ...).
* Selecting a file fills the column lists and selects a sensible default plot;
  choose any x column and one or more y columns, logarithmic axes, line / step /
  scatter style and a title, then **Plot** (or double-click the file).
* Select several files of the same kind (Ctrl/Shift-click) to overlay them, e.g.
  different snapshots, or galpost and ExtractGalaxies outputs copied in the same
  directory.
* Text files (`info`, `cgm`, `bh`, ...) are shown as text; maps as images in kpc.
* The matplotlib toolbar gives zoom, pan and axis settings; **PNG… / PDF…** save
  the figure at the chosen resolution.

Preferences (window size, last directories, SLURM settings) are stored in
`~/.config/galpost/gui.json`.

## Plots without the interface

The same plots can be produced in batch:

```bash
galpost-plot profiles.snap_134.gal0 --columns              # list the columns
galpost-plot kinematics.snap_134.gal0 --out kin.pdf        # default plot of the file kind
galpost-plot profiles.snap_1*.gal0 --x 1 --y 4 6 --logx --logy --out sigma.png
```

Columns are numbered from 1, as in the file headers.
