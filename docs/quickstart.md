# Quick start

This page walks through a complete analysis of one snapshot: halo finding,
galaxy analysis and a first look at the results.

## 1. Halo catalogue (optional)

galpost can select galaxies from

* a **FoF catalogue** (`groups_NNN/group_tab_NNN.*`), written by OpenGadget on
  the fly or by the standalone P-FoF code;
* a **SubFind catalogue** (`groups_NNN/sub_NNN.*`), written by OpenGadget in
  snapshot format;
* an explicit **list of positions** in the parameter file or in a file.

If the simulation did not produce SubFind catalogues, they can be computed in
post-processing with OpenGadget itself; see
[halo finding with OpenGadget](subfind_opengadget.md).

## 2. A parameter file

The parameter file uses the ExtractGalaxies syntax (see
[the parameter file](parameter_file.md)). A minimal one, selecting the most
massive FoF group and computing surface density profiles of stars and gas:

```
GadgetParamfile     ../paramfile.par     # units and periodicity of the run
MetalsFile          ../metals.dat
Outdir              galpost-out
Htrue
SurfaceDensity
NumberOfBins        50
RotateGalaxyFrame   0
CmType              -1                   # cold gas + stars define centre and frame
Cut                 0
Pressure            20000
DecompositionMethod 0
AlignWithAngularMomentum
AnalyzeGalaxies

ReadFofHalos
GroupDir            ../SubFindPP         # where groups_NNN is (default: snapshot directory)
FofMinMass          -1                   # the most massive group
FofRmin             0.5
FofRmax             rvir
FofRgal             0                    # rgal = rvir/10
FofRalign           8.0
FofRlog

Quantities
# Quantity Weight Type DENS/AVE DIFF/CUMU
MASS  NULL  4  DENS DIFF
MCLD  NULL  0  DENS DIFF
SFR   NULL  0  DENS DIFF
TEMP  MCLD  0  AVE  DIFF
End
```

For zoom-in simulations it is wise to add

```
MaxLowResFraction   0.01
```

which discards galaxies whose virial sphere contains more than 1% of mass in
low-resolution particles.

## 3. Run

```bash
galpost EGparams -snap ../snap_134
```

The snapshot given with `-snap` overrides the `Snapshot` tag. On a compute node
of a cluster, wrap the command in a batch script, for example with SLURM:

```bash
#!/bin/bash
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=32G --time=01:00:00
source ~/.venvs/galpost/bin/activate
for n in 120 125 130 134; do
    galpost EGparams -snap ../snap_$n
done
```

## 4. Results

`galpost-out/` now contains, for galaxy 0 of snapshot 134:

| file | content |
|---|---|
| `info.snap_134.gal0` | centre, frame, masses, spins, inertia tensors, B/T |
| `profiles.snap_134.gal0` | the requested profiles |
| `jhist.snap_134.gal0` | histogram of stellar circularities |
| `sfr.snap_134.gal0` | star formation history of the stars within rgal |
| `rotation_matrix.snap_134.gal0` | box → galaxy frame rotation |
| `list_of_selected_gals.snap_134` | one line per galaxy |
| `list_of_blocks.dat` | names of the available quantities |

All files are plain text with a commented header; see [output files](outputs.md).
Reading a profile in Python:

```python
import numpy as np
p = np.loadtxt("galpost-out/profiles.snap_134.gal0")
r, sigma_star, n_star = p[:, 0], p[:, 3], p[:, 4]
```

## 5. More analyses

Add any of the extra analyses to the parameter file, for instance

```
KinematicDecomposition
StellarPopulations
DiscStructure
KSRelation
ScalingRelations
Maps
```

and run again. Each one writes its own file; see
[extra analyses](extra_analyses.md).

## 6. Galaxy snapshots

```
WriteGalaxySnapshot
GalaxySnapshotRadius  rgal
AlignGalaxySnapshot
GalaxySnapshotInCodeUnits
```

writes `zoom.snap_134`, a GADGET format-2 snapshot with all particles within
rgal, centred on the galaxy and rotated in its frame, containing every block of
the original snapshot.
