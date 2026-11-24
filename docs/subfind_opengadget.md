# Halo finding with OpenGadget

galpost reads halo catalogues; it does not compute them. FoF groups and
SubFind subhaloes can be obtained with OpenGadget itself, either during the run
or afterwards on existing snapshots ("post-processing"). This page describes
the post-processing route, which needs nothing but the snapshots and the
OpenGadget source used for the run.

## 1. A dedicated executable

Build a separate executable from the **same source and `Config.sh` used for the
simulation**, so that the snapshot layout (number of metals, MUPPI, dust, ...)
is identical, with these changes:

```
SUBFIND                       # substructure finder
MAX_NGB_CHECK=3
KEEP_HSML_AS_GUESS
WRITE_SUB_IN_SNAP_FORMAT      # catalogue in GADGET snapshot format (read by galpost)
SO_BAR_INFO                   # optional: baryonic SO properties of groups
SO_VEL_DISPERSIONS            # optional: SO velocity dispersions
#INPUT_IN_DOUBLEPRECISION     # remove it if the snapshots are written in single precision
```

`FOF` and its linking options must already be on. Keep the build in a copy of
the source tree (or use a different `EXEC` name) so that the production
executable is untouched.

Notes:

* OpenGadget stops with *"Code compiled with INPUT_IN_DOUBLEPRECISION, but
  input files are in single precision"* if this flag does not match the
  snapshots. Snapshots are usually written in single precision (unless
  `OUTPUT_IN_DOUBLEPRECISION`), even when the run used `INPUT_IN_DOUBLEPRECISION`
  for double-precision initial conditions.
* Without `HAVE_HDF5`, `WRITE_SUB_IN_SNAP_FORMAT` produces GADGET format-2
  files `groups_NNN/sub_NNN.*` with an INFO block. groups, subhaloes and the
  particle ID list appear as "particle types" 0, 1 and 2.

## 2. A dedicated parameter file

Copy the parameter file of the run and change

```
OutputDir     ./SubFindPP/     # separate directory: groups_NNN of the run are not overwritten
ICFormat      2                # snapshots are read with ICFormat, not SnapFormat
```

`OutputDir` must end with `/`. The snapshot is read from
`OutputDir/SnapshotFileBase_NNN`, so link it there:

```bash
mkdir -p SubFindPP
ln -s ../snap_134 SubFindPP/
```

The IMF and SF files of the chemical evolution (`IMFsFileName`, `SFsFileName`)
are also read from `OutputDir`: link them too.

```bash
ln -s ../Kroupa_1993.IMF ../SFs SubFindPP/
```

Do **not** link the `.muppi` companion file: its blocks are restart variables
of MUPPI that SubFind does not need (and see the first known issue below).

## 3. Run

```bash
mpirun -np 72 ./P-Gadget3-SubFind paramfile.subfind.par 3 134
```

`RestartFlag 3` reads snapshot 134, runs FoF and SubFind and writes
`SubFindPP/groups_134/` (`group_tab_134.*` and `sub_134.*`). The memory and
the number of MPI tasks needed are similar to those of the simulation at that
redshift.

A SLURM script processing a range of snapshots, skipping those already done:

```bash
#!/bin/bash
#SBATCH --nodes=12 --ntasks-per-node=6 --cpus-per-task=6 --time=06:00:00
# usage: sbatch go_subfind.sh <first> [<last>]
first=$1; last=${2:-$1}
EXE=./P-Gadget3-SubFind; PAR=paramfile.subfind.par; OUT=SubFindPP
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
mkdir -p $OUT
for f in Kroupa_1993.IMF SFs; do [ -e $OUT/$f ] || ln -s ../$f $OUT/$f; done
for ((n=10#$first; n<=10#$last; n++)); do
    num=$(printf "%03d" $n)
    [ -e $OUT/groups_$num/sub_$num.0 ] && continue
    [ -e $OUT/snap_$num ] || ln -s ../snap_$num $OUT/snap_$num
    mpirun -np $SLURM_NTASKS $EXE $PAR 3 $num > $OUT/log.subfind.$num 2>&1
    [ -e $OUT/groups_$num/sub_$num.0 ] || echo "snapshot $num FAILED, see $OUT/log.subfind.$num"
done
```

## 4. Use in galpost

```
ReadSubFindGals        # or ReadFofHalos: group_tab files of the same run
GroupDir   SubFindPP
```

For zoom-in simulations select galaxies with `MaxLowResFraction` (FoF groups
never contain low-resolution particles, but their virial spheres may), or rely
on the SubFind selection, which excludes subhaloes containing low-resolution
particles.

## Known issues

Two bugs were found in an OpenGadget version (revision r36112) while setting
this up; check whether your version is affected.

1. **`.muppi` companion files cannot be read** (`IO/read_ic.c`, function
   `read_file`): the companion file is closed right after its first block has
   been read, while the following auxiliary blocks still read from it, so the
   second one fails with *"I/O error (fread): Bad file descriptor"*. This
   affects every read of a MUPPI snapshot together with its `.muppi` file
   (restart from a snapshot, post-processing). Fix: close the companion file
   only once, together with the main file at the end of `read_file`.
   Workaround in post-processing: do not link the `.muppi` file.
2. **Wrong Vmax of large subhaloes** (`FofSubfind/subfind_serial.c`, subhalo
   properties, `collective == 1` branch): the MPI reduction of the radius of
   maximum circular velocity stores its result in `SubVmax` instead of
   `SubVmaxRad`. In affected catalogues `VMAX` of the largest subhaloes is a
   small number (a radius in code units) and `RMAX` is 0; smaller subhaloes,
   processed serially, are correct. The same wrong value is used for the
   velocity dispersion that SubFind assigns to black holes. Fix: reduce the
   radius into `SubVmaxRad`.

Other things to be aware of:

* with `RestartFlag 3` the snapshot is read with `ICFormat`; runs started from
  format-1 initial conditions keep `ICFormat 1` in their parameter file and
  fail with *"incorrect header format"* until it is changed;
* a `.muppi` file truncated by a job killed while writing the snapshot gives a
  generic *"end of file"* error; galpost detects truncated files and reports
  the incomplete block.
