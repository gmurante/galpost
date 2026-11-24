# Main-progenitor tracking

```bash
galpost-track <parameter file> --snapshots ../snap_134 ../snap_120 ../snap_1?0 ... \
              [--tracers 1000] [--outdir DIR]
```

follows the selected galaxies back in time through any set of snapshots (not
necessarily consecutive), using particle IDs only: no FoF or SubFind
catalogues are needed except, if the parameter file selects galaxies from a
catalogue, at the latest snapshot.

## Method

1. **Reference.** The latest snapshot of the list is the reference: galaxies
   are selected there with the rules of the parameter file (FoF, SubFind,
   galaxy list or `TrackFile`).
2. **Tracers.** For each galaxy, the N (`--tracers`, default 1000) most bound
   dark-matter particles within rgal are chosen: lowest value of the `POT`
   block, or closest to the centre if the snapshot has no `POT` block. DM IDs do
   not change with time (unlike gas and star IDs).
3. **Going back.** In each earlier snapshot:
   * the tracers are found by ID;
   * the centre of the main progenitor is found with a shrinking-sphere
     iteration on the tracers (Power et al. 2003, MNRAS 338, 14), started from
     the previous centre: the sphere follows the densest group of tracers,
     i.e. the main branch, through mergers;
   * R200c and M200c (radius within which the mean density is 200 times the
     critical density at that redshift) are measured around the centre with
     all particles;
   * the fraction of the previous tracers that lies within rgal of the new
     centre is recorded;
   * the N most bound DM particles of the progenitor become the tracers for
     the next step, so the chain survives the dispersal of the original
     tracers.

## Output

`track.<Fileflag>.gal<ID>` in `Outdir` (or `--outdir`), one line per snapshot:

| columns | content |
|---|---|
| 1–3 | snapshot number, scale factor, redshift |
| 4–6 | centre (kpc, ExtractGalaxies units of that snapshot) |
| 7–9 | centre (code units) |
| 10, 11 | M200c (Msun), R200c (kpc) |
| 12–15 | rmin, rgal, ralign, rmax (kpc), with the rules of the selection mode of the parameter file (rgal = R200c/10 unless `FofRgal` is given; rmax = R200c unless `FofRmax`/`SubRmax` is given) |
| 16 | logarithmic bins (1/0) |
| 17 | fraction of the previous tracers within rgal of the new centre |
| 18 | number of tracers |

A fraction close to 1 means a smooth evolution; a low value flags a major
merger, a large gap between snapshots, or the very early assembly phase. The
header records whether `Htrue` was set: the track file must be used with the
same setting.

## Analysing the main progenitor

Add to a parameter file

```
TrackFile   galpost-out/track.gal0
```

and run galpost on any snapshot of the track: the galaxy is taken from the
track file (centre in code units, radii and log flag), with centre tracking and
alignment as for FoF galaxies. `TrackFile` overrides `ReadFofHalos`,
`ReadSubFindGals` and `ListOfGals`. With `WriteTables` and
`ScalingRelations` each run appends one line to the `table-*` files, building
the evolution of the galaxy.

Example with SLURM, one array task per snapshot (or use the Run tab of
`galpost-gui`):

```bash
#SBATCH --array=0-14
SNAPS=(134 130 120 110 100 090 080 070 060 050 040 030 020 015 010)
galpost params_track.txt -snap ../snap_${SNAPS[$SLURM_ARRAY_TASK_ID]}
```

## Validation

On an Aquila C5 MUPPI zoom (8 snapshots from z = 9.5 to 0), the tracked centre
agrees with the centre of the most massive FoF group found by ExtractGalaxies
within 0.6 kpc at z ≤ 3 and within 2 kpc up to z = 9.5, and R200c agrees with
the ExtractGalaxies virial radius within a few per cent. Starting galpost from
the track file, the refined centre coincides with the ExtractGalaxies one
within the 1 pc convergence tolerance of the centre tracking.

## Cost

Each snapshot requires reading the DM IDs and the particles around the
progenitor: about 1–2 s per snapshot for 8 × 10⁶ particles, and a few tens of
seconds for 6 × 10⁷ particles. Memory is dominated by the ID array of the DM
particles.
