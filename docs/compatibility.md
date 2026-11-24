# Compatibility with ExtractGalaxies

galpost re-implements ExtractGalaxies 0.4.3 (P. Monaco, G. Murante, L.
Tornatore and collaborators). The goal is that the same parameter file produces
the same files, so that existing analysis scripts keep working.

## What is reproduced

* parameter file syntax and all tags (see [parameters](parameters.md));
* galaxy selection from P-FoF, OpenGadget FoF, SubFind and galaxy lists, with
  the same radii (including the ExtractGalaxies virial radius formula);
* unit conversions and physical constants (`libreadgad` values);
* MUPPI phase re-assignment: multi-phase particles with a hot-phase temperature
  below 2000 K are treated as single phase; single-phase particles are cold
  below 10⁵ K, otherwise hot; the SFR of single-phase particles is set to zero;
* centre tracking, alignment, the six analysis sets within rgal and rvir,
  inertia tensors, spins, half-mass radii, jz/jcirc histograms, the four
  bulge/disc decompositions, star formation histories;
* all derived blocks and the three profile geometries, with the same binning
  and normalisation;
* numerical precision: particle data are kept in single precision exactly as in
  ExtractGalaxies and sums are done in the same order, so that results agree to
  the last printed digit;
* output files and formats, including the particle catalogues and the
  `table-*` time series.

## Validation

galpost was compared with the ExtractGalaxies outputs of an Aquila C5 MUPPI
simulation (8 snapshots from z = 9.5 to z = 0, two configurations: surface
profiles with `AnalyzeGalaxies`, and cumulative volume profiles), about 10⁵
numbers in total, comparing every number of corresponding lines of the two
sets of output files.

* 7 snapshots out of 8: all files agree to the last printed digit (differences
  ≤ 10⁻⁶, from rounding of the last digit);
* one snapshot (z = 0): the tracked centre differs by 1.5 × 10⁻⁵ kpc (2 × 10⁻¹⁰
  relative). This moves single particles across bin edges (profile differences
  up to 10⁻³ in a few bins, a few low-count bins of `jhist`) and changes the
  small DM spin component within rgal by a few per cent through one boundary
  particle. Imposing the ExtractGalaxies centre, galpost reproduces those
  outputs exactly; the origin of the tiny difference in the centre was not
  identified;
* galaxy snapshots: same particles and IDs, positions, velocities and masses
  equal within single precision.

## Bugs of ExtractGalaxies fixed in galpost

| issue | ExtractGalaxies | galpost |
|---|---|---|
| `OUTF`/`INFL` of particles with zero flux | keep the value of the previous particle | 0 |
| `GalaxySnapshotRadius rmax` with logarithmic bins | uses log₁₀(rmax) as radius | rmax |
| `SubUsePartAssignment` | compares the subhalo index of each particle with the position of the galaxy in the list (correct only if they coincide) | uses the galaxy's own subhalo |
| `.muppi` companion blocks in galaxy snapshots | written as zeros (the companion file was never read) | real values |
| galaxy snapshots | fixed list of blocks, no INFO block | all blocks, with INFO block |
| `WriteGalaxyCatalog` file names | documentation says `bl` = type 2, `ds` = type 3 | as the ExtractGalaxies code: `ds` = type 2, `bl` = type 3 |

## Intentional differences

* Galaxy snapshots contain the original snapshot values of `MHOT`, `MCLD`,
  `MSF`, `EHOT`, `SFR` and `MFST`; ExtractGalaxies wrote the values modified by
  the phase re-assignment.
* Only the particles in a sphere around the galaxies are read (large enough to
  contain all analysis volumes), instead of the whole snapshot; this is
  transparent for the results.
* Any snapshot block can be used as a profile quantity (in code units).
* The inertia-tensor frame (without `AlignWithAngularMomentum`) can differ from
  ExtractGalaxies by axis flips: eigenvector signs depend on the linear-algebra
  library (GSL vs LAPACK). Profiles are unaffected; the angular-momentum frame
  is identical.

## ExtractGalaxies behaviours that are reproduced on purpose

Kept for compatibility, so that old and new outputs can be compared:

* the virial radius uses the exponent 0.333333 instead of 1/3 (relative effect
  ≈ 5 × 10⁻⁶);
* the half-mass radius is the radius of the particle *after* the one that
  crosses half of the mass;
* particle-to-galaxy assignments accumulate across galaxies, so that during
  centre tracking the particles already assigned to previous galaxies are
  excluded;
* in `info.*` the angular momentum in the galaxy frame is labelled "box frame";
* without `AnalyzeGalaxies`, the cold gas mass in `list_of_selected_gals` uses
  the criterion "Mhot > 0 or T < 10⁵ K".
