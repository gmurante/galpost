# Output files

All files are plain text (except snapshots, `.npy` arrays and `.png` maps) and
start with a commented header. File names follow the ExtractGalaxies scheme:

```
<Outdir>/<prefix>.<Fileflag>.<snapshot>.gal<ID>      per galaxy and snapshot
<Outdir>/<prefix>.<Fileflag>.<snapshot>              per snapshot
<Outdir>/table-<name>.<Fileflag>.gal<ID>             one line appended per snapshot
```

`<Fileflag>.` is omitted when `Fileflag` is empty. `<ID>` is the galaxy ID: 0
when a single FoF group is selected, the group index for several FoF groups,
the subhalo index for SubFind selections, the line number for galaxy lists.

Units: kpc, km/s, Msun, Gyr, Msun/yr; h=1 units (kpc/h, Msun/h, ...) unless
`Htrue` is set. "Galaxy frame" means centred on the galaxy, at rest with it and
rotated so that z is the disc axis.

The first part of this page describes the files of ExtractGalaxies, which
galpost writes in exactly the same format; the second part the files of the
[extra analyses](extra_analyses.md).

---

## ExtractGalaxies files

### `profiles.*.gal<ID>` (if a `Quantities` block is given)

Header: number of bins, list of quantities, galaxy number. One line per bin:

| column | content |
|---|---|
| 1 | bin centre (for log bins: 10^(centre of log bin)); grid: cell size |
| 2, 3 | inner and outer edge; grid: x and y of the cell centre |
| 4, 6, 8, ... | profile of quantity 1, 2, 3, ... |
| 5, 7, 9, ... | number of particles of the corresponding type in the bin |

For grids the N×N cells are listed with x varying fastest.

### `info.*.gal<ID>` (with `AnalyzeGalaxies`)

Human-readable summary:

* rgal, rvir, profile range;
* centre of mass (kpc, box frame) and velocity (km/s) of the `CmType` particles
  within ralign;
* their inertia tensor (10¹⁰ Msun kpc²), eigenvalues, specific angular momentum
  (box frame), rotation matrix, angular momentum in the galaxy frame (labelled
  "box frame" as in ExtractGalaxies), angle between angular momentum and the
  major inertia axis;
* within rgal and within rvir, for six particle sets — gas, DM, cold gas, hot
  gas, stars, cold gas + stars: number of particles, mass, inertia eigenvalues,
  specific angular momentum in box and galaxy frame, its modulus and angle with
  the galaxy z axis. Cold and hot gas masses are the MUPPI phase masses;
* baryonic mass and baryon fraction within rvir, fraction of mass in the galaxy;
* specific angular momentum of the halo (rgal < r < rvir) DM, stars and hot gas;
* bulge, disc and stellar halo masses and B/T of the decomposition chosen with
  `DecompositionMethod`.

### `jhist.*.gal<ID>` (with `AnalyzeGalaxies`)

Histogram of stellar circularities ε = jz / jcirc(r), with
jcirc = r √(G M(<r)/r), for the stars within rgal.

| column | content |
|---|---|
| 1 | ε (50 bins between −2 and 2) |
| 2 | stellar mass in the bin, all radii within rgal (Msun) |
| 3–52 | stellar mass in 50 radial bins between 0 and rgal (centres in the header) |

### `sfr.*.gal<ID>` (with `AnalyzeGalaxies` and an `AGE` block)

Star formation history of the stars that are within rgal at the time of the
snapshot: column 1 cosmic time (Gyr), column 2 SFR (Msun/yr), 500 bins from the
Big Bang to z=0. Initial stellar masses (`iM`) are used when available.

### `rotation_matrix.*.gal<ID>` (with `AnalyzeGalaxies`)

The 3×3 matrix R such that x_galaxy = R (x_box − x_centre).

### `list_of_selected_gals.*` (without `WriteTables`)

One line per galaxy:

| columns | content |
|---|---|
| 1 | galaxy ID |
| 2–4 | snapshot number, scale factor, redshift |
| 5–7 | centre (kpc, box frame) |
| 8–10 | velocity (km/s) |
| 11–15 | stellar, gas, cold gas mass within rgal, FoF mass (Msun), SFR (Msun/yr) |
| 16–19 | stellar bulge, disc, halo mass within rgal, B/T |
| 20–24 | rmin, rmax, rgal, rvir, ralign (kpc) |
| 25–27 | inertia eigenvalues (10¹⁰ Msun kpc²) |
| 28–30 | specific angular momentum (kpc km/s) |
| 31 | angle between angular momentum and inertia tensor (deg) |
| 32–34 | centre (code units) |
| 35–37 | velocity (code units) |
| 38–44 | flags: Central, LogFlag, swap, ComputeCM, TrackCM, Align, Pot |

Without `AnalyzeGalaxies` only the masses within rgal are computed (bulge, disc,
B/T and SFR are 0, as in ExtractGalaxies).

### `table-*.gal<ID>` (with `WriteTables`)

One line is appended for each analysed snapshot:

* `table-info`: as `list_of_selected_gals` (without SFR, B/T and code units);
* `table-mass`: masses of the six sets within rgal and rvir, bulge/disc masses
  and B/T, half-mass radii of the six sets within rgal and rvir;
* `table-specspin`, `table-spin`: specific and total angular momenta (components
  and modulus) of the twelve sets;
* `table-inertia`: inertia tensors (6 components), eigenvalues and eigenvectors
  of the twelve sets.

### `to_find_gals.*` / `table-to_find_gals.gal<ID>` (with `ToFindGalaxy`)

Centre and velocity in code units and the rotation matrix, to locate the galaxy
in other codes.

### `list_of_blocks.dat`

Names of the known and derived quantities available for profiles.

### Particle catalogues

With `WriteGalaxyCatalog` (full catalogues) or `WriteParticleList` (ID, type
and box position only), for the particles of each galaxy used for the profiles:

| file | particles |
|---|---|
| `catalog_mp.*`, `catalog_sp.*` | multi-phase and single-phase gas |
| `catalog_st.*` | stars |
| `catalog_dm.*` | type 1 |
| `catalog_ds.*` | type 2 |
| `catalog_bl.*` | type 3 |

(`particle_list_*` for `WriteParticleList`). The columns are described in each
header: positions and velocities in the galaxy and box frames, mass, ID,
cylindrical radius; for gas also density, circularity, temperatures, internal
energy, MUPPI phase masses and energies, molecular mass, pressures, filling
factor, phase densities, SFR, outflow/inflow rates, electron fraction, density
gradient, number of multi-phase steps and element masses; for stars age,
circularity, component flag of the bulge/disc decomposition, element masses
and initial mass. `particle_list_*` files can be used as input of
`LoadParticleList`.

`WriteFullCatalog` / `WriteStarCatalog` write `fullcatalog_<type>.*` with all
particles of the snapshot (or all stars), with positions and velocities in the
box frame and in the frames of the first five galaxies.

### Galaxy snapshots: `zoom.*` (with `WriteGalaxySnapshot`)

`zoom.<Fileflag>.<snapshot>` (one galaxy) or
`zoom.<Fileflag>.gal<ID>.<snapshot>` (several galaxies): GADGET format-2
snapshots with an INFO block, containing **all blocks of the original snapshot**
(companion `.muppi` blocks included) for the particles within
`GalaxySnapshotRadius`, with the geometry of the profiles.

* header: copy of the original one with the particle numbers of the galaxy, all
  masses in the `MASS` block (mass table set to zero), `num_files = 1`;
* `AlignGalaxySnapshot`: `POS` centred and rotated, `VEL` relative to the galaxy
  and rotated, `GRAD` and `ACCE` rotated, `VBLK` treated as a velocity;
* without `GalaxySnapshotInCodeUnits`: `POS` (kpc), `VEL` (km/s), `MASS`, `MHOT`,
  `MCLD`, `MSF` (Msun), `RHO` (Msun pc⁻³), `EHOT` (10⁵⁰ erg) and `GRAD` are in
  physical units; all other blocks are always copied in code units.

---

## Files of the extra analyses

All tables have the format

```
# free-text header lines
#  1: description of column 1
#  2: description of column 2
...
<rows>
```

and can be read with `numpy.loadtxt`. Radial profiles are in annuli of the disc
plane (`ExtraBins` bins from 0 to `ExtraRmax`, linear or with `ExtraLog`) and
use particles within `ExtraZmax` of the plane, unless stated otherwise.

| file | tag | content |
|---|---|---|
| `kindec.*` | `KinematicDecomposition` | stellar mass within rgal, disc mass and D/T with ε(E) > `DiscCircularity`, spheroid = 2 × counter-rotating mass and B/T, thin (ε > 0.8) and thick (0.5 < ε ≤ 0.8) disc masses, mean ε |
| `jhistE.*` | `KinematicDecomposition` | histogram of ε(E) = jz/jcirc(E): 80 bins in [−2, 2], total and in 20 radial bins within rgal |
| `kinematics.*` | `KinematicProfiles` (and `KappaRot` in the header) | mean v_φ, σ_R, σ_φ, σ_z, V/σ and number of particles, for stars and cold gas |
| `stellarpop.*` | `StellarPopulations` | stellar mass, mass-weighted age, Z, [Fe/H], [O/Fe], [Mg/Fe] per annulus |
| `mdf.*` | `StellarPopulations` | [Fe/H] distribution of the stars within rgal |
| `ofe_feh.*` | `StellarPopulations` | stellar mass in the [O/Fe]–[Fe/H] plane (matrix; bin ranges in the header) |
| `disc.*` | `DiscStructure` | stellar surface density, median and rms height, exponential scale height per annulus; header: projected and 3D half-mass radii, exponential scale length and central surface density |
| `toomre.*` | `ToomreQ` | circular velocity, epicyclic frequency, surface densities and dispersions, Q of gas, stars and two-component |
| `bar.*` | `BarStrength` | A2(R) and phase of the m=2 mode; header: A2max, its radius, bar length, A2 within the half-mass radius |
| `flows.*` | `SphericalFlows` | per shell: outflow, inflow and net rates of gas, cold gas, hot gas, stars, DM and gas metals; mass loading |
| `ks.*` | `KSRelation` | Σ of total, cold, H₂ and atomic gas, Σ_SFR, depletion times |
| `gradients.*` | `MetallicityGradients` | 12+log(O/H) of cold gas (cold-mass and SFR weighted), stellar [Fe/H]; header: fitted slopes |
| `dust.*` | `DustAnalysis` | Σ of large and small grains, dust-to-gas and dust-to-metal ratios; header: totals within rgal |
| `cgm.*` | `HotHaloAnalysis` | one line: gas, hot and cold masses, hot temperature and metallicity of the CGM, baryon fraction |
| `bh.*` | `CentralBH` | black holes within rgal: ID, mass, accretion rate, distance from the centre; header: M_BH/M_bulge |
| `satellites.*` | `Satellites` | SubFind subhaloes within rvir: index, group, central flag, distance, stellar/gas/DM/low-resolution/total mass, Vmax, half-mass radius |
| `scaling.*` | `ScalingRelations` | one line per galaxy with the global quantities listed below |
| `track.gal<ID>` | `galpost-track` | main-progenitor track, see [tracking](tracking.md) |
| `table-scaling.gal<ID>` | `ScalingRelations` + `WriteTables` | the same line appended for each snapshot |
| `map_<comp>_<proj>.*.png/.npy/.extent` | `Maps` | surface density maps (Msun pc⁻²); `.npy` arrays have rows = y (face-on) or z (edge-on), columns = x; `.extent` gives xmin xmax ymin ymax in kpc |

### Columns of `scaling.*`

galaxy ID, snapshot number, redshift; within rgal: stellar, gas, cold gas and
H₂ mass, SFR, specific SFR; total mass within rvir, rvir, rgal; projected and
3D stellar half-mass radius, disc scale length, circular velocity at 2.2 R_d,
maximum circular velocity; B/T of the ExtractGalaxies decomposition, D/T and
B/T from ε(E); κ_rot, κ_co; mass-weighted stellar Z, [Fe/H] and age;
12+log(O/H) of cold gas at R_e and its gradient; mass and accretion rate of the
most massive black hole within rgal; A2max; hot CGM mass; baryon fraction
within rvir; low-resolution mass fraction within rvir. Quantities of analyses
that are not enabled are written as NaN.
