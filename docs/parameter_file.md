# The parameter file

galpost reads the parameter files of ExtractGalaxies unchanged. All parameters
are listed, with defaults, in the [parameter reference](parameters.md); this
page describes the syntax and the concepts.

## Syntax

* `Tag value` for numbers and strings, `Tag` alone for logical switches
  (present = on, absent = off);
* lines starting with `#` or `%` are comments; text after `#` or `;` following
  a value is ignored;
* tags are case-insensitive (`GadgetParamfile` = `GadgetParamFile`);
* unknown tags produce a warning and are ignored; a tag given twice is an error;
* the order of the lines does not matter, except inside the blocks below.

Command line:

```
galpost <parameter file> [-snap <snapshot>]
```

## Units

All outputs are in physical units: kpc, km/s, Msun, Gyr (peculiar velocities;
positions divided by (1+z)). With `Htrue` the values use the Hubble parameter
of the snapshot header; without it they are in h=1 units (kpc/h, Msun/h).
Code units are read from `GadgetParamfile` (`UnitLength_in_cm`,
`UnitMass_in_g`, `UnitVelocity_in_cm_per_s`, `PeriodicBoundariesOn`,
`ComovingIntegrationOn`).

## Galaxy selection

Exactly one mode is used:

| mode | tags | galaxies |
|---|---|---|
| FoF | `ReadFofHalos`, `FofMinMass`, `FofRmin`, `FofRmax`, `FofRgal`, `FofRalign`, `FofRlog` | groups above a total mass, or the N most massive (`FofMinMass -N`) |
| SubFind | `ReadSubFindGals`, `SubMinStellarMass`, `SubRmin`, `SubRmax`, `SubRalign`, `SubRlog` | subhaloes above a stellar mass, or the N most massive in stars; subhaloes containing low-resolution particles are excluded |
| list | `ListOfGals ... End` block or `Galfile` | explicit positions |

Each galaxy has five radii:

| radius | meaning | FoF default | SubFind default |
|---|---|---|---|
| `rmin`, `rmax` | range of the profiles | `FofRmin`; `FofRmax` or rvir | `SubRmin`; `SubRmax` or R200c |
| `rvir` | halo radius | from the FoF mass (ExtractGalaxies formula) | R200c of the parent group |
| `rgal` | galaxy radius for global quantities | `FofRgal` or rvir/10 | rvir/10 |
| `ralign` | radius used to define centre and frame | `FofRalign` or rgal | `SubRalign` or rgal |

`MaxLowResFraction` (galpost) discards galaxies contaminated by low-resolution
particles (`LowResTypes`, default types 2 and 3) within rvir. FoF groups never
contain low-resolution particles, because those are not linked, so the check is
done on the particle distribution.

`GroupDir` (galpost) is the directory containing `groups_NNN`; by default the
directory of the snapshot.

### `ListOfGals`

```
ListOfGals
# x y z rmin ralign rmax [flags]
68913.8 68920.6 69201.5  0.5 8.0 235.5  LOG TRACK_CM ALIGN
End
```

Positions and radii are in kpc (in h units unless `Htrue`). `rgal` is set to
`ralign` and `rvir` to `rmax`. Flags:

* `LOG`: logarithmic radial bins;
* `COMP_CM`: centre = centre of mass of all `CmType` particles (one galaxy only,
  no periodic boundaries);
* `POT`: centre = most bound `CmType` particle within rmax (needs the `POT` block);
* `TRACK_CM`: centre refined iteratively with the centre of mass of `CmType`
  particles within `ralign`, until it moves by less than 1 pc;
* `ALIGN`: rotate to the galaxy frame.

The same lines can be put in a separate file given with `Galfile`.

## Centre and frame

`CmType` selects the particles defining centre, velocity and frame: a GADGET
type 0–5, or (any other value) cold gas plus stars; `CmType 6` means cold gas
only. Cold gas = multi-phase MUPPI particles or single-phase gas with T < 10⁵ K.

The frame is defined within `ralign`:

* default: principal axes of the inertia tensor (x, y, z = smallest, middle,
  largest eigenvalue), z oriented with the angular momentum;
* `AlignWithAngularMomentum`: z along the angular momentum, x along y_box × z;
* `AquilaAlignment`: as above, with `ralign` = 3 × the half-mass radius of
  `CmType` particles within rgal;
* `RotateGalaxyFrame 1` / `-1` cyclically rotates the axes.

## Profiles: the `Quantities` block

```
Quantities
# Quantity Weight Type DENS/AVE DIFF/CUMU
MASS  NULL  4  DENS CUMU
TEMP  MCLD  0  AVE  DIFF
End
```

Each line defines one profile:

* **Quantity** and **Weight**: block names (weight `NULL` = no weight);
* **Type**: GADGET particle type 0–5;
* **DENS**: density, sum(q·w)/sum(w) × N / volume; with weight NULL this is the
  surface (Msun pc⁻²) or volume (Msun pc⁻³) density of q. **AVE**: weighted
  average sum(q·w)/sum(w);
* **DIFF**: per bin; **CUMU**: cumulative within the outer edge of the bin.

The geometry is set by one of

* `VolumeDensity`: spherical shells from rmin to rmax;
* `SurfaceDensity`: cylindrical shells in the disc plane, height 2 × min(Cut, rmax);
* `GridDensity`: `NumberOfBins` × `NumberOfBins` square cells of side 2 rmax / N
  in the disc plane (no cumulative profiles).

`Cut` limits the vertical extent (a number in kpc or one of
rmin/rmax/rvir/rgal/ralign); `Cut 0` means no limit. It also sets the
height of the planes used for the OUTF/INFL mass fluxes.

### Available quantities

Resolution order, as in ExtractGalaxies:

1. **snapshot blocks known to ExtractGalaxies**: `MASS`, `RHO` (Msun pc⁻³), `U`,
   `NE`, `SFR` (Msun yr⁻¹, density in Msun yr⁻¹ kpc⁻²), `MHOT`, `MCLD`, `MSF`
   (Msun), `EHOT` (10⁵⁰ erg), `AGE` are in physical units; the other known
   blocks (`NH`, `HSML`, `Z`, `TSTP`, `POT`, `NMF`, `EOUT`, ...) in code units;
2. **derived blocks**:
   * gas: `TEMP` (K), `ENTR` (T ρ^-2/3), `THOT` (hot phase temperature), `HOTN`,
     `CLDN` (particle densities of the hot/cold phase, cm⁻³), `FILH` (hot filling
     factor), `PRES` (K cm⁻³), `MMOL` (molecular mass from the pressure law with
     normalisation `Pressure`);
   * all particles, galaxy frame (km/s): `VX VY VZ`, `VR VT VP` (spherical),
     `VRp VTp VZp` (cylindrical), with prefixes `S` (square) and `M` (modulus):
     `SVX`, `MVR`, ..., `SV`, `MV`, `SVp`, `MVp`; `HZ`, `SHZ`, `MHZ` (height);
     `OUTF`, `INFL` (mass flux through the planes at Cut and 2 Cut, or through
     spheres for `VolumeDensity`);
   * chemistry (if the `Zs` block exists): `E<el>` element mass (Msun),
     `Z<el>` mass fraction (`YHe` for helium), `EMet` and `ZMet` total metals;
     element names come from `MetalsFile`. Mass fractions of stars refer to the
     initial stellar mass;
3. **any other block of the snapshot** (galpost extension), in code units.

`list_of_blocks.dat` lists the known and derived names for the current run.

## Particle restrictions

* `LoadParticleList <file>`: analyses restricted to the particles listed in a
  text file with lines `ID type` (types in ascending order). IDs of gas
  particles are compared with the star-formation generation bits masked, and
  stars not found are searched among gas particles, so the list can come from a
  later snapshot.
* `SubUsePartAssignment` (with `ReadSubFindGals`): only particles bound to the
  galaxy's subhalo are used.
