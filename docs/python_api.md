# Python API

Everything the command line does is available from Python, which is convenient
for notebooks and custom analyses.

## Running an analysis

```python
from galpost.pipeline import run

res = run("EGparams.surf", "snap_134")      # same as: galpost EGparams.surf -snap snap_134
```

`run` writes the usual files and returns a dictionary with the objects used:

| key | content |
|---|---|
| `params` | the parameters (`galpost.params.Params`) |
| `snapshot` | the snapshot reader |
| `units` | unit conversions |
| `galaxies` | list of `Galaxy` objects (centre, radii, frame, masses, ...) |
| `data` | `ParticleData`: the loaded particles in physical units |
| `analyser` | `GalaxyAnalyser`: frames, analysis results, particle assignment |
| `profiles`, `counters` | the profiles, as arrays |

```python
g = res["galaxies"][0]
g.com, g.vel, g.rotation          # centre (kpc), velocity (km/s), box -> galaxy rotation
g.rgal, g.rvir, g.stellar_mass, g.bovert

an = res["analyser"]
sets = an.analysis[0]             # 12 SetProperties: (gas, DM, cold, hot, stars, cold+stars) x (rgal, rvir)
sets[4].mass, sets[4].half_mass_radius, sets[4].specific_spin

pos, vel = an.frame(g)            # all loaded particles in the galaxy frame
```

## Particle data in the galaxy frame

`ParticleData` holds the particles of the region around the galaxies, ordered
by type as in the snapshot, in ExtractGalaxies units and precision:

```python
d = res["data"]
s = d.type_slice(4)               # stars
mstar = d.mass[s]                 # Msun
age_a = d.age                     # formation scale factor of stars
gas_T = d.temp                    # K (gas)
cold = d.cold_gas()               # boolean, gas
zs = d.metals                     # element masses (gas first, then stars)
```

For the extra analyses, `galpost.extra.common.GalaxyView` gives the particles
within a radius in the galaxy frame, with cylindrical coordinates and derived
quantities:

```python
from galpost.extra.common import GalaxyView
from galpost.profiles import QuantityResolver, read_metal_names

metals = read_metal_names(res["params"]["MetalsFile"], d.nmet)
resolver = QuantityResolver(d, an, res["params"], metals)
v = GalaxyView(an, resolver, 0, radius=30.0)
stars = v.of(4)
R, z, m = v.R[stars], v.z[stars], v.mass[stars]
h2 = v.gas_array("mmol")          # molecular mass of gas particles (0 for others)
```

## Catalogues

```python
from galpost.catalogs import FofCatalog, SubfindCatalog

fof = FofCatalog.from_dir("SubFindPP", 134)       # finds groups_134
fof.mass, fof.cm, fof.len_type, fof.extra["bh_mass"]

cat = SubfindCatalog.from_dir("SubFindPP", 134)
cat.group("MCRI"), cat.group("RCRI")              # any group block
cat.subhalo("SMST")                               # any subhalo block
c = cat.central_subhalo(0)
ids = cat.subhalo_ids(c)                          # particle IDs of the subhalo
sats = cat.subhalos_of_group(0)
```

## Galaxy selection and contamination

```python
from galpost.selection import select_from_subfind, lowres_fraction
from galpost.units import Units
from galpost.gadget_io import Snapshot

snap = Snapshot("snap_134")
u = Units.from_files("paramfile.par", snap.header, htrue=True)
gals = select_from_subfind(cat, u, snap.header, min_stellar_mass=-5,
                           rmin=0.5, rmax=0, ralign=8.0, log_bins=True,
                           do_not_align=False)
lowres_fraction(snap, gals, u)                    # sets g.lowres_fraction
```

## Parameter documentation

```bash
python -m galpost.params --markdown > docs/parameters.md
```

regenerates the [parameter reference](parameters.md) from the code.
