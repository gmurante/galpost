# galpost documentation

galpost analyses the galaxies of OpenGadget / MUPPI simulations. It is a
Python re-implementation of ExtractGalaxies, with the same parameter files and
output files, plus a set of new analyses.

1. [Installation](installation.md)
2. [Quick start](quickstart.md)
3. [The parameter file](parameter_file.md) and the [parameter reference](parameters.md)
4. [Output files](outputs.md)
5. [Extra analyses](extra_analyses.md)
6. [Halo finding with OpenGadget](subfind_opengadget.md)
7. [Snapshot reader](snapshot_reader.md)
8. [Progenitor tracking](tracking.md)
9. [Graphical interface](gui.md)
10. [Python API](python_api.md)
11. [Compatibility with ExtractGalaxies](compatibility.md)
12. [Performance](performance.md)

## Typical workflow

```
simulation snapshots ──► OpenGadget FoF/SubFind (optional, post-processing)
        │                         │
        └──────────┬──────────────┘
                   ▼
   galpost <parameter file> -snap snap_NNN
                   │
                   ▼
  profiles, galaxy properties, extra analyses, galaxy snapshots, maps
```
