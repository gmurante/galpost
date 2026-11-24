# Performance

## What galpost reads

Only the particles inside a sphere around the selected galaxies are loaded. Its
radius covers every volume used by the analysis: rvir, the profile volume
(√2 rmax for cylinders, √3 rmax for grids), the galaxy snapshot and the extra
analyses, plus a margin; if the centre moves during tracking beyond what the
margin allows, a larger region is loaded automatically. Selecting the region
requires one pass over the positions of all particles; the other blocks are
read only for the selected particles, through memory maps.

`WriteFullCatalog`, `WriteStarCatalog` and `COMP_CM` need the whole snapshot.

## Typical cost

Zoom-in simulation of a Milky-Way-sized halo with 6 × 10⁷ particles (Aquila C4
resolution), one galaxy, one snapshot, one core:

| configuration | time | memory |
|---|---|---|
| volume profiles (`VolumeDensity`, 6 quantities) | 0.5–1.5 min | 3.3 GB |
| surface profiles + `AnalyzeGalaxies` (18 quantities, metals) | 1.5–3 min | 5 GB |
| as above + all extra analyses | ~1.5 min | 10 GB |

The first read of a snapshot is limited by the file system; repeated runs on
the same snapshot are faster because the data are in the page cache. Memory
scales with the number of particles in the loaded region.

For many snapshots, run one galpost process per snapshot (e.g. a SLURM job
array or several processes on one node): processes are independent.

## Transparent huge pages

On Linux machines where *transparent huge pages* (THP) are set to `always`
(`cat /sys/kernel/mm/transparent_hugepage/enabled`) and memory is filled by the
page cache, each large allocation can force the kernel to compact memory. On a
cluster login node this made allocating 3 GB take 11 s of kernel time, and an
analysis run up to six minutes instead of 40 s.

galpost therefore disables THP for its own process at start-up
(`prctl(PR_SET_THP_DISABLE)`, no privileges needed). Set
`GALPOST_KEEP_THP=1` to keep the system setting. The same trick helps any
memory-hungry Python code on such nodes:

```python
import ctypes
ctypes.CDLL(None).prctl(41, 1, 0, 0, 0)    # PR_SET_THP_DISABLE
```

## Tips

* Run large snapshots on compute nodes rather than on login nodes.
* Keep snapshots on a local or parallel file system; reading from network
  mounts of desktop machines (or `/mnt/c` under WSL2) is slow.
* Profiles of many quantities cost little: the expensive part is reading and
  selecting the particles.
* `Verbose` prints the intermediate steps (selection, centre tracking,
  alignment, particle counts), useful to compare with ExtractGalaxies logs.
