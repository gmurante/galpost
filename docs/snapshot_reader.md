# Snapshot reader

`galpost.gadget_io` is a standalone reader for GADGET / OpenGadget binary files,
usable independently of the rest of galpost. It needs only numpy.

## Supported files

* **format 2** (`SnapFormat 2`, block labels), with or without INFO block;
* **format 1** (no labels), with a block table;
* multi-file snapshots (`snap_NNN.0`, `snap_NNN.1`, ...);
* companion files with extra blocks and the same particle order (by default
  `snap_NNN.muppi`, MUPPI restart variables); they are also found next to the
  real file when the snapshot is opened through a symbolic link;
* SubFind catalogues written with `WRITE_SUB_IN_SNAP_FORMAT` (`sub_NNN.*`);
* little- and big-endian files, blocks larger than 4 GB (whose 32-bit record
  markers overflow);
* truncated files (e.g. a job killed while writing): complete blocks remain
  readable, the incomplete one raises an explicit error.

## How the layout of the blocks is found

For each block the reader needs the data type, the number of components and the
particle types that carry it. In order of priority:

1. a **block table** given by the user;
2. the **INFO block** written by OpenGadget (`WRITE_INFO_BLOCK`);
3. **inference** from the size of the block and the particle numbers of the
   header, helped by a built-in list of common blocks (POS, VEL, ID, MASS, U,
   RHO, the MUPPI blocks, the SubFind integer fields, ...). Since the true
   layout is always among the candidates, an inference is accepted only if it
   is unique; otherwise the block is reported as ambiguous and reading it
   raises an error asking for a block table. Integer and floating-point data of
   the same size cannot be distinguished from the size alone: unknown blocks
   are assumed to be single-precision floats (double with `double_floats=True`).

On OpenGadget MUPPI snapshots and SubFind catalogues the inference reproduces
the INFO block exactly for every block.

### Block tables

A text file, one line per block (in file order, which matters for format 1):

```
# NAME  TYPE    NDIM  p0 p1 p2 p3 p4 p5
POS     FLOATN  3     1  1  1  1  1  1
ID      LONG    1     1  1  1  1  1  1
MASS    FLOAT   1     1  0  1  1  1  1
Zs      FLOATN  21    1  0  0  0  1  0
```

TYPE is an INFO type name (`FLOAT`, `FLOATN`, `DOUBLE`, `DOUBLEN`, `LONG`,
`LONGN`, `LLONG`, `LLONGN`) or a numpy code (`f4`, `f8`, `i4`, `u4`, `i8`,
`u8`); `p0`–`p5` flag the particle types that carry the block. A table
overrides the INFO block. In the parameter file use `BlockTable <file>`.

The easiest way to write a table is to export it from a file that has an INFO
block:

```python
from galpost.gadget_io import Snapshot, write_block_table
write_block_table(Snapshot("snap_134").info, "blocks.txt")
```

Examples for OpenGadget/MUPPI snapshots and SubFind catalogues are in
`examples/`.

## Usage

```python
from galpost.gadget_io import Snapshot

s = Snapshot("snap_134")              # companions=('.muppi',), block_table=None
print(s)                              # redshift, particle numbers, files
print(s.blocks)                       # all blocks (companions included)
print(s.blocks_for(4))                # blocks carried by stars
h = s.header                          # npart, mass, time, redshift, boxsize, omega0, ...

pos  = s.read("POS", 4)               # (N, 3) float32, code units
zs   = s.read("Zs", 0)                # (Ngas, 21)
mass = s.read("MASS", 1)              # synthesised from the mass table if needed

import numpy as np
idx = np.array([10, 20, 30])          # sorted indices within the type
sub = s.read("EOUT", 0, idx)          # only those particles (memory-mapped)

s.info["Zs"]                          # BlockInfo: dtype, ndim, present, source
s.ambiguous                           # blocks whose layout could not be inferred
```

SubFind catalogue:

```python
c = Snapshot("groups_134/sub_134", companions=())
mcrit = c.read("MCRI", 0)             # groups: type 0
smst  = c.read("SMST", 1)             # subhaloes: type 1, mass per particle type
ids   = c.read("PID", 2)              # particle IDs: type 2
```

See also `galpost.catalogs` ([Python API](python_api.md)), which adds
convenience functions on top of it.

## Units

The reader returns code units. `galpost.units.Units` converts them to the
ExtractGalaxies units (kpc, km/s, Msun, ...) given the GADGET parameter file.
