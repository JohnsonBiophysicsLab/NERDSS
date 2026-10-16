# Restart with an add file

An add-file restart, `nerdss -r DATA/restart.json -a add.inp`, continues a run
from its restart file while bringing in the molecules, reactions and parameter
changes that `add.inp` holds. Nothing else under `regression/` drives this path:
`add_file_molecules_check` (`make checks`) covers
`generate_coordinates_for_restart()` on its own, and `restart_continuation`
restarts without an add file. `check.sh` runs it end to end:

```bash
./check.sh <path-to>/nerdss [work-dir]
REFERENCE_BINARY=<path-to>/nerdss-before ./check.sh <path-to>/nerdss [work-dir]
```

For each case in `cases.tsv` the model is run twice in one directory, the way
an add-file restart is used:

| leg | command | parameters |
| --- | --- | --- |
| 1 | `nerdss -f parms.inp -s <seed1>` | `nItr = N`, `restartWrite = N` |
| 2 | `nerdss -r DATA/restart.json -a add.inp -s <seed2>` | `add.inp` with `nItr = 2N` |

and every file under `DATA/` and `RESTARTS/` is hashed. The case passes when

- both legs exit 0 and leg 2 reports the add file on stdout;
- the header of `DATA/copy_numbers_time.dat` written by leg 2 has every species
  the case's `adds` column names, and the header leg 1 wrote has none of them.
  For both testAdd models those are `A(a!1).B(b!1)` and `B(m!1).M(m!1)`, the
  products of the two reactions `add.inp` brings, so the reactions were
  entered;
- the final `DATA/restart.json` has `simItr = 2N - 1`: leg 2 ran to the add
  file's `nItr`, not the original one;
- a second pair of legs with the same two seeds gives the same hashes for every
  file: the whole run is reproducible;
- with `REFERENCE_BINARY` set, that build gives the same hashes as well: the
  add-file path produces the same bytes as before a change.

Exit status 0 means every case passed. The two cases take a few seconds and need
bash, awk and python3.

## The second seed is essential

Leg 2 must be given its own `-s`. Nothing writes `DATA/rng_state` unless
`RNGwrite = true`, so `read_rng_state()` finds no file at the restart, prints
"RNG state was not restored" and the stream stays the one `gsl_rng_set()` seeded
from the command line, or from the clock when `-s` is absent. Two unseeded
restarts of the same leg-1 directory differ in `final_coords.xyz`,
`trajectory.xyz`, `restart.json`, every `RESTARTS/restart*.json` written after
the restart and the `DATA/*_time.dat` series that count events. With a fixed second
seed, two runs agree on all 32 files, which is what makes the comparison between
builds meaningful.

`RNGwrite = true` is not used here on purpose: an add file changes the system,
so no uninterrupted run can be compared with, and the question is only whether
the add-file path is deterministic and unchanged.

## The restart runs in leg 1's directory

That is how the feature is used: the restart appends to `DATA/trajectory.xyz`
and the time series, and rewrites `DATA/copy_numbers_time.dat` and the other
headers with the new species. The trajectory-length check at a restart reads
`params.trajFile`, which is `trajectory.xyz` in the working directory and not
`DATA/trajectory.xyz`, so it finds no file and only warns; this harness relies
on that.

## What it does not prove

- Nothing here checks where the added molecules were placed, only that the run
  is deterministic. `add_file_molecules_check` covers the placement.
- Nothing is compared against an uninterrupted run; there is none to compare
  with.
- The add file's parameter changes are covered only as far as `nItr`.
- A change to the add-file path that is deterministic and reproduces the same
  bytes on these two models passes; a change meant to alter the output needs a
  new reference, which is the point.

## Cases

Both models are `sample_inputs/testAdd`: A with sites `a` and `m` binds an
implicit lipid M in a 470 x 470 x 752 nm box, or on a sphere of radius 341 nm
with the same lipid density, and `add.inp` brings 100 B with A + B and B + M.
N is 2000 steps of 0.1 us. With the table's seeds, leg 2 ends on the sphere
with 16 A and 8 B bound to M, so the new B + M reaction fires there; in the box
only 2 A are bound and both new columns stay at zero, so that case covers the
entry of the reactions and the propagation of the added molecules, not the new
reactions firing. A + B fires in neither. A longer N would reach them, at the
cost of the few seconds the table takes now.
