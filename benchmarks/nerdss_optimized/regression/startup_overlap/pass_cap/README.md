# A startup overlap pass that runs out of iterations

`generate_coordinates()` places every molecule at random, then runs
`fixOverlappingMolecules()` over all of them: a pass that looks for pairs of
interfaces closer than the bindRadius (sigma) of a reaction between them and
moves one molecule of each pair to a new random place. It repeats until a pass
finds nothing, at most 50 times. If `-c coords.pdb` is given, the file's
molecules are then moved onto its positions and a second run of passes looks at
those molecules only (its lines on stdout start with `|-`).

Until this check was added, a run of passes that was still finding overlaps on
its 50th iteration just stopped. The simulation started with them, exited 0, and
the only sign was fifty `Overlapping Proteins (true)` lines on stdout. Now the
last iteration prints a warning to stderr naming the run of passes (the one over
all molecules, or the one over the `-c` molecules) and the number of overlapping
pairs it still found. The run still goes on and still exits 0: the warning
does not stop it.

## The model

P and Q, each with one interface 0.5 nm from its center, and
`P(p) + Q(q)` at `sigma = 5`. `check.sh` substitutes the box edge into
`parms.inp` and runs it three ways:

| run   | box    | `-c`          | what has to happen                                        |
|-------|--------|---------------|-----------------------------------------------------------|
| full  | 1.5 nm | no            | full pass: 50 iterations, stderr warning, exit 0          |
| coord | 3 nm   | P at (0,0,0)  | `-c` pass: 50 iterations, stderr warning, exit 0          |
| clear | 200 nm | no            | full pass ends early, no warning, exit 0                  |

In the 1.5 nm box two centers are at most 2.6 nm apart, so the two interfaces
are at most 3.6 nm apart, short of sigma wherever the pass puts them. In the 3 nm
box with P held at the center by the file, Q's center is again at most 2.6 nm
from P's. Both distances also stay inside the 3 nm (`r_P + r_Q + 2`) within which
`are_molecules_in_vicinity()` lets the pair be compared at all. The harness this
model comes from used a 4 nm box for the `-c` case, where a Q center in a corner
is 3.46 nm from P's: the comparison is skipped there, and some seeds get out of
the overlap that way instead of using up the 50 iterations.

Run it with:

```bash
./check.sh <path-to>/nerdss
```

`SEED` selects the seed (default 7). The three runs take well under a second.

## What it shows

With the warning (seeds 1-20 and 12345 all pass):

```
full: 1.5 nm box, no -c
  PASS  exits 0
  PASS  the full pass ran 50 iterations
  PASS  stderr warns about the full pass
  PASS  stdout has no warning
coord: 3 nm box, P at the center from -c
  PASS  exits 0
  PASS  the -c pass ran 50 iterations
  PASS  stderr warns about the -c pass
clear: 200 nm box, no -c
  PASS  exits 0
  PASS  the full pass finished after 1 iteration(s)
  PASS  no warning

PASS
```

On the `nerdss-optimized-1.4` tip before it (4017ece8), the same seeds fail
exactly the two warning checks: both runs of passes give up after 50 iterations
and nothing is written to stderr.
