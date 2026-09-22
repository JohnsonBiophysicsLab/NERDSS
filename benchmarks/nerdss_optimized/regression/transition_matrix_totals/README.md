# Transition-matrix totals

This check guards the transition-matrix bookkeeping of molecule types that set
`countTransition = true`. `check.sh` runs
[`rev_3D`](../../../../sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D)
with transition counting switched on for both of its types, and checks every
matrix it writes against the total that matrix must have.

```bash
./check.sh <path-to>/nerdss
```

Exit status 0 means every matrix totalled exactly one entry per complex per step.
`SEED` and `NITR` override the seed (20260810) and the length (20000 steps,
about seven seconds).

## Why the total is exact

The matrices are filled lazily. When an event changes how many R a complex holds
(an association in `associate_box()`/`associate_sphere()`, a dissociation in
`break_interaction()`), every complex holding R is credited with the steps since
the previous update for R, `Parameters::lastUpdateTransition[R]`. The complexes
that took part get the stays before this step plus one transition; every other
complex gets the whole interval as stays. Either way, each complex holding R
gains one entry per step. (The one exception is a complex that splits into two
which both hold R: it is credited with both transitions.)

In `rev_3D` every molecule has a single site, so a complex holds at most one R,
nothing can split that way, and the 1000 copies of R always sit in exactly 1000
complexes. The matrix for R
must therefore total `1000 x lastUpdateTransition[R]` at every write, whatever
the trajectory, seed or compiler. `check.sh` requires every record to be a
multiple of 1000, and the final record to equal that product exactly, with
`lastUpdateTransition` read back from `DATA/restart.dat`.

The matrix for A stays empty. A is the first reactant of the only reaction, so it
is always `reactMol1`, the complex whose count is compared, and the A count of
that complex never changes. It is checked anyway (`0 = 1000 x 0`).

## What it catches

A destroyed complex keeps its slot in `complexList` until the end of the
timestep, with `numEachMol` cleared and `isEmpty` set. The loops that credit the
"unchanged complexes" walked every slot, so they read past the end of a
destroyed complex's `numEachMol`. Built with `-O3` that read returns the old
counts, and a complex that no longer existed was credited with stays. The first
such complex was always wrong: `associate_box()` destroys the second reactant's
complex just before its loop, then counted those old contents a second time, on
top of the stays it had just credited them explicitly. Under AddressSanitizer
with libc++ (Apple clang), every one of these reads is a `container-overflow`.

## Results

| build | seed | R matrix total | 1000 x lastUpdateTransition[R] | verdict |
| --- | --- | --- | --- | --- |
| `nerdss-optimized` at `01bf678` | 20260810 | 19948590 | 19931000 | FAIL, +17590 |
| fixed | 20260810 | 19931000 | 19931000 | PASS |

Across seeds 20260810, 1, 2, 3, 7 and 42, the base overshoots by 16585 to 17991
entries (0.08-0.09%) and fails the multiple-of-1000 test at 9 or 10 of its 10
non-empty records. The fixed build passes on all six, and a `g++-15` build of it
passes with the default seed.
Both builds follow the same trajectory: `lastUpdateTransition[R]` agrees for
every seed, and only the matrix differs.

The excess is roughly one interval per association-driven update, against
`copies` intervals for the whole system, so its relative size scales as
1/copies: about 0.09% here with 1000 copies, 0.6-0.7% on the 100-copy
`VALIDATE_SUITE/trimer`.
