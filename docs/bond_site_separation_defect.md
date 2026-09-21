# A bond written between interfaces that are not touching

The serial build could write a bond between two interfaces many sigma apart. The
bond then glued two complexes that had never met into one reported complex: a
6BNO actin run would report a single 75-mer made of filaments that are visibly
separate in the trajectory.

Measured 2026-09-17 on `nerdss-optimized` at 8025fd7, Apple M-series, Apple clang
`-O3 -std=c++0x`, on the 75-molecule 6BNO actin model in
[`run_code_tests/BondSiteSeparation`](../run_code_tests/BondSiteSeparation).

| iterations | seeds | runs with a far-apart bond, unfixed | fixed | worst separation, unfixed |
| ---: | ---: | ---: | ---: | ---: |
| 100 000 | 1-60 | 9 (15%) | — | 34.4 sigma |
| 1 000 000 | 1-100 | 23 (23%) | **0** | 34.4 sigma |

"Far-apart" means more than 1.15 sigma, checked on every 1000th iteration at
1 000 000 iterations and every 100th below that. Bonds live for about 290 000
iterations in this model (`offRatekb` 8.13 /s at a 0.424 us step), so a bad bond
that forms and breaks between two checks is rare enough to ignore.

## Symptom

`DATA/COMPLEXES/*.json` gives, per bonded complex, each member's centre of mass,
the quaternion taking its `.mol` template to its current orientation, and the
bond list. Reconstructing each bonded site as `com + R(q) . template_site` and
dividing by the bond's sigma gives, for a clean run of this model, a largest
separation of 1.063 sigma across every bond of every frame. An affected run
reports separations of 5 to 34 sigma, on hundreds of frames, because a bad bond
persists until it dissociates.

## Mechanism

A crossing is a pair of interfaces that `check_bimolecular_reactions()` decided
might react this timestep. The criterion depends on whether the two molecules are
in the same complex:

- **Different complexes.** `determine_{1,2,3}D_bimolecular_reaction_probability()`
  admits the pair on a diffusion-based `Rmax`, roughly
  `3 * sqrt(6 * Dtot * dt) + sigma`. That is far looser than sigma, and it is
  supposed to be: associating the pair means physically moving the two complexes
  together first, which `associate_box()`/`associate_sphere()` do in their
  `reactCom1.index != reactCom2.index` branch.
- **One complex.** `evaluate_binding_within_complex()` records a crossing only
  when the separation is below `bindRadSameCom * bindRadius`, with
  `bindRadSameCom` defaulting to 1.1. It has to: the association is a loop
  closure, and the loop-closure branch bonds the pair where it stands, rotating
  and translating nothing. Whatever separation the pair is at becomes the bond
  length.

Crossings for the whole system are recorded once per timestep, and then consumed
one molecule at a time by the reaction loop in `EXEs/nerdss.cpp` (the loop under
*"decide whether to perform reactions for each protein"*). Both steps read
complex membership, but membership changes underneath the second one, and a
crossing admitted under the first rule can be consumed under the second.

Concretely, `-s 101` at iteration 52423:

```
ITR:52423,BOND,A,8,1,A,9,0       <- different complexes; merges them
ITR:52423,BOND,A,45,3,A,69,2     <- now a loop closure, at 5.2 sigma
```

Molecule 8 associates with molecule 9, which merges the trimer [8, 43, 45] with
the 10-mer holding 9 and 69. The crossing between 45 and 69 was recorded while
they were still in separate complexes. `associate_box()` zeroes the crossings of
the two molecules that reacted, and no others, so when the loop reaches molecule
45 its crossing with 69 still fires. The two are now in one complex, so
`associate_box()` takes the loop-closure branch, moves nothing, and writes a bond
between `aa2f` and `aa2b` 5.2 sigma apart. At iteration 52500 all four molecules
sit in one 13-mer.

Instrumenting the loop-closure branch to print the separation it is handed makes
the split plain. Thirteen of the seeds below that the fix changes were run with
it, each to just past the iteration where the fix first acts: 614 loop closures.
601 are at most 1.063 sigma. The other 13 are the defect, one per seed, from 1.110
to 34.4 sigma, and in every one both molecules are already marked `propagated` —
the flag `associate_box()` sets on every member of the complexes it merges.

`propagated` alone does not identify the defect, though. 53 of the 601 legitimate
closures have both molecules `propagated` too: they close a loop inside a complex
that an earlier association in the same timestep had just merged, which is
perfectly valid, because the merge moved that complex rigidly and left the
separation inside it untouched.

## Fix

`associate()` (`src/reactions/association_dispatch.cpp`) now re-applies the
same-complex criterion at the moment of association, against the complexes as
they are then rather than as they were when the crossing was recorded:

```cpp
if (reactCom1.index == reactCom2.index) {
    double sameComSep { currRxn.bindRadSameCom * currRxn.bindRadius };
    double R1 { calc_interface_distance(...) };
    if (!(R1 < sameComSep))
        return;
}
```

`associate()` is the only caller of `associate_box()` and `associate_sphere()`,
and the serial loop and `perform_bimolecular_reactions()` are the only callers of
`associate()`, so one guard covers both association geometries and both reaction
loops. The comparison is the one admission uses — `get_distance()` records a
crossing when `R1 < Rmax` — so nothing admitted as a loop closure can be refused
at the same geometry. It is negated so that a NaN separation is refused rather
than bonded.

The refused association is dropped the way the existing `cancelAssoc` paths drop
one: nothing is written, and the pair is simply reconsidered next timestep, where
`evaluate_binding_within_complex()` will not record a crossing for it at all
while it stays out of range.

`calc_interface_distance()` is the separation computation lifted out of
`get_distance()` so that the guard and the admission check cannot drift apart.
It carries the geometry — geodesic on a sphere, along x on a fiber, in-plane on a
flat membrane, otherwise Euclidean — and `get_distance()` now calls it.

### Why not skip molecules whose complex already moved

The MPI reaction loop in `src/reactions/perform_bimolecular_reactions.cpp` skips
any molecule whose `trajStatus` is `propagated`, and the serial loop does not.
Adding that skip to the serial loop would also close this hole, but it is a
blunter instrument:

- It would refuse **legitimate** same-step loop closures. If an association
  merges complex C with complex D, every member of C is marked `propagated` — but
  a loop closure inside C is still perfectly valid, because complexes move
  rigidly and the separations inside C did not change. This is not hypothetical:
  it is the 53 legitimate closures above, out of 601.
- It is asymmetric. A crossing between a merged molecule and an untouched one
  fires or does not depending on which of the two the outer loop reaches first.
- `propagated` is also set by dissociation and by unimolecular reactions earlier
  in the timestep, so the skip would change behaviour well beyond this defect.

The guard in `associate()` refuses exactly the associations whose geometry does
not support them, and nothing else.

## What the fix changes, and what it does not

The 100 seeds at 1 000 000 iterations, unfixed against fixed, comparing the full
association/dissociation event log (`assocDissocWrite`):

- **70 seeds: identical.** Their event logs match, and so does every other file
  the sweep kept under `DATA/`, byte for byte (the per-frame `COMPLEXES` JSON is
  discarded for clean runs). The guard never fired.
- **30 seeds: diverge, each at exactly one refused association.** In every case
  the logs are identical up to an iteration that carries two `BOND`s, and the fix
  keeps the first and drops the second. From there on the trajectories part, as
  they must.

23 of the 30 are the runs the 1.15 sigma check flags, at 1.2 to 34.4 sigma. The
other seven — seeds 4, 16, 20, 74, 88, 92 and 99 — are refused just over the 1.1
limit, at **1.1104 or 1.1146 sigma**, which is why the unfixed runs looked clean
to a 1.15 sigma check. They are the same defect: a crossing admitted across two
complexes, consumed as a loop closure after an earlier association in the same
timestep merged them, both molecules `propagated`. Each of the two values recurs
to seven digits across unrelated seeds, which says they are fixed arrangements of
this model's filament geometry rather than chance.

That deserves a note for whoever builds these models. A loop closure at 1.11
sigma can never form the legitimate way, because
`evaluate_binding_within_complex()` will not record a crossing above 1.1, and the
complex is rigid, so the gap never shrinks. The unfixed build formed these
closures *only* through the defect. If they are wanted, the model should raise
`bindRadSameCom` for these reactions (1.12 would admit both), rather than rely on
the coincidence that used to bond them.

## Verification

**Bitwise.** All 13 cases of `benchmarks/nerdss_optimized/cases.tsv` and all 5
of `coverage_cases.tsv`, seed 20260810, every output file hashed: all 18
byte-identical, across 365 output files, every run exiting 0. None of these
models closes a loop out of range, so none of them should move, and none does.

**Standalone checks.** `make checks` runs two new ones alongside the existing
two:

- `loop_closure_separation_check` drives `associate()` directly on two molecules
  in one complex: it must bond them at 0.9 and 1.2 sigma and must refuse them at
  exactly the limit, at 5.2 and 27.8 sigma, and at a NaN separation. Against the
  unfixed `association_dispatch.cpp` the four refusals fail.
- `interface_distance_check` compares the refactored `get_distance()` and the new
  `calc_interface_distance()` bit for bit against a verbatim copy of the previous
  body, over 256 000 fixtures spanning every branch. No sample model sets
  `isPromoter`, so the fiber arm — the one whose promoter override of `sep` was
  rewritten from `(a && !b) || (!a && b)` to `a != b` — is unreachable by the
  suite; 32 000 of the fixtures exercise it. No mismatches.

The tree builds cleanly, and all four checks pass, under Apple clang 21, GCC 15
and GCC 16 (`make CC=g++-15 serial checks`); the repository's CI builds with GCC,
but only for `master`. Under GCC the pre-refactor `get_distance.cpp`, compiled on
its own, also agrees bit for bit with the refactored one over the same 256 000
fixtures.

**MPI.** `make mpi` builds with the change, with no new diagnostics. The MPI
build cannot say more than that here. Two-rank runs crash alike with and without
the change: the 6BNO model on both builds, and two of three hexamer seeds on both
builds, with bus errors and traps. That is the known boundary-ownership defect
that PR #21 reports on 12 of 67 inputs. `nerdss_mpi` also does not reproduce
itself run to run, so it is not compared bitwise. The guard sits in
`associate()`, which the MPI reaction loop calls exactly as the serial one does.

## Regression tests

- [`run_code_tests/loop_closure_separation_check.cpp`](../run_code_tests/loop_closure_separation_check.cpp),
  run by `make checks`, pins the criterion without depending on any trajectory.
- [`run_code_tests/BondSiteSeparation`](../run_code_tests/BondSiteSeparation)
  reproduces the defect end to end: five seeds, 40 000 iterations each, asserting
  that no bonded pair in `DATA/COMPLEXES` is further apart than 1.15 sigma. All
  five fail on the unfixed build. Its seeds stay meaningful only while their
  trajectories do, which is the reason the standalone check exists.

See also [`bond_bookkeeping_defects.md`](bond_bookkeeping_defects.md) for the
bookkeeping around the same `bnd*` pushes in `associate_box()`.
