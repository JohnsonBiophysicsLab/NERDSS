# Serial vs MPI validation on non-spherical sample inputs

Build: branch `mpi-serial-build-modes` @ `015789c` (includes the `isDissociated`
fix). Open MPI 5.0.9, Apple clang 21, `-O3`. 20 seeds per configuration, three
configurations per case: `serial`, `mpi -np 1`, `mpi -np 2`.

## How spherical systems were excluded

Spherical systems are known-broken and are excluded by reading the parameter
file, never by directory name:

```
isSphere = true      # or
sphereR  = <value>
```

Name-based exclusion would have been wrong in both directions in this tree:
`RefinedGagSphere/` is **not** a spherical system, and `testAdd/sphere/` **is**.
The check lives in both `screen2.sh` and `bench2.sh`, so a spherical case cannot
reach the benchmark by accident.

Spherical inputs found and excluded: `sphere/parms_sphere.inp`,
`gagsphere/parms.inp`, `VALIDATE_SUITE/sphere/parms_sphere.inp`,
`testAdd/sphere/parms.inp`, `testAdd/sphere/1/parms.inp`.

## What is measured

Two independent questions, because they can disagree:

1. **Species abundance** — the final row of `DATA/copy_numbers_time.dat`.
   Welch's t-test per species (unequal variances; the MPI groups are noisier),
   corrected across the whole family with Holm and Benjamini-Hochberg.

2. **Assembly pathway** — `DATA/histogram_complexes_time.dat`, which lists at
   each recorded time how many copies of each distinct *complex composition*
   exist. Two runs can agree on every species count and still get there by
   different routes (many small complexes vs few large ones), so this is tested
   separately:
   * a label-permutation test on the L1 distance between the serial and MPI mean
     histograms — one global test per case, no distributional assumption, and
     appropriate because per-composition counts are neither independent nor
     normal (mass conservation couples them);
   * per-composition Welch tests with BH correction, to say *which* complex
     types differ.

   MPI writes compositions in a different order than serial, so every
   composition is canonicalised (`IL: 1. A: 1.` and `A: 1. IL: 1.` both become
   `A:1,IL:1`) before comparison.

Scripts: `bench2.sh` (runs), `analyze.py` (species), `analyze_pathway.py`
(pathway). Data: `bench2_species.tsv`, `bench2_runs/*/hist.dat`.

### Both tests were validated before use

The statistics are hand-rolled (no scipy on this machine), so they were checked
against known answers rather than trusted:

* Student-t implementation vs published tables: `t(0.975, df=10) = 2.2281`,
  `t(0.975, df=30) = 2.0423`, `p(t=2.0, df=30) = 0.0546` — all exact.
* Pathway test **negative control**: `implicit_lipid`, which the `isDissociated`
  fix repaired, returns `p_perm = 0.67`, L1 = 0.97% of total — indistinguishable.
* Pathway test **positive control**: the `sphere` case, known broken, returns
  `p_perm = 0.0249`, L1 = 94.8% of total, with all three compositions flagged.
  So the test has power and is not merely failing to reject.

## Cases excluded, and why

| input | reason |
|---|---|
| `sphere`, `gagsphere`, `VALIDATE_SUITE/sphere`, `testAdd/sphere*` | spherical per `.inp` (user-directed exclusion) |
| `compartment/clath_compartment.inp` | **segfaults in the serial build** (nItr=100000) — not an MPI issue |
| `VALIDATE_SUITE/create_destroy/create.inp` | **segfaults in the serial build** (exit 139, nItr=20000) |
| `genetic_oscillator/clock_model.inp` | **segfaults in the serial build** (exit 139, nItr=200000) |
| `clathrin_coat/flat_clat-ap2-pip2.dir` | serial fine, but `mpi -np 1` exceeded a 240 s cap at nItr=2000 |
| `gag_coat/solution/parms.inp` | assembly cost grows super-linearly; a single nItr=20000 serial run did not finish in 10 min |

The three serial segfaults are independent of MPI and worth a separate look —
`compartment` and `create_destroy` and `genetic_oscillator` all die in the plain
serial build.

One unverified observation: the killed `gag_coat/solution` run left a
well-formed 8-column row reading `gag(homo) = -10`, a negative copy number. The
same file also contains a duplicated header, so the output was in an odd state
and I cannot separate a genuine negative count from damage done by killing the
process. Recorded as a lead, not a finding.

## Results

11 non-spherical cases, 20 seeds per configuration, 660 runs.

### np=1 is not a test of the parallel machinery

At `nprocs == 1` `is_ghosted()` cannot return true (its branches are `rank > 0`
and `rank < nprocs-1`), `prepare` computes `ratio = 1.0` and keeps the whole box,
and every send/recv sits behind `if (tempRank)`. No ghost zones, no boundary
molecules, no messages. np=1 answers "does compiling with -Dmpi_ change the
physics", which is a real question, but it says nothing about whether the
decomposition works. The np>=2 rows are the ones that matter.

### np=1: 10 of 11 cases indistinguishable

Only `auto_phos` differs: 5/5 species shifted, pathway p_perm = 0.00005,
L1 = 36.5% of total. The pathway says how -- MPI over-forms A-A dimers
(A:2, 3.65 -> 14.10) and A-Phos complexes (5.25 -> 7.95) while depleting free A
(95.45 -> 71.85) and free Phos (3.75 -> 1.05). The `isDissociated` fix improved
this case but did not repair it; it is a genuine physics difference in the MPI
code path, independent of decomposition.

Every other case is clean on both metrics: `implicit_lipid`, `clathrin_flat`,
`clathrin_pucker`, `rev_3D`, `rev_2D`, `rev_3Dto2D`, `uni_state_rev`,
`michaelis`, `mem_loc_IL`, `mem_loc_EL`.

### np=2: mostly cannot complete, and wrong when it does

| case | np=2 completion | mass conservation | impossible complexes | pathway |
|---|---|---|---|---|
| `implicit_lipid` | 20/20 | conserved | none | **differs**, p=0.0166 |
| `michaelis` | 16/20 | **VIOLATED** (S 108->109) | none | **differs**, p=0.0413 |
| `mem_loc_EL` | 14/20 | **VIOLATED** (A, B, M) | **`A:2`** | indistinguishable |
| `rev_2D` | 9/20 | **VIOLATED** (A, R: 792-801 vs 800) | **`A:2`, `R:2`** | indistinguishable |
| `mem_loc_IL` | 5/20 | conserved | none | indistinguishable |
| `clathrin_flat` | 1/20 | **VIOLATED** (clat 100->54) | none | n too small |
| `clathrin_pucker` | 0/20 | - | - | - |
| `auto_phos` | 0/20 | - | - | - |
| `rev_3D` | 0/20 | - | - | - |
| `rev_3Dto2D` | 0/20 | - | - | - |
| `uni_state_rev` | 0/20 | - | - | - |

Five of eleven cases cannot complete a single multi-rank run in twenty attempts.
Of the six that complete at least once, four violate mass conservation.

Serial and np=1 conserve mass **exactly** -- every run, every species,
min = max = the declared input count. np=2 does not.

`rev_2D`'s only reaction is `A(a) + R(r) <-> A(a!1).R(r!1)`, so a complex of two
A molecules cannot exist. The np=2 histograms contain them literally
(`1  A: 2.`), 14 occurrences over 9 surviving runs; serial and np=1 produce zero
in 20 runs each. Molecule totals drift both below and above the fixed 800
(792-801), so material is lost *and* duplicated.

Caveat on the novel-composition check: it also flags `clathrin_pucker` np=1 with
`clat:18/19/20/31`. Clathrin assembles to arbitrary size, so those are legitimate
compositions that serial simply never sampled in 20 runs. The check is only
conclusive where the reaction network forbids the composition, as in `rev_2D`
and `mem_loc_EL`.

### The two tests disagree in both directions, usefully

`implicit_lipid` np=2: pathway differs (p=0.0166) while **no individual
composition** is flagged -- the global permutation test pools the distribution
and is more sensitive than per-bin tests.

`rev_2D` np=2: global test says indistinguishable (p=0.60) while the `A:2`
composition is flagged (BH p=0.029) -- an impossible-but-rare complex is
significant per-bin and negligible in total L1.

Keep both tests; neither subsumes the other.

