# Whole-suite serial profile of 1.4, and the next optimization candidate

Date 2026-11-04. `nerdss-optimized-1.4` at `9200537` (after PR #100).
This is the first profile of the branch taken over **every** runnable sample
input at once, and the first on Linux x86-64; the earlier ones (RESULTS.md §9,
§13, §16, `profile_1_test_2026-10-14`, BENCH_2026-11-03) were `sample(1)` on an
Apple M5 and covered one or a few cases each.

## Setup

| | |
| --- | --- |
| Host | Linux 6.18 x86_64, Intel Xeon @ 2.80 GHz, 4 vCPU, shared cloud VM |
| Compiler | g++ 13.3.0, `make serial -j4 CFLAGS='-O3 -g -fno-omit-frame-pointer'` (`-std=c++0x` from the Makefile) |
| GSL | 2.7.1 (Ubuntu `libgsl-dev`), stripped shared library |
| Profiler | gperftools 2.15 `libprofiler.so` by `LD_PRELOAD`, `CPUPROFILE_FREQUENCY=1000`; the VM delivers ~250 samples/s |
| Cases | `cases.tsv` (14), `coverage_cases.tsv` (5), `scale_2k` and `scale_10k` from `make_scale_inputs.sh`, `enzyme` from `scale_cases.tsv`: 22 cases, nItr as in the tables |
| Seed | 20260810 |
| Command | `tools/profile_cases.sh bin/nerdss <out> cases.tsv coverage_cases.tsv` (and the same with `INPUT_ROOT` for the scale cases) |

Files here:

| File | Content |
| --- | --- |
| [`timings.tsv`](timings.tsv) | wall time and sample count per profiled case (under the profiler) |
| [`profiles/`](profiles/) | `<case>.flat.txt` (self time) and `<case>.cum.txt` (inclusive) from `google-pprof --text` |
| [`AGGREGATE.md`](AGGREGATE.md) | per-function mean over the 22 cases, per-case self and inclusive tables |
| [`PAIRWALK_LINES.md`](PAIRWALK_LINES.md) | line-level share of the pair walk's hot lines, per case |
| [`occupancy_prefilter.patch`](occupancy_prefilter.patch) | the prototype measured below (not applied to the branch) |
| [`ab_timing.tsv`](ab_timing.tsv) | raw interleaved A/B wall times of that prototype |
| [`tools/`](tools/) | `profile_cases.sh`, `aggregate.py`, `pairwalk_lines.py`, `ab_timing.sh`, `ab_summary.py` |

Caveats. (1) pprof reports inlined frames (`std::vector::operator[]`,
`push_back`, `QuatRotation::rotate`, ...) as functions of their own, so their
"self" time belongs to whichever function inlined them; read the inclusive
columns for attribution. (2) libgsl is stripped, so its static functions are
billed to the nearest exported symbol: `gsl_ran_wishart` in every table is
really the generator (`gsl_ran_wishart` is at `0x142050` and the next exported
symbol, `gsl_rng_set`, at `0x1495f0`; the static `mt_get` lives in that gap).
(3) The VM is shared, so absolute times are only comparable within one
interleaved A/B pass, not with the M5 numbers elsewhere in this directory.

## 1. Where the time goes

Mean over the 22 cases, every case weighted equally (full table in
[`AGGREGATE.md`](AGGREGATE.md)):

| frame | mean inclusive | top-5 by self time in |
| --- | ---: | ---: |
| `measure_separations_to_identify_possible_reactions` (the pair walk) | **37.1%** | **20 of 22** cases |
| of which `check_bimolecular_reactions` | 20.8% | 10 |
| `Complex::propagate` | 15.5% | 0 |
| `create_complex_propagation_vectors` | 10.8% | 7 |
| `determine_3D_bimolecular_reaction_probability` | 8.4% | 1 |
| `Complex::update_properties` | 6.5% | 10 |
| `SimulVolume::update_memberMolLists` (cell-list rebuild) | 6.2% | 0 |
| `gsl_ran_gaussian_ziggurat` + generator (`gsl_ran_wishart`, see caveat 2) | 6.0% + 2.9% | 8 |

The cases split into two regimes, and the next candidates differ between them:

**Sparse / large boxes** - `rev_3D`, `scale_2k`, `scale_10k`, `rev_3Dto2D`:
about one molecule per sub-cell (`rev_3D`: 2 000 molecules in 27 000 cells).
Here `check_bimolecular_reactions` is only 1-2%; the time is the *bookkeeping
of the grid itself*:

| case | pair walk incl. | neighbour scan (lines 303-307) | first target load (line 199) | cell-list rebuild incl. | of which `add_member` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `rev_3D` | 33.2% | **18.7%** | 6.4% | **19.5%** | 13.2% |
| `scale_2k` | 34.4% | **19.7%** | 6.2% | **20.1%** | 13.8% |
| `scale_10k` | 33.1% | **19.1%** | 7.5% | **28.2%** | 22.4% |
| `rev_3Dto2D` | 22.5% | 10.9% | 6.2% | 10.6% | 6.0% |
| `hexamer` | 35.3% | 8.0% | 1.8% | 4.8% | 2.4% |
| `gagsphere` | 28.8% | 8.0% | 0.7% | 4.5% | 2.5% |

Lines are those of
[`measure_separations_to_identify_possible_reactions.cpp`](../../../src/reactions/measure_separations_to_identify_possible_reactions.cpp)
at `9200537`; the per-case numbers for all 22 cases are in
[`PAIRWALK_LINES.md`](PAIRWALK_LINES.md).

**Dense, assembling models** - `trimer`, `hetTrimer`, `homoTrimer`,
`closed_homoTrimer`, `cluster_gagsphere`, `unimol_state_change`:
`check_bimolecular_reactions` is 36-66% inclusive and
`determine_3D_bimolecular_reaction_probability` 18-32%, spent in
`pirr_pfree_ratio_psF` (`exp`, `erf`, `erfc`), `passocF` and the overlap sweep
`sweep_separation_box` (22.8% in `homoTrimer`). That is the physics of pairs
that really are in range; RESULTS.md §11, §16 and §22-26 have already taken the
cheap layout and libm wins there.

`Complex::propagate` and its callees (propagation vectors, Euler rotation,
ziggurat Gaussians, reflection) are the third block, 15% on average and 39% in
`clathrin` and `sphere`; RESULTS.md §9, §12 and §16 already worked this path.

## 2. Why the neighbour scan is expensive

For each target molecule the walk reads the `typeMask` of each of its 13
half-shell neighbour sub-cells before it looks at any partner:

```cpp
for (int neighCellItr : neighborList)
  if ((simulVolume.subCellList[neighCellItr].typeMask & targMask) != 0)
    liveCells[numLiveCells++] = neighCellItr;
```

`SubVolume` is about 80 bytes (four indices, two `std::vector` headers and the
mask), so 27 000 of them are ~2 MB, and the 13 neighbours of a cell are spread
over three z-planes of that array. In a sparse box nearly every one of those
13 loads misses cache *and finds an empty cell*: `rev_3D` occupies ~1 800 of
27 000 cells (6.7%). The answer it needs - "is anything there?" - is already
held in `SimulVolume::occupancyMask`, one bit per cell (3.4 KB for 27 000
cells, resident in L1), which `add_member()` sets together with `typeMask`.

## 3. Candidate 1 (recommended next): test the occupancy bit before the SubBox

[`occupancy_prefilter.patch`](occupancy_prefilter.patch) adds one test in
front of the `typeMask` load:

```cpp
if (!((occupancy[neighCellItr >> 6] >> (neighCellItr & 63)) & 1))
  continue;
```

**Why it is exact.** Every writer that sets a `typeMask` bit (`add_member()`)
also sets the cell's occupancy bit, and every writer that clears
`occupancyMask` (`clear_member_lists()`, the `MpiContext` overload of
`update_memberMolLists()`) clears `typeMask` in the same pass; the MPI
deserializer rebuilds the mask from non-empty member lists. So
`typeMask != 0` implies the bit is set, a cell skipped by the new test would
have failed the old one, and the accepted cells and their order are unchanged.
No arithmetic, no random draw and no iteration order moves.

**Bitwise check.** `run_suite.sh` with the baseline and the prototype (same
flags), `compare_suites.sh`: **byte-identical on all 19 cases** of
`cases.tsv` (14) and `coverage_cases.tsv` (5) - every `DATA/`, `RESTARTS/`
and `PDB/` file.

**Speed** - interleaved A/B (`tools/ab_timing.sh`, A and B alternated per
repetition, 3 repetitions, median wall time, raw data in
[`ab_timing.tsv`](ab_timing.tsv)):

| case | reps | baseline median (s) | prototype median (s) | speed-up | per-rep ratios |
| --- | ---: | ---: | ---: | ---: | --- |
| rev_3D | 3 | 18.30 | 16.65 | 1.099x | 1.099, 1.127, 1.142 |
| rev_3Dto2D | 3 | 14.99 | 14.55 | 1.030x | 0.932, 1.043, 0.954 |
| homoTrimer | 3 | 13.70 | 13.55 | 1.011x | 1.058, 1.044, 0.913 |
| hetTrimer | 3 | 13.21 | 13.27 | 0.995x | 1.007, 1.023, 0.992 |
| hexamer | 3 | 14.59 | 14.47 | 1.008x | 1.008, 0.981, 1.098 |
| clathrin | 3 | 6.45 | 6.63 | 0.973x | 1.092, 0.948, 0.931 |
| sphere | 3 | 15.56 | 15.82 | 0.984x | 0.994, 0.977, 0.977 |
| gagsphere | 3 | 7.63 | 7.15 | 1.066x | 1.091, 1.100, 0.859 |
| scale_10k | 3 | 27.09 | 25.10 | 1.079x | 1.005, 1.082, 1.105 |
| enzyme | 3 | 9.64 | 9.54 | 1.011x | 1.025, 0.918, 1.092 |
| **sum of medians** | | 141.15 | 136.74 | **1.032x** | |

Case nItr as in the suite tables except `scale_10k` (2 000) and `enzyme` (4 000), to keep the
pass short. On this shared VM a single repetition moves by up to ±9% (see the per-rep ratios of
the dense cases, where the change cannot matter), so only a consistent shift means anything:

* **Where the profile predicted a gain, every repetition is faster:** `rev_3D`
  **1.10x** (1.10, 1.13, 1.14) and `scale_10k` **1.08x** (1.01, 1.08, 1.11), the cases where
  the neighbour scan is ~19% of the run. `gagsphere` (8% scan) and `rev_3Dto2D` (11%) lean
  positive but are within the noise.
* **Where the scan is small, nothing measurable changes** (`hetTrimer` 0.995x, `hexamer`
  1.008x, `homoTrimer` 1.011x, `enzyme` 1.011x, `clathrin` 0.973x with reps on both sides of 1).
  **`sphere` is the one to watch:** 0.984x, and all three repetitions below 1 (0.99, 0.98,
  0.98). Its scan is only 1.5% of the run, so the extra bit test cannot save much there, and a
  ~2% cost from it is plausible. Re-time `sphere` on the reference host before landing; if it
  holds, gate the test on grid sparsity (e.g. occupied cells < 25% of all cells, computed once
  per step in `refresh_occupied_cells()`), which keeps the change exact.
* Sum of medians over the ten cases: 1.032x.

Verdict: a byte-identical, two-statement change worth ~10% on sparse and large boxes,
the regime the branch's scale work targets, neutral on dense cases, and possibly ~2% slower on
`sphere`. It is the recommended
next optimization.

**Cost and risk.** Two statements in one file; the serial and MPI builds share the
code path. To land it: apply the patch, `make serial`, `make serial legacy`,
`make mpi`, run `cases.tsv` + `coverage_cases.tsv` with `bin/nerdss` and
`bin/nerdss_legacy` against a baseline run, and re-time with
`interleaved_timing.sh` on the reference M5 host.

## 4. Candidate 2: an incremental cell-list update

The cell-list rebuild, `SimulVolume::update_memberMolLists()`, empties every
occupied cell and re-bins every molecule each step. On the sparse/large cases
it is the single biggest frame after the pair walk: **19.5% `rev_3D`, 20.1%
`scale_2k`, 28.2% `scale_10k`**, of which `add_member()`'s `push_back` into a
randomly placed `SubVolume` is 13-22%. Its share grows with N, so it matters
most for the large systems the branch is meant to reach.

Most molecules do not change cell in a step: the sub-cell edge is `rMaxLimit`
(tens of nm) while the per-step displacement is `sqrt(6 D dt)` (about a
nanometre for the suite's parameters). An incremental update would compute
each molecule's new bin (the arithmetic it does today), and only for a molecule
whose bin changed erase it from the old list and insert it into the new one at
its sorted position.

What keeps it bitwise: the full rebuild produces every member list in
ascending molecule index (it loops over `moleculeList` in order), and the pair
walk visits members in list order, so the incremental version must keep each
list sorted (sorted insert, ordered erase), clear `typeMask`/occupancy bits
when a cell empties, and recompute `typeMask` for a cell that lost a member.
Creation (`create_molecule_and_complex_from_rxn` already calls
`add_member()`), destruction, restarts and the `simItr % 1000` bounds-check
branch must all route through the same two operations. Expected payoff is
most of the 10-28% above on sparse and large boxes, close to nothing on the
dense assembling cases. Larger than candidate 1 and touches the invariants
documented in `class_SimulVolume.hpp`, so it should follow it, with the
regression folders that create and destroy molecules
(`creation_placement_crowding`, `restart_continuation`, `add_file_restart`)
added to the check. The `MpiContext` overload has the same structure and the
full-grid sweep BENCH_2026-11-03 identified; doing the serial overload first
gives a template for it.

## 5. Other frames, and why they rank lower

| Frame | Share | Assessment |
| --- | --- | --- |
| line 199, `moleculeList[target].isImplicitLipid` | 6-7.5% on sparse cases | The first touch of the target `Molecule`; the miss would move to line 235 (`freelist`) rather than disappear. Implicit lipids are never binned by `update_memberMolLists()`, so the test itself is redundant in the serial build, but removing it saves only the branch. |
| `check_bimolecular_reactions` rejection cascade | 44-66% incl. on dense cases; its prologue alone 4.6% of `hetTrimer` | A gate-field mirror was built and reverted (RESULTS.md §22.7). A caller-side pre-test is possible but must reproduce the exemptions for the volume-exclusion, surface and fiber paths; needs its own study. |
| `determine_3D_bimolecular_reaction_probability` / `pirr_pfree_ratio_psF` | 18-32% on trimers | `exp`/`erf`/`erfc` of a per-pair argument; nothing to cache across steps. Already examined in §16. |
| `sweep_separation_box` | 22.8% `homoTrimer`, 6% mean | Overlap sweeps of assembled complexes; worth a dedicated profile with `clusterOverlapCheck` cases. |
| ziggurat + Mersenne Twister | ~9% | Stream-changing to replace; out of scope for a bitwise branch. |
| `write_json_restart` | 3.7% `homoTrimer` | An artefact of the shortened runs: with `checkPoint` unset the checkpoint interval defaults to `nItr / 10`, so a 6 000-step run writes ten restarts; the model as shipped (`nItr = 1e9`) writes one per 10^8 steps. |

## 6. Reproducing

```bash
make serial -j4 CFLAGS='-O3 -g -fno-omit-frame-pointer'
cd benchmarks/nerdss_optimized
P=profile_suite_2026-11-04/tools
$P/profile_cases.sh ../../bin/nerdss /tmp/prof/main cases.tsv coverage_cases.tsv
bash make_scale_inputs.sh /tmp/scale
INPUT_ROOT=/tmp/scale/sample_inputs $P/profile_cases.sh ../../bin/nerdss /tmp/prof/scale <(grep -E '^scale_(2k|10k)' scale_cases.tsv)
python3 $P/aggregate.py AGGREGATE.md /tmp/prof/*/profiles
python3 $P/pairwalk_lines.py ../../bin/nerdss /tmp/prof/*/profiles/*.prof
```

For the A/B: build the prototype in a separate worktree (`git worktree add`,
`git apply occupancy_prefilter.patch`, same `make` line), then
`run_suite.sh` + `compare_suites.sh` for the bitwise check and
`$P/ab_timing.sh 3 <base> <proto> ab_timing.tsv <cases>` +
`$P/ab_summary.py ab_timing.tsv` for the timing.
