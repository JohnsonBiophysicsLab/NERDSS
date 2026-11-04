# Benchmarks and profiling history - index

This is the catalog of every benchmark, profile, validation run and regression
reproducer kept on `nerdss-optimized-1.4`. Nothing has been moved: each entry
links to the file where it already lives, so the links in
[`docs/nerdss_optimized.md`](../docs/nerdss_optimized.md), the scripts and the
pull requests keep working.

- [1. Timeline](#1-timeline) - every measurement, oldest first
- [2. By topic](#2-by-topic) - serial speed, profiles, MPI, statistics, regressions
- [3. Directory map](#3-directory-map) - what each folder holds
- [4. Harness and case tables](#4-harness-and-case-tables) - the scripts that produce the numbers
- [5. Related records outside `benchmarks/`](#5-related-records-outside-benchmarks)
- [6. Conventions for adding a new record](#6-conventions-for-adding-a-new-record)

Dates are the measurement dates stated inside each report. Many files were
imported to this branch together in `22828ec` (2026-10-14), so `git log` dates
are not measurement dates. Unless a row says otherwise, the host is an Apple M5
(macOS, Apple clang 21, GSL 2.8) and the seed is `20260810`.

## 1. Timeline

| Date | Record | Kind | Build compared | Headline |
| --- | --- | --- | --- | --- |
| 2026-08-10 | [RESULTS.md §1-7](nerdss_optimized/RESULTS.md#1-bitwise-identity-opta-vs-master) | bitwise + timing + statistics | `master` 260f6e2 vs optA / optB (issues #8-#12) | optA byte-identical on 13 cases; optB statistically equivalent |
| 2026-08-10 | [RESULTS.md §8](nerdss_optimized/RESULTS.md#8-reverse-unimolecular-state-change) | correctness | | reverse unimolecular state change |
| 2026-08-11 | [RESULTS.md §9](nerdss_optimized/RESULTS.md#9-profile-driven-follow-up-batch-three-result-preserving-optimizations) | profile-driven optimization | | `invCbrtDr` cache, occupied-cell clearing, third item; all bitwise |
| 2026-08-11 | [RESULTS.md §10](nerdss_optimized/RESULTS.md#10-pairwise-serial-benchmark-of-current-branch-against-mainline) | pairwise timing | branch vs mainline | |
| 2026-08 | [RESULTS.md §11](nerdss_optimized/RESULTS.md#11-find_which_reaction-attempted-measured-reverted) | optimization, reverted | | `find_which_reaction()` rewrite did not pay |
| 2026-08 | [RESULTS.md §12](nerdss_optimized/RESULTS.md#12-how-much-of-section-9s-speedup-is-an-artifact-of-small-test-cases) | scaling check | | how much of §9 survives at larger N |
| 2026-08 | [RESULTS.md §13](nerdss_optimized/RESULTS.md#13-is-any-profile-frame-worth-parallelizing-with-openmp) | profile + OpenMP feasibility | | no frame worth OpenMP |
| 2026-08 | [RESULTS.md §14](nerdss_optimized/RESULTS.md#14-merging-coord-and-vector-into-one-vec3d) | refactor | | `Coord` + `Vector` -> `Vec3D` |
| 2026-08 | [RESULTS.md §15](nerdss_optimized/RESULTS.md#15-coverage-of-the-validation-suite-and-what-it-does-not-reach) | coverage | | what the suite does not reach (`known_uncovered.tsv`) |
| 2026-08-14 | [RESULTS.md §16](nerdss_optimized/RESULTS.md#16-the-allocator-and-the-libm-calls-one-large-win-one-small-one-already-done) | profile-driven optimization | | allocator win; libm `cos`/`exp` ~2% |
| 2026-08 | [RESULTS.md §17](nerdss_optimized/RESULTS.md#17-functions_for_spherical_systemcpp-the-arithmetic-the-copies-and-the-sqrt-that-had-to-stay) | optimization | | spherical-system arithmetic |
| 2026-08 | [RESULTS.md §18](nerdss_optimized/RESULTS.md#18-the-pole-guard-the-caller-never-had) | correctness | | pole guard |
| 2026-08-24 | [RESULTS.md §19](nerdss_optimized/RESULTS.md#19-the-pairwise-cell-grid-seven-corrections-to-how-it-is-sized-and-walked) | correctness + speed | | seven cell-grid corrections |
| 2026-08-25 | [RESULTS.md §20](nerdss_optimized/RESULTS.md#20-a-latitude-longitude-index-for-the-spherical-membrane) | optimization | | spherical shell index (`shell_index_test.cpp`) |
| 2026-08-26 | [RESULTS.md §21](nerdss_optimized/RESULTS.md#21-re-audit-cumulative-drift-and-three-residual-defects) | re-audit | | cumulative drift, three residual defects |
| 2026-08-27 | [RESULTS.md §22](nerdss_optimized/RESULTS.md#22-flattening-molecule-to-structure-of-arrays-built-measured-reverted) | layout experiment, reverted | | structure-of-arrays did not pay |
| 2026-08-27/28 | [RESULTS.md §23-26](nerdss_optimized/RESULTS.md#23-narrowing-molecule-stage-1-the-reweighting-vectors) | layout, staged | | narrowing `Molecule`; stage 4 not worth it |
| 2026-10-09 | [BRANCH_REPORT.md](nerdss_optimized/BRANCH_REPORT.md) | branch vs `master` | frozen reference `3505932d` | change notes + pairwise benchmark |
| 2026-10-09 | [PROSPECTIVE_CHANEG_NOTES_v1.2.md](nerdss_optimized/PROSPECTIVE_CHANEG_NOTES_v1.2.md) | change notes | v1.2 | every change classified stream-changing or not |
| 2026-10-14 | [profile_1_test_2026-10-14/REPORT.md](nerdss_optimized/profile_1_test_2026-10-14/REPORT.md) | profile (user model) | v1.2 reference | `trajWrite = -1` wrote every step: 16x slowdown; MPI not feasible for that box |
| 2026-10-14 | [fiber_b8_3/README.md](nerdss_optimized/fiber_b8_3/README.md) | A/B on fiber model | `master` vs 1.3 | B8.3 fiber branch, rMaxLimit sized by 1D reach |
| 2026-10 (pre-1.4) | [mpi_validation/FINDINGS.md](mpi_validation/FINDINGS.md) | serial vs MPI statistics | branch `mpi-serial-build-modes` 015789c | superseded by the 1.4 MPI records below |
| 2026-10-24 | [nerdss_optimized/mpi_validation/README.md](nerdss_optimized/mpi_validation/README.md) | MPI validation | 1.4 tip `bf0eccd4` | what was still wrong in `nerdss_mpi` |
| 2026-10-27/28 | [AUDIT_2026-10-27.md](nerdss_optimized/mpi_validation/AUDIT_2026-10-27.md) | MPI audit | 1.4 at `4017ece8` | cluster sweeps across region boundary, ghost-ghost pairs, restarts at np>1 |
| 2026-10-29 | [PROPOSAL_2026-10-29.md](nerdss_optimized/mpi_validation/PROPOSAL_2026-10-29.md) | proposal | | 13 numbered MPI items |
| 2026-10-30 | [DEFERRAL_2026-10-30.md](nerdss_optimized/mpi_validation/DEFERRAL_2026-10-30.md) | cost measurement | | cost of the deferred first-order pass (item 10) |
| 2026-10-30 | [SCALING_2026-10-30.md](nerdss_optimized/mpi_validation/SCALING_2026-10-30.md) + [results_scaling_2026-10-30/](nerdss_optimized/mpi_validation/results_scaling_2026-10-30/) | MPI scaling | | wide-box scaling vs serial |
| 2026-11-03 | [BENCH_2026-11-03.md](nerdss_optimized/mpi_validation/BENCH_2026-11-03.md) + [results_bench_2026-11-03/](nerdss_optimized/mpi_validation/results_bench_2026-11-03/) | MPI benchmark + np=1 profile | 1.4 `b4d7438e` vs standalone `f1dd0590` | `nerdss_mpi` beats serial from np=2 on a 42-column box; np=1 gap is the full-grid sweep + `remove_empty_slots` |
| 2026-11-04 | [profile_suite_2026-11-04/REPORT.md](nerdss_optimized/profile_suite_2026-11-04/REPORT.md) | serial profile, whole suite, Linux x86-64 | 1.4 `9200537` | pair walk 37% of the mean case; on sparse boxes 19% is the 13-neighbour `typeMask` scan and 20-28% the cell-list rebuild. An occupancy-bit prefilter is byte-identical on 19/19 cases and 1.10x on `rev_3D` (prototype, not applied) |

## 2. By topic

### Serial speed and bitwise identity

| Record | What it answers |
| --- | --- |
| [RESULTS.md §1-2](nerdss_optimized/RESULTS.md#1-bitwise-identity-opta-vs-master), [§10](nerdss_optimized/RESULTS.md#10-pairwise-serial-benchmark-of-current-branch-against-mainline) | Is the branch byte-identical / how much faster than `master`? |
| [BRANCH_REPORT.md §3](nerdss_optimized/BRANCH_REPORT.md#3-pairwise-benchmark-against-master) | Same question for the frozen 2026-10-09 reference |
| [RESULTS.md §12](nerdss_optimized/RESULTS.md#12-how-much-of-section-9s-speedup-is-an-artifact-of-small-test-cases), `scale_cases.tsv` | Does the speed-up hold at 2k / 10k / 40k molecules? |
| [fiber_b8_3/](nerdss_optimized/fiber_b8_3/README.md) | Fiber model A/B, with runs, timings and figures committed |

### Profiles (where the time goes)

| Date | Record | Workload | Profiler |
| --- | --- | --- | --- |
| 2026-08-11 | [RESULTS.md §9](nerdss_optimized/RESULTS.md#9-profile-driven-follow-up-batch-three-result-preserving-optimizations) | suite | `sample(1)` |
| 2026-08 | [RESULTS.md §13.1](nerdss_optimized/RESULTS.md#13-is-any-profile-frame-worth-parallelizing-with-openmp) | `clathrin` 150k steps | `sample(1)` |
| 2026-08-14 | [RESULTS.md §16](nerdss_optimized/RESULTS.md#16-the-allocator-and-the-libm-calls-one-large-win-one-small-one-already-done) | suite | `sample(1)` + `atos` |
| 2026-10-14 | [profile_1_test_2026-10-14/](nerdss_optimized/profile_1_test_2026-10-14/REPORT.md) | user model `1_test`, 895 molecules | `sample(1)`, `tools/aggsample.py` |
| 2026-11-03 | [BENCH_2026-11-03.md, "What is left at np=1"](nerdss_optimized/mpi_validation/BENCH_2026-11-03.md#what-is-left-at-np1) | `rev_3D`, serial vs MPI np=1 | `sample(1)` |
| 2026-11-04 | [profile_suite_2026-11-04/](nerdss_optimized/profile_suite_2026-11-04/REPORT.md) | all of `cases.tsv`, `coverage_cases.tsv`, `scale_2k/10k`, `enzyme` | gperftools, Linux |

### Optimizations tried and reverted (do not retry without new evidence)

| Record | Why it was reverted |
| --- | --- |
| [RESULTS.md §11](nerdss_optimized/RESULTS.md#11-find_which_reaction-attempted-measured-reverted) | `find_which_reaction()` |
| [RESULTS.md §13](nerdss_optimized/RESULTS.md#13-is-any-profile-frame-worth-parallelizing-with-openmp) | OpenMP: ceiling too low |
| [RESULTS.md §22](nerdss_optimized/RESULTS.md#22-flattening-molecule-to-structure-of-arrays-built-measured-reverted) | structure-of-arrays `Molecule` |
| [RESULTS.md §26](nerdss_optimized/RESULTS.md#26-narrowing-molecule-stage-4-the-ceiling-measured-before-the-work) | `Molecule` narrowing stage 4 |

### Random numbers and statistics

| Record | What it checks |
| --- | --- |
| [RESULTS.md §3-5](nerdss_optimized/RESULTS.md#3-issue-10-random-orientation-uniformity) | quaternion uniformity (#10), `GaussV` (#12), statistical equivalence |
| [rng_quality/](nerdss_optimized/rng_quality/) | `rng_quality.cpp` + `run.sh` behind §3-4 |
| `statistical_check.sh`, `compare_observables.py` | seed-averaged comparison for stream-changing changes |

### MPI

| Date | Record |
| --- | --- |
| pre-1.4 | [mpi_validation/FINDINGS.md](mpi_validation/FINDINGS.md) + [harness](mpi_validation/README.md) (superseded) |
| 2026-10-24 | [nerdss_optimized/mpi_validation/README.md](nerdss_optimized/mpi_validation/README.md) |
| 2026-10-27 | [AUDIT_2026-10-27.md](nerdss_optimized/mpi_validation/AUDIT_2026-10-27.md) |
| 2026-10-29 | [PROPOSAL_2026-10-29.md](nerdss_optimized/mpi_validation/PROPOSAL_2026-10-29.md) |
| 2026-10-30 | [DEFERRAL_2026-10-30.md](nerdss_optimized/mpi_validation/DEFERRAL_2026-10-30.md), [SCALING_2026-10-30.md](nerdss_optimized/mpi_validation/SCALING_2026-10-30.md) |
| 2026-11-03 | [BENCH_2026-11-03.md](nerdss_optimized/mpi_validation/BENCH_2026-11-03.md) |

### Regression reproducers

Each folder under [`nerdss_optimized/regression/`](nerdss_optimized/regression/)
holds a minimal model and a README with the defect, the fix and how to re-run it.

| Folder | Defect it pins |
| --- | --- |
| [add_file_restart](nerdss_optimized/regression/add_file_restart/README.md) | restart with an add file |
| [bimol_state_change](nerdss_optimized/regression/bimol_state_change/README.md) | bimolecular state change |
| [bimol_state_change_rates](nerdss_optimized/regression/bimol_state_change_rates/README.md) | reversible bimolecular state change: equilibrium vs rates |
| [coordinate_file_formats](nerdss_optimized/regression/coordinate_file_formats/README.md) | `-c` XYZ / CSV / PDB |
| [creation_placement_crowding](nerdss_optimized/regression/creation_placement_crowding/README.md) | creating molecules into a crowded box |
| [histogram_mol_states](nerdss_optimized/regression/histogram_mol_states/README.md) | `histogramMolStates` |
| [restart_continuation](nerdss_optimized/regression/restart_continuation/README.md) | exact continuation from a checkpoint |
| [restart_format](nerdss_optimized/regression/restart_format/README.md) | restart format round trip |
| [sphere_no_implicit_lipid](nerdss_optimized/regression/sphere_no_implicit_lipid/README.md) | spherical boundary without implicit lipid |
| [startup_overlap/](nerdss_optimized/regression/startup_overlap/) | six startup-overlap cases: [file_direction](nerdss_optimized/regression/startup_overlap/file_direction/README.md), [pass_cap](nerdss_optimized/regression/startup_overlap/pass_cap/README.md), [per_pair](nerdss_optimized/regression/startup_overlap/per_pair/README.md), [replaced_recheck](nerdss_optimized/regression/startup_overlap/replaced_recheck/README.md), [starting_states](nerdss_optimized/regression/startup_overlap/starting_states/README.md), [vicinity_gate](nerdss_optimized/regression/startup_overlap/vicinity_gate/README.md) |
| [transition_matrix_totals](nerdss_optimized/regression/transition_matrix_totals/README.md) | transition-matrix totals |
| [unimol_state_change_reversible](nerdss_optimized/regression/unimol_state_change_reversible/README.md) | reversible unimolecular state change |

## 3. Directory map

```
benchmarks/
├── README.md                     this index
├── shell_index_test.cpp          unit test for the spherical shell index (RESULTS.md §20)
├── mpi_validation/               pre-1.4 serial-vs-MPI harness and FINDINGS.md (superseded)
│   └── headers/                  species column names per case
└── nerdss_optimized/
    ├── RESULTS.md                26 numbered sections, 2026-08-10 .. 08-28
    ├── BRANCH_REPORT.md          branch vs master at the 2026-10-09 reference
    ├── PROSPECTIVE_CHANEG_NOTES_v1.2.md
    ├── *.tsv, *.sh, *.py         harness and case tables (section 4)
    ├── fiber_b8_3/               fiber model A/B: model, runs, timing, figs, scripts
    ├── mpi_validation/           1.4 MPI: README, AUDIT, PROPOSAL, DEFERRAL, SCALING, BENCH
    │   ├── models/               MPI-specific models (wide_box, deferral, ...)
    │   ├── results_scaling_2026-10-30/
    │   └── results_bench_2026-11-03/
    ├── profile_1_test_2026-10-14/  user-model profile: REPORT, raw profiles, tools/
    ├── profile_suite_2026-11-04/   whole-suite Linux profile: REPORT, tables, tools/
    ├── regression/               minimal reproducers, one folder per defect
    └── rng_quality/              RNG quality comparison (#10, #12)
```

Not committed (see `.gitignore`): `nerdss_optimized/results/` (output of
`run_suite.sh`), `nerdss_optimized/reference*` (reference outputs), and the
per-run `*.counters.txt`, `*.overlaps.txt`, `*.summary.txt`,
`*.manifest.sha256` files.

## 4. Harness and case tables

| File | Role |
| --- | --- |
| [`run_suite.sh`](nerdss_optimized/run_suite.sh) | run a case table with one build; hash every output, record timings |
| [`compare_suites.sh`](nerdss_optimized/compare_suites.sh) | per-case bitwise verdict and speed-up between two `run_suite.sh` results |
| [`interleaved_timing.sh`](nerdss_optimized/interleaved_timing.sh) | timing with builds interleaved per case (what the reported speed-ups use) |
| [`run_mpi_suite.sh`](nerdss_optimized/run_mpi_suite.sh) | same as `run_suite.sh` for `nerdss_mpi` at one rank count |
| [`statistical_check.sh`](nerdss_optimized/statistical_check.sh), [`compare_observables.py`](nerdss_optimized/compare_observables.py) | seed-averaged equivalence for stream-changing changes |
| [`make_scale_inputs.sh`](nerdss_optimized/make_scale_inputs.sh) | generate `scale_2k/10k/40k` from `rev_3D` |
| [`profile_suite_2026-11-04/tools/profile_cases.sh`](nerdss_optimized/profile_suite_2026-11-04/tools/profile_cases.sh) | gperftools profile of every case in a table (Linux) |
| [`profile_1_test_2026-10-14/tools/`](nerdss_optimized/profile_1_test_2026-10-14/tools/) | `sample(1)` profiling and aggregation (macOS) |
| [`mpi_validation/bench_mpi.sh`](nerdss_optimized/mpi_validation/bench_mpi.sh), `*_check.sh`, `*.py` | MPI benchmark and per-defect checks (see that README) |

| Case table | Used for |
| --- | --- |
| [`cases.tsv`](nerdss_optimized/cases.tsv) | the 14-case bitwise and timing suite (also np>1) |
| [`coverage_cases.tsv`](nerdss_optimized/coverage_cases.tsv) | serial-only cases that reach code `cases.tsv` does not |
| [`scale_cases.tsv`](nerdss_optimized/scale_cases.tsv) | size sweep, plus `enzyme` |
| [`pairwise_cases.tsv`](nerdss_optimized/pairwise_cases.tsv) | paired runtime / stochastic-output comparison with `master` |
| [`debugremoval_cases.tsv`](nerdss_optimized/debugremoval_cases.tsv) | debug-function removal measurement |
| [`known_broken.tsv`](nerdss_optimized/known_broken.tsv) | sample inputs that cannot run, and why |
| [`known_uncovered.tsv`](nerdss_optimized/known_uncovered.tsv) | code no case executes |
| [`mpi_validation/cases_short.tsv`](nerdss_optimized/mpi_validation/cases_short.tsv), [`scaling_cases.tsv`](nerdss_optimized/mpi_validation/scaling_cases.tsv), [`coverage_cases_short.tsv`](nerdss_optimized/mpi_validation/coverage_cases_short.tsv) | MPI runs |

## 5. Related records outside `benchmarks/`

| Record | Content |
| --- | --- |
| [`docs/nerdss_optimized.md`](../docs/nerdss_optimized.md) | what changed and why, per issue; links into RESULTS.md |
| [`docs/nerdss_cell_grid_audit.html`](../docs/nerdss_cell_grid_audit.html) | cell-grid audit behind RESULTS.md §19 |
| [`docs/molecule_layout_plan.md`](../docs/molecule_layout_plan.md) | the plan RESULTS.md §22-26 measured |
| [`docs/NERDSS_Parallel_Developer_Guide.md`](../docs/NERDSS_Parallel_Developer_Guide.md) | MPI design notes |
| [`sample_inputs/TIMESTEP_AUDIT.md`](../sample_inputs/TIMESTEP_AUDIT.md) | timestep audit of all sample inputs (`python_scripts/audit_sample_timesteps.py`) |
| [`docs/changelog.md`](../docs/changelog.md) | user-facing changelog |

## 6. Conventions for adding a new record

1. One folder or one Markdown file per measurement, named with the measurement
   date: `<topic>_YYYY-MM-DD/REPORT.md` or `<TOPIC>_YYYY-MM-DD.md`.
2. State at the top: date, commit, host, compiler, flags, seed, nItr, and the
   exact command. A timing means nothing without the iteration count.
3. Commit the scripts that produced the numbers under the record's `tools/`.
   Commit raw tables when they are small (< ~1 MB); otherwise say how to
   regenerate them.
4. Add one row to [§1 Timeline](#1-timeline) and to the matching table in
   [§2](#2-by-topic) in the same commit.
