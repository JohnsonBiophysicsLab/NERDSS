# MPI validation, 2026-10-24 (nerdss-optimized-1.4)

What was still wrong with `nerdss_mpi` at the 1.4 tip (`bf0eccd4`, = tag
`nerdss-optimized-1.3-ref-2026-10-24`), what this branch fixes, and what is
still open.  The earlier write-up in `benchmarks/mpi_validation/FINDINGS.md`
describes the state before PRs #47, #52, #53 and #61; most of what it lists
was fixed by those.  Its one np=1 physics difference (auto_phos) is gone: with
32 seeds each, np=1 and serial agree on all five auto_phos species (smallest
p = 0.11).

## How the audit was run

* Builds, all from one tree: `-O3` (`make mpi`); a checker build
  (`-O3 -DCHECK_STRUCTURE=true -DCHECK_OWNERSHIP=true -DCHECK_OWNERSHIP_EVERY=1`,
  possible from the command line since 20fc7303); AddressSanitizer with both
  checkers (`-g -O1 -fsanitize=address -fsanitize-recover=address`,
  `ASAN_OPTIONS=detect_leaks=0:halt_on_error=0`).  The tip was audited with the
  exchange-conflict report (e2135edb) patched in and nothing else.
* `run_mpi_suite.sh` over `cases.tsv` and `coverage_cases.tsv` at np=2, 3 and
  4: the checker build with seeds 1-3, `-O3` with seed 4, ASan with seed 1 on
  the quarter-length tables `cases_short.tsv` / `coverage_cases_short.tsv`.
  `sphere` and `cluster_gagsphere` are spherical and skipped by the suite.
* Serial reference: `bin/nerdss` at the tip, seeds 1-8 (`cases.tsv`) and 1-4
  (`coverage_cases.tsv`), and more seeds for the follow-ups below (up to 168
  for homoTrimer), run through the same suite by `serial_as_mpirun.sh`.
* `schedule.py` runs the suites under a core budget, strictly in order.
  `analyze_audit.py` tabulates failures, checker / ASan / conflict counts, mass
  conservation (molecules of each type summed over the merged complex
  histogram, every block, against serial's constant total) and copy numbers:
  each column of the merged `copy_numbers_time.dat` averaged over the second
  half of the run (equilibrium) and over the first half after time 0
  (relaxation, where a rate error shows), per seed, as Welch t of the MPI
  ensemble against the serial one.  `audit_tables.py` prints the tables below.
* Probes for what the suites cannot see: `pair_count_check.sh`,
  `reweight_history_check.sh`, `coordinate_file_check.sh`, and
  `conflict_control.sh` with the fast hetTrimer models in `models/`.

The machine was shared with other sessions throughout (load 5-140 on 10
cores), so no timing here means anything.

## Defects, by severity

### Wrong physics in every np>1 run

| # | defect | evidence | fix |
|---|---|---|---|
| 1 | Most candidate pairs across each rank boundary were never listed.  `measure_separations_to_identify_possible_reactions()` skipped ghosts as the outer molecule of its half-shell cell walk, so a ghost/owned pair was met only from the owned side; across a boundary that is usually the ghost's side, and the left neighbour drops the pair as its right ghost's, so nobody drew it | step 1, identical coordinates at every np: homoTrimer 2071 pairs at np=1, 1994/2016/1883 at np=2/3/4; hetTrimer 929 vs 895/887/859 | b699efe2: ghosts pair with owned molecules (ghost-ghost pairs still skipped); pairs inside the left ghost column are dropped (the left neighbour's, and a straddling complex let both ranks list them); after: 2071 and 929 at every np (`pair_count_check.sh`) |
| 2 | The pair reweighting history was lost or corrupted at the edges.  `Molecule::serialize()` does not carry `prevReweight`, both deserializers overwrite a received molecule whole, and `remove_empty_slots()` (every MPI step) remapped every molecule index except the history's | homoTrimer np=2, 3000 steps: history found for 1.7% of pairs held by edge molecules vs 77.4% in the interior (np=1: 78.3%); about a fifth of all pairs were drawn unreweighted | 825a8c3a: keep the local history on overwrite, remap and prune it on compaction and departure; after: 76.7% vs 77.8% (`reweight_history_check.sh`) |
| 3 | Every rank seeded its generator with the same -s value, so the ranks drew overlapping windows of one stream | by reading `main()`: at step 1 the first molecule on every rank got the same displacement | fa765010: ranks above 0 reseed with seed + rank after set-up; np=1 unchanged |

Defects 1 and 2 pull in opposite directions -- fewer pairs, but the edge pairs
that were listed were drawn without reweighting -- which is part of why the
tip's copy numbers showed them only in homoTrimer, the smallest box with the
most binding (8 bonds short of serial at np>1, p = 0.02; below).

### Wrong for a class of model, silently

| # | defect | evidence | fix |
|---|---|---|---|
| 4 | Compartments do not exist under MPI: `Membrane::serialize()` drops `hasCompartment`/`compartmentR`, and the MPI step has no transmission phase | `compartment`: serial transmits 3 complexes by step 30000, np=2/3/4 none; clatComp(cd3) t = 9-13 at every np | 22a5ae8c refuses compartment models |
| 5 | Sphere: each rank kept the whole `totalSA` while `nSites` was split, so the implicit lipid's density was 1/np | `VALIDATE_SUITE/sphere`, 1000 steps, bound A, seeds 1-3: serial 177/179/195, np=1 187/183/182, np=2 113/117/119 | 3357fbaa (exact: a sphere's area between two planes is proportional to their distance); np=2 after: 168/186/198 |
| 6 | Sphere: zeroth-order creation drew with the whole `sphereVol` over the whole sphere on every rank: np times the rate, products in other ranks' slabs | by reading `check_for_zeroth_order_creation()` / `create_random_coords()` (a box already uses its slab) | 22a5ae8c refuses it at np>1 |
| 7 | `nerdss_mpi -c` was parsed and then dropped | rev_3D with every A at x = -400 nm, np=2: rank 1 still held 490 A at time 0 | 84e8c773 (`coordinate_file_check.sh`) |

### Broken state at some rank counts

| # | defect | evidence | fix |
|---|---|---|---|
| 8 | More ranks than x columns: the extra ranks own none (startCell > endCell) | trimer (3 columns) at np=4: every run fails, 6759 ownership and 11 structure reports and 9 conflicts over 5 runs, ASan container-overflow in `send_data_to_left_neighboring_ranks()` (send_receive.cpp:86) | 3998d2f1, folded into 29319315 |
| 9 | Narrow ranks act, in the same half-step, on the column the right neighbour binds into: a middle rank owning 1 or 2 columns, rank 0 owning 2 | trimer np=3 at the tip: 1 exchange conflict; trimer with onRate 1000 / offRatekb 1e5 at np=3: 3 conflicts, then ids owned by nobody (itr 494, 848) | 29319315 |
| 10 | Step-top first-order reactions (dissociation, state change, destruction) on a complex the right neighbour binds into during its first half, through a ghost copy that predates them | fast hetTrimer controls at np=4, after #9: 4 conflicts in 2 of 12 runs | 0363a361; 0 of 12, and 0 of 12 on the final build |

### Restarts

| # | defect | evidence | fix |
|---|---|---|---|
| 12 | `nerdss_mpi` wrote no restart file at all (`write_output()`'s restart block was `if (0)`), at any rank count | by reading | e447989c: `write_gathered_restart()`, below |
| 13 | `nerdss_mpi -r <file>` opened `restart.dat<rank>` whatever the argument, by testing `Parameters::rank`, which nothing initializes (so which file either build opened depended on that int) | `nerdss_mpi -np 1 -r restart500.json`: "could not find restart file" | 5d5f0866: `-r` reads the file it names |

### Latent

| # | defect | evidence | fix |
|---|---|---|---|
| 11 | A complex slot recycled from `Complex::emptyComList` kept the destroyed complex's `ownerRank`, which `claim_unowned_complexes()` then never overrides | by reading; unreachable today (the list is emptied each MPI step before the next dissociation) | 2fc60f87 (`run_code_tests/complex_destroy_reset_check.cpp`) |

Tooling: 20fc7303 (checker macros from the compiler line, suite timeout,
ASan column), e2135edb (EXCHANGE CONFLICT report, exchange_conflicts
column), b699efe2 and 825a8c3a (CANDIDATE PAIRS and REWEIGHT HISTORY
diagnostics in CHECK_STRUCTURE builds), `run_code_tests/mpi_column_split_check.cpp`.

### The column rule (29319315)

A step acts on each rank's x columns in two halves with an exchange after each
(`region_mid_index()`).  The right neighbour binds into this rank's last owned
column (bin n-2, its bin 0) in its first half, so this rank's first half must
not reach bin n-2; an initiator pairs one column further than itself, so the
left region has to end two columns short of it.  A middle rank cannot move its
split (bin 1 has to be finished before the first exchange), so it needs three
columns; rank 0 has no left constraint and does everything in its second half
when it owns two or fewer; the last rank has no right neighbour.
`MpiContext::column_share()` keeps the even split wherever every middle rank
already gets three, gives each middle rank exactly three otherwise, and
`prepare` refuses with the largest usable -np when even that is impossible.

So at np=3/4 some cases.tsv models are now refused that the tip ran: trimer
(3 columns) and michaelis_menten (4) run at np<=2 only, unimol_state_change (7)
at np<=3.  The tip ran michaelis_menten and unimol_state_change at np=4 with no
conflict detected in 3 seeds each; the refusal rests on the mechanism and on
the positive controls, not on a failure seen in those models.  unimol (7) at
np=3 now splits 2,3,2, hetTrimer (8) 1,3,3,1 and homoTrimer (10) 2,3,3,2 at
np=4; every other case keeps its decomposition.

### The deferred first-order pass (0363a361)

Owned complexes with a member in bin n-2 or beyond carry the MPI-only
`Molecule::firstOrderHeld` through the step-top pass and get their
first-order reactions right after the first exchange, on the state that holds
the neighbour's first-half bonds.  A complex the neighbour bound into sits
the pass out, as a reactant of an association does in a serial step; what the
late pass reacts keeps its candidate pairs for overlap avoidance only, as a
serial `isDissociated` molecule does.  (`isGhosted` is not usable for the
mark: under `mpi_`, `check_dissociation_implicitlipid()` counts this rank's
free proteins by it.)

### The known leftovers from earlier sessions

* Two owners after a cross-rank association: real, but only between a
  bimolecular pass and the exchange that follows it, where no checkpoint can
  look (bb998c75 removed the one there).  Its only possible consequence is the
  other rank acting on the molecule in the same half-step: #9 and #10, which
  the conflict report counts.
* Two ranks binding the same stripe molecule in one half-step: #9 and #10;
  now reported, and gone in every run and control since.
* `ownerRank` on recycled slots: #11.
* coverage_cases.tsv at np>1: in the audit at np=2/3/4, checkers and ASan.
* np=4: #8 and #9.

## Audit results

Tip = `bf0eccd4` with the conflict report; branch = `825a8c3a`.  Runs per
case and np: checker seeds 1-3, `-O3` seed 4, ASan seed 1 (short tables);
coverage cases: checker seed 1 and ASan seed 1 on the branch, seeds 1-3 and
ASan on the tip.

### Failures and reports

Every failure on the branch is a refusal at set-up: `compartment` (#4),
trimer at np=3/4, michaelis_menten at np=3/4 and unimol_state_change at np=4
(too few columns, #8/#9).  No ownership, structure, ASan, exchange-conflict or
mass-conservation report anywhere on the branch, in any model, at any np.

| case | np | tip: runs/failed, reports | branch: runs/failed, reports |
|---|---|---|---|
| clathrin | 2 | 5/0 | 5/0 |
| clathrin | 3 | 5/0 | 5/0 |
| clathrin | 4 | 5/0 | 5/0 |
| closed_homoTrimer | 2 | 5/0 | 5/0 |
| closed_homoTrimer | 3 | 5/0 | 5/0 |
| closed_homoTrimer | 4 | 5/0 | 5/0 |
| cluster_mem_loc | 2 | 4/0 | 2/0 |
| cluster_mem_loc | 3 | 4/0 | 2/0 |
| cluster_mem_loc | 4 | 4/0 | 2/0 |
| compartment | 2 | 4/0 | 2/2 |
| compartment | 3 | 4/0 | 2/2 |
| compartment | 4 | 4/0 | 2/2 |
| gagsphere | 2 | 4/0 | 2/0 |
| gagsphere | 3 | 4/0 | 2/0 |
| gagsphere | 4 | 4/0 | 2/0 |
| hetTrimer | 2 | 5/0 | 5/0 |
| hetTrimer | 3 | 5/0 | 5/0 |
| hetTrimer | 4 | 5/0 | 5/0 |
| hexamer | 2 | 5/0 | 5/0 |
| hexamer | 3 | 5/0 | 5/0 |
| hexamer | 4 | 5/0 | 5/0 |
| hist_mol_states | 2 | 4/0 | 2/0 |
| hist_mol_states | 3 | 4/0 | 2/0 |
| hist_mol_states | 4 | 4/0 | 2/0 |
| homoTrimer | 2 | 5/0 | 5/0 |
| homoTrimer | 3 | 5/0 | 5/0 |
| homoTrimer | 4 | 5/0 | 5/0 |
| implicit_lipid | 2 | 5/0 | 5/0 |
| implicit_lipid | 3 | 5/0 | 5/0 |
| implicit_lipid | 4 | 5/0 | 5/0 |
| mem_localization | 2 | 5/0 | 5/0 |
| mem_localization | 3 | 5/0 | 5/0 |
| mem_localization | 4 | 5/0 | 5/0 |
| michaelis_menten | 2 | 5/0 | 5/0 |
| michaelis_menten | 3 | 5/0 | 5/5 |
| michaelis_menten | 4 | 5/0 | 5/5 |
| rev_2D | 2 | 5/0 | 5/0 |
| rev_2D | 3 | 5/0 | 5/0 |
| rev_2D | 4 | 5/0 | 5/0 |
| rev_3D | 2 | 5/0 | 5/0 |
| rev_3D | 3 | 5/0 | 5/0 |
| rev_3D | 4 | 5/0 | 5/0 |
| rev_3Dto2D | 2 | 5/0 | 5/0 |
| rev_3Dto2D | 3 | 5/0 | 5/0 |
| rev_3Dto2D | 4 | 5/0 | 5/0 |
| trimer | 2 | 5/0 | 5/0 |
| trimer | 3 | 5/0 1 confl | 5/5 |
| trimer | 4 | 5/5 6759 own, 11 struct, 2 asan, 9 confl | 5/5 |
| unimol_state_change | 2 | 5/0 | 5/0 |
| unimol_state_change | 3 | 5/0 | 5/0 |
| unimol_state_change | 4 | 5/0 | 5/5 |

### Copy numbers against serial

Largest |t| over a case's copy-number columns against the 8-seed (cases.tsv)
or 4-seed (coverage) serial ensemble, uncorrected for the number of columns
or of cases (about 150 tests per table, so a few |t| near 3 are expected by
chance); "-" where a case writes too few rows for a first half, or was
refused.  Spherical cases are not run.

| case | np | tip: max\|t\| 2nd half (p) / 1st half | branch: max\|t\| 2nd half (p) / 1st half |
|---|---|---|---|
| clathrin | 2 | 2.63 (0.04) / 1.47 | 1.37 (0.23) / 3.06 |
| clathrin | 3 | 2.47 (0.052) / 1.91 | 1.96 (0.078) / 3.14 |
| clathrin | 4 | 2.76 (0.032) / 1.85 | 3.49 (0.0063) / 1.77 |
| closed_homoTrimer | 2 | 2.45 (0.041) / - | 2.32 (0.046) / - |
| closed_homoTrimer | 3 | 2.78 (0.022) / - | 2.60 (0.027) / - |
| closed_homoTrimer | 4 | 1.96 (0.1) / - | 1.39 (0.19) / - |
| cluster_mem_loc | 2 | 1.00 (0.36) / - | - |
| cluster_mem_loc | 3 | 0.79 (0.47) / - | - |
| cluster_mem_loc | 4 | 1.57 (0.22) / - | - |
| compartment | 2 | 9.32 (0.00073) / 3.09 | - |
| compartment | 3 | 12.92 (5.1e-05) / 3.95 | - |
| compartment | 4 | 11.77 (8.6e-05) / 3.95 | - |
| hetTrimer | 2 | 1.07 (0.31) / - | 2.04 (0.12) / - |
| hetTrimer | 3 | 0.78 (0.47) / - | 0.49 (0.65) / - |
| hetTrimer | 4 | 0.81 (0.44) / - | 2.45 (0.048) / - |
| hexamer | 2 | 0.81 (0.45) / - | 1.90 (0.11) / - |
| hexamer | 3 | 0.73 (0.5) / - | 2.56 (0.029) / - |
| hexamer | 4 | 0.61 (0.56) / - | 1.36 (0.22) / - |
| hist_mol_states | 2 | 2.18 (0.13) / 1.00 | - |
| hist_mol_states | 3 | 1.19 (0.32) / 1.00 | - |
| hist_mol_states | 4 | 0.98 (0.39) / 0.39 | - |
| homoTrimer | 2 | 2.45 (0.041) / - | 2.32 (0.046) / - |
| homoTrimer | 3 | 2.78 (0.022) / - | 2.60 (0.027) / - |
| homoTrimer | 4 | 1.96 (0.1) / - | 1.39 (0.19) / - |
| implicit_lipid | 2 | 1.12 (0.32) / 0.59 | 0.23 (0.83) / 0.37 |
| implicit_lipid | 3 | 0.41 (0.69) / 1.24 | 1.03 (0.34) / 1.54 |
| implicit_lipid | 4 | 2.32 (0.065) / 0.40 | 0.02 (0.99) / 0.19 |
| mem_localization | 2 | 2.57 (0.057) / - | 0.65 (0.53) / - |
| mem_localization | 3 | 1.64 (0.17) / - | 0.65 (0.53) / - |
| mem_localization | 4 | 0.92 (0.38) / - | 1.82 (0.099) / - |
| michaelis_menten | 2 | 0.45 (0.67) / 0.67 | 0.75 (0.49) / 1.27 |
| michaelis_menten | 3 | 1.44 (0.18) / 1.53 | - |
| michaelis_menten | 4 | 0.73 (0.49) / 1.53 | - |
| rev_2D | 2 | 0.19 (0.86) / 0.19 | 3.71 (0.0059) / 0.39 |
| rev_2D | 3 | 0.35 (0.75) / 0.79 | 2.48 (0.055) / 1.49 |
| rev_2D | 4 | 0.19 (0.85) / 1.16 | 0.97 (0.37) / 0.11 |
| rev_3D | 2 | 1.90 (0.095) / 0.10 | 0.08 (0.94) / 0.40 |
| rev_3D | 3 | 0.79 (0.46) / 0.60 | 0.61 (0.56) / 0.77 |
| rev_3D | 4 | 1.38 (0.22) / 1.60 | 0.53 (0.62) / 1.36 |
| rev_3Dto2D | 2 | 0.91 (0.39) / 1.16 | 1.22 (0.27) / 1.71 |
| rev_3Dto2D | 3 | 0.75 (0.49) / 0.59 | 2.26 (0.055) / 1.21 |
| rev_3Dto2D | 4 | 0.12 (0.9) / 0.42 | 2.30 (0.053) / 1.34 |
| trimer | 2 | 0.99 (0.39) / 1.15 | 0.92 (0.42) / 1.34 |
| trimer | 3 | 1.01 (0.36) / 1.16 | - |
| unimol_state_change | 2 | 1.21 (0.29) / 0.90 | 2.57 (0.052) / 2.28 |
| unimol_state_change | 3 | 0.93 (0.37) / 0.89 | 0.72 (0.49) / 1.94 |
| unimol_state_change | 4 | 1.27 (0.25) / 1.02 | - |

The two kinds of value that stand out were followed with more seeds:

* **rev_2D** (branch np=2, |t| 3.71 with 4 seeds): with 24 seeds each and 32
  serial ones, bound pairs serial 115.1, branch 115.7 (np=2) / 114.5 (np=3),
  tip 115.3 / 116.4, every p > 0.5.  Chance.
* **homoTrimer and clathrin** (|t| 2-3.5 at the tip and on the branch): with
  24 serial seeds instead of 8, clathrin agrees everywhere (12 branch seeds at
  np=2 and 4, every column p > 0.13).  homoTrimer bonds at t = 1.5 ms: serial
  75.8 (168 seeds), np=1 73.6 (156, p = 0.08 against serial), branch np=2-4
  72.0 (28; p = 0.46 against np=1), tip np=2-4 67.8 (12, p = 0.02 against
  serial).  The tip's further deficit at np>1 is what defect #1 predicts; the branch adds nothing
  to np=1.  Whether np=1 itself sits 3% low is not settled at this power.
* **auto_phos** (`hist_mol_states`): 32 seeds of serial and np=1 agree on all
  five species (smallest p = 0.11).

### Probes on the final build

| check | result |
|---|---|
| `pair_count_check.sh` | homoTrimer 2071 and hetTrimer 929 candidate pairs at step 1 for np=1-4 |
| `reweight_history_check.sh` | history found for 76.7% of edge pairs, 77.8% interior |
| `coordinate_file_check.sh` | rank 0 1000 A, rank 1 0 at time 0 |
| `conflict_control.sh` | het_fast np=3 and np=4, het_faster np=4, seeds 1-4 each: 0 conflicts, 0 structure, 0 ownership |
| `run_code_tests` (`make checks`) | all pass, including complex_destroy_reset_check and mpi_column_split_check |
| serial byte identity | `run_suite.sh` on cases.tsv and coverage_cases.tsv: all 273 + 92 output hashes equal to the tip's, after every commit that touches a file the serial build compiles |
| GCC 15 clean copy | 196 serial sources and the MPI build with g++-15, no warnings; `make checks` passes |


## Checkpoints from any rank count

`write_gathered_restart()` (`src/io_mpi/write_gathered_restart.cpp`) writes
the serial build's restart files, on its schedule and under its names:
`DATA/restart.json` every `restartWrite` steps and at the end of the run,
`RESTARTS/restart<step>.json` every `checkPoint` steps, the `.dat` format when
`legacyRestartFormat` is set, and with `rngwrite` rank 0's generator state
(`DATA/rng_state`, `RESTARTS/rng_state<step>`).  Every rank sends rank 0 what
it owns, with every reference written as an id; rank 0 numbers the whole
system once, maps the ids back, and writes it with `write_restart()`.  Per-rank
tallies are summed (free lipid sites per state, counters, event histograms,
transition matrices, each observable's change from the common start) and
cluster lifetimes concatenated; what the split rescaled (water box bounds and
volume, membrane area, sites) comes from the system as parsed; what depends on
the numbering is rebuilt (`bndpartner` in `bndlist`'s order, `monomerList`,
`numEachMolType`, `mySubVolIndex`, and `bindPairList` by `associate_box()`'s
rule, which the serial population-based dissociation follows).  It reads
state only: with it, np=2/3 runs write the same outputs byte for byte.

Such a file restarts with serial `nerdss -r` or `nerdss_mpi -np 1 -r`
(np>1 restarts stay refused; see below).  Verified:

* every cases.tsv / coverage model at np=2, 3 and 4 (quarter-length runs,
  checkpoint at mid-run; compartment and the too-narrow configurations are
  refused): both continuations run with no checker report, and where the MPI
  run wrote a row at the checkpoint step both start from exactly that row;
  under ASan at np=3, no report in the writer or the np=1 restart;
* exact continuation at np=1 (the gathered numbering is then rank 0's own):
  all 14 models write the same copy-number rows from the checkpoint on as the
  uninterrupted run;
* homoTrimer at np=4 to step 2500, serial from the checkpoint to 6000, 24
  seeds: bonds at step 5000 72.8 +- 1.8 against serial's 71.0 +- 0.8 (168
  seeds), p = 0.39;
* `bindPairList` as written equals, bond for bond, the list the serial build
  keeps itself (checked on a serial checkpoint by the same rule);
* implicit lipid: free sites plus bound interfaces equal nSites in every
  checkpoint, with the whole box's area; the `.dat` format restarts too.

`checkpoint_check.sh` repeats the np=1 exact continuation, the np>1 restart
with both builds, and the `bindPairList` rule.  Not carried, as in the serial
format: the reaction-candidate state of the step in progress (a checkpoint is
taken between steps), and a model's `currSimTime` differs in its last digit
after a restart (time arithmetic relative to the restart step).

## Still open

* **Restarting at np>1.**  Checkpoints are written at any np (above) and
  continue serially or at np=1, but `-r` is still refused at np>1.  What
  reading them on several ranks would take:
  1. a split for a bonded system: `prepare_rank_data()` copies molecules
     without remapping `partnerIndex`/`bndpartner`, trims each complex to the
     members the rank sees (a complex has to travel whole), and seeds
     ownership by "has a member this rank owns", which gives a straddling
     complex two owners -- a new run gets away with all three only because
     nothing is bound at step 0;
  2. the per-rank shares of what the writer sums: the implicit lipid's free
     sites per state (`finish_rank_setup()` sets them to the rank's nSites
     only for a new run), and per-rank RNG streams (seed + rank, as a new run
     does, since a single saved state cannot be split).
* **Fiber (promoter) molecules across ranks**: `Molecule::serialize()` does
  not carry `isPromoter` nor `Complex::serialize()` `onFiber`, so a 1D fiber
  molecule a neighbour sends back would arrive as an ordinary 3D one.  Found
  while writing the checkpoint (which carries both itself); not yet measured
  on a fiber model under MPI.
* **Ghost-ghost pairs a rank should decide.**  A pair of two molecules that
  are both ghosts of the left neighbour's complexes but sit in columns this
  rank decides (bin 1, or bins 0 and 1) has no initiator anywhere: this rank
  skips it, the left neighbour drops it as its right ghost's.  Needs a complex
  straddling the boundary; not measured.
* **What the deferral costs**: a complex the neighbour bound into takes no
  first-order draw that step (its older bonds lose one dissociation chance),
  and a molecule that changes state in the late pass cannot also bind that
  step, which in serial it could, with its new state.  Both second order in
  the rates and confined to the edge column.
* **Spherical models** stay outside the validation set by user decision; #5
  and #6 were found by reading and one probe, and nothing else about the
  sphere's decomposition was checked.
* **Output merging** covers copy numbers, the complex histogram and PDB/CIF
  frames only; bound-pair, mono/dimer, event-counter, observable and
  transition-matrix files stay per rank, as partial tallies to be summed.
* **`-c` itself** expects coordinates in the centred frame (-L/2..L/2) while
  NERDSS writes PDB frames in 0..L, so a NERDSS PDB frame is not accepted as
  it is by either build (serial exits 1).  Not MPI-specific.
* **Performance**: `find_molecule()`/`find_complex()` are linear scans called
  per received molecule and per partner, so each exchange is quadratic in the
  stripe population; np=2 runs were roughly 3-5x slower than serial here.
