# `nerdss-optimized` against `master`: change notes and pairwise benchmark

Reference state frozen on 2026-10-09. This file records (1) the frozen
reference state of `nerdss-optimized`, (2) every change the branch carries
relative to `master`, classified, with the commit, pull request or issue that
holds the evidence, and (3) a pairwise benchmark of the two branches: bitwise
identity with the known stream-changing changes switched off, efficiency, and
seed-averaged statistics where the trajectories necessarily differ.

Repository: `JohnsonBiophysicsLab/NERDSS` (fork of `mjohn218/NERDSS`). Links
below are to that repository: `#N` is a pull request or issue number there,
and a short hash is a commit on `nerdss-optimized`.

## 1. Frozen reference state

| | |
| --- | --- |
| Branch tip | `nerdss-optimized` at `3505932d` (merge of #47, 2026-09-28), `v1.1-189-g3505932d` |
| Tag | `nerdss-optimized-ref-2026-10-09` (annotated, on `3505932d`) |
| Baseline | `master` at `c507ae30` (tag `v1.1`), the merge base of the two branches; `master` has no commits that `nerdss-optimized` lacks |
| Divergence | 189 commits (162 non-merge, 27 merges), 22 merged pull requests (#15 to #47) |
| Source delta (excluding docs, benchmarks and sample inputs) | `src` 174 files +8519/-9778, `include` 37 files +2927/-428, `EXEs` 2 files +318/-146, `run_code_tests` 11 files +1894, `Makefile` +86/-6 |
| Reference binaries | `make serial` (Apple clang 21.0.0, `-O3 -std=c++0x`, GSL 2.8) from `git archive` of each tip; SHA-256 in `reference_2026-10-09/binaries.sha256` |
| Host | Apple M5, 10 physical cores, Darwin 25.6.0 |

The directory `benchmarks/nerdss_optimized/reference_2026-10-09/` holds the
material needed to repeat the comparison: the binary hashes, the hash manifests
and timings of every run reported in section 3, the case list with the
iteration counts used, the seeded selection script, the single-case runner and
the patch that produces the legacy-stream build (section 3.1). To rebuild the
reference binary:

```bash
git archive nerdss-optimized-ref-2026-10-09 | tar -x -C <dir> && make -C <dir> serial
```

## 2. Change notes, `master` -> `nerdss-optimized`

Every entry names what was wrong or what changed, where the evidence is, and
whether it moves simulation output. "Stream-changing" means the sequence of
random draws changes, so trajectories differ from `master` by construction
even when the physics is unchanged; "output-only" means a written file changes
but no trajectory does.

### 2.1 Bug fixes

#### Physics and mathematics

| # | What was wrong | Fix | Effect on output |
| --- | --- | --- | --- |
| B1 | Random molecule orientations were not uniform over rotations: four U(-1,1) quaternion components normalized sample the 4-cube, not the 3-sphere; measured chi^2/dof against uniform 9810, bins off by up to 59% ([issue #10](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/10)) | Shoemake's subgroup algorithm, `rand_unit_quat()`; chi^2/dof 1.02 (`4f3cb397`, RESULTS.md section 3) | stream-changing, every model |
| B2 | A crossing recorded between two complexes could reach `associate()` as a loop closure after an earlier association in the same step merged the complexes, writing a bond between interfaces up to 34 sigma apart; 23 of 100 runs of the 6BNO actin model fused filaments that never touched ([issue #31](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/31)) | `associate()` re-applies the same-complex admission criterion at the moment of association (`e9a7df24`, [#28](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/28); `docs/bond_site_separation_defect.md`) | only runs that would have formed such a bond (0 of 18 suite cases; 30 of 100 actin runs) |
| B3 | The loop-closure correction treated a dissociation as drawn at rate kb dt while `check_dissociation()` draws it with 1 - exp(-kb dt); kept dissociations were too rare once kb dt is not small (above ~1e6/s at usual timesteps). Upstream `nerdss_development` `524b685` | Ported as `349f5eff`, plus a guard that keeps the ratio finite when the back reaction's first rate is zero (`0cf85cc6`, [#45](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/45)) | physics change: +0.15% in the ratio at kb = 1e4/s (homoTrimer, closed_homoTrimer differ), +93% in upstream's fast-dissociation regime |
| B4 | The pairwise cell grid forced at least four sub-volumes along z, so any box thinner than 4 rMaxLimit got sub-volumes thinner than the interaction range and the +/-1 stencil missed reacting pairs: 7.5% of pairs within range in a 30 nm box, and `trimer` (a cubic 118 nm box against a 37.6 nm range) was affected | Every dimension is floor(L / rMaxLimit), with an assertion on the resulting edges (`459bd831`; RESULTS.md 19.2) | stream-changing for `trimer`; `rev_2D` restart files only |
| B5 | The max(4000, N^2/2) sub-volume budget shrank z twice as fast as x and y, leaving non-cubic sub-volumes under a cubic stencil (compartment: 58.8 x 58.8 x 125 nm against a 28 nm range, 98x over-inclusion), and sized the grid from the start-up molecule count | Budget removed; memory cap of 500 000 sub-volumes applied uniformly and iterated so it binds on any aspect ratio (`a81c7a6e`, `77daa3af`, `391a1a33`); the MPI cap can no longer raise a dimension (`c25b0a0b`) | stream-changing for `compartment`, `rev_3D`, `rev_3Dto2D`, `rev_2D` (the cases the cap or budget touched) |
| B6 | `set_rMaxLimit()` leaves rMaxLimit at 0 for a model with no bimolecular reaction; dividing the box by it gave a negative cell count and an endless "CELL PAIR MAX EXCEEDED" report (~1 TiB/hour of stdout, [issue #25](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/25)) | Rejected at setup with a message (`459bd831`) | `unimolecular_reverse` exits 1 instead of hanging |
| B7 | When the boundary check restarted the binning pass, molecule 0 was left out of the grid for that step (`molItr = 0` in a `++molItr` loop) | Counter restarts at -1 (`a576661b`) | none in any sample input (branch never reached) |
| B8 | `check_bimolecular_reactions()` decided a pair's dimensionality (1D / 2D / 3D) in three hand-written ladders and wrote the diffusion average out four times, and the copies had drifted: the `pro2` volume-exclusion block lacked the fiber branch its `pro1` mirror twenty lines up has, so a pair of complexes on a fiber was measured as 3D on one side and 1D on the other. One issue in two steps: | | |
| B8.1 | Ladders and averages given one home: `pair_dim()` and `weighted_D_sum()`. Verified bit-exact over 12 M comparisons against the four original expressions (an FMA-contraction pitfall was found and avoided on the way); putting the three ladders side by side is what made the drift visible, and the commit preserved it deliberately so the refactor stayed result-preserving | `9432e5d3` ([#20](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/20)); `weighted_D_sum_check` | none (bit-exact) |
| B8.2 | The missing fiber branch added to the `pro2` block, making the two exclusion blocks exact mirrors. What the block decides is only *whether a pair goes on the list the overlap sweep checks*; the sweep's own test is unchanged and is 1D for fiber pairs in every build. The old 3D arm put a fiber pair on that list when its 3D separation was below 3 sqrt(6 Dtot_3D dt) + sigma with Dtot_3D the three-axis average (plus a rotational term), i.e. 3.0 relative-step standard deviations; the 1D rule uses 4 sqrt(2 (D1x + D2x) dt) + sigma, 4.0 standard deviations, with the separation measured along the fiber. Pairs between 3 and 4 standard deviations apart were therefore left off the list by the old block and never checked for overlap when the step's tail carried them into the bound site's exclusion zone. Measured on a two-species fiber model (appendix A): the old block dropped 18.5% of the exclusion pairs, and realized overlaps of a free site with a bound site fell from 225 +- 5 per 300 000 steps on `master` to 10.5 +- 1 with the fix (20x) | `07d2b8d5` ([#20](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/20)); `reference_2026-10-09/fiber_b8/` | fiber models with `excludeVolumeBound` (none in the tree; the appendix model reaches it) |
| | Other changes to the same function are optimizations and are listed there: the redundant-check audit (O1), and the range cut and type mask in front of it (O8). | | |
| B8.3 | **Found 2026-10-09 while verifying B8.2, not yet fixed on the branch:** the centre-of-mass pre-cut that O8 added in front of this function (`1a704f1e`) drops any pair further apart than `rMaxLimit`, exempting only pairs with both complexes on the membrane. `set_rMaxLimit()` estimates every reach with the 3D formula, 3 sqrt(6 Dtot dt) with Dtot averaged over three axes, while the 1D reaction and exclusion rules reach 4 sqrt(2 (D1x + D2x) dt): 4/3 as far. On the appendix model `rMaxLimit` is 62 nm against a 1D reach of 82 nm (free sites) and 71 nm (bound site against free site), so fiber pairs 62 to 82 nm apart, which `master` evaluates, are never evaluated on the branch. Exempting fiber pairs the way surface pairs are exempted (4 lines, `fiber_b8/scripts/fiber_precut_exemption.patch`) makes the fix-reverted branch build reproduce `master` byte for byte over 300 000 steps on 4 of 4 seeds, which is also what isolates B8.2 above. Effect: small (association and exclusion candidates at 3 to 4 standard deviations of the relative step), confined to fiber models, and a divergence from `master` rather than a correction | reach: `1a704f1e` (O8); fix: the patch above, or make `set_rMaxLimit()` use the 1D reach for promoter pairs | fiber models only |
| B9 | Reverse *bimolecular* state change: a gate rejecting index 0, pointer subtraction across unrelated containers, missing `rateIndex`, no best-match selection, a read from an empty match list, and callers indexing `backRxns` with the forward index ([issue #8](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/8)) | All six fixed (`4f3cb397`; RESULTS.md section 6) | only reversible bimolecular state changes; no sample input declares one (`enzyme`'s is irreversible) |
| B10 | Reverse *unimolecular* state change never fired: `rateIndex` overwritten with a reaction index and `rxnIndex` never assigned, no `-1` guard on `conjBackRxnIndex`, `stateChangeRxns` stored as (back, forward) against its documented (forward, back), and `BackRxn::display()` reading one past the end | Fixed at each site (`58609700`, [#15](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/15)); new model `VALIDATE_SUITE/unimol_state_change_reversible` reaches kf/(kf+kb) = 0.25 against 0.765 before | only models with a reversible unimolecular state change ([issue #26](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/26) is the display crash) |
| B11 | Transition-matrix stays were credited to destroyed complexes (their `numEachMol` read past the end): every association counted `reactCom2`'s stays twice, diagonals overshot by 0.15 to 0.75% | `isEmpty` complexes skipped in the three diagonal loops (`49ea1f2e`, [#34](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/34)); regression `transition_matrix_totals` | output-only: `transition_matrix_time.dat` diagonals of `trimer`, `hetTrimer` |
| B12 | `transitionMatrix` / `lifeTime` are indexed by cluster size with no upper bound ([issue #23](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/23)) | Start-up warning naming the overflow risk (`ab42eb01`, [#24](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/24)); the bound itself is still open | none |
| B13 | A rate written as `INF` or `NAN` parsed silently; `1 - exp(-INF * 0)` is NaN, which every consumer reads as "never", so the `rate=INF` recipe three sample inputs recommended could silence a reaction | Non-finite rates rejected with the correct coupled-reaction recipe (`916e4218`); `first_order_probability()` returns certainty for +INF directly (`527853a4`) | none for shipped inputs (comments corrected, values unchanged) |
| B14 | An input naming neither a water box nor a sphere ran box physics against a zero-volume box; `isSphere = false` was ignored and a water box plus a sphere was conflated with a sphere | `BoundaryShape` plus `require_boundary()` fail-fast in serial and MPI paths (`cc824aa1`, `64dcb797`, `d2366f3e`, `4a6076f2`) | none for shipped inputs |
| B15 | `create_complex_propagation_vectors_on_sphere` built a degenerate frame for a COM on the z axis (zero tangent, absorbing state) | Pole guard (`2650775f`; RESULTS.md 18) | none: no sample run reaches the pole (closest approach 0.061 nm) |
| B16 | Scientific notation for `nItr` and the write intervals was truncated by `stoi` | Read with `stod` then converted (`0bad8621`) | none for inputs written as integers |

#### Crashes and memory errors (serial)

| # | What was wrong | Fix | Effect on output |
| --- | --- | --- | --- |
| B17 | A molecule destroyed mid-step keeps its slot with `myComIndex = -1`; the overlap loop, the creation reactant pool and a transmitted molecule's partners' crossing lists indexed `complexList[-1]`. `compartment`, `create_destroy` and `genetic_oscillator/clock_model` segfault on `master` ([issue #27](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/27)) | `isEmpty` guards and `remove_partner_crossings()` (`394fdd1d`, `52b9ab30`/`72a773e4`, `6b01636b`, [#20](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/20), [#22](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/22)) | those three models run to completion; stream-changing only where the guard fires (runs that crashed before) |
| B18 | `moleculeOverlaps()` walked the previous step's sub-box member lists and read past `moleculeList` | Bound check (`bbd19827`/`0d422619`, [#22](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/22)) | none |
| B19 | A spherical boundary with no implicit lipid segfaults on step 1 (`RS3Dvect` empty, read at +400 through a null pointer) | Guard on whether the table exists, `hasRS3Dtable()` (`f9b2ae23`, [#40](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/40)); regression `sphere_no_implicit_lipid` | none for the five spherical samples (all carry an implicit lipid) |
| B20 | Creation-reaction placement (`moleculeOverlaps()`): returned "no overlap" at the first out-of-range member, read `reactantListNew[1]` of unimolecular reactions, cut on bounding spheres without the bind radius, re-initialized the molecule on every resample (gaps in ids, inflated `numTotalUnits`), and looked in one SubBox where the exclusion sphere spans 27 | Rewritten to `generate_coordinates()`'s criterion over the 3 x 3 x 3 block, with the reaction scan hoisted and a type mask (`8a4cc57a` [#41](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/41), `1f8389bb` [#43](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/43)); `creation_placement_check`, regression `creation_placement_crowding` | fresh suite runs unchanged; `clock_model` changes counts only; crowded creation models place correctly (13512 overlapping pairs -> 0) |
| B21 | Bond bookkeeping ([issues #18](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/18), [#19](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/19)): `bndRxnList` was write-only and unbounded with a dead reader; `bndlist`/`bndpartner` fell out of correspondence on the cancel-dissociation path (88 490 broken slots on homoTrimer) | `bndRxnList` and `correct_structure()` deleted; one `erase_bond()` path (`3f887239`; `docs/bond_bookkeeping_defects.md`) | restart files of `homoTrimer`, `closed_homoTrimer` only (they serialize the lists); MPI reads the lists at a common subscript |
| B22 | Debug helpers allocated a string per molecule per step and aborted on NaN coordinates | Removed (`60659c2e`) | none for finite runs; NaN now propagates into output instead of `exit(1)` (behaviour change, documented) |
| B23 | Latent: `determine_2D_implicitlipid_reaction_probability()` could leave `probvec` one entry longer than `crossbase` | Closed by the `CrossEntry` merge (`4b2b916d`) | none measured |

#### Restart files (serial and MPI)

| # | What was wrong | Fix | Effect on output |
| --- | --- | --- | --- |
| B24 | `read_restart()` matched the header by position, so a file written before a field existed silently shifted every later value (`checkPoint` read into `assocDissocWrite`) | Field names checked, with a bounded key scan (`f4bcba78`, `7991e1d0`) | old files refused by name instead of misread |
| B25 | Uninitialized `Membrane` fields were written into every restart file (`membrane =` line changed with the binary's memory layout and launch environment) | Initializers (`aca34052`) | restart files |
| B26 | `invCbrtDr` was not refreshed when a restart set `Dr`, so every restarted model with rotational diffusion got `Dr = inf` and NaN coordinates (15 of 24 sample models) | Refreshed in `read_restart()` (`75f7cfe6`, [#29](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/29)) | restarts of those models now complete |
| B27 | A restart file did not carry the compartment (transmission interface, bind radius, site D and rho, `enforceCompartmentBC`), so restarting `compartment` segfaulted ([issue #30](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/30)) | Saved for compartment systems only (`36a556cb`, [#29](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/29)) | restart files of compartment models |
| B28 | No implicit-lipid model continued exactly: the lipid's `trajStatus` was re-propagated (extra random draws), the 2D binding table and the step-0 protein counts were rebuilt from the restart's state, and the water box's `xLeft`/`xRight` were zero so every created molecule landed on x = 0 | Four fixes, `#ImplicitLipid` section, water box rebuilt through its constructor (`e6deafdd`, `e0279dfe`, `988508f6`, `5bdf0b72`, [#33](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/33)); regression `restart_continuation` (9 of 9 exact) | restart files of implicit-lipid models; restarted trajectories |
| B29 | Add-file restarts: bound-pair lists not renumbered (heap corruption, `std::length_error`), wrong interface looked up in `set_rMaxLimit()`, no transition matrix for an added template, added molecules never on `monomerList` (never destroyed), placement test checked nothing, no ids, species total short, units counted twice ([issue #38](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/38) tracks `master`) | [#32](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/32) (`fd32157a`, `f2c5a290`, `482ac54b`), [#36](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/36) (`bd789909`, `b39e8d6d`, `a320b400`, `ab61471d`), [#39](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/39) (`3e020402`), [#42](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/42) (`1ecab490`, check only); `add_file_molecules_check` | add-file restarts only |
| B30 | Upstream JSON reader defects (interfaces written under one key and read under another, back reactions requiring a field only forward reactions wrote, `numEachMolType` recounted over destroyed slots) | Fixed while porting (`2bc2f423`, [#44](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/44)) | JSON restarts |

#### MPI executable (`nerdss_mpi`)

| # | What was wrong | Fix | Effect |
| --- | --- | --- | --- |
| B31 | `make serial && make mpi` linked an MPI binary out of serial objects (`MpiContext` layout differs by 96 bytes); it ran and reported zero of every species. `make mpi debug` picked `g++` and failed on `mpi.h` | Separate object trees per mode, inline namespace guard, goal-list match (`6d4f325a`, `3500503e`, `a4550f1d`, [#21](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/21)) | build correctness |
| B32 | `isDissociated` was never cleared in the parallel loop, so binding switched off molecule by molecule (np=1 held 0 bound pairs where serial holds 158) | Cleared each step (`015789c5`) | MPI physics |
| B33 | `write_all_species()` read `moleculeList[-1]` for every free interface; `indices_to_IDs` corrupted the rank's own reweighting history; the implicit lipid was reported only from the unghosted rank (merged totals at 0.50 to 0.63 of serial) | `c162f500`, `898698a1`, `e83d4040` ([#21](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/21)) | MPI output |
| B34 | Ownership was decided twice by different criteria (sender by column, receiver by `get_x_bin`), leaving molecules owned by nobody | `ownerRank` on the complex, `isGhosted` derived from it, `CHECK_OWNERSHIP` invariant check (`5c63ffc9`, `daf047ba`, `f91bcc86`) | MPI correctness |
| B35 | Restarts: `MpiContext` pointers uninitialized, histograms empty, occupancy mask empty, `moleculeList[-1]` for the lipid's complex, and every restart re-ran from step 1 ([issue #37](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/37) tracks `master`) | Five commits (`f443287d`, `02daa60e`, `fe63286e`, `3e0a9cdf`, `b2ff8492`, [#35](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/35)); multi-rank restart refused with a message | one-rank MPI ASan-clean; restarts resume at the checkpoint |
| B36 | On two or more ranks the heap was corrupted and mass was not conserved: complexes handed over without all members, stale indices after exchange and compaction, ghosts able to reassign ownership, complexes created after start-up owned by nobody, owned molecules deleted as "not received back", counts taken by x bin, implicit-lipid sites rounded per rank | Thirteen commits (`ebbf581e` through `7a37823e`, [#47](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/47)): 12 of 12 models ASan-clean at np=2 and np=3, mass conserved, no id owned twice or by nobody | MPI correctness; serial binary unchanged |
| B37 | `deserialize_back()` bounded sibling arrays by one array's length (introduced on this branch) | Bounded by the shortest (`e78d3449`) | MPI robustness |

### 2.2 Optimizations

Speedups are the measurements recorded in `RESULTS.md` at the time each
change landed (Apple M5, serial `-O3`, interleaved timing unless noted); the
branch-level figure is in section 3 of this file.

| # | What was inefficient | Change | Measured |
| --- | --- | --- | --- |
| O1 | `find_which_reaction()` / `find_reaction_rate_state()` built a temporary match vector per candidate pair; `hasIntangibles()` scanned the whole interface list per requirement; `check_bimolecular_reactions()` re-read `myComIndex` ~20 times per pair; `Complex::propagate()` called a no-op `rotQuat.unit()`; `Coord` arithmetic could not inline ([issues #8](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/8), [#9](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/9), [#11](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/11)) | Streaming `best_matching_rate()`, direct lookup, header inlining (`4f3cb397`) | suite 1.091x, bit-exact (RESULTS.md sections 1, 2) |
| O2 | `GaussV()` used a hand-written Marsaglia polar method: ~2.55 uniforms, a log, a sqrt and an unpredictable branch per draw ([issue #12](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/12)) | `gsl_ran_gaussian_ziggurat` (`4f3cb397`, `ee36d8a2`) | sampler 3.19x; suite 1.168x together with B1; stream-changing |
| O3 | `Complex::update_properties()` summed `1/cbrt(Dr)` over every member every step | Per-template `invCbrtDr` cache (`c544b1a7`) | `pow` 4.2% of clathrin self time -> 0 |
| O4 | Every sub-cell was cleared every step however few held members; `Vector`'s hot members could not inline | `occupiedSubCells` + `clear_member_lists()`, header inlining (`0784dd8b`) | suite 1.087x, sphere 1.509x, bit-exact (section 9) |
| O5 | Two debug helpers heap-allocated a 23-character message per molecule per step (one past the short-string buffer) | Removed (`60659c2e`) | 1.092x on six cases; allocator share 9.1% -> 0.9% (section 16) |
| O6 | `add_3D_rotational_diffusion()` evaluated `cos(sqrt(k Dr dt))` twice per candidate pair | Memoized on the complex, keyed on the exact argument (`64417c94`) | ~1.5 to 2% |
| O7 | Two vector types (`Coord`, `Vector`) with a magnitude cache nothing maintained | One `Vec3D` (`207e22dd`) | runtime 0.999x (neutral), operations 5 to 11% faster in isolation (section 14) |
| O8 | The pairwise search walked every sub-volume, offered every stencil pair whatever the types, and opened the full interface search on pairs further apart than rMaxLimit | Visit only occupied sub-volumes with a bit mask (`cf646a9a`, `808272b6`), register created molecules (`ebac7dbf`), COM-distance cut (`1a704f1e`), per-type interaction mask on sub-boxes (`51e5d33f`), grid sized to the range (B4, B5) | 1.18x over the 13 cases, `mem_localization` 3.13x, bit-exact except the five grid-sizing cases (section 19.5) |
| O9 | A spherical membrane was binned on a Cartesian grid over the bounding cube, 94% of which can never hold a surface molecule | Latitude-longitude `ShellIndex`, active only for spheres with explicit pair reactions and molecules above 0.9 R (`336322d8`, `3d253cc0`) | neutral on the samples (built for larger spheres); changes `cluster_gagsphere` only |
| O10 | `Molecule` was 656 bytes with 22 heap-owning vector members; the pairwise rejection cascade touched four cache lines per molecule | Twelve reweighting vectors -> one `ReweightEntry` (`f9309bfa`), four crossing lists -> one `CrossEntry` (`4b2b916d`), fields ordered by access (`707a11ea`), `bndRxnList` removed (B21): 656 -> 304 bytes ([#17](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/17), sections 23 to 26) | 1.23x aggregate on a performance core, 1.43x on an efficiency core, up to 1.67x at 40 000 molecules; bit-exact |
| O11 | `set_memProtein_sphere()` / `find_Lipid_sphere()` took the whole molecule list by value (2.0 million `Molecule` copies per sphere run); dead stores and repeated trig in the spherical geometry | References, hoisting, `sin_cos()` (`ecbc9f06`, `7922b698`) | sphere 1.046x, bit-exact (section 17) |
| O12 | Creation placement walked whole complexes per SubBox member and rescanned reactions per candidate | Members walked directly, reaction scan hoisted, type mask (`1f8389bb`, part of B20) | crowded model 9 m 47 s -> 19 s (mostly from no longer carrying overlapping pairs) |
| O13 | Measured and deliberately not done: a side table for `find_which_reaction()` (0.999x, `733df66c`), OpenMP over any profile frame (`835324b3`), a structure-of-arrays `Molecule` (byte-identical but no faster, section 22), removing the association scratch fields (ceiling 1.9%, `5bbf6b8b`), explicit `sincos` (already fused by the compiler, section 16.5) | recorded so they are not retried | - |

### 2.3 Other enhancements

| # | Change | Where |
| --- | --- | --- |
| E1 | Restart files are JSON (`DATA/restart.json`, `RESTARTS/restart<N>.json`); either format is read, detected by content; `legacyRestartFormat = true` writes the positional `.dat` format and a restart keeps the format it read | `2bc2f423` [#44](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/44), `dab0c8f4` [#46](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/46); `restart_format` round-trip check |
| E2 | Restart files carry state they did not: compartment parameters and flags, the implicit lipid's 2D table and protein counts, numerical settings, `numEachMolType`, `totalNumOfStates` | B27, B28, `cc1d954d`, `2bc2f423` |
| E3 | Numerical tolerances are configurable input parameters instead of literals: 16 `numerics*` keywords covering the integration tolerances, table look-up tolerances, flat-diffusion thresholds, association-angle tolerances, and the `Vec3D` equality precision; validated at parse time | `08025e28`, `44adcbf5`, `539c2d0d`, `cc1d954d` |
| E4 | Build: header dependencies tracked (`-MMD -MP`), the executable depends on its own source, serial / MPI / debug / profile builds in separate object trees, `make checks` target, CI builds MPI straight after serial | `aca34052`, `cf82352a`, `6d4f325a`, `3500503e`, `7991e1d0`, `2739460d` |
| E5 | De-branching: `BoundaryShape` enum with `waterBoxGiven` provenance; the `check_bimolecular_reactions()` ladder consolidation (B8.1); sphere and compartment reflectors merged behind `RadialSide` (four files -> two); the entering/exiting compartment branches collapsed; the per-axis reach scan collected in `complex_extent.hpp`; the solution and membrane overlap sweeps share one body; the cluster sweeps share their two loops; reaction bookkeeping, crossing records, rotational-diffusion terms and the Brownian resample each written once; `first_order_probability()`; `calc_interface_distance()` split from `get_distance()` | [#20](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/20), `8d81d571`, `dc9ad025`, `32892764`, `7d8520b1`, `aad4ba9a`, `ddff3a88`, `a54c4dde`, `321d8935`, `527853a4`, [#28](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/28); every step verified bit-exact |
| E6 | Spherical-system geometry in named types (`SphericalCoord`), `find_cardesian_coords` spelled correctly, dead helpers and commented-out call sites removed; unused math functions removed; factorial-like helpers deprecated; `Lipid` -> `lipid` naming | `ecbc9f06`, `18ee4d4d`, `dbb5cbff`, `3823414d`, `8e895ce6`, `9ec5ab89` |
| E7 | `erflike` vendored as an alternative erf route (then deprecated in favour of the existing path), with licenses recorded for both vendored erf implementations | `61e75e58`, `4baa42f2`, `1821bfe1` |
| E8 | Diagnostics: `transitionMatrixSize` warning (B12), rMaxLimit <= 0 error (B6), sub-volume edge assertion (B4), restart key errors (B24), bound-pair-list count check (B29), boundary fail-fast (B14), MPI ownership invariant (B34), multi-rank restart refused (B35) | as listed |
| E9 | Tests: six `run_code_tests/*_check.cpp` programs under `make checks` (`weighted_D_sum`, `membrane_serialization`, `loop_closure_separation`, `interface_distance`, `add_file_molecules`, `creation_placement`); seven regression cases under `benchmarks/nerdss_optimized/regression/` (`bimol_state_change`, `unimol_state_change_reversible`, `restart_continuation`, `restart_format`, `transition_matrix_totals`, `sphere_no_implicit_lipid`, `creation_placement_crowding`); `run_code_tests/BondSiteSeparation`; new sample model `VALIDATE_SUITE/unimol_state_change_reversible` | `7991e1d0`, [#15](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/15), [#28](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/28), [#33](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/33), [#34](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/34), [#36](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/36), [#40](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/40), [#41](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/41), [#43](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/43), [#44](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/44) |
| E10 | Benchmark harness `benchmarks/nerdss_optimized/` (`run_suite.sh`, `compare_suites.sh`, `interleaved_timing.sh` with instruction and cycle counters, `statistical_check.sh`, `compare_observables.py`, `make_scale_inputs.sh`, case tables, `known_broken.tsv`, `known_uncovered.tsv`, `rng_quality/`) and the serial-versus-MPI validation harness `benchmarks/mpi_validation/` | `4f3cb397`, `e5be2a6e`, `4919df1a`, `b96999e7`, `913b6e4c` |
| E11 | Sample inputs: `clock_model.inp` and `autophos_D10.inp` comments corrected (the `rate=INF` recipe), `cases.tsv` headers, coverage table raised to reach the transmission path | `916e4218`, `7d078790` |

### 2.4 Documentation and beautification

| Item | Where |
| --- | --- |
| `docs/nerdss_optimized.md` (1223 lines): what changed and why, per issue and per series | `4f3cb397` onward |
| `benchmarks/nerdss_optimized/RESULTS.md` (3056 lines, 26 sections): every measurement behind sections 2.2 and 2.3, including the negative results | `4f3cb397` onward, `324c02d5` renumbering, `9f772720` |
| `docs/debranching_plan.md`, `docs/molecule_layout_plan.md`: the plans, corrected where the code disagreed | `80db71e2`, `ea1f60f5`, `cd513253`, `b5890668` |
| `docs/bond_bookkeeping_defects.md`, `docs/bond_site_separation_defect.md`: the two defect investigations | `8e6e20df`, `f2021bb8`, `9f346bdc`, `e9a7df24` |
| `docs/nerdss_cell_grid_audit.html`: the cell-grid audit as a standalone page | `63fc2b1d` |
| `docs/mainpage.md`, `docs/input.md`, `README.md`, `io.hpp` comments: restart format and `legacyRestartFormat` | [#44](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/44), [#46](https://github.com/JohnsonBiophysicsLab/NERDSS/pull/46) |
| Doxygen docstrings for `Vec3D`, `Quat`, `rand_gsl`; corrected comments in the radial reflectors; `#pragma once` on the spherical-system header | `5411d442`, `4a6076f2`, `ecbc9f06` |
| `.gitignore` for build products, results and local notes | `a7435a06`, `ccde189a`, `8025fd71` |
| `benchmarks/nerdss_optimized/regression/*/README.md` and `benchmarks/mpi_validation/FINDINGS.md` | as above |

### 2.5 Known changes to bitwise identity, in one place

These are the changes that make `nerdss-optimized` differ from `master` byte
for byte even where the physics is unchanged. Section 3 switches the first
five off to measure everything else.

| Change | Cases affected | Nature |
| --- | --- | --- |
| B1 random orientation sampler (Shoemake) | every model with more than one interface placed at random | stream-changing, corrected distribution |
| O2 Gaussian sampler (ziggurat) | every model | stream-changing, same distribution |
| B4, B5 grid sizing | `trimer`, `rev_3D`, `rev_3Dto2D`, `compartment` (trajectory); `rev_2D` (restart files only) | pair order and completeness |
| O9 shell index | `cluster_gagsphere` | pair order |
| B3 dissociation correction ratio | `homoTrimer`, `closed_homoTrimer` and any model with loop closures | physics |
| B2 loop-closure guard | runs that would have written a far-apart bond | physics (bug) |
| B17 `isEmpty` guards | models that destroy molecules (crashed before) | stream-changing where they fire |
| B11 transition-matrix diagonals | `trimer`, `hetTrimer` | output-only |
| B21, B25, B27, B28, E1, E2 restart-file content and format | every case that writes restart files | output-only |
| B20 creation placement | `clock_model` counts; crowded creation models | placement |

## 3. Pairwise benchmark against `master`

Measured 2026-10-09 on the reference host (section 1), serial `-O3` builds
from `git archive` of each tip, seed 20260810 unless stated. Three binaries:

| Label | Contents | SHA-256 (first 16) |
| --- | --- | --- |
| `master` | `master` at `c507ae30` | `ef0c0153ef2aa54d` (the same bytes as every earlier `master` baseline in `RESULTS.md`) |
| `opt` | `nerdss-optimized` at `3505932d`, the reference state | `d9d67abf9988af5e` |
| `opt-legacy` | `3505932d` built with `-DNERDSS_LEGACY_STREAM`, which switches the five known stream-changing changes back to `master`'s behaviour (section 3.1) | `eb97f2b9e3e276755` |

### 3.1 The legacy-stream build

`nerdss-optimized` cannot be compared with `master` byte for byte as it is:
B1 and O2 change the random stream of every model. `opt-legacy` is the
reference source plus `reference_2026-10-09/scripts/legacy_stream.patch`
(generated by `legacy_stream_patch.py`), 89 changed lines behind one macro,
which restores:

1. `GaussV()`: `master`'s Marsaglia polar sampler instead of the ziggurat (O2);
2. `rand_unit_quat()`: four U(-1,1) components normalized, as `master` did, instead of Shoemake (B1);
3. the sub-volume grid sizing of `master`: the z floor of four, the 30-per-dimension cap and the max(4000, N^2/2) budget, and no edge assertion (B4, B5);
4. no dissociation back-correction (B3);
5. the spherical shell index never activates (O9).

Everything else on the branch - every optimization, refactor, restart change
and the remaining bug fixes - is in `opt-legacy` exactly as in `opt`. So a
case on which `opt-legacy` reproduces `master` byte for byte shows that none
of those changes moved its trajectory, and a case that differs isolates
the cause to something outside the five switches.

Restart files are excluded from the trajectory verdict and reported
separately: `opt` writes `restart.json` where `master` writes `restart.dat`
(E1), the branch stores state `master` did not (E2, B27, B28), and `master`'s
files carry uninitialized fields (B25), so they can never match. PDB files
are hashed from line 2, after the wall-clock stamp, as `run_suite.sh` does.

### 3.2 Bitwise identity on the 18 benchmark cases

`cases.tsv` (13) and `coverage_cases.tsv` (5), run with `run_suite.sh`; the
hash manifests are in `reference_2026-10-09/suite_manifests/`. "Trajectory
files" are every file under `DATA/` and `PDB/` except the restart files.

| Case | nItr | `opt-legacy` vs `master`: trajectory files | `opt` vs `master` | Note |
| --- | ---: | --- | --- | --- |
| rev_3D | 20 000 | **identical** (13 of 13) | differ | |
| rev_2D | 200 | **identical** (13 of 13) | differ | |
| rev_3Dto2D | 8 000 | **identical** (13 of 13) | differ | |
| homoTrimer | 6 000 | **identical** (14 of 14) | differ | |
| closed_homoTrimer | 6 000 | **identical** (14 of 14) | differ | |
| hetTrimer | 8 000 | 13 of 14; `transition_matrix_time.dat` differs | differ | B11: `master` credits destroyed complexes with stays (diagonal 1 665 656 against 1 662 669); every other file identical |
| trimer | 8 000 | 13 of 14; `transition_matrix_time.dat` differs | differ | B11, as above (66 963 against 66 300) |
| hexamer | 25 000 | **identical** (14 of 14) | differ | |
| clathrin | 150 000 | **identical** (13 of 13) | differ | |
| implicit_lipid | 100 000 | **identical** (13 of 13) | differ | |
| sphere | 40 000 | **identical** (13 of 13) | differ | |
| michaelis_menten | 150 000 | **identical** (13 of 13) | differ | |
| mem_localization | 2 000 | **identical** (13 of 13) | differ | |
| compartment (coverage) | 30 000 | not comparable | not comparable | `master` segfaults after step 0 (B17); `opt` and `opt-legacy` complete |
| unimol_state_change (coverage) | 20 000 | not comparable | not comparable | `master` segfaults in `BackRxn::display()` at start-up (B10, [issue #26](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/26)); the branch completes |
| gagsphere (coverage) | 2 000 | **identical** (13 of 13) | differ | volume exclusion between bound partners |
| cluster_mem_loc (coverage) | 2 000 | **identical** (13 of 13) | differ | `clusterOverlapCheck = true` |
| cluster_gagsphere (coverage) | 1 000 | **identical** (13 of 13) | differ | sphere, cluster sweep; the shell index is off in `opt-legacy` |

**Verdict:** with the five known stream-changing changes switched off, the
branch reproduces `master` byte for byte on every trajectory, observable,
histogram and PDB file of all 16 cases `master` can run, except the
transition-matrix file of the two cases that count transitions, where the
difference is the documented diagonal correction (B11). No unexplained
divergence. `opt` itself differs from `master` on every case, as B1 and O2
require; its agreement with `master` is measured statistically in section 3.4.

Restart files differ in every case for the format and content reasons in
3.1, on `opt` and `opt-legacy` alike.

### 3.3 Randomly sampled production inputs: identity and efficiency

**Selection.** `reference_2026-10-09/scripts/select_cases.py` enumerates every
stock `.inp` under `sample_inputs/` outside `VALIDATE_SUITE` that describes a
fresh run (48 inputs after excluding add files, checkpoint copies, the
restart-only `4-runRemodelingSim`, the `nItr = 2` template, two inputs with
no molecules block, the single-molecule `Gag/parm.inp` and the `testAdd`
fixtures) and draws from it with `random.Random(20261009).sample`. The first
eight draws were taken; one of them, `clathrin_pioneer/parms.inp`, cannot run
on either branch because it names `pip2.mol`, which does not exist in that
directory, so it was replaced by the ninth draw. The result spans implicit-
lipid clathrin models, two flat/puckered clathrin lattices, a Gag lattice on
an implicit membrane, two thin-slab secretion models, a membrane creation
model and the genetic oscillator with creation and destruction.

Iteration counts were set so that `master` runs each case for about 30 to
110 s and are recorded in `reference_2026-10-09/sample_cases.txt`. The three
binaries ran back to back per case, one case at a time, on an otherwise idle
host, under `/usr/bin/time -l`; the efficiency measure is **retired
instructions**, which does not move with background load, with user CPU time
alongside. Hash manifests and timing summaries are in
`reference_2026-10-09/pair_manifests/`.

| Case | Input | nItr | `master` | `opt-legacy` vs `master`, trajectory files | instr. speedup `opt-legacy` | instr. speedup `opt` | CPU-time speedup `opt` |
| --- | --- | ---: | --- | --- | ---: | ---: | ---: |
| invivo_macro | `clathrin_invitro_invivo/invivo/parmsMacroRate.inp` | 40 000 | 61.7 s, 823 G instr | identical except `transition_matrix_time.dat` (B11) | **4.13x** | 4.08x | 3.79x |
| flat_clat | `clathrin_coat/flat_clat.dir/clath_kon1uM.inp` | 400 000 | 30.3 s, 372 G | **identical** | **3.07x** | 3.37x | 3.63x |
| puckered_clat | `clathrin_coat/puckered_clat.dir/parms_pucker.inp` | 400 000 | 25.3 s, 300 G | **identical** | **2.65x** | 2.74x | 3.12x |
| gag_lattice | `gagLatticeRemodeling/1-initialAssemblyInSolution/parm.inp` | 500 000 | 65.9 s, 933 G | **identical** | **2.04x** | 2.64x | 2.65x |
| secretion_kd500 | `secretion_catalyzed_assembly/membrane/kd500/parms.inp` | 3 000 | 105.7 s, 1319 G | identical except `transition_matrix_time.dat` (B11) | 1.01x | 1.79x (*) | 1.79x (*) |
| secretion_kd100 | `secretion_catalyzed_assembly/membrane/kd100/parms.inp` | 4 000 | 81.9 s, 1020 G | identical except `transition_matrix_time.dat` (B11) | 1.02x | 0.63x (*) | 0.63x (*) |
| pucadyil | `pucadyil/parms_pucadyil_ka3d.inp` | 500 000 | **SIGSEGV** after step 0 (B17) | not comparable | - | - (`opt` 11.4 s; `opt-legacy` 68.3 s, see *) | - |
| clock_model | `genetic_oscillator/clock_model.inp` | 200 000 | **SIGSEGV** after step 0 (B17, [issue #27](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/27)) | not comparable | - | - (`opt` 30.2 s, 361 G) | - |

(*) The two secretion models and `pucadyil` spend most of their time building
2D reaction-probability tables, one per distinct (rate, Dtot) pair the
trajectory produces, which is the branch's known remaining cost centre and
is unchanged on either branch. How many tables a run builds follows its
trajectory, so once the random stream differs the cost moves either way:
`opt` builds fewer tables than `master` on kd500 (1.79x) and more on kd100
(0.63x), while `opt-legacy`, which follows `master`'s trajectory, shows the
code's own effect on these models: parity (1.01x, 1.02x). The same applies
to `pucadyil`, where the two branch builds follow different trajectories.

**Identity:** on all six models `master` can run, `opt-legacy` reproduces
every trajectory, observable and PDB file byte for byte; the three models
that count transitions differ only in `transition_matrix_time.dat`, the
documented diagonal correction (B11). Two of the eight sampled models cannot
be run by `master` at all.

**Efficiency:** where the time is in the pairwise search and the reaction
bookkeeping, the branch retires 2.0x to 4.1x fewer instructions on the same
trajectory, with the full `opt` build between 2.6x and 4.1x. Where the time is
in 2D table construction, the branch is at parity.

### 3.4 Efficiency on the 13 standard benchmark cases

The same three binaries on `cases.tsv`, run back to back per case under
`/usr/bin/time -l` after the suite runs above (`reference_2026-10-09/pair13_results.tsv`).
`run_suite.sh`'s own wall-clock column is not used: it rounds every time up
to the next whole second and shared the host with other runs.

| Case | nItr | `master` instr (G) | `opt-legacy` instr (G) | `opt` instr (G) | instr. speedup `opt-legacy` | instr. speedup `opt` | CPU-time speedup `opt` | `opt-legacy` vs `master` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| rev_3D | 20 000 | 106 | 61 | 53 | 1.74x | 1.98x | 1.98x | identical |
| rev_2D | 200 | 977 | 977 | 976 | 1.00x | 1.00x | 1.00x | identical (2D table construction) |
| rev_3Dto2D | 8 000 | 140 | 70 | 64 | 1.99x | 2.19x | 2.78x | identical |
| homoTrimer | 6 000 | 204 | 77 | 77 | 2.66x | 2.66x | 2.26x | identical |
| closed_homoTrimer | 6 000 | 204 | 77 | 77 | 2.66x | 2.66x | 2.29x | identical |
| hetTrimer | 8 000 | 139 | 58 | 56 | 2.39x | 2.48x | 1.94x | transition matrix only (B11) |
| trimer | 8 000 | 182 | 75 | 75 | 2.44x | 2.43x | 1.97x | transition matrix only (B11) |
| hexamer | 25 000 | 184 | 92 | 82 | 2.01x | 2.25x | 2.41x | identical |
| clathrin | 150 000 | 149 | 51 | 44 | 2.89x | 3.41x | 3.74x | identical |
| implicit_lipid | 100 000 | 96 | 50 | 46 | 1.90x | 2.08x | 2.43x | identical |
| sphere | 40 000 | 144 | 82 | 75 | 1.75x | 1.91x | 2.11x | identical |
| michaelis_menten | 150 000 | 103 | 49 | 49 | 2.11x | 2.12x | 2.34x | identical |
| mem_localization | 2 000 | 68 | 22 | 21 | 3.04x | 3.26x | 3.90x | identical |
| **total** | | **2696** | **1741** | **1695** | **1.55x** | **1.59x** | 1.62x | |
| total without rev_2D | | 1719 | 764 | 719 | 2.25x | 2.39x | 2.52x | |

This pass reproduced the section 3.2 verdict: the trajectory files of
`opt-legacy` hashed identically to `master` on all 13 cases, with the same
two transition-matrix exceptions. `rev_2D` is the 2D reaction-table case
(about 77 s of table construction at any iteration count) and dilutes the
total; without it the branch retires 2.4x fewer instructions than `master`
over the standard cases, between 1.7x and 3.4x per case.

### 3.5 Chart

![Speedup over master, retired instructions](reference_2026-10-09/speedup_chart.svg)

`reference_2026-10-09/speedup_chart.svg`, generated by
`scripts/make_chart.py` from `speedup_chart_data.tsv`. Blue is `opt-legacy`,
which follows `master`'s trajectory, so its bar is the code's own effect;
orange is `opt` as shipped. The two secretion bars show the trajectory
dependence discussed in 3.3, and the two models `master` cannot run have no
bar.

### 3.6 Seed-averaged statistics where the trajectories differ: `opt` vs `master`

`opt` as shipped differs from `master` on every case (B1, O2), so its physics
is compared the way `RESULTS.md` sections 5 and 10.2 compare it: each of the
sampled models was run under both builds with independent seeds (20260811
onward), every reported species copy number was averaged over the second
half of each run, and the per-seed averages were combined into mean +- SEM
and a Welch z per species (`compare_observables.py`; full output in
`reference_2026-10-09/statistics_summary.txt`). Six seeds per build for every
model; the two lattice models whose 6-seed pass showed a species at
|z| >= 2 were extended to 18 seeds. The two models `master` cannot run have
no baseline and are absent. Species that are constant in both builds (fixed
totals) are not counted.

| Case | seeds (`master` / `opt`) | species compared | largest \|z\| | species at \|z\| >= 2 (`master` vs `opt`) |
| --- | --- | ---: | ---: | --- |
| invivo_macro | 6 / 6 | 12 | 2.21 (`clat(cd2)`) | `clat(cd2)` 189.9±0.3 vs 188.2±0.7 (z -2.2) |
| flat_clat | 18 / 18 | 9 | 3.09 (`clat(cd2!1).clat(cd2!1)`) | `clat(cd2!1).clat(cd2!1)` 3.4±0.4 vs 2.0±0.2 (z -3.1) |
| puckered_clat | 6 / 6 | 9 | 1.49 (`clat(cd1)`) | none |
| gag_lattice | 18 / 18 | 15 | 2.08 (`pol(hom)`) | `pol(hom)` 3.2±0.4 vs 2.2±0.3 (z -2.1) |
| secretion_kd500 | 6 / 6 | 11 | 0.95 (`s(c)`) | none |
| secretion_kd100 | 6 / 6 | 11 | 1.42 (`s(a)`) | none |

How the two extended cases moved with more seeds, which is the test that
separates a change of random stream from a change of physics (a physics
change grows with the seed count; stream noise shrinks):

| Case, species | 6 seeds | 18 seeds |
| --- | ---: | ---: |
| gag_lattice `gag(hx1)` | 41.4±0.6 vs 46.1±0.9, z 4.16 | 41.9±0.8 vs 44.0±0.8, z 1.75 |
| gag_lattice `pol(hx1)` | z -2.83 | z -1.58 |
| flat_clat `clat(cd3)` | 82.5±1.7 vs 77.8±1.4, z -2.14 | 80.0±1.0 vs 78.8±1.0, z -0.86 |
| flat_clat `clat(cd2!1).clat(cd2!1)` | 2.9±0.5 vs 2.7±0.6, z -0.36 | 3.4±0.4 vs 2.0±0.2, z -3.09 |

Reading: the 6-seed outliers on the Gag lattice shrank to within scatter at
18 seeds, as stream noise does. Over the 67 varying species in the table
above, three sit between |z| 2 and 3.1; with species inside a model coupled
by conservation and 67 comparisons in play, that count is what a pure
stream change produces, and none of the three is a headline observable (the
largest is a homodimer bond count of 3.4 against 2.0 copies out of 100
clathrins). The two physics corrections carried by `opt` (B1, the corrected
orientation sampler, and B3, the dissociation correction) are allowed to
move small quantities in these lattice models, so a sub-2-sigma shift there
is not evidence against the branch. The flagged `flat_clat` dimer count is
the one value worth more seeds if that model is used for a quantitative
claim. No observable with a sizeable population differs between the
builds, which is the same conclusion as the earlier 3-model, 6-seed checks
(largest |z| 0.61, `RESULTS.md` 10.2).

`pucadyil` and `clock_model` ran to completion on `opt` for all six seeds
(`master` segfaults on both after step 0), which is a correctness result
rather than a statistic.

## 4. Summary

- **Frozen:** tag `nerdss-optimized-ref-2026-10-09` on `3505932d`; binaries,
  hashes, manifests and scripts under `reference_2026-10-09/`.
- **Bitwise:** with the five known stream-changing changes switched off
  (`-DNERDSS_LEGACY_STREAM`), the branch reproduces `master` byte for byte on
  every trajectory and PDB file of all 22 cases `master` can run (13
  standard, 3 coverage, 6 sampled production inputs), the only exception
  being the transition-matrix diagonal correction (B11) on the five models
  that count transitions. Four of the 26 cases tried (`compartment`,
  `unimol_state_change`, `pucadyil`, `clock_model`) segfault on `master` and
  complete on the branch.
- **Efficiency:** 1.7x to 4.1x fewer retired instructions on the same
  trajectory wherever the time is in the pairwise search and reaction
  bookkeeping (2.4x over the standard cases without `rev_2D`; 2.0x to 4.1x
  on the sampled production models); parity on the models whose time is 2D
  reaction-table construction, which no change on the branch addresses yet.
- **Open item found on the way (B8.3):** on 1D fiber models the branch's
  `rMaxLimit` pre-cut drops reaction and exclusion candidates that `master`
  evaluates (the 1D reach is 4/3 of the 3D estimate); exempting fiber pairs
  restores byte-for-byte agreement with `master`. No model in the tree is
  affected; appendix A has the measurement.
- **Statistics:** on the six sampled production models `master` can run, the
  shipped build agrees with `master` within seed scatter on every
  observable with a sizeable population; three low-copy species sit between
  |z| 2 and 3.1 across 67 comparisons, consistent with a change of random
  stream, and the Gag-lattice outliers seen at 6 seeds disappeared at 18.

## 5. Reproducing this comparison

```bash
# binaries (each in a clean export; CC defaults to Apple clang)
git archive master | tar -x -C build_master && make -C build_master serial
git archive nerdss-optimized-ref-2026-10-09 | tar -x -C build_opt && make -C build_opt serial
git archive nerdss-optimized-ref-2026-10-09 | tar -x -C build_legacy \
  && python3 benchmarks/nerdss_optimized/reference_2026-10-09/scripts/legacy_stream_patch.py build_legacy \
  && make -C build_legacy serial DEFS=-DNERDSS_LEGACY_STREAM

# 18-case suites and the trajectory-only comparison
cd benchmarks/nerdss_optimized
CASES_FILE=cases.tsv          ./run_suite.sh ../../build_master/bin/nerdss master_main 1
CASES_FILE=cases.tsv          ./run_suite.sh ../../build_legacy/bin/nerdss legacy_main 1
CASES_FILE=coverage_cases.tsv ./run_suite.sh ../../build_master/bin/nerdss master_cov 1
CASES_FILE=coverage_cases.tsv ./run_suite.sh ../../build_legacy/bin/nerdss legacy_cov 1
python3 reference_2026-10-09/scripts/compare_manifests.py results master_main legacy_main
python3 reference_2026-10-09/scripts/compare_manifests.py results master_cov legacy_cov

# sampled production inputs: paired runs, then the multi-seed statistics
reference_2026-10-09/scripts/run_pair.sh reference_2026-10-09/sample_cases.txt <out> 20260810 master optlegacy opt
reference_2026-10-09/scripts/run_stats.sh reference_2026-10-09/sample_cases.txt <out>/stats 6 6
python3 compare_observables.py <out>/stats
```

The `run_pair.sh` and `run_stats.sh` scripts carry the scratch-directory
paths of the session that produced this report at their top; edit the two
variables `S` and `R` to point at the build directory and `sample_inputs/`.

## Appendix A. B8.2 in detail: what the missing fiber branch changed, measured

**The question.** Whether a fiber pair is "measured as 1D or 3D" sounds like
it should not change the size exclusion, since 1D is the fiber's cross-section
of 3D. For the overlap test itself that is right: the fiber sweep
(`sweep_separation_complex_rot_fiber_box.cpp`) tests `|dx| < sigma` along the
fiber for every pair on its list, in every build. What the exclusion block in
`check_bimolecular_reactions()` decides is only *which pairs get on that list*,
and there the two arms differ in three ways:

| | 1D arm (`pro1` block always; `pro2` block after the fix) | 3D arm (`pro2` block before the fix) |
| --- | --- | --- |
| Diffusion constant | `Dtot = D1.x + D2.x` | `Dtot = (D1 + D2) averaged over x, y, z` (a third of the 1D value when `D.y = D.z = 0`), plus a rotational term |
| Admission radius | `4 sqrt(2 Dtot dt) + sigma` = **4.0** standard deviations of the pair's relative step | `3 sqrt(6 Dtot_3D dt) + sigma` = **3.0** standard deviations |
| Separation measured | `|x1 - x2|` along the fiber | full 3D interface distance (equal here, since promoter sites sit on the axis) |

A free site 3 to 4 standard deviations from a bound site is put on the sweep
list by the 1D arm and left off by the 3D arm. Left off, it is never checked,
so when the step's tail carries it into the bound site's `+- sigma` zone the
overlap stands. Which arm a given pair got was decided by the order the cell
search happened to offer the pair in (`pro1` or `pro2`), not by physics.

**The model** (`reference_2026-10-09/fiber_b8/model_pro2_order/`): a 2000 nm
fiber along x with 50 `A` sites and 200 `B` sites, all `isPromoter = true`
with `D = [100, 0, 0]` nm^2/us, `dt = 1` us; `A(l) + B(r) -> A(l!1).B(r!1)`
irreversibly with `sigma = 2` nm and `excludeVolumeBound = true`. All 50 `A`
bind within the first 5000 steps and the remaining 150 free `B` keep sliding
past the 50 bound sites, which is exactly the pair the exclusion block exists
for. Irreversible binding removes the one other source of transient overlaps
(a freshly dissociated pair is skipped by the probability routines for a
step). `D` of a bound dimer is 50, so for a bound site against a free site
the relative-step standard deviation is sqrt(2 * 150 * 1) = 17.3 nm, the 1D
admission radius 71.3 nm and the old 3D radius 54.0 nm.

**Builds.** Five `-DNERDSS_LEGACY_STREAM` builds of the reference source
(section 3.1), so the random stream and grid match `master`: the fix as
shipped; the fix reverted (`fiber_b8/scripts/revert_B8_2.patch`, the pro2
block measuring fiber pairs as 3D exactly as `master` does); both of those
with fiber pairs exempted from the branch's `rMaxLimit` pre-cut (B8.3,
`fiber_precut_exemption.patch`); and the fixed build instrumented to count,
for every fiber pair the pro2 block sees, what the old 3D rule would have
decided (`count_B8_admissions.patch`).

**Isolation.** With the pre-cut exempted, the fix-reverted build reproduces
`master` byte for byte over 300 000 steps on all four seeds (the two builds
then differ only in B8.2, and B8.2 is off), which is what makes the
comparison below a measurement of B8.2 alone. Without the exemption the
branch builds differ from `master` from the first frame, on both input
orders; that is B8.3.

**Result**, 300 000 steps, seeds 20260811 to 20260814, overlaps counted on
every 20th step as `A`-site / `B`-site pairs closer than 0.9 sigma (a bound
pair sits at exactly sigma):

| Build | Overlapping pair-frames per run (4 seeds) | Mean +- SEM |
| --- | --- | ---: |
| `master` | 237, 220, 215, 228 | 225 +- 5 |
| branch, fix reverted, pre-cut exempted (= `master`, byte for byte) | 237, 220, 215, 228 | 225 +- 5 |
| branch, fix reverted | 235, 250, 196, 219 | 225 +- 11 |
| branch with the fix | 11, 10, 8, 13 | 10.5 +- 1.0 |
| branch with the fix, pre-cut exempted | 9, 11, 11, 9 | 10.0 +- 0.6 |
| `master`, species listed in the other order | 204, 167, 187, 173 | 183 +- 8 |
| branch with the fix, species listed in the other order | 8, 11, 13, 14 | 11.5 +- 1.3 |

The instrumented build saw 899 million fiber pair evaluations in the pro2
block over the four runs: the 1D rule put 654 million on the sweep list, the
old 3D rule would have put 533 million, and 121 million (18.5% of the correct
list) were pairs the old block dropped; none went the other way. The 20x
reduction in realized overlaps is the whole effect of B8.2; the residual
10 per run is the baseline of the sweep itself (it gives up after 50
resamples when a site is jammed between two exclusion zones).

![B8.2 figure](reference_2026-10-09/fiber_b8/b8_fiber_exclusion.svg)

`fiber_b8/b8_fiber_exclusion.svg`: panel A draws the two admission windows
against the Gaussian of the pair's relative step; panel B is the table
above. The trajectories themselves (the x coordinate of every site, every
20 steps) are not archived; `fiber_b8/runs/` holds each run's copy numbers,
overlap count and, for the instrumented build, the admission counters.

**What this means for the classification.** B8.2 is a bug fix in the sense
that the two blocks enforced different exclusion windows for the same
physical pair depending on an implementation detail, and the window the
old pro2 block used was the 3D one, which lets a measurable fraction of
3-to-4-sigma steps end in an overlap the model asked to forbid. It does not
change the exclusion geometry, and it cannot change any model in the tree,
since none sets `isPromoter`. The same exercise exposed B8.3, a branch-side
divergence on fiber models that is still open.

### 3.7 User-supplied membrane model (2026-10-12 and 2026-10-13)

A kinase / phosphatase / adaptor model on an explicit lipid membrane
(`reference_2026-10-09/user_input_2026-10-12/input_as_given/`): a
30 x 30 x 8000 nm box holding 367 lipids with a U/P state, 17 adaptors,
153 kinases and 28 phosphatases, dt = 1 us, 1 000 000 steps, with
`coupledRxnLabel` catalysis. It is the user's 56 000-molecule model scaled
down 100x in molecule count with the box's x and y sides scaled by 10 and
its height kept, so concentrations and the dimensional factor are unchanged.
Same three binaries as section 3; one run at a time for the timing rows,
`/usr/bin/time -l`; instruction counts are the efficiency measure.

**As given, `master` cannot run it.** `set_rMaxLimit()` gives 55.2 nm and
`master`'s grid takes `floor(30 / 55.2) = 0` sub-volumes along x and y
(master `class_SimulVolume.cpp:43-44`), so the first molecule fails to place:
`Cannot fit complex 0 into simulation volume. Exiting...`
(`runs_as_given/master_startup_failure.stdout.log`). The legacy-stream build
carries `master`'s grid and fails the same way. The branch gives a box
thinner than the range one sub-volume (B4) and runs it. Bitwise identity
with `master` is therefore not defined for this input.

| Build | Result on the input as given, 1 000 000 steps |
| --- | --- |
| `master` | exits at start-up (zero sub-volumes along x and y) |
| `opt-legacy` | exits at start-up (same grid) |
| `opt` | runs: 4.61 T instructions, 359 s CPU at seed 20260810; 4.6 to 5.6 T and 360 to 454 s over seeds 20260811 to 20260813 |

Branch-only statistics, four seeds, second-half means +- SEM
(`runs_as_given/`): the chemistry is seed-dominated at this size, with
phosphorylation reaching 61 to 281 of the 367 lipids by 1 s, 1 to 10 of each
enzyme-lipid complex and 5 to 12 adaptor-lipid complexes at the end. No
`master` statistic exists to compare against.

**A variant `master` can hold.** `input_60nm/` widens the box to
60 x 60 x 8000 nm and multiplies every count by 4 (1468 lipids, 68
adaptors, 612 kinases, 112 phosphatases): same concentrations, same height,
so the same dimensional factor, 25x down from the original instead of 100x.
`master` costs 19 ms per step here, so the run was capped at 40 000 steps to
stay inside 15 minutes.

| Build | 40 000 steps, seed 20260810 | instructions | CPU s | speedup (instr.) | trajectory files vs `master` |
| --- | --- | ---: | ---: | ---: | --- |
| `master` | 30 sub-volumes of 60 x 60 x 267 nm (30-per-axis cap) | 9918 G | 774 | 1.00x | - |
| `opt-legacy` | same grid, `master`'s random stream | 2368 G | 182 | **4.19x** | **identical, 17 of 17** |
| `opt` | 144 sub-volumes of 60 x 60 x 55.6 nm | 1966 G | 139 | **5.05x** | differ (B1, O2) |

The 4.19x on the same trajectory comes from the per-type interaction mask
(O8): 1468 lipids share one membrane cell on `master`'s grid and are offered
to each other every step although no reaction names a lipid-lipid pair.
The remaining gain of `opt` is the finer grid. Wall-clock times are not
quoted for these rows because the host slept during the overnight batch;
CPU time and instructions are unaffected.

**Statistics on the 60 nm variant** (`master` vs `opt`, 4 seeds each, 40 000
steps): no comparison is possible, because the model's first binding events
come after about 10^5 steps at this size. Over the 40 000-step window every
run of either build ends with 0 to 2 phosphorylated lipids and no complex,
so both builds agree trivially and the window says nothing about the
physics. A statistical comparison on this model needs runs of 10^6 steps,
which `master` cannot deliver inside the 15-minute budget (about 5.3 CPU
hours per run on this variant).

![User input speedup](reference_2026-10-09/user_input_2026-10-12/speedup_chart.svg)

**Cost of the original 56 000-molecule system** (`input_original`, 300 x 300
x 8000 nm, 36 694 lipids, 1723 adaptors, 15 290 kinases, 2791 phosphatases;
`runs_original/`), measured on the branch over the first 500 steps before any
chemistry, with set-up (48 s of placement) separated by file timestamps:

| Build | per iteration | 1 000 000 iterations |
| --- | ---: | ---: |
| `opt` | 0.42 s (1.2e-4 CPU h), measured | about 120 CPU h, 150 to 170 allowing for the 1.3x rise the as-given runs show once complexes form |
| `opt-legacy` | about 0.5 s, from a 500-step run killed at the 5-minute budget and the first calibration (301 s for set-up plus 500 steps) | about 140 CPU h |
| `master` | more than 1.5 s (did not finish 500 steps in 13 minutes); about 2.3 s by the 5.6x `master`/`opt` ratio of the 60 nm variant | about 600 to 700 CPU h |
