# Merge ("condense") candidates on `nerdss-optimized-1.4`

Survey of `nerdss-optimized-1.4` at `b4d7438` (2026-11-04). Nothing in the source was
changed: this is a proposal. The scan scripts are in [`tools/`](tools/) and are
described under [Method](#method); run them from the repository root, e.g.
`python3 docs/merge_candidates/tools/extract.py && python3 docs/merge_candidates/tools/similar.py`
(intermediate files go to `docs/merge_candidates/tools/out/`, or `$DUPSCAN_OUT`).

Every file:line below is relative to the repository root and refers to `b4d7438`.

**Scope rule assumed throughout:** a merge is only proposed as "bitwise" when
every `DATA/` and `RESTARTS/` byte of `benchmarks/nerdss_optimized/cases.tsv` +
`coverage_cases.tsv` must stay identical for `bin/nerdss` *and*
`bin/nerdss_legacy`, i.e. same floating-point expressions (same statement
boundaries - see the FMA note in `docs/debranching_plan.md` "Axis 2"), same RNG
draw order and count, same tolerances, same container insertion order.

---

## 0. What is already consolidated (not re-proposed)

| Done | Where | Notes |
| --- | --- | --- |
| `Coord`+`Vector` -> `Vec3D` | `include/classes/class_Vec3D.hpp` | docs/nerdss_optimized.md |
| Reflector dispatch, box/sphere/compartment choice in one place | `src/boundary_conditions/reflect_dispatch.cpp` | |
| Sphere + compartment reflectors merged behind `RadialSide` | `reflect_traj_complex_radial.cpp`, `reflect_traj_tmp_crds_radial.cpp` | debranching step 3 |
| Shared complex-extent scan (`scan_axis_extent*`, `for_each_complex_point`, `ExtremePoint`) | `include/boundary_conditions/complex_extent.hpp` | box reflectors deliberately left as box-specific algorithms |
| Association dispatch | `src/reactions/association_dispatch.cpp` | |
| Sweep dispatch; `sweep_separation_complex_rot_box` + `_memtest_box` merged into `sweep_separation_box(..., ignoreZ)` | `sweep_separation_dispatch.cpp`, `sweep_separation_complex_rot_box.cpp:36` | |
| Cluster-sweep partner collection / resample loops | `src/trajectory_functions/cluster_sweep_common.cpp` | the rest of cluster box/sphere **explicitly declined** ("five flags") - see G below |
| Five association angles for state changes | `src/reactions/apply_state_change_rotations.cpp` | associate_* explicitly not folded - see D2 below, which shows the fold is a pure argument mapping for `associate_box` |
| `resample_complex_trajectory` (6-draw order in one place) | `src/trajectory_functions/resample_complex_trajectory.cpp` | the *surface* variant is still inlined 7x - see C |
| Shared bookkeeping (`record_crossing_pair`, `zero_partner_probvec`, `remove/add_to_monomer_list`, `erase_bond`, `update_state_change_observable`, ...) | `src/reactions/reaction_bookkeeping.cpp` | |
| `check_compartment_reaction` entering/exiting collapsed; `determine_compartment_probability` | `src/reactions/determine_compartment_probability.cpp` | |
| `finish_implicit_lipid_binding`, `association_rejection_checks`, `association_predicates` | `src/reactions/` | |
| `populate_reaction_lists_for_add` folded | `e6e3eca` | |
| `RS3D_of_molType` cached lookup | `include/classes/class_Membrane.hpp:208` | 2 inline copies remain in `EXEs/nerdss.cpp` - see F |
| 17 unreachable files + dead functions inside live files deleted | `e7be6f7`, `deda0b3`, `974a047` | serial side is now clean; MPI side is not - see A |

---

## 1. Ranked summary

Score = payoff (1-5, LOC + maintainability) x safety (1-5, 5 = bitwise by
construction and covered by the suite). "Dup LOC" = duplicated significant
lines (comments/blank/brace-only lines excluded); "saves" is the rough net.

| Rank | ID | Family | Dup LOC (saves) | Sites | Bitwise? | Payoff | Safety | Score |
| ---: | --- | --- | --- | ---: | --- | ---: | ---: | ---: |
| 1 | B | `parse_input` / `parse_input_for_add` | 177 common of 189/198 (~170) | 2 defs, 5 callers | yes (parse-time, no FP/RNG) | 4 | 5 | 20 |
| 2 | A | Dead MPI functions + unbuilt `mpi_proto/` | ~330 in src + 1345 in mpi_proto | 0 callers | n/a (deletion) | 3 | 5 | 15 |
| 3 | D | Association/state-change shared tails: `finalize_bond`, `resolve_state_change_rxn`, `commit_state_change`, `associate_box` rotations -> `apply_state_change_rotations`, `commit_tmp_association_coords` | ~85 + 4x30 + 2x45 + 50 + 14x4 (~300) | 2+4+4+1+14 | yes | 4 | 4-5 | 18 |
| 4 | C | Surface-aware resample helper (7 inline copies) + `resample_traj` k1/k2 twin blocks | ~75 (~60) | 7 | yes (same draw order) | 3 | 5 | 15 |
| 5 | H | Molecule/complex creation: `claim_molecule_slot`/`claim_complex_slot`, molecule "shell" + registration | ~90 (~70) | 3 mol + 2 cplx slots, 5 init fns | yes, if differences are parameters | 3 | 4 | 12 |
| 6 | K | `CreateDestructMol` == `TransmissionMol`; templated restart read/write + ctors | ~150 (~110) | 6 blocks | yes (RESTARTS compared) | 3 | 4 | 12 |
| 7 | I | `remove_from_subcell` (6x), zeroth-order `creation_volume` + Poisson draw (2x each), coupled-destruction block (2x) | ~110 (~80) | 6+2+2 | yes | 3 | 4 | 12 |
| 8 | F | RS3D-by-reaction lookup (5x) -> `Membrane::RS3D_of_rxn`; 2 inline `RS3D_of_molType` scans in serial main | ~45 + 20 | 7 | yes, keep the D-average as one statement; keep each site's guard | 2 | 4 | 8 (high maintainability: 300-char lines) |
| 9 | G | `sweep_separation_complex_rot_memtest_cluster_box` vs `_sphere` via compile-time geometry traits | 126 common of 148/158 (~120) | 2 | yes | 3 | 4 | 12 (reverses a documented decision) |
| 10 | M | Serial-vs-MPI setup drift: delete `initialize_paramters_for_implicitlipid_model`, make serial main call `parse_input_for_*` / `check_overlap` | ~150-300 | 2-4 | serial yes, MPI changes by design | 4 | 3 | 12 |
| 11 | J | Dissociation product-iface resolution (`break_interaction` / `_implicitlipid`) | 2x22 (~20) | 2 | yes | 1 | 5 | 5 |
| 12 | L | `print_complex_hist` serial vs `io_mpi` | 28 of 30/51 (~28) | 2 | serial yes | 1 | 5 | 5 |
| 13 | N | `Quat::from_axis_angle` (19 sites); `orient_crds_to_template` vs json `molecule_orientation` | ~40 + 43 | 19 + 2 | yes (parameterise the tmpComCoord oddity) | 2 | 4 | 8 |
| 14 | S | Small single-definition helpers: `dd_table_step`, row-1 table interpolation (`get_prevNorm`/`get_prevSurv`), `sim_time_seconds` (23x), box lipid-dip loop (5x), `write_traj`/`write_xyz` frame writer, `inner_coord_frames` | ~120 total | many | yes if each expression is kept verbatim as one statement | 2 | 4 | 8 |
| 15 | F2 | Bimolecular probability skeleton (1D/2D/3D, IL 2D/3D): canonical pair order, previous-reweight lookup, sep<0 clamp, currReweight emplace | ~25 x5 (~80) | 5 | yes for the shared parts; 1D uncovered | 2 | 3 | 6 |
| 16 | G2 | `sweep_separation_complex_rot_sphere` + `_memtest_sphere` (like the box merge) | 59 common of 121/150 | 2 | yes, but `_memtest_sphere` is uncovered | 3 | 2 | 6 |
| 17 | P | MPI left/right twins, `serialize`/`deserialize` mirrors | 48+31+... (~250) | many | MPI not deterministic - cannot be shown | 3 | 2 | 6 |

**Not recommended** (looked at, merge would change bits or merge nothing): see [section Q](#q-families-examined-and-not-recommended).

---

## 2. Dead code

Found by transitive reachability from `main` (both EXEs), every class method and
every `run_code_tests/*.cpp` function (script `reach.py`), then confirmed by
`grep -w` over `src/ EXEs/ include/ run_code_tests/ mpi_proto/` (only commented
call sites exist). **Serial build: no unreachable free functions remain** (the
`e7be6f7`/`deda0b3` cleanups did their job). Everything below is MPI-only
(`src/mpi`, `src/debug` are compiled only by `make mpi`).

| Function | Location | Lines | Evidence |
| --- | --- | ---: | --- |
| `check_molecule_coordinates` | `src/mpi/check_molecule_coordinates.cpp:40` | 44 | only commented calls `EXEs/nerdss_mpi.cpp:847`, `:1052` |
| `exit_ghosted_zone`, `enter_ghosted_zone` | `src/mpi/check_molecule_coordinates.cpp:96`, `:108` | 12 | no callers at all |
| `update_copyCounters_enter_ghosted_zone` / `_exit_` (98% identical to each other) | `src/mpi/check_molecule_coordinates.cpp:120`, `:147` | 40 | only callers are the two dead functions above; commented calls in `send_receive.cpp:373,378,514,519` |
| `move_complexes_based_on_propagation` | `src/mpi/check_molecule_coordinates.cpp:176` | empty body | no callers |
| `check_enter_ghosted_zone`, `check_exit_ghosted_zone` | `src/mpi/deserialize.cpp:850`, `:864` | 22 | only commented calls `deserialize.cpp:562,571,753,762` |
| `serialize_molecules_partial` | `src/mpi/serialize.cpp:95` | 30 | no callers |
| `delete_disappeared_complexes` | `src/mpi/delete_disappeared.cpp:98` | 24 | "which nothing calls" (`serialize.cpp:28` comment); `_partial` is the live one |
| `disconnect_molecule_partners` | `src/mpi/delete_disappeared.cpp:308` | 21 | no callers |
| `debug_bndpartner_interface` | `src/debug/debug.cpp:138` | 44 | only commented calls `remove_empty_slots.cpp:124`, `nerdss_mpi.cpp:1186` |
| `write_debug_information` | `src/debug/write_debug_information.cpp:10` | empty body | no callers (whole file) |
| whole directory `mpi_proto/` (`complex.cpp`, `full_par_comm.cpp`, own `main`) | `mpi_proto/` | 1345 | not in `Makefile` `INCLUDE_FOLDERS` nor `CMakeLists.txt` `GLOB src/*/*.cpp`; referenced only by `REPLACEMENT_REPORT.md` |

Plus the matching declarations in `include/mpi/mpi_function.hpp` (lines 29, 85,
91, 97, 103, 110, 133, 138, 149, 309, 355) and `include/debug/debug.hpp:34,71`.
`debug_print`/`debug_print_complex` look unreachable to a call graph but are
used through `include/macro.hpp:36-46` macros - **keep**. `LEGACY_read_restart`
/ `LEGACY_write_restart` are live (the `.dat` restart format).

Not dead, but ~2,450 lines of commented-out code (`// ...;`/`{`/`}` lines) remain,
concentrated in exactly the files the merges below touch:
`associate_box.cpp` (127), `EXEs/nerdss.cpp` (126),
`check_for_unimolecular_reactions_population.cpp` (95),
`perform_bimolecular_state_change_sphere.cpp` (90, including a 100-line block
at `:159-252`), `class_Molecule_Complex.cpp` (80), `break_interaction.cpp` (69),
`associate_ImplicitLipid_{sphere,box}.cpp` (68/65), `io_mpi/write_pdb.cpp`
(an 80-line commented copy of the function at `:9-88`). Removing them first
(stage 0 below) shrinks every later diff and is verifiable at object level.

---

## 3. Per-family detail

### B. `parse_input` vs `parse_input_for_add` (rank 1)

- `src/parser/parse_input.cpp:227` `parse_input` (225 body lines) and
  `src/parser/parse_input.cpp:519` `parse_input_for_add` (235).
  Token similarity 0.95 normalized / 0.94 raw / 0.93 line; 177 of 189/198
  significant lines identical (`diff -w`).
- Callers: `parse_input` - `EXEs/nerdss.cpp:159`, `src/parser/parse_input_for_a_new_simulation.cpp` (MPI),
  run_code_tests; `parse_input_for_add` - `EXEs/nerdss.cpp:333`,
  `src/parser/parse_input_for_add_file.cpp:88`, `run_code_tests/add_file_molecules_check.cpp:350`.
- What actually differs (complete list from the diff):
  1. offsets captured on entry (`addForwardRxnNum`, `addBackRxnNum`,
     `addCreateDestructRxnNum`, `addMolTemplateListNum`) and `molTypeIndex += addMolTemplateListNum`;
  2. error/log strings ("Add reaction file cannot be opened", "Add input file parsing complete");
  3. new-sim sets `hasParsedMol` and exits if reactions precede molecules; add has a TODO instead;
  4. new-sim does `MolTemplate::numEachMolType = std::vector<int>(numMolTypes)`; add does
     `numEachMolType.emplace_back(0)` per template and sizes `transitionMatrix`/`lifeTime`;
  5. `RxnBase::totRxnSpecies`: new-sim `= totSpecies + 1` inside the reactions block; add
     recomputes from `Interface::State::totalNumOfStates` + bimolecular count after the file;
  6. `populate_reaction_lists(..., addForwardRxnNum)` start index (0 for new).
- The ~50-line "TEMPORARY COUPLED RXN" resolution block (`parse_input.cpp:430-480` and its twin)
  is identical in both and is itself two copies (match by `absRxnIndex` vs by `label`) - fold to
  `resolve_coupled_rxn(fwd, back, cd, matchByIndex)`.
- **Bitwise:** yes. Parse-time only, no FP beyond what is copied, no RNG. Only stdout strings
  differ, and they can be kept exactly with a `const char*` from the mode struct.
- **Strategy:** one `static void parse_input_impl(..., const AddContext* add)` (`add == nullptr`
  = new simulation) with the six differences as `if (add)` branches in place; the two public
  functions become 3-line wrappers. Keep difference 3 and 5 gated exactly (do not "fix" the add
  path's missing order check in the same commit).
- **Validation:** suite + `make checks` (add_file_molecules_check) + any restart-with-add case
  in `benchmarks/nerdss_optimized/regression/`.

### A. Dead MPI code + `mpi_proto/` (rank 2)

See [section 2](#2-dead-code). Strategy: delete functions + declarations in one
commit, `mpi_proto/` in a second. **Verification:** `make serial` and
`make serial legacy` cannot change (the files are not in the serial object list;
binary hashes should be identical - `run_suite.sh` prints `binary_sha256`);
`make mpi` must link; `run_mpi_suite.sh` smoke at np=1,2.

### D. Association / state-change shared blocks (rank 3)

Big bodies (`associate_box` 566 lines vs `associate_sphere` 359; normalized
similarity 0.72 but only 170 common lines of 499/312) are *not* mergeable as a
whole (debranching plan step 6, "least mergeable pair"). Five sub-blocks are
verbatim and geometry-free:

| Sub-ID | Block | Copies | Lines | Helper |
| --- | --- | --- | ---: | --- |
| D1 | bond bookkeeping tail: `remove_from_monomer_list` x2, partner/iface fields, `assocDissocFile` BOND line, `isBound`/`index`, `bndlist`/`bndpartner`, freelist swap-pop x2, `zero_partner_probvec` x2, `ncross=-1`, crossings clear, `update_Nboundpairs`, `copyNumSpecies` +-, `bindPairList`, observable | `src/reactions/associate_box.cpp:915-1015`, `src/reactions/associate_sphere.cpp:464-561` | ~85 each | `finalize_bond(iter, ifaceIndex1, ifaceIndex2, reactMol1, reactMol2, reactCom1, currRxn, ...)` in `reaction_bookkeeping.cpp` (diff -w: only comments/formatting differ) |
| D2 | five angle rotations | `src/reactions/associate_box.cpp:310-372` | ~50 | call existing `apply_state_change_rotations(reactIface1, reactIface2, /*stateChangeIface=*/ifaceIndex2, /*facilitatorIface=*/ifaceIndex1, /*stateChangeMol=*/reactMol2, /*facilitatorMol=*/reactMol1, reactCom2, reactCom1, currRxn.assocAngles, currRxn, ...)`. Checked argument by argument (theta1/theta2/omega/phi1/phi2 and their isnan guards) - identical calls. **`associate_sphere.cpp:147-197` is NOT identical: it has no `isnan` guard on theta1/theta2** (box has); leave sphere alone or fix in a separate, declared behaviour-change commit. |
| D3 | state-change prologue: resolve forward/back index, `bindRadius`, `newState`, `assocAngles`, `relStateIndex` (+ the two `set_tmp_association_coords` loops) | `perform_bimolecular_state_change_box.cpp:15-53`, `_sphere.cpp:16-54`, `perform_implicitlipid_state_change_box.cpp:13-42` (+`:87-91`), `_sphere.cpp:14-43` (+`:87-91`) | ~30 x4 | `StateChangeRxn resolve_state_change_rxn(rxnItr, stateChangeMol, stateChangeIface, forwardRxns, backRxns, molTemplateList)` returning a struct; note bimol_sphere binds `ForwardRxn&` (non-const, it writes `bindRadius2D`), the others `const ForwardRxn&` - return indices, not references |
| D4 | state-change epilogue ("change the state" .. `crossings.clear()`) | bimol: `perform_bimolecular_state_change_box.cpp:301-351` == `_sphere.cpp:387-437`; IL: `perform_implicitlipid_state_change_box.cpp:299-387` == `_sphere.cpp:208-296` | 45 + 45 (two pairs) | `commit_bimolecular_state_change(...)`, `commit_implicitlipid_state_change(...)`; box/sphere stays a call-site difference |
| D5 | "write tmp coords to real coords" loop (`comCoord = tmpComCoord; iface.coord = tmpICoords[i]; clear_tmp_association_coords()`) | 14 loops in 9 files (`associate_box.cpp:808,817`, `associate_sphere.cpp`, `finish_implicit_lipid_binding.cpp:26`, `check_dissociation_implicitlipid.cpp:259`, all 4 state-change files x2) | 4 x14 | `Molecule::commit_tmp_association_coords()`; callers keep their own `trajStatus`/`need_to_send`/memberList extras (these differ: state-change loops do **not** set `need_to_send`, which matters only for MPI) |

Also in the box files: the "lipid below the box bottom" `zchg` loop is written
5 times (`perform_bimolecular_state_change_box.cpp:206-241`,
`perform_implicitlipid_state_change_box.cpp:204-239`, `associate_ImplicitLipid_box.cpp:188-210`),
each with a dead "lipid ABOVE" branch that computes and discards `ztmp`.
Helper `max_lipid_dip_below_box(com, moleculeList, membrane, zchg)` keeping the
`if (ztmp > zchg)` comparison and the facilitator-then-stateChange order.

- **Bitwise:** yes for D1, D3, D4, D5 (integer bookkeeping, copies, same order).
  D2 is a call with the same arguments. Coverage: association in box and sphere
  cases, state changes in `regression/bimol_state_change`, `unimol_state_change`;
  **confirm IL state change coverage** with the llvm-cov recipe in
  `known_uncovered.tsv` before D3/D4 on the IL files.
- Call sites: each block 2-4 definitions, each called from 1 dispatcher.

### C. Surface-aware trajectory resample (rank 4)

The 7-line block "if on sphere surface: `create_complex_propagation_vectors_on_sphere` -> trajTrans,
then 3 `GaussV()` for trajRot" is inlined 7 times:
`src/trajectory_functions/sweep_separation_complex_rot_sphere.cpp:185`,
`sweep_separation_complex_rot_memtest_sphere.cpp:184` and `:219`,
`cluster_sweep_common.cpp:131`, `src/classes/class_Cluster.cpp:135` and `:167`,
`create_complex_propagation_vectors.cpp:14`.
`class_Cluster.cpp:108` `resample_traj` (69 lines) is two identical halves for
`k1`/`k2` (and its else-branches are an inline `resample_complex_trajectory`);
it is 0.86 containment-similar to `resample_partner_trajectories`
(`cluster_sweep_common.cpp:105`).

- **Strategy:** `void resample_complex_trajectory_on_sphere(Complex&, const Parameters&)`
  next to `resample_complex_trajectory`, and
  `resample_complex_trajectory_for(Complex&, params, bool onSphereSurface)`; `resample_traj`
  becomes a loop over `{(p1,k1),(p2,k2)}` calling a shared `resample_once(k, didMove, ...)`
  also used by `resample_partner_trajectories`.
- **Bitwise:** yes - same calls, same draw order (on_sphere draws 2, then 3 for rotation).
  **Do not** route the *box* branch of `create_complex_propagation_vectors`
  (`:33-40`) through `resample_complex_trajectory`: it skips the `GaussV()` draw when
  `D == 0` (a stream-relevant difference), the others do not.
- Coverage: `sphere`, `gagsphere`, `cluster_gagsphere`.

### H. Molecule / complex creation (rank 5)

- Slot claiming, "reuse an empty slot or append": molecules at
  `src/reactions/create_molecule_and_complex_from_rxn.cpp:199-217`,
  `create_molecule_and_complex_from_transmission_rxn.cpp:12-28`,
  `src/system_setup/generate_coordinates_for_restart.cpp:31-47`; complexes at
  `create_molecule_and_complex_from_transmission_rxn.cpp:30-47`,
  `generate_coordinates_for_restart.cpp:50-66` (pair similarity 0.78-0.86).
  `create_molecule_and_complex_from_rxn` has the complex reuse **commented out**
  and always appends - that must stay a parameter, since complex indices reach
  the outputs.
- Molecule initialisers (0.83-1.00 pairwise; `initialize_molecule_after_zeroth_reaction`
  vs `_uni_reaction` normalized 1.00): `src/system_setup/initialize_molecule.cpp:5`,
  `src/reactions/initialize_molecule_after_reaction.cpp:72`, `:111`,
  `initialize_molecule_after_transmission_reaction.cpp:63`,
  `generate_coordinates_for_restart.cpp:86`. Common: type/mass/isLipid/isPromoter,
  `freelist` iota, `interfaceList` sizing, `isEmpty=false`, `index`, counters.
  Differences to keep: coordinate drawer; `numTotalUnits` touched (zeroth/uni/transmission)
  or not (restart, initialize_molecule); `index` = arg vs `numberOfMolecules`;
  **transmission sets neither `isPromoter` nor `id`** (see defects).
- **Strategy:** `int claim_molecule_slot(moleculeList)`, `int claim_complex_slot(complexList)`,
  `Molecule molecule_shell(const MolTemplate&)` + `void register_molecule(Molecule&, index,
  params*, bool assignId)`; each initialiser keeps its own drawer call between the two.
- **Bitwise:** yes (no FP; `maxID`/counters are not on the RNG stream). Note the
  `try/catch (std::out_of_range)` around `emptyMolList.back()` can never fire; keep the
  semantics, do not "fix" it in the merge commit.

### K. `CreateDestructRxn` vs `TransmissionRxn` (rank 6)

- `CreateDestructRxn::CreateDestructMol` (`include/classes/class_Rxns.hpp:629`) and
  `TransmissionRxn::TransmissionMol` (`:701`) are the same struct (molTypeIndex, molName,
  interfaceList) minus `serialize`.
- Duplicated readers/writers: `src/io/read_restart.cpp:318` vs `:342` (json, 0.89 sim,
  22/24 lines); legacy reader `read_restart.cpp:1532-1654` vs `:1660-1815` (runs of 38, 25,
  18, 18 identical lines); legacy writer `src/io/write_restart.cpp:961-1023` vs `:1030-1091`
  (36 identical lines); ctors `src/classes/class_Rxns.cpp:304` vs `:402` (cont 0.94).
- **Strategy:** one `struct RxnMolSpec` with both nested names as aliases
  (`using CreateDestructMol = RxnMolSpec;`), then `template <class Rxn>` read/write helpers
  with a traits accessor for the one field that differs (`creationRadius` vs `bindRadius`).
- **Bitwise:** yes - the suite compares `RESTARTS/`. Legacy `.dat` coverage needs
  a `.dat` restart case (e.g. `sample_inputs/gagLatticeRemodeling/.../restart1.dat`)
  run before/after; the json path needs a case with a transmission reaction
  (`compartment`).

### I. Destruction / creation helpers (rank 7)

- "erase molecule from its sub-cell member list" loop, 6 copies:
  `check_for_unimolecular_reactions_population.cpp:101`, `:508`,
  `check_dissociation_implicitlipid.cpp:318`, `check_for_unimolecular_reactions.cpp:217`,
  `perform_transmission_reaction.cpp:68`, `check_dissociation.cpp:166` ->
  `SimulVolume::remove_member(int subCell, int molIndex)`.
- Coupled destruction after dissociation: `check_dissociation.cpp:~140-175` vs
  `check_dissociation_implicitlipid.cpp:~217-262` (IL adds IL-state counting and a
  monomerList erase the non-IL copy lacks - see defects).
- `check_for_zeroth_order_creation.cpp:25-69` vs `:102-148`: creation volume
  (sphere / inside / outside compartment / box) and the Poisson inversion loop are
  each written twice in one function (`numEvents` is `unsigned` in one copy and `int`
  in the other; `lambda`/`prob` are `long double`). Helpers
  `double creation_volume(membrane, tmpl)` and
  `unsigned draw_poisson_events(long double lambda)` - keep `long double`, keep the
  `factor > bigNum -> gsl_ran_poisson` fallback (it re-draws; RNG count must match).
- **Bitwise:** yes if types are preserved.

### F. RS3D lookups (rank 8)

- "RS3D for this reaction" scan - match `RS3Dvect[i]`=bindRadius, `[i+100]`=rate,
  `[i+200]`= the 3-term D average, within 1E-15 - five copies:
  `src/reactions/associate_ImplicitLipid_box.cpp:17-23`,
  `associate_ImplicitLipid_sphere.cpp:23-29`,
  `perform_implicitlipid_state_change_box.cpp:43-69`,
  `perform_implicitlipid_state_change_sphere.cpp:44-70`,
  `check_dissociation_implicitlipid.cpp:181-187`.
  -> `double Membrane::RS3D_of_rxn(const ForwardRxn&, const std::vector<MolTemplate>&, double fallback = -1.0) const`.
- Inline `RS3D_of_molType` scans in the serial main: `EXEs/nerdss.cpp:1093-1102` and
  `:1157-1165` (MPI's `check_overlap.cpp` already uses the method).
- **Bitwise:** yes, with two cautions. (1) The D-average must be computed as **one
  statement** textually identical to `initialize_parameters_for_implicitlipid_and_compartment_model.cpp:136-138`
  (the stored value it is compared to with 1E-15) - FMA contraction on arm64/clang is decided
  per statement. (2) The two serial-main sites use **different guards**:
  `:1094` tests `membraneObject.implicitLipid`, `:1159` tests `hasRS3Dtable()`. For a
  compartment-only model these differ (table exists, implicitLipid false), so each site must
  keep its own guard around the call.

### G. Cluster sweep box vs sphere (rank 9)

- `src/trajectory_functions/sweep_separation_complex_rot_memtest_cluster_box.cpp:7` (172 lines)
  vs `_sphere.cpp:7` (182): normalized 0.87, raw 0.92, 126 common of 148/158 lines.
- Exactly five live differences (diff): `nStep` 10 vs 5 (`:70`); interface placement on the
  surface via `calculate_update_position_interface` (sphere `:142-174`); the
  "can't solve overlap" warning printed only by sphere (`:201-203`); fallback
  `sweep_separation_complex_rot_memtest_box/_sphere` (`:185,188` / `:208,211`); box zeroes
  `trajTrans`/`trajRot` at the end, sphere does not (`:258-274` vs `:281-286`).
- `cluster_sweep_common.cpp` header records the decision *not* to merge ("a body carrying five
  flags would be harder to read"). A compile-time traits struct avoids runtime flags:
  `struct BoxCluster { static constexpr int nStep = 10; static constexpr bool surfacePlacement = false; ... static void fallback(...); };`
  and `template <class G> void sweep_cluster(...)`. Bitwise by construction; covered by
  `cluster_mem_loc` and `cluster_gagsphere`. Worth it only if the owners accept reversing
  the recorded decision.

### M. Serial vs MPI setup drift (rank 10)

- `initialize_paramters_for_implicitlipid_and_compartment_model`
  (`src/system_setup/initialize_parameters_for_implicitlipid_and_compartment_model.cpp:7`,
  called by `EXEs/nerdss.cpp:229`, `:453`) vs `initialize_paramters_for_implicitlipid_model`
  (same file `:169`, called only by `parse_input_for_a_new_simulation.cpp:146` and
  `parse_input_for_a_restart_simulation.cpp:137`, i.e. MPI). Normalized 0.89, containment 0.98:
  the MPI one is the serial one minus transmission (compartment) RS3D rows, plus it fills
  `RS3Dvect` with 500 x -1 even without IL/compartment, and **never calls
  `forget_RS3D_matches()`**. Strategy: delete `_model`, call the compartment version from
  MPI with its `transmissionRxns`. Serial unchanged (bitwise); MPI behaviour changes
  (compartment RS3D rows appear; table empty for plain models, which `hasRS3Dtable()` readers
  already handle) - document as an MPI fix.
- `parse_input_for_a_new_simulation` (MPI) vs the inline block `EXEs/nerdss.cpp:155-283`:
  line ratio 0.65 and they have drifted (serial: `require_compartment_fits`,
  `_and_compartment_model`; MPI: `check_transition_matrix_size`, `_model`). Same pattern for
  restart (`parse_input_for_a_restart_simulation.cpp:152-169` == `nerdss.cpp:466-483`) and add
  (`parse_input_for_add_file.cpp:74-86` == `nerdss.cpp:319-331`).
- `check_overlap` (`src/reactions/check_overlap.cpp:37`, MPI only) vs the serial inline loop
  `EXEs/nerdss.cpp:1090-1180`: serial tests `OnSurface` and has a fiber branch; MPI tests
  `std::abs(D.z) < 1E-10` and has none.
- Strategy: make the serial code the reference, move it verbatim into the shared functions,
  call them from both mains. Serial bitwise by construction; MPI changes and is validated
  statistically (`statistical_check.sh`, `run_mpi_suite.sh`) - two commits each.

### J. Dissociation product interfaces (rank 11)

`src/reactions/break_interaction.cpp:221-260` and
`break_interaction_implicitlipid.cpp:10-32`: identical symmetric/state-matched
product resolution (`absIface1/2`), differing only in how much is printed before
`exit(1)`. -> `std::pair<int,int> resolve_dissociation_product_ifaces(...)`. Bitwise.

### L. `print_complex_hist` serial vs MPI (rank 12)

`src/io/print_complex_hist.cpp:78` vs `src/io_mpi/print_complex_hist.cpp:11`:
28 of 30 serial lines shared; MPI adds "skip complexes with no non-IL member" and
"count only on `ownerRank`". One body with a `const std::function<bool(const Complex&)>&
include` (serial passes nullptr). `unordered_map` iteration order is unchanged because
the insertion sequence is unchanged. Bitwise for serial.

### N. Rotations (rank 13)

- Axis-angle quaternion `Quat(cos(a/2), sin(a/2)*ax.x, ...)` + `normalize()` written at 19
  sites (`omega_rotation.cpp` 4, `phi_rotation.cpp` 4, `orient_crds_to_template.cpp` 3,
  `write_bonded_complex_json.cpp` 3, `theta_rotation.cpp` 2, `association_predicates.cpp` 2,
  `associate_box.cpp` 1) -> `Quat::from_axis_angle(const Vec3D& unitAxis, double angle)`.
  Check per site whether `normalize()` follows (do not add one where it did not).
- `src/reactions/orient_crds_to_template.cpp:5` vs `src/io/write_bonded_complex_json.cpp:226`
  `molecule_orientation` (normalized 0.90, 43 of 59/56 lines; the JSON copy is "slightly modified
  ... does not make weird state modifications"). Make the core a pure function on
  `(com, std::vector<Vec3D>& ifaces, comForSecondStep)` and call it on `tmpICoords` or on a copy.
  Both copies subtract `tmpComCoord` from already-relative vectors in the second step
  (`orient_crds_to_template.cpp:96-97`, json `:85-86`), so `comForSecondStep` must be a
  parameter to stay bitwise.
- `omega_rotation` vs `phi_rotation` (0.87 normalized, 0.22 line): same rotate / check /
  reverse / retry skeleton but different angle function, axis, rotation pivots and exit
  conditions - only the quaternion helper is worth sharing.

### S. Small single-definition helpers (rank 14)

| Helper | Copies | Note |
| --- | --- | --- |
| `dd_table_step(Dtot, dt) = sqrt(Dtot*dt)/50` | `DDpirr_pfree_ratio_ps.cpp:15`, `create_DDMatrices.cpp:7`, `get_prevSurv.cpp:8` | comment in two of them: "must be the SAME definition" |
| row-1 linear interpolation of a 2xN gsl table | `get_prevNorm.cpp:3`, `get_prevSurv.cpp:3` (0.89 normalized) | `get_prevSurv` = step + `get_prevNorm`; keep the return expression as one statement |
| `sim_time_seconds(itr, params)` | 23 occurrences of `(simItr - params.itrRestartFrom) * params.timeStep * {1E-6|usToSeconds} + params.timeRestartFrom` | `usToSeconds` is `1E-6`, identical bits; one statement |
| periodic-output block | `EXEs/nerdss.cpp:1489-1509` == `:1590-1610` | note drift: `"system_information.dat"` vs `"DATA/system_information.dat"` |
| xyz frame writer | `src/io/write_traj.cpp:6` vs `write_xyz.cpp:6` (0.78) | differ in name truncation, `isEmpty` skip, EMTY padding, header line - 4 parameters for ~15 lines; marginal |
| `inner_coord_frames(com, comNew)` | the pair `inner_coord_set` + `inner_coord_set_new` called together at `class_Molecule_Complex.cpp:906`, `:1009`, `calculate_update_position_interface.cpp:26` | the two functions themselves differ (not mergeable) |
| `check_for_structure_overlap_system` reactCom1/reactCom2 loops | `src/reactions/check_for_structure_overlap_system.cpp:63-87` vs `:112-137` | asymmetric (see defects); loop over `{com, checker}` pairs |
| `rename io_mpi write_pdb` | `src/io_mpi/write_pdb.cpp:90` writes `.cif`, not PDB | rename + delete commented copy `:9-88`; not mergeable with `src/io/write_pdb.cpp:6` (different format) |

### F2. Bimolecular probability skeleton (rank 15)

`determine_{1D,2D,3D}_bimolecular_reaction_probability.cpp` (97/99/71 lines) and
`determine_{2D,3D}_implicitlipid_reaction_probability.cpp` (66/47): pairwise
normalized similarity 0.50-0.60, only 22-29 common lines per pair. Shared:
`isDissociated` guard + `selected_rate`, `sep < 0` clamp, the proA/proB canonical
ordering (`3D:93-111` == `2D:115-131` == `1D:128-137`), the previous-reweight lookup
loop, `crossings.back().prob`, `currReweight.emplace_back`. Kernels differ by design
(`passocF`/`pirr_pfree_ratio_psF`, 2D tables + `DDpirr_pfree_ratio_ps`, `passocF_1D`, IL tables).
Recommend small helpers (`canonical_pair()`, `find_prev_reweight()`, `store_reweight()`), not a
template method. Bitwise for the shared parts (copies/compares only). The 1D file and
`passocF_1D`/`pirr_pfree_ratio_psF_1D` are in `known_uncovered.tsv` - test by a
differential harness, not the suite.

### G2. Sphere pair sweeps (rank 16)

`sweep_separation_complex_rot_sphere.cpp:6` vs `_memtest_sphere.cpp:18`: 0.79
normalized but only 59 common lines. Differences: com1 placed by
`calculate_update_position_interface` when on surface; com1 surface resample; `p2 >
startProIndex` guard (memtest only); live warnings (memtest) vs commented
(sphere). Also: memtest_sphere fills `memCheckList` (`:51-72`) but never reads it;
sphere has an empty `if (sqrt(...) > sphereR) {}` (`:133-137`) and a
diagnostic loop that writes nothing (`:208-240`). Merge would mirror the box merge
(`sweep_separation_box(..., flags)`), but `_memtest_sphere` is in
`known_uncovered.tsv` (only gagsphere at nItr 20000 reaches it, 11 min), so the suite
cannot vouch for it. Do it after C, with the slow case run once.

### P. MPI twins (rank 17)

`deserialize_complexes` vs `_right` (`src/mpi/deserialize.cpp:132` / `:295`, 0.96 normalized,
48 of 54/61 lines), `deserialize_molecules` vs `_right` (`:430` / `:623`, 0.79),
`receive_right/left_neighborhood_zones` (`src/mpi/send_receive.cpp:280` / `:420`, 0.96
normalized but 31 common lines), `send_data_to_left/right_neighboring_ranks` (`:77` / `:173`, 0.91),
`is_ghosted` vs `is_owned_by_processor` (`src/mpi/check_molecule_coordinates.cpp:6`/`:23`, 0.99),
`find_complex`/`find_molecule` (`src/mpi/id_index.cpp:131`/`:144`, 1.00 normalized),
`increase_recv/send_buffers` (`include/mpi_context.hpp:158`/`:168`),
`prepare_data_structures_for_parallel_execution`/`_restart` (`src/mpi/prepare.cpp:282`/`:454`, 0.57),
and the `serialize`/`deserialize` mirror pairs in
`include/classes/class_{Parameters,MolTemplate,Molecule_Complex,Membrane,Rxns,copyCounters}.hpp`
(0.94-1.00 normalized; ~250 lines that a single `template <class Ar> void transfer(Ar&)` visitor
would halve). The left/right twins have **drifted** (e.g. the left copy re-parents
`justBoundThisStep` members, `deserialize.cpp:191`, the right copy does not; snapshot types
differ: `vector<MolSnapshot>` vs `unordered_map<int, MolSnapshot>`), and
`docs/debranching_plan.md` records that `nerdss_mpi` is not reproducible run to run, so no
bitwise check exists. Lowest safety; do only with the MPI owners and the invariant checkers
(`check_structure_invariant`, `check_ownership_invariant`).

### Q. Families examined and not recommended

| Family | Why not |
| --- | --- |
| `reflect_traj_*_box` vs `_radial`, `reflect_complex_rad_rot_box/_sphere`, `reflect_traj_check_span_box/_sphere`, `_nocheck_box/_sphere` | Already reduced to the shared scan; box = per-axis wall algorithm, sphere = farthest radial point (debranching plan "What survived contact"). `_nocheck_*` and `check_span_*` are uncovered. |
| `check_if_spans_box` vs `_sphere` | Different algorithms (per-axis translation vs radius+radial shift); already refactored. |
| `create_complex_propagation_vectors` vs `_on_sphere` | Different levels: the latter is the tangent-plane step primitive called by the former; only the inline caller pattern (C) is shared. |
| `associate_box` vs `associate_sphere` bodies, `associate_implicitlipid_box/_sphere`, `perform_*_state_change_box/_sphere` geometry parts | Geometry genuinely differs (plane/z vs geodesic arcs, `find_position_after_association`, `update_association_coords_sphere`). Take the shared blocks (D) instead. |
| `passocF` vs `passocF_1D`, `pirr_pfree_ratio_psF` vs `_1D` vs `DDpirr_pfree_ratio_ps` | Different Green's functions. Only shared fragment is the `ep1*erfc` / `exp(-a^2)*erfcx` fallback (2 copies, `passocF.cpp:20-31`, `pirr_pfree_ratio_psF.cpp:28-38`), which could become one inline helper (bitwise). |
| `prob_entering_compartment` vs `prob_exiting_compartment` (`functions_implicitlipid.cpp:312`/`:349`, 0.99 normalized) | Besides `r = +-dr + R`, the operand order differs: `(r+R-sigma)` vs `(R-sigma+r)`, `(r-R-sigma)` vs `(R-sigma-r)`. A sign-parameterised merge changes bits on the covered (exiting) path. |
| `get_geodesic_distance` vs the sphere arm of `calc_interface_distance` | `acos(dot/(r1*r2))` vs `acos(dot/r1/r2)` - not bit-identical. |
| `create_survMatrix` / `create_normMatrix` / `create_pirMatrix` / `create_DDMatrices` | Different integrands and GSL routines; only file co-location would help. (`create_normMatrix` allocates a 10M-entry workspace per R step - a perf note, not a merge.) |
| `inner_coord_set` vs `inner_coord_set_new` | Different frame construction. |
| `theta_rotation` / `phi_rotation` / `omega_rotation` / `reverse_rotation` / `rotate` / `transform` / `calc_angular_displacement` vs `calc_one_angular_displacement` | Distinct geometry; only N (quaternion helper) applies. |
| `check_dissociation` vs `check_dissociation_implicitlipid` | 63 common of 103/222 lines; IL path has `#ifdef mpi_` free-protein counting (`:80`) and IL state accounting. Extract only the coupled-destruction block (I). |
| `check_for_unimolecular_reactions` vs `_population` | 0.32 normalized; different event models. |
| `SimulVolume::update_memberMolLists` serial vs MPI overload (`class_SimulVolume.cpp:385`/`:570`) | 67 common of 131/134 lines, but the serial one re-bins with a restart loop and the MPI one has rank-edge rules; folding needs `xOffset=0, rank=0, nprocs=1` to be proven equal - possible later, low payoff. |
| `check_bimolecular_reactions` excludeVolumeBound half | Already noted in `known_uncovered.tsv` as the largest remaining duplication and the least safe one (0 hits in every case). |

---

## 4. Defects / drift found while comparing (do NOT fix inside a merge commit)

1. `initialize_molecule_after_transmission_reaction` (`src/reactions/initialize_molecule_after_transmission_reaction.cpp:63`)
   assigns neither `tmp.id = Molecule::maxID++` nor `isPromoter`; `Molecule::id` has no initializer
   (`include/classes/class_Molecule_Complex.hpp:379`), so transmission-created molecules carry an
   indeterminate id. Every other initialiser assigns one.
2. `associate_sphere` theta rotations lack the `isnan` guard `associate_box` has (`associate_sphere.cpp:147,155`).
3. `initialize_paramters_for_implicitlipid_model` (MPI) never calls `forget_RS3D_matches()` and omits transmission RS3D rows.
4. `check_for_structure_overlap_system.cpp`: reactCom1 members are tested with
   `measure_overlap_protein_interfaces` (`:91`), reactCom2 members with
   `measure_overlap_free_protein_interfaces` (`:140`).
5. Coupled destruction after dissociation: the IL copy removes the destroyed molecule from
   `monomerList` (`check_dissociation_implicitlipid.cpp:325`); the plain copy (`check_dissociation.cpp:~150-170`) does not, although `break_interaction` may just have added it. Verify.
6. Serial main: RS3D guard `implicitLipid` at `EXEs/nerdss.cpp:1094` vs `hasRS3Dtable()` at `:1159`
   (differs for compartment-only models); system-info file `system_information.dat` (`:1489`) vs
   `DATA/system_information.dat` (`:1590`).
7. Serial overlap loop tests `OnSurface` and handles fibers; MPI `check_overlap` tests `|D.z| < 1E-10` and has no fiber branch.
8. `sweep_separation_complex_rot_memtest_sphere` builds `memCheckList` but measures full 3D
   separations (box memtest drops z for surface partners).
9. MPI left/right deserializers have drifted (`justBoundThisStep` handling on one side only).
10. `try { ... emptyMolList.back() ... } catch (std::out_of_range&)` (3 copies) cannot catch; an
    all-stale `emptyMolList` would call `back()` on an empty vector.

---

## 5. Safe staged procedure

**Where:** a separate worktree, so the baseline build stays untouched, e.g.
`git worktree add ../NERDSS-merge -b claude/condense-1 nerdss-optimized-1.4`.

**Baseline once (at the branch point):**
```bash
make serial && make serial legacy && make mpi
./benchmarks/nerdss_optimized/run_suite.sh bin/nerdss        base   1   # cases.tsv
./benchmarks/nerdss_optimized/run_suite.sh bin/nerdss_legacy base-L 1
# + the same for coverage_cases.tsv (as run_suite.sh/compare_suites.sh are driven there)
make checks                                                     # run_code_tests/*_check.cpp
```

**Stage 0 - comment-only cleanup** of the files a family will touch. Verifiable without running:
compile each touched TU before/after with the same flags and compare `objdump -d
--no-show-raw-insn` of `.text` (raw `.o` can differ only where `__LINE__`/`assert`
strings move). Identical disassembly = done.

**Then one family (or sub-family) per commit**, in this order:
A (dead, 2 commits) -> B -> D1, D2, D3, D4, D5 (separate commits) -> C -> J -> H -> I -> K
-> F -> S (one helper per commit) -> L -> G (if accepted) -> N -> F2 -> M (serial-preserving
commit, then the MPI behaviour commit) -> G2 -> P.

Per commit:
1. `make serial`, `make serial legacy`, `make mpi` (all three must build; `make mpi` exercises the
   `#ifdef mpi_` blocks inside shared files such as `check_dissociation_implicitlipid.cpp:80`).
2. `make checks`.
3. `run_suite.sh` for `cases.tsv` and `coverage_cases.tsv` with `bin/nerdss` and `bin/nerdss_legacy`;
   `compare_suites.sh base <label>` must report every case identical (DATA/ and RESTARTS/).
4. Coverage gate: if any function the commit changes is in `known_uncovered.tsv` (1D probability,
   fiber sweep, check_span, `_nocheck_*`, memtest_sphere, compartment-entering) or is not reached
   by a case (re-run the llvm-cov recipe in that file's header when unsure), the suite proves nothing:
   add a randomized differential harness that calls the old (copied verbatim into the harness) and
   new function on >=10^6 random states and compares bit for bit - the method the reflector merge used.
5. FP gate for families that move arithmetic across a function boundary (F, S `sim_time_seconds`,
   S interpolation, F2, N): rebuild with `-O0` and with `-O3 -ffp-contract=off` and re-run 3 cases;
   and run once on an FMA-capable toolchain (clang/arm64, or x86 `-march=native`), since baseline
   x86-64 GCC never emits FMA and so cannot reveal contraction differences
   (`docs/debranching_plan.md`, "The floating-point constraint is FMA contraction").
6. MPI-touching commits (A, L, M, P, D5 if `need_to_send` is involved): `run_mpi_suite.sh` /
   `mpi_validation` at np=1 and 2 as smoke, invariant checkers on, and `statistical_check.sh`
   against the baseline MPI build for M/P (bitwise is impossible there).
7. Commit message: the differences found, how each is preserved (parameter, guard kept at call
   site, statement kept verbatim), coverage evidence, and "byte-identical on N/N cases".

Defects in section 4 go in their own later commits, each declared as a behaviour change and
validated statistically, per the branch's two-group rule.

---

## Method

- `extract.py`: comment-stripping brace matcher that extracts 871 function bodies (28,391
  significant lines) from `src/ EXEs/ include/ mpi_proto/` (excluding `json.hpp`, `Faddeeva.cpp`)
  into `funcs.json`.
- `similar.py`: all pairs with >=10 lines; 6-token shingle Jaccard/containment prefilter, then
  `difflib` ratio on identifier/number-normalized tokens ("norm", type-2 clones), raw tokens
  ("raw") and whitespace-normalized lines ("lines") -> `sim.txt` (476 pairs).
- `windows.py 8`: verbatim runs of >=8 significant lines shared across files/functions ->
  `win8.txt` (166 runs).
- `reach.py`: transitive name-reachability from both `main`s, class methods and
  `run_code_tests` -> unreachable functions; `dead.py`: zero-reference functions; both confirmed with `grep -w`.
- `common.py A B`: common significant lines between two files/ranges (numbers quoted above).

Top pairs (norm / raw / lines / containment, body lines):

| norm | raw | lines | cont | nA | nB | A | B |
| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 0.95 | 0.94 | 0.93 | 0.94 | 225 | 235 | `src/parser/parse_input.cpp:227` parse_input | `src/parser/parse_input.cpp:519` parse_input_for_add |
| 0.87 | 0.92 | 0.84 | 0.92 | 172 | 182 | `sweep_separation_complex_rot_memtest_cluster_box.cpp:7` | `..._cluster_sphere.cpp:7` |
| 0.89 | 0.89 | 0.52 | 0.98 | 125 | 153 | `initialize_parameters_for_implicitlipid_and_compartment_model.cpp:7` | same file `:169` (`_model`) |
| 0.80 | 0.80 | 0.75 | 0.86 | 215 | 166 | `perform_implicitlipid_state_change_box.cpp:12` | `..._sphere.cpp:13` |
| 0.79 | 0.81 | 0.51 | 0.89 | 170 | 138 | `sweep_separation_complex_rot_memtest_sphere.cpp:18` | `sweep_separation_complex_rot_sphere.cpp:6` |
| 0.79 | 0.78 | 0.56 | 0.81 | 175 | 133 | `associate_ImplicitLipid_box.cpp:11` | `associate_ImplicitLipid_sphere.cpp:14` |
| 0.74 | 0.68 | 0.58 | 0.77 | 216 | 239 | `perform_bimolecular_state_change_box.cpp:9` | `..._sphere.cpp:10` |
| 0.73 | 0.72 | 0.66 | 0.70 | 216 | 215 | `perform_bimolecular_state_change_box.cpp:9` | `perform_implicitlipid_state_change_box.cpp:12` |
| 0.72 | 0.69 | 0.46 | 0.82 | 566 | 359 | `associate_box.cpp:12` | `associate_sphere.cpp:10` |
| 0.96 | 0.96 | 0.91 | 0.96 | 88 | 85 | `src/mpi/deserialize.cpp:132` | `:295` (`_right`) |
| 0.96 | 0.89 | 0.56 | 0.92 | 95 | 92 | `src/mpi/send_receive.cpp:280` receive_right | `:420` receive_left |
| 0.91 | 0.88 | 0.68 | 0.90 | 70 | 82 | `send_receive.cpp:77` send_left | `:173` send_right |
| 0.90 | 0.88 | 0.80 | 0.89 | 69 | 69 | `orient_crds_to_template.cpp:5` | `write_bonded_complex_json.cpp:226` |
| 0.79 | 0.76 | 0.65 | 0.79 | 104 | 96 | `deserialize.cpp:430` | `:623` (`_right`) |
| 0.78 | 0.77 | 0.70 | 0.92 | 51 | 52 | `create_molecule_and_complex_from_rxn.cpp:194` | `create_molecule_and_complex_from_transmission_rxn.cpp:4` |
| 0.86 | 0.85 | 0.72 | 0.81 | 40 | 52 | `generate_coordinates_for_restart.cpp:26` | `create_molecule_and_complex_from_transmission_rxn.cpp:4` |
| 1.00 | 0.98 | 0.94 | 1.00 | 18 | 18 | `initialize_molecule_after_reaction.cpp:72` zeroth | `:111` uni |
| 0.89 | 0.88 | 0.83 | 1.00 | 22 | 24 | `read_restart.cpp:318` createdestruct_from_json | `:342` transmission_from_json |
| 0.81 | 0.81 | 0.79 | 1.00 | 33 | 48 | `src/io/print_complex_hist.cpp:78` | `src/io_mpi/print_complex_hist.cpp:11` |
| 0.89 | 0.84 | 0.38 | 0.70 | 10 | 11 | `get_prevNorm.cpp:3` | `get_prevSurv.cpp:3` |
| 0.99 | 0.98 | 0.83 | 0.91 | 35 | 35 | `functions_implicitlipid.cpp:312` entering | `:349` exiting (not bitwise-mergeable) |
| 0.87 | 0.74 | 0.22 | 0.68 | 57 | 68 | `omega_rotation.cpp:3` | `phi_rotation.cpp:5` (skeleton only) |
| 0.60 | 0.57 | 0.44 | 0.67 | 99 | 71 | `determine_2D_bimolecular_reaction_probability.cpp:5` | `determine_3D_...:4` |
| 0.71 | 0.70 | 0.54 | 0.89 | 136 | 144 | `class_SimulVolume.cpp:385` | `:570` (MPI overload) |
