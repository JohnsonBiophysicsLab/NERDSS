# `nerdss-optimized` test and benchmark results

All numbers below were measured on 2026-08-10. See
[`docs/nerdss_optimized.md`](../../docs/nerdss_optimized.md) for what changed and why.

## Setup

| | |
| --- | --- |
| Host | Apple M5, 10 physical cores, Darwin 25.5.0 arm64 |
| Compiler | Apple clang 21.0.0 (clang-2100.1.1.101), `-O3 -std=c++0x` |
| GSL | 2.8 (Homebrew), generator `mt19937` |
| Build | `make serial` (serial only, no MPI, no OpenMP) |
| Seed | 20260810 for every run unless stated otherwise |

Three builds were compared:

| Label | Contents | Binary SHA-256 (first 16) |
| --- | --- | --- |
| `baseline` | `master` at `260f6e2` | `ef0c0153ef2aa54d` |
| `optA` | issues #8, #9, #11 - the result-preserving subset | `acfc15286bde59f9` |
| `optB` | all five issues (#8, #9, #10, #11, #12) - the branch as delivered | `0051be593d6ea02e` |

`optA` exists only to separate the two kinds of change. Splitting them is what
makes the bitwise claim below possible: because `optA` is byte-identical to
`master`, every difference `optB` shows is attributable to issues #10 and #12
alone.

**The host was not idle.** Steam, WeChat and WindowServer were active throughout,
and one repetition of `homoTrimer` under `optB` took 535 s instead of 6.4 s (the
simulation's own reported wall time confirms the process really was starved, and
two repetitions of the identical run took 6.37 s and 6.40 s). All reported times
are medians over three repetitions, which absorbs that outlier, and the headline
timing table below interleaves the three builds per case so all of them meet the
same machine conditions.

## 1. Bitwise identity: `optA` vs `master`

`run_suite.sh` hashes every file under `DATA/`, `RESTARTS/` and `PDB/`. PDB files
are hashed from line 2 onward because line 1 records the wall-clock creation
time; nothing else in the outputs depends on the clock, which was verified by
running `master` twice and diffing.

**All 13 cases produced byte-identical output.** 8 to 40 output files per case,
including `DATA/COMPLEXES`, all restart files, and every PDB frame body.

| Case | Input | nItr | Bitwise identical |
| --- | --- | --- | --- |
| rev_3D | `VALIDATE_SUITE/bimolecular_reversible/rev_3D/parms3d.inp` | 20,000 | yes |
| rev_2D | `VALIDATE_SUITE/bimolecular_reversible/rev_2D/parms2D.inp` | 200 | yes |
| rev_3Dto2D | `VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/parms3dto2d.inp` | 8,000 | yes |
| homoTrimer | `VALIDATE_SUITE/homoTrimer/parmTri6.inp` | 6,000 | yes |
| closed_homoTrimer | `VALIDATE_SUITE/closed_homoTrimer/parmTri6.inp` | 6,000 | yes |
| hetTrimer | `VALIDATE_SUITE/hetTrimer/parm_autodiff_hetTri.inp` | 8,000 | yes |
| trimer | `VALIDATE_SUITE/trimer/parm_hetTri.inp` | 8,000 | yes |
| hexamer | `VALIDATE_SUITE/hexamer/parms_phex.inp` | 25,000 | yes |
| clathrin | `VALIDATE_SUITE/clathrin/parms_clath_kon1uM.inp` | 150,000 | yes |
| implicit_lipid | `VALIDATE_SUITE/implicit_lipid/parms.inp` | 100,000 | yes |
| sphere | `VALIDATE_SUITE/sphere/parms_sphere.inp` | 40,000 | yes |
| michaelis_menten | `VALIDATE_SUITE/michaelis_menten/michaelis.inp` | 150,000 | yes |
| mem_localization | `VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/parms.inp` | 2,000 | yes |

`optB` differs from `master` on all 13 cases, as it must: issues #10 and #12
change the random stream.

## 2. Timing

Interleaved pass, `interleaved_timing.sh 3`: for each case the three builds run
back to back, then the pattern repeats, so no build gets a quiet or a busy
stretch to itself. Medians of three repetitions, seconds.

| Case | nItr | baseline | optA | optB | optA speedup | optB speedup |
| --- | --- | --- | --- | --- | --- | --- |
| clathrin | 150,000 | 6.285 | 4.920 | 4.032 | 1.277x | **1.559x** |
| hexamer | 25,000 | 8.610 | 7.752 | 6.300 | 1.111x | 1.367x |
| implicit_lipid | 100,000 | 4.280 | 3.864 | 3.220 | 1.108x | 1.329x |
| closed_homoTrimer | 6,000 | 8.614 | 6.647 | 6.490 | 1.296x | 1.327x |
| hetTrimer | 8,000 | 7.685 | 6.285 | 5.949 | 1.223x | 1.292x |
| trimer | 8,000 | 7.677 | 6.436 | 6.054 | 1.193x | 1.268x |
| rev_3D | 20,000 | 5.878 | 6.012 | 4.690 | 0.978x | 1.253x |
| sphere | 40,000 | 11.057 | 9.768 | 8.927 | 1.132x | 1.239x |
| homoTrimer | 6,000 | 8.452 | 7.285 | 7.156 | 1.160x | 1.181x |
| michaelis_menten | 150,000 | 5.234 | 4.788 | 4.476 | 1.093x | 1.169x |
| rev_3Dto2D | 8,000 | 6.841 | 6.417 | 5.851 | 1.066x | 1.169x |
| mem_localization | 2,000 | 2.923 | 2.769 | 2.592 | 1.056x | 1.128x |
| rev_2D | 200 | 39.919 | 40.172 | 39.984 | 0.994x | 0.998x |
| **suite total** | | **123.455** | **113.115** | **105.721** | **1.091x** | **1.168x** |

Reading the spread:

- The result-preserving changes alone (`optA`) buy 1.09x over the suite while
  producing identical bytes. The gain concentrates in models with several
  interfaces and ancillary-interface requirements per reaction - clathrin
  (1.28x) and the trimers (1.19-1.30x) - which is where `find_which_reaction()`,
  `hasIntangibles()` and the `Coord` arithmetic are called most often.
- `rev_3D` at 0.978x and `rev_2D` at 0.994x are the two models with a single
  interface and a single rate state per reaction, so they have nothing for the
  reaction-matching work to save. Both sit within the run-to-run scatter on this
  host; treat them as parity, not as a regression.
- `rev_2D` is dominated by fixed 2D lookup-table construction, roughly 40 s
  regardless of iteration count (1,000 iterations also took ~40 s during
  calibration). Its ratio is close to 1 for all builds because almost none of
  its time is in the code that changed.
- The separate non-interleaved run gives the same picture for `optA`:
  1.02x-1.30x per case, and per-repetition spread under 2% on nearly every case.

## 3. Issue #10: random orientation uniformity

`rng_quality/run.sh`, 4,000,000 samples per sampler. The statistic is the z
component of `(0,0,1)` rotated by the sampled quaternion. For uniformly
distributed rotations that is exactly uniform on [-1,1], and NERDSS depends on
that when it assigns initial molecular orientations. 50 bins.

| Sampler | chi^2 | dof | chi^2/dof | max bin deviation |
| --- | --- | --- | --- | --- |
| old: normalized U(-1,1)^4 | 480,707.6 | 49 | **9,810.36** | **58.80%** |
| new: Shoemake subgroup | 50.2 | 49 | **1.02** | **0.70%** |

The old sampler is not marginally biased. Individual orientation bins were off
by up to 59%, and chi^2/dof of 9,810 rules out uniformity by any margin. The
replacement is statistically indistinguishable from uniform.

## 4. Issue #12: `GaussV()` distribution and speed

Same harness, 20,000,000 samples per sampler, 100 bins over [-5,5] compared
against the standard normal CDF.

| Sampler | mean | variance | skewness | kurtosis | chi^2/dof | outside +-5 | seconds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| old: Marsaglia polar | -0.00016 | 0.99980 | -0.00036 | 3.00084 | 1.01 | 14 | 0.287 |
| new: GSL ziggurat | 0.00018 | 0.99976 | -0.00064 | 2.99849 | 0.68 | 11 | **0.090** |

Both samplers are correct - a standard normal has mean 0, variance 1, skewness 0
and kurtosis 3, and both match to four or five digits. The ziggurat is **3.19x
faster** in isolation. Whole-program, issues #10 and #12 together take the suite
from 1.091x to 1.168x, i.e. they contribute about 7% on top of the
result-preserving work.

## 5. Statistical equivalence of `optB` with `master`

`statistical_check.sh`, 6 independent seeds per build per model. Each run's
species copy numbers are averaged over the second half of the trajectory, those
per-seed averages are combined into a mean and standard error, and the two builds
are compared with a Welch z-score. This is the test that replaces bitwise
comparison for the stream-changing group.

| Model | nItr | Species | master mean +- SEM | optB mean +- SEM | z |
| --- | --- | --- | --- | --- | --- |
| implicit_lipid | 100,000 | `IL(m)` | 342.302 +- 1.119 | 341.642 +- 0.848 | -0.47 |
| | | `B(b)` | 200.000 +- 0.000 | 200.000 +- 0.000 | 0.00 |
| | | `B(m)` | 42.302 +- 1.119 | 41.642 +- 0.848 | -0.47 |
| | | `B(m!1).IL(m!1)` | 157.698 +- 1.119 | 158.358 +- 0.848 | 0.47 |
| michaelis_menten | 150,000 | `S(ser~U)` | 103.954 +- 0.276 | 104.024 +- 0.425 | 0.14 |
| | | `S(ser~P)` | 0.859 +- 0.300 | 0.934 +- 0.196 | 0.21 |
| | | `E(kin)` | 5.812 +- 0.370 | 5.958 +- 0.515 | 0.23 |
| | | `E(kin!1).S(ser~U!1)` | 3.188 +- 0.370 | 3.042 +- 0.515 | -0.23 |
| rev_3D | 60,000 | `A(a)` | 337.396 +- 5.037 | 333.604 +- 3.675 | -0.61 |
| | | `R(r)` | 337.396 +- 5.037 | 333.604 +- 3.675 | -0.61 |
| | | `A(a!1).R(r!1)` | 662.604 +- 5.037 | 666.396 +- 3.675 | 0.61 |

Largest `|z|` across all species and models: **0.61**. The two builds agree well
inside the seed-to-seed scatter, so replacing the two samplers moved the
simulated physics by nothing measurable at this sample size.

## 6. Reverse bimolecular state change (issue #8 correctness part)

This path is reached by no input under `sample_inputs/`, so
[`regression/bimol_state_change`](regression/bimol_state_change) was written for
it: an enzyme-facilitated state change laid out so that
`forwardRxns[1].conjBackRxnIndex == 0`, the arrangement master's
`conjBackRxnIndex > 0` gate rejects outright.

**It could not be validated end to end, because the machinery around it fails on
master for reasons outside issues #8-#12.** Three separate pre-existing defects,
each verified by isolation:

1. `set_rMaxLimit()` inspects only `ReactionType::bimolecular`, so a model whose
   only bimolecular reaction is a state change keeps `rMaxLimit == 0`;
   `SimulVolume::Dimensions` then divides the box length by zero, producing a
   negative cell count and an endless `CELL PAIR MAX EXCEEDED` rescale. Worked
   around in the model by adding a scaffolding association reaction, which
   brings `rMaxLimit` to 19.4339.
2. Executing the state change then produces NaN coordinates within a few hundred
   iterations. Setting the state-change rate to zero removes the NaN and the run
   completes, which localizes the fault to the state-change execution and not to
   the scaffolding.
3. Supplying concrete `assocAngles` avoids the NaN but corrupts the copy counters
   instead: `K(k)` reaches -54 and `K(kd)` reaches 114 with only 30 K molecules
   present, while no reaction is recorded as having fired.

Outcome per build, same model and seed:

| Build | Outcome |
| --- | --- |
| `baseline` | NaN, molecule 22, iteration 139 |
| `optA` | NaN, molecule 22, iteration 139 - identical |
| `optB` | NaN, molecule 0, iteration 210 - same fault, different random stream |

`optA` failing at exactly the same molecule and iteration as master is the useful
result: the issue #8 changes are inert here, consistent with the corrected branch
being unreachable on master. The corrections stand on removing undefined
behavior and on the index convention being provably right, and on the 13
bitwise-identical suite cases showing they change nothing for models that work
today. Making the bimolecular state change usable needs the three defects above
addressed first, each with its own validation.

## 7. Pre-existing failures, unchanged

Three `sample_inputs` cases fail on master; details and reproduction in
[`known_broken.tsv`](known_broken.tsv). Verified across all three builds at
nItr 1000 with a 30 s cap:

| Case | master | optA | optB |
| --- | --- | --- | --- |
| create_destroy | SIGSEGV, 11 output blocks | SIGSEGV, 11 output blocks | SIGSEGV, 22 output blocks |
| unimolecular_reverse | no progress, ~110M `CELL PAIR MAX EXCEEDED` lines in 30 s | same | same |
| clock_model | SIGSEGV at seed 20260810 | SIGSEGV at seed 20260810 | completes at seed 20260810 |

`clock_model` is seed-dependent on every build: it completes on both master and
`optB` with seeds 1 through 5, and crashes on master with seed 20260810. `optB`
surviving that particular seed is a consequence of different initial
orientations, not a fix. Nothing here is caused or repaired by this branch.

## 8. Reverse unimolecular state change

Measured 2026-08-10, same host and compiler as above, seed 20260810 unless
stated. This is a follow-on correction, not part of issues #8-#12: it was
recorded as a known defect and left alone at the time because the path was
believed reachable from existing models.

Unlike sections 1-7, which compare against `master`, this section's baseline is
`nerdss-optimized` at `dbb5cbf`, the branch tip these changes apply to. That is
three commits past the `optB` of section 1, so the binary hashes differ from the
ones there. Two builds:

| Label | Contents | Binary SHA-256 (first 16) |
| --- | --- | --- |
| `prechange` | `nerdss-optimized` at `dbb5cbf`, unmodified | `2b23b2da88e15e2a` |
| `fixed` | `prechange` plus the four corrections in this section | `1c044a13c233de94` |

A third build isolates one of them, because the display crash otherwise stops the
new model before the simulation begins:

| Label | Contents | Binary SHA-256 (first 16) |
| --- | --- | --- |
| `displayfix` | `prechange` plus the `BackRxn::display()` fix only | `20c79d910a2e9c88` |

The measurements below were first taken against `4f3cb39` and then repeated
against `dbb5cbf` after rebasing. Every number was unchanged, including the
per-species equilibrium fractions, so the intervening math-function commits do
not touch this path or its random stream.

### 8.1 No existing model changes

**All 13 suite cases are byte-identical between `prechange` and `fixed`.** 268
files hashed per build under `DATA/`, `RESTARTS/` and `PDB/`; the two
`manifest.sha256` files are byte-identical as whole files, and every case exited
0. `michaelis_menten` is one of the 13 and does use a unimolecular state change.

`auto_phos` is not in the suite. Run separately at nItr 200,000, all 14 `DATA/`
files and all `RESTARTS/` files are byte-identical between the two builds.

That is expected rather than lucky, and the reason is checkable: **no input under
`sample_inputs/` declares a unimolecular state change with `<->`.**
`michaelis_menten`, `unimolecular_reverse` and `auto_phos` write theirs as separate
`->` reactions, so `conjBackRxnIndex` is `-1` and the reverse branch is never
entered. (`enzyme` has no unimolecular state change at all; its one state change
is bimolecular - see section 8.4.) The earlier note that this path was "reachable
from existing models" was too pessimistic: the *function* is reached constantly,
but its reverse branch is not reached by anything in the repository.

No timing claim is made here. The suite ran while other measurements shared the
machine, and the per-case numbers from that pass are load noise, not signal.

### 8.2 The new model, and what it separates

[`sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible`](../../sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible)
holds two reversible unimolecular state changes with `kf = 1000` and
`kb = 3000` s-1, so each must equilibrate at `kf/(kf+kb) = 0.25` in the P state.
An irreversible rate-0 reaction at forward index 1 offsets the two index spaces,
putting the state changes at forward 2 and 3 but back 1 and 2 - an arrangement a
single reversible pair cannot produce.

Time-averaged over the equilibrated second half of the run:

| build | A(ser~P) | B(thr~P) | A(tag~P) | expected | verdict |
| --- | --- | --- | --- | --- | --- |
| `prechange` | - | - | - | 0.25 | SIGSEGV during reaction display |
| `displayfix` | 0.765 | 0.753 | 0 | 0.25 | FAIL, P is absorbing |
| `fixed` | 0.229 | 0.269 | 0 | 0.25 | PASS |

The middle row is the informative one. With the crash removed but the index
assignments untouched, the model runs and the reverse direction simply never
fires: the populations decay one way from 200 U toward 200 P, and 0.765 is where
that decay has got to after about 8 time constants, not an equilibrium. The
corrected build instead settles at the predicted fraction, a factor of four away.

`A(tag~P)` staying at 0 on every build confirms the rate-0 offset reaction never
fires, so it perturbs only the index spaces and not the dynamics.

Across seeds 20260810, 1, 2, 3 and 7 the corrected build gives P fractions of
0.229/0.269, 0.229/0.260, 0.224/0.249, 0.257/0.241 and 0.254/0.233 against a
predicted 0.25, with a binomial spread of 0.031 on 200 copies. The agreement is
not seed-specific. The mean of those ten values is 0.244.

Reproduce with:

```bash
./benchmarks/nerdss_optimized/regression/unimol_state_change_reversible/check.sh bin/nerdss
```

which exits 0 only if both state changes reach `kf/(kf+kb)` and the rate-0
reaction stayed silent. It exits 1 on `displayfix` and on `prechange`.

### 8.3 The display crash was flaky, which is worth recording

`BackRxn::display()` read one element past the end of `rate.otherIfaceLists`.
On the new model that crashed 4 runs in 5 with an identical seed and identical
inputs, and completed the fifth, because whether the read lands on mapped memory
depends on heap layout. It also completed every time under `lldb`. A defect that
disappears under a debugger and passes one run in five is the kind that gets
attributed to the model rather than the code, so it is called out here
explicitly.

### 8.4 Two models outside the suite, and a correction to section 6

`auto_phos` and `enzyme` both use interface states and neither is in
`cases.tsv`, so both were run separately against `prechange` and `fixed`:

| model | nItr | outputs compared | result |
| --- | --- | --- | --- |
| `sample_inputs/auto_phos/autophos_D10.inp` | 200,000 | 14 `DATA/` files + `RESTARTS/` | byte-identical |
| `sample_inputs/enzyme/parms_clat_enzyme.inp` | 20,000 | all `DATA/` + `RESTARTS/` | byte-identical |

`enzyme` also corrects section 6 and the corresponding claim in
`docs/nerdss_optimized.md`. Both said no input under `sample_inputs/` contains a
bimolecular state change. `parms_clat_enzyme.inp` line 116,
`syn(pi) + pip2(head~U) -> syn(pi) + pip2(head~P)`, is one; the build's own
reaction dump reports `Type: Bimolecular state change` for it, alongside 22
bimolecular associations.

It does not change section 6's conclusion. The reaction is irreversible, so
`conjBackRxnIndex == -1`, the corrected reverse gate rejects it, and the
`perform_bimolecular_state_change_*` change is inside the `isStateChangeBackRxn`
branch, which that reaction cannot reach. The issue #8 corrections stay inert for
it, which the byte-identical run above confirms directly. What was wrong was the
justification, not the verdict.

The reason it went unnoticed is worth recording: `enzyme` appears in neither
`cases.tsv` nor `known_broken.tsv`, so the one sample input that exercises the
forward bimolecular state-change path is not covered by the suite. It runs clean
for 20,000 iterations, so adding it to `cases.tsv` is cheap.

## 9. Profile-driven follow-up batch: three result-preserving optimizations

A CPU profile of the branch (not of `master`) put most of the remaining time in
places none of issues #8-#12 touch. Three of them were addressed as one batch.
All three are result-preserving, and the suite confirms that directly.

### 9.1 What changed

| # | Change | Files |
| --- | --- | --- |
| 1 | Cache `1/cbrt(Dr)` per axis on `MolTemplate` as `invCbrtDr`, instead of calling `pow(Dr, 1.0/3.0)` for every member molecule of every complex on every timestep. `MolTemplate::Dr` is written only by `set_value()` at parse time, so the quantity is a per-template constant. `cache_diffusion_derivatives()` refreshes it, and is called from both `set_value()` and `deserialize()`; the value is derived, so it is deliberately not added to the MPI wire format. The cached expression is written exactly as the inline one was, so the stored double is bit-for-bit what was computed before. | `class_MolTemplate.{hpp,cpp}`, `class_Molecule_Complex.cpp` |
| 2 | `SimulVolume::update_memberMolLists()` emptied every sub-cell on every timestep. Added `occupiedSubCells` plus `clear_member_lists()`, so only sub-cells that actually hold members are cleared. For `clathrin` that is 100 molecules against 2744 cells: a 494 nm box divided by a 33.7 nm interaction limit gives 14^3 cells, so the old sweep did 27 cell-clears per molecule. Serial overload only -- the `MpiContext` overload keeps its full sweep, because the MPI ranks also mutate `memberMolList` directly from `prepare.cpp` and `deserialize.cpp` and that path is not exercised here. | `class_SimulVolume.{hpp,cpp}` |
| 3 | Issue #11's treatment applied to `Vector`, which never received it: trivial constructors, `dot`, `cross`, `calc_magnitude` and `normalize` moved into the header so they inline; `Vector(Coord)` takes `const Coord&` rather than a by-value copy; read-only operators marked `const`; `operator/=` now modifies `*this` and returns a reference. | `class_Vector.{hpp,cpp}` |

`operator+(const Coord&)` was deliberately left non-`const`. The free
`operator+(const Vector&, const Coord&)` is an equally exact match for a `const
Vector` but returns `Coord` rather than `Vector`, so which overload a call site
selects -- and therefore the type of the result -- is currently decided by
whether the left operand is `const`. Marking the member `const` makes the two
ambiguous, which is how the pair was found. Untangling it changes result types at
call sites, so it does not belong in a result-preserving change.

### 9.2 Two defects that had to be fixed first

Neither is part of the batch; both blocked measuring it.

**The Makefile did not track header dependencies.** `make` compared each `.o`
only against its `.cpp`, so editing a header rebuilt nothing. Because candidates
1 and 2 add a member to a struct declared in a header, an incremental build
produced an executable in which some translation units used the new layout and
the rest still used the old one. It linked, and the first test binary built that
way was discarded. Fixed by adding `-MMD -MP` and `-include $(OBJS:.o=.d)`. Any
header change made before this fix was silently at risk.

**`Membrane` had uninitialized members that `write_restart()` prints.**
`nSites`, `No_free_lipids`, `No_protein`, `totalSA`, `Dx`-`Drz` and `offset` had
no initializers, and line 31 of `write_restart.cpp` prints four of them
unconditionally. Measured behaviour: the indeterminate values are *deterministic
for a given binary* -- the same executable run twice produces byte-identical
`restart.dat` -- but they change whenever a code change perturbs memory layout.

That is why the first isolated comparison of this batch reported
`bitwise_identical NO` for all 13 cases while every other output file was
byte-identical: the only file that differed was `DATA/restart.dat`, in the
`membrane = ` line, e.g. `-1 -492228480 0 -492309832 -208766624 244036` against
`-1 -490136896 0 -492229096 -2023371160 244036`. The three differing fields are
exactly the three uninitialized ints.

The consequence for the harness is worth stating plainly: **before this fix,
`compare_suites.sh` would report a false `NO` for any change that shifts memory
layout, including a purely result-preserving one.** Fixed with `{ 0 }`
initializers. Both fixes are applied to the reference and the candidate build
alike, so they cancel out of the comparison below.

### 9.3 Bitwise identity

Reference: `4baa42f` + the two section 9.2 fixes. Candidate: the same, plus
candidates 1-3. Both built with Apple clang 21, `-O3`, seed 20260810, in a
dedicated `git worktree` so that unrelated concurrent edits in the main working
tree could not contaminate either binary.

`run_suite.sh` + `compare_suites.sh`, all 13 cases in `cases.tsv`, comparing
every file under `DATA/`, `PDB/` and `RESTARTS/` (PDB from line 2, since line 1
carries a wall-clock timestamp):

| case | nItr | bitwise identical |
| --- | --- | --- |
| `clathrin` | 150,000 | yes |
| `closed_homoTrimer` | 6,000 | yes |
| `hetTrimer` | 8,000 | yes |
| `hexamer` | 25,000 | yes |
| `homoTrimer` | 6,000 | yes |
| `implicit_lipid` | 100,000 | yes |
| `mem_localization` | 2,000 | yes |
| `michaelis_menten` | 150,000 | yes |
| `rev_2D` | 200 | yes |
| `rev_3D` | 20,000 | yes |
| `rev_3Dto2D` | 8,000 | yes |
| `sphere` | 40,000 | yes |
| `trimer` | 8,000 | yes |

13 of 13 byte-identical.

### 9.4 Timing

`interleaved_timing.sh`, 5 repetitions, median per case. The interleaved harness
is used rather than `run_suite.sh`'s own timings, which are quantized to roughly
a second on this host and produced a non-credible 1.998x for `clathrin`.

| case | nItr | ref (s) | candidate (s) | speedup |
| --- | --- | --- | --- | --- |
| `clathrin` | 150,000 | 3.979 | 3.538 | 1.125 |
| `closed_homoTrimer` | 6,000 | 6.581 | 6.468 | 1.017 |
| `hetTrimer` | 8,000 | 5.870 | 5.760 | 1.019 |
| `hexamer` | 25,000 | 6.025 | 5.437 | 1.108 |
| `homoTrimer` | 6,000 | 6.403 | 6.089 | 1.052 |
| `implicit_lipid` | 100,000 | 3.261 | 2.797 | 1.166 |
| `mem_localization` | 2,000 | 2.569 | 2.458 | 1.045 |
| `michaelis_menten` | 150,000 | 4.407 | 4.390 | 1.004 |
| `rev_2D` | 200 | 43.057 | 40.658 | 1.059 |
| `rev_3D` | 20,000 | 5.514 | 4.990 | 1.105 |
| `rev_3Dto2D` | 8,000 | 5.934 | 5.632 | 1.054 |
| `sphere` | 40,000 | 9.161 | 6.070 | 1.509 |
| `trimer` | 8,000 | 6.121 | 5.879 | 1.041 |
| **total** | | **108.882** | **100.166** | **1.087** |

`sphere` gains most (1.509x) because it has the largest cell-count-to-molecule
ratio, which is what candidate 2 addresses. `michaelis_menten` is flat (1.004x),
as expected: it is dominated by reaction bookkeeping rather than by cell
maintenance or diffusion-constant updates.

### 9.5 The speedup is where the profile said it would be

`sample(1)` on `clathrin`, 8 s window, 1 ms period, same seed and input for both
builds. Self-time share of total samples:

| frame | ref (6344 samples) | candidate (5493 samples) |
| --- | --- | --- |
| `pow` | 4.2% | absent |
| `SimulVolume::update_memberMolLists` | 6.3% | 2.1% |
| `Vector::Vector(Coord)` | 3.1% | absent |
| `Complex::update_properties` | 3.4% | 3.3% |

Candidate 1 removes `pow` entirely; `update_properties` barely moves, so its
remaining cost is the rest of the function rather than the cube root. Candidate 2
cuts cell maintenance by three times. Candidate 3 inlines `Vector`'s constructor
out of existence.

The three eliminated shares sum to about 11.5 points of the `clathrin` profile,
against a measured 11.1% wall-clock gain on the same case. The profile and the
clock agree.

### 9.6 One hypothesis that was tested and rejected

The profile showed roughly 7% of runtime in `malloc`/`free`. The obvious suspect
was `clear_reweight_vecs()`, which copy-assigns six `std::vector`s per molecule
per timestep and then clears the source -- the canonical case for swapping
instead of copying, and provably result-preserving.

Implemented and benchmarked across all 13 cases: **0.990x**, i.e. no gain, within
noise. The vectors are almost always empty, so the copy-assignment never
reallocates and the swap saves nothing. Reverted.

The allocator traffic is therefore still unexplained and still worth about 5-6%.
It sits in two adjacent call sites inlined into `main`, which a `-g` build would
resolve to source lines. That is the next thing to scope, not to assume.

### 9.7 The committed state reproduces the measured state

Section 9.3's comparison was made in a worktree, so it is worth confirming that
what landed on the branch is what was measured. Built from committed
`nerdss-optimized` (`0784dd8`) and re-run over `cases.tsv`: **byte-identical to
the validated candidate binary on all 13 cases.**

That comparison also crosses a second variable, incidentally. The validated
candidate was built from `4baa42f` and contains none of the `Quat` rewrite in
`aea2003`/`656a2d9`/`cd22362`; the committed build contains all of it. Their
outputs match byte-for-byte on all 13 cases, so that rewrite is result-preserving
as well, at least at this seed and these iteration counts. That was not expected:
replacing `inverse()` with the conjugate for a unit quaternion drops a division
by a `norm()` that is only 1 to within rounding, which would normally perturb low
bits. Worth a separate look rather than being taken as general.
