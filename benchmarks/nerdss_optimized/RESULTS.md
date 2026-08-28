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

## 10. Pairwise serial benchmark of current branch against mainline

Measured on 2026-08-11 on the same Apple M5 host, with Apple clang 21.0.0,
GSL 2.8 and the Makefile's serial `-O3 -std=c++0x` build. This repository calls
its mainline branch `master`; there is no `main` ref, so `master` is the baseline
for this comparison.

| Label | Revision | Binary SHA-256 (first 16) |
| --- | --- | --- |
| `master` | `260f6e2a861b` | `ef0c0153ef2aa54d` |
| `optimized` | `a7435a0690ab` | `4f981984ef6833a9` |

Both revisions were exported into isolated temporary source trees and compiled
from committed files. The unrelated uncommitted reaction-table work in the main
working tree was therefore excluded.

### 10.1 Paired total runtime

Three inputs were selected to cover ordinary 3-D reversible binding,
mixed-dimensional implicit-lipid binding, and state-change chemistry. Each pair
used the same input, iteration count and seed (`20260810`). Runtime is external
wall time for the entire serial process, including initialization, simulation
and output.

There were 10 adjacent pairs per case. Five pairs ran `master` first and five ran
`optimized` first, balancing order effects. Speedup is `master / optimized` for
each pair. The table reports the median paired ratio as the robust point estimate
and, separately, the geometric mean with a two-sided 95% Student-t interval on
the 10 log ratios.

| Case | nItr | master median (s) | optimized median (s) | median paired speedup | geometric mean speedup (95% CI) | optimized faster |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `rev_3D` | 60,000 | 19.200 | 13.429 | **1.403x** | 1.417x (1.375x-1.460x) | 10/10 |
| `implicit_lipid` | 100,000 | 5.716 | 3.695 | **1.557x** | 1.560x (1.544x-1.576x) | 10/10 |
| `michaelis_menten` | 150,000 | 6.389 | 5.463 | **1.186x** | 1.220x (1.107x-1.346x) | 10/10 |
| **three-case total** | | **31.087** | **22.426** | **1.379x** | **1.396x (1.355x-1.438x)** | **10/10** |

For each model, 10/10 paired wins has an exact two-sided sign-test `p = 0.00195`.
The total row sums the three model times within each repetition before forming
its paired ratio; it also has 10/10 wins and the same sign-test result.

One `master` Michaelis-Menten repetition took 10.062 s while its other nine took
5.446-6.770 s. No run failed and its adjacent optimized run was normal. It was
retained as scheduling noise rather than deleted. The median estimate is robust
to it, while the wider geometric-mean interval transparently includes its
effect.

The six independent-seed runs used for the output comparison below also report
NERDSS's internal whole-program wall time, measured from before input parsing to
the end of the simulation. Extracting those values gives a second timing result
in which every pair follows a different stochastic trajectory. The interval is
again computed on paired log ratios, now with six seed pairs.

| Case | Seeds | master median (s) | optimized median (s) | median paired speedup | geometric mean speedup (95% CI) | optimized faster |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `rev_3D` | 6 | 17.655 | 12.722 | **1.389x** | 1.383x (1.345x-1.422x) | 6/6 |
| `implicit_lipid` | 6 | 5.135 | 3.289 | **1.549x** | 1.553x (1.511x-1.596x) | 6/6 |
| `michaelis_menten` | 6 | 6.839 | 5.173 | **1.301x** | 1.282x (1.170x-1.405x) | 6/6 |
| **three-case total** | **6** | **29.538** | **21.195** | **1.381x** | **1.384x (1.360x-1.409x)** | **6/6** |

Six wins out of six gives an exact two-sided sign-test `p = 0.03125` per model
and for the total. More importantly, the independent-seed total (1.381x median)
agrees with the balanced fixed-seed timing result (1.379x median). The speedup is
therefore not an artifact of one unusually cheap trajectory from the changed
random stream.

### 10.2 Stochastic-output comparison

Bitwise output comparison is not appropriate because the orientation and
Gaussian samplers intentionally change the random stream. Instead, each build
ran each model with six independent seeds (`20260811` through `20260816`). For
each run, every reported species copy number was averaged over the second half
of the trajectory. The table combines those six per-seed averages as mean +-
SEM and reports `(optimized - master) / sqrt(SEM_master^2 + SEM_optimized^2)`.

| Model | Species | master mean +- SEM | optimized mean +- SEM | standardized difference |
| --- | --- | ---: | ---: | ---: |
| `implicit_lipid` | `IL(m)` | 342.302 +- 1.119 | 341.642 +- 0.848 | -0.47 |
| | `B(b)` | 200.000 +- 0.000 | 200.000 +- 0.000 | 0.00 |
| | `B(m)` | 42.302 +- 1.119 | 41.642 +- 0.848 | -0.47 |
| | `B(m!1).IL(m!1)` | 157.698 +- 1.119 | 158.358 +- 0.848 | 0.47 |
| `michaelis_menten` | `S(ser~U)` | 103.954 +- 0.276 | 104.024 +- 0.425 | 0.14 |
| | `S(ser~P)` | 0.859 +- 0.300 | 0.934 +- 0.196 | 0.21 |
| | `E(kin)` | 5.812 +- 0.370 | 5.958 +- 0.515 | 0.23 |
| | `E(kin!1).S(ser~U!1)` | 3.188 +- 0.370 | 3.042 +- 0.515 | -0.23 |
| `rev_3D` | `A(a)` | 337.396 +- 5.037 | 333.604 +- 3.675 | -0.61 |
| | `R(r)` | 337.396 +- 5.037 | 333.604 +- 3.675 | -0.61 |
| | `A(a!1).R(r!1)` | 662.604 +- 5.037 | 666.396 +- 3.675 | 0.61 |

The largest absolute standardized difference is **0.61**. No observable is near
the conventional magnitude-2 screening threshold, so this sample finds no
measurable shift in the simulated species distributions. Six seeds per build is
a useful regression check, not proof that the full stochastic processes are
identical; more seeds would narrow the SEMs if tighter equivalence bounds are
needed.

The focused timing list is in `pairwise_cases.tsv`. Reproduce either execution
order by setting `CASES_FILE` when invoking `interleaved_timing.sh`; reproduce
the ensemble comparison with `statistical_check.sh`.

## 11. `find_which_reaction()`: attempted, measured, reverted

`find_which_reaction()` was the largest remaining self-time entry not addressed by
sections 1-9. An attempt to optimize it produced correct code with no speedup,
and was reverted. Recorded because the reason it failed is more useful than the
attempt, and because the same reasoning error had already occurred once in
section 9.6.

### 11.1 Where the time actually goes

Built with `-O3 -g` and profiled with `sample(1)` on `clathrin`, then mapped
offsets to lines with `atos`. `find_which_reaction()` held 582 of 5436 samples,
10.7% of self time, distributed as:

| line | samples | statement |
| --- | --- | --- |
| 18 | 67 | `const ForwardRxn& oneRxn = forwardRxns[rxnItr];` |
| 30 | 66 | `if (absIface1 == oneReactant.absIfaceIndex)` |
| 17 | 45 | `for (auto rxnItr : currState.myForwardRxns)` |
| 22 | 43 | `reactItr < oneRxn.reactantListNew.size()` |
| 24 | 36 | `molTemplateList[oneReactant.molTypeIndex].isImplicitLipid` |
| 121 | 61+51 | `best_matching_rate(...)`, inlined |
| 78-79 | 52+21 | reverse state-change gate, never taken for this model |
| 120, 140, 144 | 33+8+52 | match test and return |
| prologue | 55+24 | |

Grouped: the reactant-matching scan is about 55% of the function, the rate
selection about 33%, and the reverse state-change gate about 12% -- the last
being pure cost for `clathrin`, which has no bimolecular state change but still
tests for one on every call.

### 11.2 The attempt

The three hottest lines each dereference something: into the large `ForwardRxn`,
out of it into `reactantListNew`'s own heap block, and into the much larger
`MolTemplate` for one bool. That reads like a memory-latency problem, so the scan
was given a compact side table: `CompactReactant` holding just
`{absIfaceIndex, isImplicitLipid}` at 8 bytes, eight to a cache line, with all of
a state's reactant slots stored end to end, plus per-candidate copies of the
static fields the gates read (`conjBackRxnIndex`, `rxnType`, `relRxnIndex`, and
the symmetric-bimolecular test). `ForwardRxn` was then touched only when a
reaction had already matched, or when the reverse branch was entered. The scan's
comparisons were left textually unchanged so it selects the same reactant slots.

The table was keyed on `Interface::State::index`, built on first use, and
rebuilt if a state's `myForwardRxns` length ever changed. Reaction lists are
final before the timestep loop -- the add-file path runs at `nerdss.cpp:400`, the
loop at line 811 -- so a table built lazily cannot go stale mid-run.

### 11.3 Result: correct, and no faster

Bitwise: **13 of 13 cases byte-identical**, so the rewrite preserved behavior.

Timing, CPU time (user+sys) rather than wall clock, interleaved, 7 repetitions,
median per case:

| case | nItr | base (s) | table (s) | ratio | base sd | table sd |
| --- | --- | --- | --- | --- | --- | --- |
| `clathrin` | 150,000 | 4.440 | 4.500 | 0.987 | 0.050 | 0.033 |
| `homoTrimer` | 6,000 | 8.720 | 8.790 | 0.992 | 0.298 | 0.144 |
| `michaelis_menten` | 150,000 | 4.930 | 4.880 | 1.010 | 0.208 | 0.233 |
| `rev_3D` | 20,000 | 5.510 | 5.450 | 1.011 | 0.254 | 0.188 |
| **total** | | **23.600** | **23.620** | **0.999** | | |

Every ratio is inside one standard deviation. The decisive measurement is the
function's own share: `find_which_reaction()` went from 582/5436 samples (10.7%)
to 622/5532 (11.2%). The table made the function it was meant to speed up
slightly *slower*.

An intermediate version was genuinely wrong and was fixed before this
measurement: the table accessor lived in its own translation unit and neither
build file enables link-time optimization, so it compiled to an out-of-line call
added to the hottest path in the program. Splitting it into an inlinable fast
path in the header and an out-of-line builder changed nothing measurable, which
is what settled the question.

### 11.4 Why it failed, and what that says about reading profiles

Line-level sample attribution shows which instructions **retire**, not which ones
**stall**. This model has one or two candidate reactions per state, each with two
reactant slots, and `find_which_reaction()` is called millions of times per
timestep against that same handful of bytes. The reaction data is therefore
L1-resident in steady state: the loads on lines 18, 24 and 30 were already hits,
so making them contiguous removed no latency while the table lookup added work.

This is the same error as section 9.6's `clear_reweight_vecs()` hypothesis --
inferring a memory bottleneck from the fact that the hot lines contain
dereferences. Both times the fix was cheap to build and the measurement took
longer than the implementation. The pattern worth keeping: a line-level profile
localizes cost, it does not explain it, and on this codebase the explanation has
been wrong twice.

### 11.5 What would actually help

`find_which_reaction()` is **call-count-bound, not per-call-cost-bound**. There is
no per-call fat left worth removing; the lever is to call it less.

The opening is structural. `check_bimolecular_reactions.cpp:96` already tests
`absIface2 == statePartner` before calling, so the caller knows *which* partner
matched -- and then `find_which_reaction()` re-derives which reaction that partner
belongs to by scanning. If `State::rxnPartners` were parallel to
`State::myForwardRxns`, the caller could pass the reaction index and the scan
would disappear rather than get cheaper, moving the work to parse time with no
table and no lifetime question.

They are currently **not** parallel, which is what makes this a real change rather
than a cleanup: `populate_reaction_lists.cpp` pushes to `myForwardRxns` alone at
lines 17 and 24, pushes to `rxnPartners` alone at lines 39 and 48, and pushes to
both together at lines 64-65, 68-69, 76-77 and 80-81. Establishing and maintaining
the parallel invariant is the work, and it belongs in its own issue with its own
verification, because getting it wrong silently mis-selects reactions.

### 11.6 A measurement caveat worth recording

The first timing pass for this attempt reported 0.821x aggregate with `rev_2D` at
133 s against 43 s in section 9.4. That was contamination, not a regression: load
average 7.9 on a ten-core host, with an unrelated application at 96% CPU, a
`python3` process at 96%, and a second concurrent NERDSS benchmark at 86%. It was
discarded and re-measured on CPU time, which is far less sensitive to competing
load than wall clock.

Separately, `interleaved_timing.sh` aborted mid-run with
`line 69: 5: command not found` because the script was being edited in another
working session while it ran; it has since gained a `CASES_FILE` override. The
medians were recovered from the raw `timings.tsv` rather than by re-running. Both
incidents argue for checking `uptime` before trusting a wall-clock number on a
shared machine.

## 12. How much of section 9's speedup is an artifact of small test cases

Section 9.4 reports 1.087x aggregate over `cases.tsv`. That number is partly a
property of the suite rather than of the code, and this section quantifies which
part.

### 12.1 The suite's cases are dilute, and that is what candidate 2 rewards

Candidate 2 replaced an O(sub-cells) sweep with an O(occupied sub-cells) one, so
its benefit is governed by sub-cells per molecule. Measured per case, against the
section 9.4 speedups:

| case | reactions | molecules | sub-cells | cells/mol | speedup |
| --- | --- | --- | --- | --- | --- |
| `clathrin` | 7 | 100 | 2744 | 27.44 | 1.125 |
| `implicit_lipid` | 1 | 201 | 3375 | 16.79 | 1.166 |
| `sphere` | 1 | 501 | 8000 | 15.97 | 1.509 |
| `rev_3D` | 1 | 2000 | 27000 | 13.50 | 1.105 |
| `hexamer` | 1 | 1000 | 10648 | 10.65 | 1.108 |
| `rev_3Dto2D` | 1 | 3850 | 27000 | 7.01 | 1.054 |
| `mem_localization` | 3 | 3955 | 10830 | 2.74 | 1.045 |
| `rev_2D` | 1 | 1600 | 3600 | 2.25 | 1.059 |
| `homoTrimer` | 1 | 1000 | 1000 | 1.00 | 1.052 |
| `closed_homoTrimer` | 1 | 1000 | 1000 | 1.00 | 1.017 |
| `michaelis_menten` | 2 | 117 | 64 | 0.55 | 1.004 |
| `hetTrimer` | 3 | 999 | 512 | 0.51 | 1.019 |
| `trimer` | 3 | 300 | 36 | 0.12 | 1.041 |

The six cases at cells/mol >= 7 span 1.054-1.509x; the seven at <= 2.74 span
1.004-1.059x. Six of thirteen cases hold 1000 molecules or fewer.

### 12.2 Direct test: same model, same box, twenty times the molecules

`clathrin` with the copy number raised from 100 to 2000, everything else
unchanged, so the sub-cell count stays at 2744 and only density moves. 4000
iterations, CPU time, 5 repetitions, median:

| copies | molecules | sub-cells | cells/mol | ref (s) | candidates 1-3 (s) | speedup |
| --- | --- | --- | --- | --- | --- | --- |
| 100 | 100 | 2744 | 27.44 | 0.130 | 0.110 | **1.182** |
| 2000 | 2000 | 2744 | 1.37 | 24.650 | 24.030 | **1.026** |

Same code, same model, same box: the gain falls from 18% to 2.6% on density
alone. The 1.182x at 100 copies also cross-checks section 9.4's 1.125x for the
same model at 150,000 iterations, so it is not a start-up artifact, and the
2000-copy runs are 24 s of CPU, far past any parse overhead.

This is structural, not incidental. `class_SimulVolume.cpp:80-85` caps each
dimension at 30, so the sub-cell count cannot exceed 27,000 however large the
system gets, while molecule count scales with the model. Sub-cells per molecule is
therefore bounded by 27000/N and falls as N grows. A production run of 10^4 to
10^5 molecules sits near or below 1 cell per molecule, i.e. in the ~1.02x regime,
not the ~1.2x one.

### 12.3 What does carry over

The residual 1.026x at 2000 copies is candidates 1 and 3 still working. Both cost
scale per molecule per timestep rather than per sub-cell: candidate 1 removes one
`pow` per member molecule of every moved complex, candidate 3 removes an
out-of-line call per `Vector` construction. Their share of runtime is roughly
independent of N, so they should hold at production scale. Candidate 2 is the one
that should be expected to fade.

The honest summary of section 9.4 is therefore: **about 1.02-1.03x is robust to
system size, and the remainder is a dilute-system effect that the suite
over-represents.**

### 12.4 The same question applied to section 11's failure

For `find_which_reaction()` the relevant notion of "small" is the reaction
network, not the molecule count. The scan's per-call cost depends on
`myForwardRxns.size()` and `reactantListNew.size()` (always 2), and neither
depends on molecule count or iteration count -- a larger system makes more calls
at the same cost each. So section 11's verdict is unaffected by N.

It is affected by network size, and the suite is small there too: nine of thirteen
cases define a single reaction, `clathrin` defines 7, and the richest input in the
entire `sample_inputs` tree, `enzyme`, defines 12 across 4 molecule types. The
working set the scan touches is a few hundred bytes. Pushing it out of L1 would
take on the order of hundreds of reactions per interface state, which is a
different class of model than anything in this repository.

One caveat is worth stating rather than glossing: at high density the pair loop
streams through a much larger `moleculeList`, and `hasIntangibles()` reads each
molecule's `interfaceList`, so the reaction data could plausibly be evicted
between calls in a way it is not at 100 molecules. That is the one condition under
which the reverted table might pay off, and it was **not** tested -- the density
experiment in 12.2 was run against the reverted tree. Anyone revisiting section 11
should start there, with a 2000-copy case rather than the suite defaults.

## 13. Is any profile frame worth parallelizing with OpenMP?

Asked of the section 11.1 profile, and answered by measuring both sides of the
trade rather than by inspection. Conclusion: **no frame on this workload**, and for
the one frame with clean inner parallelism the ceiling is too low to matter even
in the limit. This branch contains no OpenMP and none was added; the investigation
ran in a detached worktree and no source file changed.

### 13.1 Self time on the current branch

`clathrin`, 150,000 iterations, `sample(1)`, 5436 samples:

| frame | self | parallel over | verdict |
| --- | --- | --- | --- |
| `check_bimolecular_reactions` | 19.6% | candidate pairs | already built on `openmp-production-lane`; 0.72x at 8 workers, 0.54x at 10, on this host |
| `main` (self) | 14.7% | per-molecule bookkeeping | independent, but one region per timestep x 150,000 |
| `find_which_reaction` | 10.7% | none inside (1-7 reactions) | pair level only, i.e. frame 1 |
| `Complex::propagate` | 8.2% | complex members | analysed below |
| allocator | ~8.0% | none | threads increase contention |
| RNG (`ziggurat` + `mt_get`) | 5.7% | none | changes the stream; ends reproducibility |
| `get_distance` | 5.0% | per pair | frame 1 |
| libm `cos`/`exp` | ~4.7% | per complex | SIMD, not OpenMP |
| `determine_3D_bimolecular_reaction_probability` | 3.2% | per pair | frame 1 |
| `Complex::update_properties` | 3.0% | per complex | inside propagate |
| `clear_reweight_vecs` | 2.2% | per molecule | as `main` |

Frames 1, 3, 7 and 9 are one opportunity, not four: all sit inside the same pair
loop, which needs conflict-free scheduling to stay deterministic because pairs
share molecules and complexes. That is what `openmp-production-lane` implements,
and it measured a regression on this host. `Complex::propagate`'s member loop is
the only frame whose inner iterations are independent with no shared writes, no
RNG and no ordering constraint.

### 13.2 Cost of a parallel region on this host

`#pragma omp parallel for reduction`, Apple M5, Apple clang 21, Homebrew `libomp`,
200,000 repetitions per row, overhead measured against the identical serial loop:

| threads | wait policy | width 8 | width 128 | width 2048 |
| --- | --- | --- | --- | --- |
| 10 | `passive` (default) | 55,866 ns | 55,694 ns | 80,114 ns |
| 10 | `active` | 6,777 ns | 10,166 ns | 4,891 ns |
| 4 | `active` | 624 ns | 630 ns | **-211 ns** |
| 2 | `active` | 388 ns | 396 ns | **-59 ns** |

The negative entries are real wins: at width 2048 two or four threads beat serial.
The crossover exists, it is simply far above the widths this program produces.

### 13.3 Cost of one member iteration

`Complex::propagate` is 8.19% of `clathrin`'s 4.44 s of CPU, so 0.364 s, spread
over roughly 100 molecules x 150,000 timesteps = 1.5e7 member iterations, giving
**about 24 ns per member**. This is an order-of-magnitude figure: it assumes every
complex propagates on every timestep, and if fewer do then the true per-member
cost is higher, which raises the break-even width and strengthens the conclusion.

Solving `24W = overhead + 24W/T` for the break-even width `W`:

| threads | wait policy | break-even width |
| --- | --- | --- |
| 2 | `active` | ~32 members |
| 4 | `active` | ~34 members |
| 10 | `active` | ~231 members |
| 10 | `passive` | ~2600 members |

### 13.4 How wide the loop actually gets

`gagsphere` -- 2500 gag molecules, 436 nm box, timestep 0.1, the most
assembly-heavy input in `sample_inputs` -- run for 300,000 iterations. Percentages
are of total propagate work, i.e. weighted by member count:

| sim time (s) | complexes | max size | mean width | %work >=8 | >=32 | >=128 |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0000 | 2500 | 1 | 1.00 | 0.0 | 0.0 | 0.0 |
| 0.0072 | 1288 | 23 | 1.94 | 2.2 | 0.0 | 0.0 |
| 0.0120 | 1020 | 30 | 2.45 | 12.2 | 0.0 | 0.0 |
| 0.0168 | 818 | 36 | 3.06 | 26.8 | 1.4 | 0.0 |
| 0.0240 | 646 | 39 | 3.87 | 41.2 | 2.8 | 0.0 |
| 0.0300 | 552 | 43 | 4.53 | 51.4 | 3.0 | 0.0 |

Complexes grow steadily and had not plateaued at 300,000 iterations, but the mean
loop is 4.53 members wide and the largest complex holds 43. Only 3.0% of propagate
work is at or above the 4-thread break-even and none is above the 10-thread one.
Against an 8.19% frame that caps the achievable gain at roughly **0.18% of total
runtime**.

For contrast, `clathrin` at 150,000 iterations ends with its largest complex at 5
members, and the 2000-copy variant from section 12.2 at 4.

### 13.5 The bound that does not depend on complex size

`Complex::propagate` is 8.19% of runtime, so even with perfect scaling and zero
overhead the whole idea is worth at most 8.19% x (1 - 1/4) = **6.1% at four
threads**, 7.4% at ten -- and only if nearly all of that time sat in very wide
loops. The asymptotic case, all 2500 gag molecules coalesced into a single sphere,
is also the case where every other source of parallelism has disappeared, because
there is then one complex to propagate.

### 13.6 What generalizes: the runtime configuration dominates

The most useful result here is not about complex size. Default `passive` waiting
costs 56-80 us per parallel region, **90 times** the 620 ns of four threads with
`active` waiting, and ten threads costs eight times more per region than four
because `OMP_PROC_BIND=spread` places workers on the M5's efficiency cores.

That shape -- per-region cost rising with thread count -- matches the
`openmp-production-lane` measurements on this host exactly: 1.05x at two workers,
0.72x at eight, 0.54x at ten. Some or all of that regression may be a runtime
configuration artifact rather than a limit of the wave scheduling. Before writing
further OpenMP anywhere in NERDSS, that lane is worth re-measuring with
`OMP_WAIT_POLICY=active` and a two-to-four worker cap. Note that lane's own
recorded diagnosis, ~2459 waves of ~41 pairs per timestep, is a wave-occupancy
problem that a cheaper barrier would mitigate but not remove.

## 14. Merging `Coord` and `Vector` into one `Vec3D`

`Coord` and `Vector` were the same three doubles twice over: `Vector` derived
from `Coord` and added a cached `magnitude`. They are now one type,
[`Vec3D`](../../include/classes/class_Vec3D.hpp), with the two halves of the API
reconciled onto one set of names and one set of operators.

Reference: `324c02d`. Candidate: the merge. Both Apple clang 21, `-O3`, serial,
seed 20260810, Apple M5.

### 14.1 What the two types disagreed about

Merging is not a rename, because the two halves behaved differently in ways call
sites had come to depend on:

| | `Coord` | `Vector` | `Vec3D` |
| --- | --- | --- | --- |
| size | 24 B | 32 B | 24 B |
| length | `get_magnitude()`, `magnitude_squared()` | `calc_magnitude()` + `magnitude` field | `length()`, `length_squared()` |
| `cross()` | absent | returns the *normalized* cross product | `cross()` is the cross product, `unit_cross()` normalizes |
| angle | absent | `dot_theta()`, off two cached magnitudes | `angle_between()`, lengths measured or passed |
| projection | absent | `vector_projection()`, returns the *rejection* | `rejection_from()` |
| `operator<<` | `%12.6g` columns | `[xi + yj + zk]` | the columns; the other is `write_ijk()` |
| `a + b` | `Coord` | `Vector` or `Coord`, decided by whether the left operand was const | `Vec3D` |
| serialize | 3 doubles | 4 doubles | 3 doubles |

### 14.2 The magnitude cache, and the four call sites that depended on it

`Vector::magnitude` was written by `calc_magnitude()` and then maintained by
nothing: no operator that changed x, y or z updated it. Of the 78
`calc_magnitude()` calls in the tree, nearly all sat one or two lines above the
read they served, and for those `length()` returns the same bits. Four sites
depended on the cache holding something *other* than the current length, and
each decided a real branch:

- **`create_arbitrary_vector.cpp`** took an angle against a fresh
  `Vector(1,0,0)` whose cache was never written. `dot_theta` returns 0 for any
  operand with a cached magnitude under 1E-8, so the guard was `0 != 0` and the
  x-axis arm of the function has never executed. Kept as the y-axis arm, written
  out.
- **`requiresSignFlip.cpp`** measured `axis` on entry, rotated it, and then took
  another angle against the pre-rotation length. Now an explicit
  `const double axisLength { axis.length() }` at the top.
- **`transform.cpp`** read the caller's cached magnitude, and `check_bases.cpp`
  reaches it through `calculate_phi()` with a vector it builds inline and never
  measures - so `transform()` sees an angle of 0, concludes the axis is already
  aligned, and returns without transforming anything. `transform()` and
  `calculate_phi()` now take the length as an argument and `check_bases.cpp`
  passes `0.0` with a comment saying why.
- **`functions_for_spherical_system.cpp::rotate_on_sphere()`** rebuilt `targjk`
  and then re-tested the cache the rebuild had reset to zero, which made its
  NaN guard always take the first arm and its `exit(1)` arm dead.

Reinstating the "correct" behaviour at any of these would change results, so
none of them was reinstated. They are now visible in the source instead of
implied by an object's history.

### 14.3 Bitwise identity

`run_suite.sh` + `compare_suites.sh`, all 13 cases in `cases.tsv`, every file
under `DATA/`, `PDB/` and `RESTARTS/` (PDB from line 2).

**13 of 13 byte-identical**: `clathrin`, `closed_homoTrimer`, `hetTrimer`,
`hexamer`, `homoTrimer`, `implicit_lipid`, `mem_localization`,
`michaelis_menten`, `rev_2D`, `rev_3D`, `rev_3Dto2D`, `sphere`, `trimer`.

Two things had to be held fixed to get there, both because the host is arm64,
where the compiler will fuse a multiply and an add into one `fmadd` within a
single source expression:

- every expression in `class_Vec3D.hpp` is character for character what the
  corresponding `Coord` or `Vector` member computed, in the same order;
- `rejection_from()` keeps `normal * coefficient` and the subtraction as two
  statements. Written as one expression they contract into an `fmsub` and the
  low bits move.

### 14.4 Whole-simulation timing

`interleaved_timing.sh`, 3 repetitions, median per case.

| case | nItr | before (s) | after (s) | ratio |
| --- | --- | --- | --- | --- |
| `clathrin` | 150,000 | 3.548 | 3.586 | 0.989 |
| `closed_homoTrimer` | 6,000 | 6.561 | 6.508 | 1.008 |
| `hetTrimer` | 8,000 | 6.130 | 6.252 | 0.980 |
| `hexamer` | 25,000 | 5.553 | 5.490 | 1.011 |
| `homoTrimer` | 6,000 | 6.437 | 6.399 | 1.006 |
| `implicit_lipid` | 100,000 | 2.883 | 2.876 | 1.002 |
| `mem_localization` | 2,000 | 2.567 | 2.582 | 0.994 |
| `michaelis_menten` | 150,000 | 4.596 | 4.513 | 1.018 |
| `rev_2D` | 200 | 40.962 | 40.961 | 1.000 |
| `rev_3D` | 20,000 | 4.391 | 4.482 | 0.980 |
| `rev_3Dto2D` | 8,000 | 6.075 | 6.077 | 1.000 |
| `sphere` | 40,000 | 6.252 | 6.404 | 0.976 |
| `trimer` | 8,000 | 6.094 | 6.014 | 1.013 |
| **total** | | **102.049** | **102.144** | **0.999** |

Per-case geometric mean 0.998, range 0.976-1.018. Repetition spread *within* a
single build reaches 8% (`implicit_lipid`) and exceeds 3% in five cases, so the
whole ±2% band here is noise. Excluding `rev_2D`, which is 40% of the total and
lands on exactly 1.000, the ratio is 0.998.

**The change is runtime-neutral.** That is the result, not a hedge: no case
moved outside its own run-to-run spread in either direction.

### 14.5 Where the operations themselves did move

`benchmarks/vec3d_benchmark.cpp` measures the three operations the merge
altered, against the pre-merge `Coord`/`Vector` kept verbatim in the same binary.
It also checks, over 500,000 random vectors, that both give bit-identical
results: **0 mismatches** for `normalize()`'s x/y/z and **0** for the angle.

| operation | legacy | `Vec3D` | ratio |
| --- | --- | --- | --- |
| `normalize()`, length not read after (the common call site) | 0.62 ns | 0.57 ns | 1.096 |
| `normalize()`, length read after | 0.99 ns | 0.94 ns | 1.057 |
| angle between two vectors, both lengths unknown | 4.75 ns | 4.28 ns | 1.109 |
| streaming length sum, 4,096 elements (L1) | | | 1.075 |
| streaming length sum, 131,072 elements (1 MiB) | | | 1.077 |
| streaming length sum, 4,194,304 elements (128 MiB) | | | 1.098 |

Two sources. `Vector::normalize()` took two square roots - one to measure, one
to refresh the cache afterwards - and most call sites feed x, y and z straight
into a quaternion and never ask for the length; that second root is now gone.
And 32 bytes became 24, which shows up when an array of them is streamed.

5-11% on these operations and 0% on the simulation is consistent, not
contradictory: `Complex::propagate` is 8.19% of runtime (section 13.5) and the
association angle code is a fraction of a percent, so a tenth off a tenth of a
percent is not measurable at the suite level.

### 14.6 What is not covered

- The MPI build was not compiled: `mpicxx` is not installed on this host.
  `Vector::serialize()` wrote a fourth double for the cache and now writes
  three, which shrinks every serialized `Complex` by 8 bytes. Both ends of every
  message are built from that one definition, the layout is written and read
  purely sequentially through `PUSH`/`POP`, and there is no `sizeof`-based or
  literal byte count anywhere in `src/mpi` or `include/mpi` that could go stale
  - the buffers are fixed 100 MB and 50 MB slabs from `macro.hpp`. So the
  shrink is safe by construction, but that is reasoned, not measured.
- Two console lines changed. The reaction parser's normal echo keeps its
  `[xi + yj + zk]` format through `write_ijk()`, but
  `create_arbitrary_vector()` no longer prints a spurious "angle between vectors
  with at least one of magnitude 0" warning per call, because it no longer takes
  that angle. Neither is in a hashed output file.

## 15. Coverage of the validation suite, and what it does not reach

`cases.tsv` was built to be fast and to cover the common reaction paths. It is
not a coverage suite, and treating it as one is a mistake this section exists to
prevent: a bitwise suite that never executes a function reports "identical" for
any change to it, correct or not, in exactly the same words it uses for code it
does execute.

Measured by building with

```
make serial CFLAGS="-O1 -fprofile-instr-generate -fcoverage-mapping"
```

running every case, merging with `llvm-profdata merge -sparse` and reading
`llvm-cov report -show-functions` over the three directories that the
consolidation work touched.

### 15.1 What cases.tsv misses

Whole families of functions sat at 0.00% region coverage: every compartment
routine (three reflectors, two transmission probabilities,
`check_compartment_reaction`), the cluster overlap sweep in both geometries, the
unimolecular observable counter, all three span-check reflectors and the three
`_nocheck` reflectors they call, and the entire `excludeVolumeBound` half of
`check_bimolecular_reactions.cpp` - which is, by line count, the largest
remaining duplication target in `reactions`, and therefore the least safe.

`coverage_cases.tsv` closes the part of that gap the existing sample inputs can
reach. `known_uncovered.tsv` records the rest, with what was tried, so the next
person does not have to rediscover it.

### 15.2 A flag that parses is not a path that runs

`clusterOverlapCheck` defaults to false, so nothing in `cases.tsv` reaches
`sweep_separation_complex_rot_memtest_cluster_*`. Turning it on in four models
changed the output of all four, which looks like proof that the code now runs.
It is not. `write_restart.cpp` echoes the flag into every restart file, so three
of the four differed **only** in `restart*.dat`:

| model | differing files | differing outside `restart*.dat` |
| --- | --- | --- |
| `implicit_lipid` | 20 | 0 |
| `sphere` | 20 | 0 |
| `clathrin` | 2 | 0 |
| `mem_localization/SmallBox/FastDsol` | 34 | **14** |

Only the last reaches the cluster code, and `llvm-cov` agrees: with all four
cases merged, `_cluster_box` is covered and `_cluster_sphere` is still 0.00%.
The lesson is that "the output changed" answers a different question than "the
code ran", and only the second one licenses a refactor.

### 15.3 Retroactive verification

`coverage_cases.tsv` run against the pre-consolidation binary
(`nerdss-optimized` at `207e22d`) and against the branch tip, `compare_suites.sh`
between them:

| case | nItr | bitwise identical |
| --- | --- | --- |
| `compartment` | 4,000 | yes |
| `unimol_state_change` | 20,000 | yes |
| `cluster_mem_loc` | 2,000 | yes |

That is new evidence, not a restatement: these three cases are what first
executed `count_unimolecular_observable`, `determine_compartment_probability`,
`collect_cluster_partners` and `resample_partner_trajectories`. Until this table
existed, those four had been refactored and shipped without a single test having
run them.

`cases.tsv` needs no re-run here: the tip binary hashes to
`4492514f0ec3ec63`, the same bytes already compared over all 13 cases and 268
output files.

## 16. The allocator and the libm calls: one large win, one small, one already done

Measured on 2026-08-14, Apple M5, Apple clang 21.0.0, GSL 2.8, `make serial`
with `-O3 -std=c++0x`, seed 20260810. Profiles are `sample(1)` on a `-O3 -g`
build with offsets mapped by `atos`.

This section answers the two items section 13.1 left as the only non-OpenMP
targets it had identified: the allocator at roughly 8% and libm `cos`/`exp` at
roughly 4.7%, the latter marked "SIMD, not OpenMP". One turned out to be much
larger than 8% and trivially removable, one is worth about 2%, and one was
already being done by the compiler.

### 16.1 Where the allocator time actually was

`clathrin`, 150,000 iterations, 5,436 samples. Grouping leaf symbols by image:

| group | share of leaf samples |
| --- | ---: |
| `libsystem_malloc` (`_xzm_free`, `_xzm_xzone_malloc_tiny`, `_free`, `malloc_type_malloc`, `operator new`/`delete`) | 9.10% |
| `libsystem_m` plus the `cos`/`__sincos_stret` stubs | 4.08% |
| `_platform_memset`, `__bzero`, `_platform_memmove` | 2.43% |

Attributing the allocator group to callers put **450 of 5,436 samples, 8.3%, on
a single line**: `EXEs/nerdss.cpp:1284`, which was

```cpp
for (auto &moli : moleculeList) {
  debug_check_nan_Mol(moli, simItr, "Before overlap checking");
}
```

`debug_check_nan_Mol()` took a `const std::string&`. Every caller passed a
string literal, so every call constructed a temporary, and libc++ heap-allocates
a `std::string` whose contents exceed its 22-character short-string buffer.
`"Before overlap checking"` is 23 characters. `"After overlap checking"`, the
otherwise identical call after the overlap sweep, is 22, fits in SSO, and does
not appear anywhere in the allocator profile. One character separated a free
call from 8.3% of runtime.

### 16.2 Removing both debug helpers

`debug_check_nan_Mol` and `debug_print_wrong_Mol` were deleted rather than
repaired, together with the two whole-`moleculeList` loops that drove them, the
four call sites around association, and the commented-out `monitor_iter`
scaffolding. This gives up the guard that called `exit(1)` on NaN coordinates:
NaN now reaches trajectory and restart output instead of aborting the run.

All 13 cases in `cases.tsv` are bitwise identical, which is the expected result
for code that only read coordinates and allocated.

Interleaved paired timing, `interleaved_timing.sh 5`, medians:

| case | nItr | before (s) | after (s) | speedup |
| --- | ---: | ---: | ---: | ---: |
| `rev_3D` | 20,000 | 4.610 | 3.898 | **1.183x** |
| `implicit_lipid` | 100,000 | 2.781 | 2.480 | **1.121x** |
| `clathrin` | 150,000 | 3.506 | 3.253 | 1.078x |
| `hexamer` | 25,000 | 5.278 | 4.942 | 1.068x |
| `michaelis_menten` | 150,000 | 4.366 | 4.127 | 1.058x |
| `mem_localization` | 2,000 | 2.789 | 2.667 | 1.046x |
| **total** | | **23.330** | **21.367** | **1.092x** |

The profile after the change confirms the mechanism rather than merely the
result: `malloc`/`free` fell from 9.10% to **0.94%** of leaf samples and
memset/bzero from 2.43% to 1.13%. An 8.2-point drop in allocator share against a
measured 1.092x is as close to accounting for itself as this kind of measurement
gets.

The case table is `debugremoval_cases.tsv`. `rev_2D` is excluded because roughly
40 s of its runtime is fixed 2D lookup-table construction independent of
iteration count, which dilutes a per-timestep effect and costs 80 s per
repetition pair to measure.

### 16.3 A build bug found while measuring this

The first build after editing `EXEs/nerdss.cpp` reported success and produced a
**byte-identical binary with the change absent**. `bin/nerdss` listed only
`$(OBJS)` as prerequisites; the executable source appeared on the recipe line but
never as a dependency, so make compared the binary against the objects alone and
concluded there was nothing to do.

This is the same defect the header-dependency block at the bottom of the
`Makefile` already documented for objects, left in place for the one translation
unit that holds the entire timestep loop. It was caught only because the binary
hashes were being compared. **Any earlier measurement in this repository that
touched `EXEs/nerdss.cpp` and nothing else may have timed a stale binary.**

Fixed by making `$(EDIR)/$(_EXEC).cpp` a prerequisite and running `DEPFLAGS` on
the link step so header edits relink too. Verified: touching the executable
source now triggers 1 invocation where it triggered 0, touching a header it
includes triggers 52, and an unchanged tree triggers 0.

### 16.4 The libm calls: `add_3D_rotational_diffusion`

After the allocator work, the whole remaining libm cost sat in
`add_3D_rotational_diffusion()`, which evaluated
`cos(sqrt(k * Dr.z * timeStep))` twice per candidate pair, once per partner.
Disassembling `determine_3D_bimolecular_reaction_probability` confirmed exactly
two calls to the `_cos` stub and none to `__sincos_stret`.

That value depends on nothing but `Dr.z` and the timestep, so it is constant for
as long as the complex's composition is. Caching it on the complex converts an
O(candidate pairs) count of `cos()` calls into O(complexes whose `Dr` changed),
which is why the gain tracks pair density rather than molecule count.

The cache is two `mutable` slots on `Complex`, one for `k = 2` and one for
`k = 4`, keyed on the exact argument compared bit-for-bit. Keying on the argument
rather than invalidating from `update_properties()` means it cannot go stale: any
change to `Dr.z` or the timestep is a miss. The argument is built from the same
operands in the same order as the expression it replaced, and `-O3` without
`-ffast-math` may not reassociate it, so a miss recomputes the identical double
and a hit returns the one computed from it. The fields are not serialized,
following `deleteIfNotReceivedBack`, so restart files are unchanged.

All 13 cases bitwise identical.

`interleaved_timing.sh 9`, **minimum** per case rather than median:

| case | before (s) | after (s) | speedup |
| --- | ---: | ---: | ---: |
| `rev_3D` | 4.574 | 4.403 | 1.039x |
| `hexamer` | 5.227 | 5.128 | 1.019x |
| `michaelis_menten` | 4.313 | 4.237 | 1.018x |
| `mem_localization` | 2.612 | 2.565 | 1.018x |
| `implicit_lipid` | 2.527 | 2.530 | 0.999x |
| **total** | **22.795** | **22.152** | **1.029x** |

Minimums, not medians, and the reason is worth recording. This host was at load
average 4 to 11 throughout. On a first 5-repetition pass `clathrin`'s five
`before` timings were 3.678, 8.770, 6.666, 9.153 and 4.332 s -- a 2.5x spread on
identical work -- which put its median ratio at 0.808x and dragged the aggregate
to 0.945x, i.e. reported a regression. A second 9-repetition pass put the same
case at 1.577x by median. Neither number means anything. The cases whose
per-timestep work is largest, `rev_3D` and `hexamer`, held to within 2% across
repetitions in both passes and are the only ones worth reading. `clathrin` is
omitted from the table above at 1.077x by minimum for the same reason.

Taking the five listed cases at face value the change is worth **about 1.5-2%**,
the same class as the size-robust residual in section 12.3 rather than the
dilution-dependent gains the suite over-represents. A quiet host would be needed
to tighten it further, and that is the honest limit of this measurement.

### 16.5 `sincos` fusion: already done by the compiler

The remaining idea from section 13.1 was to halve the libm call count where
`sin(x)` and `cos(x)` are both needed for the same `x`:
`create_euler_rotation_matrix()` makes six such calls for three angles, and
`Complex::propagate()` makes six for the quaternion's three half-angles.

Apple's `__sincos` was first checked for bit-exactness against separate `sin` and
`cos`, since a fused call that returns different doubles could not be a
result-preserving change. Over 3,200,000 arguments across eight scales from 1e-8
to 1e3: **zero mismatches in both sin and cos**. So the fusion would have been
safe.

It would also have been pointless. Disassembling the two functions in the
committed binary:

| function | `__sincos_stret` calls | separate `sin`/`cos` calls |
| --- | ---: | ---: |
| `create_euler_rotation_matrix(double,double,double)` | 3 | 0 |
| `Complex::propagate` | 3 | 0 |

Three calls for three angles, and no scalar `sin` or `cos` at all. Clang at
`-O3` already recognises adjacent `sin`/`cos` of a common argument and lowers the
pair to one `__sincos_stret`. Writing the call explicitly would have emitted the
same three instructions.

This is the third time on this branch that a profile frame has been read as an
opportunity when the cost was somewhere else -- sections 9.6 and 11 are the
other two -- and the first where the optimization had already been applied by
something other than this repository. The general lesson is narrow but real: on
a frame attributed to a libm symbol, check what the compiler already emitted
before assuming the source spells out what runs. `DYLD-STUB$$cos` at 99 samples
and `cos` at 42 in the section 16.1 profile were not the euler matrix or the
quaternion at all; every one of them came from
`add_3D_rotational_diffusion()`'s cosine, which has no matching sine and which
16.4 then removed.

### 16.6 Combined result

The committed tip was rebuilt from a clean tree and compared over all 13 cases
and every file under `DATA/`, `RESTARTS/` and `PDB/` against a build of
`e5be2a6`, the branch state this work started from: **bitwise identical**, so the
committed state reproduces the measured state in the sense section 9.7 uses.

That span is wider than the two changes described here. Four other commits landed
between `e5be2a6` and this section's first one -- `321d893`, `a54c4dd`,
`ddff3a8`, `1821bfe` -- so the identity result covers those as well, and is
evidence for them too. It is not evidence about anything before `e5be2a6`.

On timing, the two changes were each measured against their own immediate
predecessor -- 1.092x in 16.2 and about 1.02x in 16.4 -- so the combined figure
of roughly **1.11x** on the six-case table is the product of two separately
measured ratios, **not** a directly measured one. A direct pre-versus-post run
was attempted and abandoned: this host would not stay quiet long enough, and two
attempts collided with each other's result directory before the contention was
noticed. Given that a single disturbed repetition already moved `clathrin` from
0.808x to 1.577x in 16.4, a combined number produced under those conditions
would be worth less than the product of two clean ones. Anyone wanting the
direct measurement should take `debugremoval_cases.tsv` to an idle machine and
compare `1821bfe` against this tip.

What is not estimated is the correctness claim. The committed tip removes about 8
points of allocator share and the last two scalar `cos()` calls in the hot path,
and alters no output byte in any of the 13 cases.

## 17. `functions_for_spherical_system.cpp`: the arithmetic, the copies, and the `sqrt` that had to stay

Seven issues were raised against this file: an unnecessary `sqrt`, `pow(x, 2)`,
`sin` and `cos` computed separately for the same angle, functions kept past
their last caller, `Vec3D` doubling as a spherical coordinate with nothing
marking which convention a value held, no statement anywhere that a "spherical
system" here is a system *on* a sphere rather than *inside* one, and variable
names that do not follow the project's camelCase.

Two of the seven were already handled by the compiler. One of them was wrong in
three of the eight places it applied. And the largest measured win in the file
was not on the list at all.

### 17.1 The trig and the `pow`: what `-O3` already emitted

Section 16.5 found that clang fuses adjacent `sin`/`cos` of a common argument
into `__sincos_stret` on its own. The same check on this file, reading
relocations out of the object built from `HEAD`:

| function (baseline) | `__sincos_stret` | `_sin` | `_cos` | `_acos` | `_asin` | `_pow` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `find_cardesian_coords` | 2 | 0 | 0 | 0 | 0 | 0 |
| `rotate_on_sphere` | 1 | 0 | 0 | 1 | 0 | 0 |
| `find_spherical_coords` | 0 | 1 | 0 | 2 | 0 | 0 |
| `find_position_after_association` | 0 | 0 | 1 | 0 | 0 | 0 |
| `calc_bindRadius2D` | 0 | 0 | 0 | 0 | 1 | 0 |

Three things that source reading would have flagged were already gone. The four
`pow(x, 2.0)` calls emitted **no `_pow` at all** -- the constant exponent was
folded to a multiply. `find_cardesian_coords` wrote `sin(theta)` twice, `cos(phi)`
and `sin(phi)` once each, and got two `__sincos_stret` calls: both the fusion and
the common-subexpression elimination of the repeated `sin(theta)` had happened.
`find_position_after_association` wrote `cos(arc1 / R)` twice and got one call.

The candidate object has exactly the same counts. So none of issues 2 and 3 is a
speedup on this compiler at this optimization level, and the file was changed
anyway for the reason the issue gave: the folding does not survive
`-fno-builtin`, it never happened on older MSVC, and it does not apply when the
exponent is not a literal. `pow` is written out as multiplication;
`include/math/sincos.hpp` wraps `__builtin_sincos`, which both compilers lower to
whichever sincos entry point the target has, and which was checked to return the
same doubles as separate `sin` and `cos`.

### 17.2 `sqrt`: valid against a threshold, invalid between two lengths

Eight comparisons in the file took a square root. Five are guards against a
*fixed tolerance*, all of them `1e-8` -- a step that did not move, a target on
top of its reference point, a tangential component that vanishes. Those became
comparisons against a `...Squared` constant formed from the tolerance by
multiplication. `sqrt` is correctly rounded and monotonic, so the two disagree
only for an input within about one part in 10^16 of the tolerance exactly, which
is not a value any step in this simulation lands on.

The other three compare **two measured lengths**, and for those the same
reasoning fails. `set_memProtein_sphere` and `find_Lipid_sphere` search for the
molecule whose temporary COM is furthest from the centre of the sphere -- over
molecules that are all *on* that sphere. The lengths being ranked are therefore
equal up to rounding by construction, and `sqrt` collapsing two of them onto one
double is the normal case rather than a remote one. Measured over 2,000,000 pairs
of points placed on a sphere of radius 100 by the construction the callers use:

| | |
| --- | ---: |
| `\|a\|` and `\|b\|` bit-identical after `sqrt` | 34.42% |
| `(\|a\| > \|b\|)` disagrees with `(\|a\|^2 > \|b\|^2)` | 1.98% |

A disagreement there does not shift a coordinate by an ulp. It selects a
different molecule as the membrane reference protein, which then sets the
orientation the whole complex is rotated into. A differential test against the
previous implementation over 80,000 fixtures reported **2,628 flips in each of
the two searches** before the change was reverted, and **0 after**.

None of the 18 cases in `cases.tsv` and `coverage_cases.tsv` shows this: in every
call they make, the search has exactly one candidate to rank, so no comparison
between two lengths ever happens. A bitwise-identical verdict from the suite
would have been reported in the same words as a verified one -- the situation
section 15 was written about, arrived at from the opposite direction.

The third, the nearer of two candidate positions in
`find_position_after_association`, keeps its roots too. Its two candidates are
the roots of a quadratic and are not near-equal by construction, and 600,000
randomized comparisons found no disagreement, but the saving is two `sqrt` per
association event against roughly 900 events in the `sphere` case -- far too
small to be worth a comparison whose failure mode is a different position rather
than a rounder one.

### 17.3 The copies nobody asked about

```
void set_memProtein_sphere(Complex reactCom, Molecule& memProtein,
    std::vector<Molecule> moleculeList, const Membrane membraneObject)
```

`Complex` by value, `std::vector<Molecule>` **by value**, `Membrane` by value.
Neither this function nor `find_Lipid_sphere` writes to any of the three. Every
call deep-copied the entire molecule list to read a handful of `isLipid` flags and
coordinates off it. `Molecule` is 656 bytes and holds **22 non-static
`std::vector` members**, so one copy is 656 bytes plus an allocation for each
member that is not empty.

`associate_sphere`, `associate_ImplicitLipid_sphere` and
`perform_implicitlipid_state_change_sphere` each call both functions once per
event, so a sphere association copied the list twice. The `sphere` case holds
1,101 molecules and records 914 successful associations in 40,000 iterations:
about **2.0 million `Molecule` copy constructions**, each with up to 22
allocations behind it, to pass data nothing writes to.

The three call sites shrank accordingly, and this is where the reading of the
change stopped being a measurement and started being a guess -- 17.6 tests it and
finds the guess wrong. `.text` instructions per object:

| caller | baseline | candidate | |
| --- | ---: | ---: | ---: |
| `associate_sphere` | 6,032 | 5,046 | -986 |
| `associate_ImplicitLipid_sphere` | 4,344 | 3,654 | -690 |
| `perform_implicitlipid_state_change_sphere` | 3,885 | 3,195 | -690 |
| `perform_bimolecular_state_change_sphere` | 1,076 | 1,076 | 0 |
| total | 15,337 | 12,971 | **-2,366** |

The last row is the control: that file calls `calc_bindRadius2D` and
`find_position_after_association` but neither of the two, and does not move.

Two smaller structural changes are in the same class. The nine-double inner
coordinate frame was passed by value -- 72 bytes copied per call, once per
molecule and once per interface per timestep -- and is now a reference. And
`calculate_inner_coord_coefficients`, declared in the header but called only by
`translate_on_sphere` in the same file, became file-local and inlined into it:

| | baseline | candidate |
| --- | ---: | ---: |
| `translate_on_sphere` | 57 | 177 |
| `calculate_inner_coord_coefficients` | 147 | inlined |
| the pair | 204 | **177** |

### 17.4 `Vec3D` stays by value, and that direction was measured

Changing the `std::array` parameters to references invites doing the same to the
`Vec3D` ones. That was tried first, and it is slower. `Vec3D` is three doubles in
a trivially copyable aggregate, so it is a homogeneous floating-point aggregate
and travels in three FP registers; a reference forces the caller to spill it to
the stack for an address to point at and makes every read in the callee a load.

| function | `const Vec3D&` | `Vec3D` |
| --- | ---: | ---: |
| `inner_coord_set` | 146 | 142 |
| `inner_coord_set_new` | 169 | 166 |
| `translate_on_sphere` | 182 | 177 |
| `rotate_on_sphere` | 124 | 118 |
| `find_position_after_association` | 132 | 126 |

The reason the two conventions differ is now written on the header, next to the
declarations that differ, so the inconsistency reads as deliberate.

### 17.5 Bitwise identity

Three independent checks, because the suite alone was shown above to be blind to
the one real behaviour change in this work.

**A differential test on the geometry.** `HEAD`'s implementation is textually
embedded in a `baseline` namespace beside the rewritten one and both are called
on the same inputs, comparing raw bit patterns: 40,000 random points per sphere
plus both exact poles, the three axes, radii from `1e-9` to `1e12`; steps that
are real, exactly zero, `1e-12` and `1e-9`; targets at the COM, purely radial,
and along each basis vector; rotation angles including 0, `1e-15` and both signs
of pi. **11,322,143 comparisons, 0 mismatches.** Perturbing the baseline copy by
one part in 10^15 makes it report 1,908,832 -- the harness detects what it claims
to detect.

**A differential test on fixtures.** `set_memProtein_sphere` and
`find_Lipid_sphere` need `Molecule`/`Complex` inputs the first test cannot build,
and they are the two whose signatures changed most. 80,000 randomized complexes
covering both the explicit-lipid path and the implicit-lipid path, comparing
every coordinate of the returned molecule: **0 mismatches.** This is the test that
caught 17.2.

**The suites.** All 13 cases in `cases.tsv` and all 5 in `coverage_cases.tsv`,
every file under `DATA/`, `RESTARTS/` and `PDB/`, against a build of the parent
commit: **18/18 bitwise identical.**

**And the committed tip is the measured tip.** Section 16.6 had to argue that its
committed state reproduced its measured state by re-running the suite and
comparing output bytes. Here the check is tighter: `git archive` of the tip
extracted to an empty directory and built with `make serial` produces a binary
whose SHA-256 is `554ca4b9a82690a96a5f1626...`, which is the same binary every
number in 17.5 and 17.6 was taken on -- byte for byte, not output for output. The
verdicts above therefore apply to the commit directly, with nothing to re-run.

### 17.6 Timing, and where the speedup is not

Measured with a purpose-built interleaved harness rather than `run_suite.sh`, for
the reason 17.7 gives. Seven repetitions, three builds per repetition per case so
all three meet the same machine conditions, medians and minimums both reported.
`ablate` is `HEAD` with **only** the by-value to const-reference change of 17.3
applied and nothing else, so its column isolates that change from the rest.

| case | base | ablate | cand | ablate/base | cand/base |
| --- | ---: | ---: | ---: | ---: | ---: |
| `sphere` (median) | 11.253 | 11.137 | 10.754 | 1.010 | **1.046** |
| `sphere` (minimum) | 11.155 | 11.043 | 10.482 | 1.010 | **1.064** |
| `gagsphere` (median) | 3.798 | 3.769 | 3.798 | 1.008 | 1.000 |

The `sphere` distributions do not overlap: across seven paired repetitions the
slowest candidate run was faster than the fastest baseline run.

**The copies are not where the time was.** 17.3 counted 2.0 million eliminated
`Molecule` copy constructions and 2,366 eliminated instructions at the call
sites, and reading those numbers it is natural to call it the largest win in the
file. The ablation says it is worth **1.010x** -- about a fifth of the 1.046x, and
consistent with 1.32 GB of `memcpy` at roughly 0.13s of an 11.2s run. The other
3.6% comes from the propagation path.

The reason is a ratio nothing in the source makes visible. `translate_on_sphere`
and `rotate_on_sphere` run once per molecule and once per interface per timestep;
`set_memProtein_sphere` and `find_Lipid_sphere` run twice per association event.
In the `sphere` case that is 44,040,000 molecule propagations against 1,828 list
copies -- **24,000 to 1**. An expensive thing done rarely lost to a cheap thing
done constantly, which is the same shape of error as sections 9.6, 11 and 16.5,
and the fourth time on this branch that a plausible cost has been mis-attributed
without a measurement to separate it.

So the credit for the 3.6% belongs to the changes on the per-molecule path: the
frame passed by reference instead of as 72 bytes of copy, the coefficient solve
inlined into its only caller, the threshold `sqrt` calls removed from the guards
every call runs through, `radius()` replaced by an inlinable `Vec3D::length()`
where it had been an out-of-line call in another translation unit, and eighteen
dead stores removed from the two frame builders. Those were not separated from
each other; the ablation was spent on the claim that looked largest.

`gagsphere` is a negative control, not a second sphere case. Its parameter file
`gagOriginalModelSolution.inp` declares `waterBox = [436, 436, 436]` and no
`isSphere`, so no code in this file executes, and 1.000x is the correct answer for
it. A useful control, but it means **`sphere` is the only case here that measures
anything**, and the 1.046x rests on one model.

A second on-sphere case was attempted and abandoned. `gagsphere/parms.inp` does
declare `isSphere = true, sphereR = 70`, and `cluster_gagsphere` runs it with
`clusterOverlapCheck = true`, which is the configuration that drives
`calculate_update_position_interface` and therefore this file hardest. But its
cost appears only once clusters have grown: at `nItr = 160` the whole run is
0.162s and has formed 11 bound pairs, so every build reads the same startup time
and nothing is measured, while `nItr = 1000` forms 64 pairs and takes 62s per
run -- 15 minutes for five interleaved repetitions of three builds. Anyone wanting
a second sphere data point should spend that on an idle machine; the harness
takes an injected parameter line, so
`gagsphere<TAB>gagsphere<TAB>parms.inp<TAB>1000<TAB>clusterOverlapCheck = true`
is the row to use.

### 17.7 A measurement bug found on the way

`run_suite.sh` cannot time this change, or any change of a few per cent. Its
`run_with_timeout` ends with `wait "$watchdog"`, and the watchdog polls with
`sleep 1`, so it can only exit after its in-flight sleep finishes -- and the
caller reads the clock after that. Every timing it reports is rounded up to the
next whole second. In one 13-case run every value was an exact integer multiple
of ~1.0095s: `rev_3D` 7.069, `rev_3Dto2D` 10.095, `homoTrimer` 11.105, `hexamer`
10.092, `clathrin` 6.060, `sphere` 11.104, `rev_2D` 77.649. The comment inside
`run_with_timeout` asserts the opposite, which is true of the `wait "$child"`
above it and not of the wait that follows.

The hash manifest is unaffected, so every bitwise verdict the suite has ever
produced still stands. Only the `seconds` column is compromised.
`interleaved_timing.sh` has no watchdog and is unaffected.

### 17.8 Found and deliberately not changed

**A missing pole guard in the caller.** Fixed in section 18, which also retracts
the reason given here for leaving it: the change does alter output at the poles,
but nothing reaches a pole, so it is bitwise-identical on both suites after all.

**A length compared against unity.** `rotate_on_sphere`'s second early-out is
`std::abs(targi.length() - 1.0) < 1e-8`, where `targi` is a projection in nm. It
reads like a unit-vector test applied to the wrong quantity, and it fires only
for a target exactly one nanometre out along the rotation axis. Preserved
verbatim; the question needs the model's author, not a refactor.

**A south pole at -pi.** `find_spherical_coords` returns `theta = -M_PI` where
`acos` would give `+M_PI` and where the documented range is `[0, pi]`. Both map
back to the same cartesian point -- `find_cartesian_coords` differs only in the
sign of a `sin(theta)` of magnitude 1.2e-16 -- and the one caller that reads
`theta` back treats them identically. Documented rather than corrected.

**Two redundant normalizations.** `inner_coord_set` ends by normalizing all three
basis vectors, but `j` and `k` come from `unit_cross`, which already normalized
them. Dropping those two would save two `sqrt` and six divides per call, and
would perturb the last bits wherever the re-measured length is not exactly 1.0.
Not taken in a commit whose claim is bitwise identity.

## 18. The pole guard the caller never had

17.8 listed a missing guard in `create_complex_propagation_vectors_on_sphere` and
left it alone, on the grounds that fixing it changes output at the poles and so
could not belong to a bitwise-identical commit. The first half of that is true.
The second does not follow, and this section retracts it: nothing in any
configuration NERDSS runs ever reaches a pole, so the guard is bitwise-neutral
everywhere the code actually goes. It is now in.

### 18.1 What the degenerate frame does

The caller seeds its tangent from `(0, 0, 1)`. On the z axis that is parallel to
the radial direction `i`, so `temp.unit_cross(i)` is the zero vector, and
`Vec3D::normalize()` deliberately leaves a zero vector alone rather than making
NaNs. `j` and `k` come back as `(0, 0, 0)`.

Driven through the real `rotate_on_sphere` at `COM = (0, 0, 100)` with a step of
`dl = 1.414`, the result is `COMnew = (0, 0, 99.99)`: **arc moved 0 against 1.414
requested, radius 100 -> 99.99.** The COM does not collapse -- `targi` survives
because `i` is still well defined -- it loses the tangential half of its step and
keeps the small inward radial half. Iterating shows the shape of it:

```
  step 1: comCoord=(0, 0, 99.9900001667)  |r|=99.9900001667
  step 2: comCoord=(0, 0, 99.9800003333)  |r|=99.9800003333
  ...
```

`x` and `y` stay exactly `0.0`, so the axis is an absorbing state: a complex that
arrives never diffuses again, and sinks `dl^2/2R` per step. Nothing corrects it.
`reflect_traj_complex_rad_rot_nocheck_sphere` returns early for an `OnSurface`
complex because "the movement only involves theta and phi, and R doesn't change" --
precisely the invariant the degenerate frame breaks. With the guard, the same
five steps leave the axis on the first one, move the full `1.414` every time, and
hold `|r| = 100` exactly.

### 18.2 Whether anything reaches it

Four routes to a COM with `x` and `y` both exactly `0.0`, none of which arrives.

**Random placement.** `Molecule::create_random_coords` writes
`comCoord.x = (2R) * rand_gsl() - R`. `gsl_rng_default` is mt19937, whose
`gsl_rng_uniform` returns `k / 2^32`, so `x` is exactly zero only for a draw of
exactly 0.5: `2^-32` per axis, **`2^-64` for both.**

**Propagation landing there.** An instrumented build counting exact hits,
guard-window entries and closest horizontal approach, on the three sphere cases:

| case | calls | exact poles | inside the guard window | closest rho |
|---|---|---|---|---|
| `sphere`, 40000 itr | 16,426,980 | 0 | 0 | 0.061 nm on R = 100 |
| `cluster_gagsphere` | 37,027 | 0 | 0 | 5.54 |
| `gagsphere` | 0 | -- | -- | never enters the surface path |

The guard window is `|z| - |COM| < 1e-8`, which on R = 100 is a cap of radius
~0.0014 nm; the closest any complex came was 43x wider than that, over 16.4M
steps.

**Restart files.** `write_restart` writes coordinates `std::fixed` under
`precision(20)`. Nothing rounds to zero, so a restart can only carry an exact
zero that already existed.

**`-c/--coordinate`.** The one route that *can* write an exact zero:
`std::stod(line.substr(30, 8))`, applied as an exact translation. It still does
not reach the bug. A molecule placed on the axis is not `OnSurface` until it
binds, and its first free 3D diffusion step gives it nonzero `x` and `y` first --
run with an `A` seeded at `(0, 0, 99.5)`, 0 exact poles. Being `OnSurface` at
t = 0 requires an explicit `isLipid` molecule, and NERDSS does not accept one on
a sphere: `class_SimulVolume.cpp` demands a lipid sit at `z = -waterBox.z/2` and
exits with "is off the membrane" otherwise. Placed exactly at the south pole it
satisfies that test and then segfaults in the main loop, before the first
propagation.

So the defect is real and latent. What the guard buys is not a corrected
trajectory today but the removal of a silent, unrecoverable failure from a path
that a future model -- an explicit lipid on a sphere, a hand-built initial
condition -- could walk into with no diagnostic at all.

**Still open, in the same family.** `Complex::propagate` calls
`inner_coord_set(COM, COMnew)`. At a pole that takes the *motion* branch, where
`com` and `comNew` are parallel and `i.unit_cross(v)` is zero as well; the guard
exists only in the no-motion branch. With the caller fixed, `COMnew` is off-axis
and the path is safe, but the hole is still there for any caller that passes two
parallel arguments.

### 18.3 Why the guard, and not a call to `inner_coord_set`

The obvious repair is to delete the hand-rolled block and call
`inner_coord_set(COM, COM)`, which takes the no-motion branch and carries the
guard. The block is that branch character for character. It is not that branch
instruction for instruction: it normalizes `j` a second time *before* building
`k`, while `inner_coord_set` builds `k` from the singly-normalized `j` and
normalizes `i` again at the tail. Two redundant normalizations in different
places, which do not cancel -- the same trailing-normalization arithmetic 17.8's
last entry describes.

Over 20,000,000 uniform points on an R = 100 sphere, comparing raw bit patterns:

| frame vs the current block | differing | max abs delta |
|---|---|---|
| `inner_coord_set(COM, COM)` | 8,819,992 (**44.10%**) | 3.33e-16 (1-2 ulp) |
| the same block, guard added | **0** | 0 |

Delegating would perturb the propagation frame at 44% of ordinary surface
positions and diverge every sphere trajectory, in exchange for correct behaviour
at a state nothing reaches. The guard in place is a no-op everywhere the code
goes. The duplication therefore stays, with a comment naming `onAxisTolerance` as
the constant it is tracking and recording why the two cannot simply be merged.

### 18.4 Bitwise identity

All 13 cases in `cases.tsv` and all 5 in `coverage_cases.tsv`, every file under
`DATA/`, `RESTARTS/` and `PDB/`, against a build of the parent commit:
**18/18 bitwise identical**, `sphere`, `gagsphere` and `cluster_gagsphere`
included, 23 output files hashed apiece.

The baseline binary is SHA-256 `554ca4b9a82690a96a5f1626...`, which 17.5 already
identified as the build of the committed tip, so the comparison is against the
recorded parent exactly rather than against a rebuild that merely ought to match
it. The 20,000,000-point frame comparison in 18.3 is the reason to expect that
verdict rather than merely to observe it: the guard cannot change a frame at any
position off the axis, and no run puts a complex on the axis.

## 19. The pairwise cell grid: seven corrections to how it is sized and walked

Measured on 2026-08-24, same host and build settings as section 1: Apple M5,
10 physical cores, Apple clang 21, `-O3 -std=c++0x`, GSL 2.8, `make serial`,
seed 20260810. Baseline is the branch tip before this work, SHA-256
`47752a9dbade23ac`; the tip after it is `6406a0a97a92c454`.

The serial bimolecular search bins molecules into a uniform Cartesian grid over
`Membrane::waterBox` and walks a 13-neighbour half-stencil. An audit of how the
grid is sized found one silent correctness bug and four ways the sizing
hyperparameters cost work; each is one commit below, each verified against its
own immediate parent.

### 19.1 What the grid was doing, per geometry

Over-inclusion is candidate pairs the stencil offers per pair actually within
`rMaxLimit`. The floor for *any* cell list with cells exactly `rMaxLimit` wide
is 27 / (4π/3) ≈ 6.4×; cases below it are simply clustered.

| case | box (nm) | rMaxLimit | grid | edge ÷ rMax | occupied | offered/step | in range | over-incl. |
|---|---|---|---|---|---|---|---|---|
| clathrin | 494³ | 33.7 | 14³ | 1.05 | 3.5% | 49 | 10.0 | 4.96× |
| michaelis_menten | 144³ | 33.9 | 4³ | 1.06 | 84% | 1 647 | 278 | 5.93× |
| rev_3D | 940³ | 16.7 | 30³ | 1.88 | 6.8% | 1 972 | 167 | 11.8× |
| rev_3Dto2D | 1000³ | 16.9 | 30³ | 1.97 | 6.1% | 45 741 | 4 148 | 11.0× |
| mem_localization | 470·470·752 | 24.4 | 19·19·30 | 1.01 | 5.1% | 165 084 | 57 433 | 2.87× |
| rev_2D | 1000·1000·10 | 5.79 | 30·30·4 | 5.75 / 0.43 | 21% | 12 198 | 233 | 52.3× |
| implicit_lipid | 200³ | 12.8 | 15³ | 1.04 | 4.8% | 397 | 122 | 3.25× |
| sphere | R=100 | 9.77 | 20³ | 1.02 | 5.6% | 760 | 185 | 4.11× |
| compartment | 1000³ | 28.4 | 17·17·8 | 2.07 / 4.41 | 4.0% | 177 | 1.8 | 97.6× |

### 19.2 The correctness bug: sub-volumes thinner than the interaction range

`SimulVolume::Dimensions` gave z a `std::max(4, ...)` floor, so any box got at
least four sub-volumes along z. Below `waterBox.z = 4 * rMaxLimit` that made
them thinner than the interaction range and the ±1 stencil stopped spanning it:
molecules two sub-volumes apart in z could be within `rMaxLimit` and were never
offered to `check_bimolecular_reactions()`. `check_dimensions()` opened with a
repair for exactly this, which re-applied the same floor and so could not fix
it.

Brute-force count of every pair within `rMaxLimit`, rev_3D model, 2000
molecules, `rMaxLimit` 16.70 nm, 19 snapshots, box thickness varied:

| box z | cells z | edge z | pairs in range | reached by stencil | missed |
|---|---|---|---|---|---|
| 30 nm | 4 | 7.50 | 29 123 | 26 951 | 2 172 (7.5%) |
| 50 nm | 4 | 12.50 | 20 475 | 20 391 | 84 (0.41%) |
| 100 nm | 5 | 20.00 | 11 350 | 11 350 | 0 |
| 940 nm | 30 | 31.33 | 1 328 | 1 328 | 0 |

Zero missed whenever the edge is at or above `rMaxLimit`, non-zero the moment
the floor forces it below. The code is byte-identical on `master`, so this is
pre-existing rather than something the branch introduced.

It is not confined to thin boxes. `trimer` is a **cubic** 118.41 nm box against
a 37.58 nm range: `floor(118.41 / 37.58) = 3`, so the floor forced four z
sub-volumes 29.60 nm thick. Reconstructing both grids over 101 trajectory
frames, 300 molecules, 454 619 pairs within `rMaxLimit`:

| grid | unreachable by the stencil |
|---|---|
| old 3 × 3 × 4 (39.47, 39.47, 29.60 nm) | 1 808 — 0.398% ± 0.009% |
| new 3 × 3 × 3 (39.47 nm) | 0 |

0.4% of reacting pairs restored. That is below what an equilibrium comparison
resolves: over 28 seeds at nItr 30000, `trimer`'s total bound pairs move
46.194 ± 1.250 → 48.025 ± 1.278, z = +1.02. The direction is right — all three
bound species up, all six free species down — but the seed-to-seed scatter is
2.7% against a 0.4% mechanism, so the +3.96% figure is scatter, not signal. The
pair count above is the measurement that carries the claim.

### 19.3 The commits, and what each one changed

| # | commit | cases not bitwise identical | why |
|---|---|---|---|
| 7 | Register created molecules in `occupiedSubCells` | 0 / 18 | |
| 1 | Size sub-volumes from the interaction range in every dimension | `trimer`, `rev_2D` | both were subject to 19.2 |
| 4 | Visit only the occupied sub-volumes | 0 / 18 | |
| 2 | Drop the max(4000, N²/2) budget | `compartment` | the only grid the budget touched |
| 3 | Bound sub-volumes by a memory budget, not 30 per dimension | `rev_3D`, `rev_3Dto2D`, `rev_2D`, `compartment` | the only grids the cap touched |
| 5 | Reject out-of-range pairs before the interface search | 0 / 18 | |
| 6 | Skip candidate pairs whose molecule types cannot react | 0 / 18 | |
| — | Track occupied sub-volumes with a bit per sub-volume | 0 / 18 | reimplements 4 |
| — | Stop the MPI sub-volume cap from raising a dimension | 0 / 18 | `#ifdef mpi_`; serial binary byte-identical |

`rev_2D` under commit 1 does not change trajectory at all: every file under
`DATA/` and `PDB/` is byte-identical and the restart files differ only in each
molecule's recorded `mySubVolIndex`, all 1600 of them by the same 2700 = 3 ×
900, the z-row offset. Every molecule sits on the membrane, so both grids put
all of them in one z-row behind the same in-plane stencil.

Each grid change lands on the finest grid that still covers its range:

| case | before | after | edge | rMaxLimit |
|---|---|---|---|---|
| rev_3D | 30³ | 56³ | 16.79 | 16.70 |
| rev_3Dto2D | 30³ | 59³ | 16.95 | 16.94 |
| rev_2D | 30·30·4 | 172·172·1 | 5.81 | 5.79 |
| compartment | 17·17·8 | 35³ | 28.57 | 28.35 |
| trimer | 3·3·4 | 3³ | 39.47 | 37.58 |

Statistical agreement for the stream-changing commits, plateau-averaged copy
numbers via `compare_observables.py`:

| commit | case | seeds | nItr | largest \|z\| |
|---|---|---|---|---|
| 1 | trimer | 28 | 30 000 | 1.11 |
| 2 | compartment | 16 | 4 000 | 1.00 |
| 3 | rev_3D | 12 | 60 000 | 0.58 |
| 3 | rev_3Dto2D | 12 | 8 000 | 1.24 |
| 3 | rev_2D | 12 | 200 | 0.58 |

### 19.4 `rMaxLimit` is not a bound for 2D reactions

Commit 5 rests on `rMaxLimit` bounding how far apart two molecules can be and
still react, which is the same bound the cell list assumes when it sizes
sub-volumes by it. Instrumenting every reacting pair shows the bound holds
except in one place:

| model | reacting pairs | beyond rMaxLimit | worst |
|---|---|---|---|
| rev_2D | 20 220 | 1 550 (7.7%) | 1.109×, a 2D pair |
| rev_3D | 1 338 617 | 0 | |
| rev_3Dto2D | 728 096 | 0 | |
| clathrin | 13 362 366 | 0 | |
| mem_localization | 1 255 726 | 0 | |

`set_rMaxLimit()` estimates a 2D reaction's reach as `3 sqrt(6 Dtot dt)`, while
`determine_2D_bimolecular_reaction_probability()` uses `3.5 sqrt(4 Dtot dt)`
over a Dtot that `add_2D_rotational_diffusion()` and `discretize_2D_Dtot()`
have both revised. `rev_2D` is the one model whose `rMaxLimit` is itself set by
a 2D reaction; everywhere else a 3D reaction sets it and covers the 2D one.

Commit 5 therefore exempts pairs with both complexes on the surface. The
underlying gap is wider than that commit and pre-existing: with sub-volumes
5.81 nm wide and a reach to 6.43 nm, `rev_2D`'s cell list already cannot
guarantee it offers every reacting pair. Closing it means changing
`set_rMaxLimit()`, hence `rMaxLimit`, hence the grid — a separate change.

### 19.5 Timing

`interleaved_timing.sh`, median of 3 repetitions, builds interleaved per case,
machine otherwise idle.

| case | base (s) | after (s) | speedup |
|---|---|---|---|
| mem_localization | 2.457 | 0.784 | **3.13×** |
| rev_3Dto2D | 5.967 | 3.437 | 1.74× |
| michaelis_menten | 4.365 | 2.597 | 1.68× |
| clathrin | 3.362 | 2.082 | 1.62× |
| closed_homoTrimer | 6.364 | 4.409 | 1.44× |
| homoTrimer | 6.389 | 4.441 | 1.44× |
| hetTrimer | 6.152 | 4.488 | 1.37× |
| trimer | 5.757 | 4.206 | 1.37× |
| implicit_lipid | 2.631 | 2.108 | 1.25× |
| hexamer | 5.112 | 4.183 | 1.22× |
| rev_3D | 3.572 | 3.386 | 1.06× |
| sphere | 5.957 | 5.865 | 1.02× |
| rev_2D | 41.332 | 42.105 | **0.98×** |
| | 99.417 | 84.091 | 1.18× |

`mem_localization` is the case commit 6 was aimed at: 3755 of its 3955
molecules are lipids, and 99.6% of the pairs the stencil offered were
lipid-lipid pairs no reaction names.

`rev_2D` is 2% slower and is the one case that gets nothing back. It is pure
2D, so commit 5's filter is exempted for almost all of its pairs, and its two
molecule types do react, so commit 6's mask never fires; its runtime is
dominated by the 2D reaction-probability tables rather than by the search. What
it does get is the finer grid, whose cost it pays without the benefit.

`sphere` is flat for the same structural reason section 18's audit gave: the
grid is laid over the (2R)³ bounding cube, 94.4% of it can never hold a
surface-bound molecule, and none of these seven commits changes that. A
surface decomposition for `OnSurface` complexes is the outstanding item.

## 20. A latitude-longitude index for the spherical membrane

Measured 2026-08-25, same host and settings as section 19. Parent is
`6111749`, SHA-256 `6406a0a97a92c454`; the tip after this work is
`cb83c595accd2df5`.

Section 19 left one item open: the Cartesian grid covers the `(2R)^3` box that
bounds a spherical system, so for molecules pinned to the shell it spends its
resolution on interior volume that is empty by construction, and it bins on
Cartesian coordinates while `get_distance()` measures those pairs by geodesic
arc.

### 20.1 What was left to win, measured first

After the seven changes in section 19, the pairwise search on the two
spherical samples looks like this:

| sample | pairs offered / step | of which surface-surface |
| --- | --- | --- |
| `sphere`, R=100 | 0 | 0 |
| `gagsphere`, R=70 | 1915 | 70 (3.6%) |

`sphere` offers nothing at all: its only reaction binds A to an implicit lipid,
so the molecule-type mask from `51e5d33` already discards every pair. On the
R=70 sphere, 96.4% of what is offered has at least one molecule off the
surface, and a shell index does not touch those.

Where that model's time actually goes, from 11662 samples of `sample`:

| samples | frame |
| --- | --- |
| 11656 | `determine_2D_bimolecular_reaction_probability` |
| 11655 | `integrator` |
| 8433 | `create_survMatrix` |
| 3222 | `create_pirMatrix` |
| 824 / 611 / 478 / 321 | `j0` / `y0` / `j1` / `y1` |

Essentially all of it is 2D reaction-table construction, which happens only for
pairs already within `Rmax` and which no neighbour structure can reduce. The
conclusion going in was therefore that a shell index cannot pay on any current
sample, and the measurements below bear that out. It was built for spheres
large enough that the bounding cube exceeds the sub-volume budget and the
interior starts taking resolution from the shell.

### 20.2 Sizing, and why the stencil is complete

Bands one cutoff wide in colatitude, each divided into cells one cutoff wide
along its own narrowest parallel. Both bounds come from the haversine identity

    sin^2(g/2) = sin^2((tA - tB)/2) + sin(tA) sin(tB) sin^2(dp/2)

whose terms are separately non-negative, so a pair at central angle `g` obeys
`|tA - tB| <= g` and `sin(|dp|/2) <= sin(g/2) / sqrt(sin tA sin tB)`.

Checked against an `O(N^2)` sweep by
[`shell_index_test.cpp`](../shell_index_test.cpp):

| R (nm) | cutoff (nm) | bands | cells | pairs in range | missed | double counted |
| --- | --- | --- | --- | --- | --- | --- |
| 70 | 19.45 | 10 | 106 | 188 539 | 0 | 0 |
| 100 | 9.77 | 28 | 962 | 23 389 | 0 | 0 |
| 1000 | 30 | 94 | 11 040 | 4 882 | 0 | 0 |
| 50 | 19.45 | 7 | 50 | 368 244 | 0 | 0 |
| 30 | 19.45 | 4 | 14 | 559 316 | 0 | 0 |
| 500 | 5 | 282 | 100 818 | 545 | 0 | 0 |
| 70 | 2 | 98 | 12 104 | 4 582 | 0 | 0 |

Over-inclusion runs 3.5x to 4.5x, against a floor of `9/pi` = 2.86x for a
nine-cell stencil in two dimensions.

### 20.3 The activation gate, and what it was worth

The first version activated for any spherical system. That made `sphere` **9%
slower** while staying bitwise identical: 448 surface molecules rebinned every
step, an `acos` and an `atan2` apiece, to fill an index whose every pair the
type mask then discarded.

The index now also requires that some two explicit molecules can react with
each other at all. With the gate:

| | `sphere` | `cluster_gagsphere` |
| --- | --- | --- |
| index | inactive | 10 bands, 106 cells |
| before gate | 0.909x | - |
| after gate | 0.997x, bitwise identical | see below |

### 20.4 Bitwise and statistical verdict

17 of the 18 benchmark cases are bitwise identical, timing 0.99x to 1.00x.
`cluster_gagsphere` is the only case where the index activates and the only one
that moves. Over 12 seeds at nItr 1000 its plateau-averaged copy numbers agree
within seed scatter, largest `|z|` **0.29**, including all three gag-gag
surface reactions:

| species | before | after | z |
| --- | --- | --- | --- |
| `gag(mem!1).IL(m!1)` | 65.500 +- 1.520 | 66.083 +- 1.510 | +0.27 |
| `gag(homo!1).gag(homo!1)` | 19.583 +- 1.264 | 19.083 +- 1.138 | -0.29 |
| `gag(het2!1).gag(het1!1)` | 1.833 +- 0.458 | 1.833 +- 0.423 | 0.00 |

### 20.5 Timing on the one case that uses it

A single seed reads 0.842x, which is not the overhead. This model's cost
follows its trajectory, because it is dominated by how many distinct 2D tables
get built, and the trajectory differs by construction once the pair order
changes. Across five seeds:

| build | median | range |
| --- | --- | --- |
| before | 44.89 s | 32.08 - 71.85 |
| after | 45.11 s | 32.28 - 64.94 |

A 2.2x spread between seeds against a 0.5% difference in medians: neutral, and
the single-seed figure is noise.

## 21. Re-audit: cumulative drift and three residual defects

Measured 2026-08-26, same host and settings as sections 19 and 20. Parent is
`63fc2b1`, SHA-256 `cb83c595accd2df5`; the tip after this section is
`fe381b55964209a7`.

Sections 19 and 20 verified every commit against its own parent. This section
asks the two questions that per-commit checking does not answer: does the whole
sequence drift anywhere it was not meant to, and is anything left in the
decomposition that should be fixed.

### 21.1 Cumulative drift, original baseline to tip

Running the suite at the tip against the manifests from `base` -- the build
before any of this work -- six of the eighteen cases differ, and they are
exactly the union of the cases each stream-changing commit was expected to
move:

| case | changed by | why |
| --- | --- | --- |
| trimer | `459bd83` | z floor forced SubBoxes below the interaction range |
| rev_2D | `459bd83`, `77daa3a` | same floor, then the 30 cap |
| compartment | `a81c7a6`, `77daa3a` | the maxPairs budget, then the 30 cap |
| rev_3D | `77daa3a` | the 30 cap |
| rev_3Dto2D | `77daa3a` | the 30 cap |
| cluster_gagsphere | `336322d` | the shell index activates |

The other twelve are byte-identical from the original baseline all the way to
the tip. **No unexpected drift.**

### 21.2 One timing outlier, and what it was

That run recorded `homoTrimer` at 1698.693 s against 6.060 s for the baseline,
a 280x slowdown, while still producing byte-identical output. It is not the
decomposition: the grid for that model is 10^3 SubBoxes of 32.1 nm against a
29.7 nm range, which is correct and small.

The same binary, input and seed re-run three times gives 3.97, 3.95 and 3.91 s
wall at 3.86, 3.88 and 3.85 s user. The same binary measured 4.057 s and
4.062 s in the two preceding suite runs. The 1698 s is the host, and this
specific case has a history of it -- the preamble to these results records a
535 s repetition of `homoTrimer` under the same conditions.

The lesson for reading section 19's and 20's tables: the bitwise verdicts are
load-independent and the interleaved timings are load-tolerant, but a
single-shot suite timing on this host is not evidence.

### 21.3 Three residual defects, all fixed and all result-preserving

| commit | defect | reachable today |
| --- | --- | --- |
| `3d253cc` | shell cutoff angle used the small-angle limit of the chord bound | yes, by margin only |
| `391a1a3` | the sub-volume budget could fail to bind | no sample reaches it |
| `a576661` | the re-bin restart skipped molecule 0 | no sample reaches it |

**`3d253cc`.** `ShellIndex` sized cells from `rMaxLimit / radiusFloor`, the
small-angle limit of `2*asin(rMaxLimit / 2*radiusFloor)`, and smaller than it
since `asin(u) >= u`. rMaxLimit is a chord bound, and a chord subtends its
widest angle at the smallest radius, so the floor is the worst case and the
exact form is the one to use. The old form was safe only by the gap between
where surface complexes sit (0.971 and 0.995 of sphereR) and the 0.9 floor.
Cost: gammaCut moves 0.4% on the R=70 sample, leaving band and cell counts
unchanged.

`shell_index_test.cpp` now checks the property the simulation relies on --
a pair within `rMaxLimit` by straight-line distance must be offered -- with
points at the radius floor, keeping the angular check as a separate count:

| R (nm) | cutoff | bands | cells | pairs within cutoff | missed | angle-missed | double counted |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 70 | 19.45 | 10 | 106 | 190 048 | 0 | 0 | 0 |
| 100 | 9.77 | 28 | 962 | 23 416 | 0 | 0 | 0 |
| 1000 | 30 | 94 | 11 040 | 4 882 | 0 | 0 | 0 |
| 50 | 19.45 | 7 | 50 | 374 103 | 0 | 0 | 0 |
| 30 | 19.45 | 4 | 12 | 584 039 | 0 | 0 | 0 |
| 500 | 5 | 282 | 100 818 | 545 | 0 | 0 | 0 |
| 70 | 2 | 98 | 12 102 | 4 583 | 0 | 0 | 0 |

**`391a1a3`.** The budget scaled all three axes by one cube-root factor and
clamped each at one SubBox. The factor assumes every axis absorbs its share and
an axis at one cannot, so the product could stay above the cap:

| grid requested | old result | new result | cap |
| --- | --- | --- | --- |
| 1 x 1 x 1 197 604 | 1 x 1 x 895 000 | 1 x 1 x 500 000 | 500 000 |
| 7 x 7 x 8 399 (from 11 x 11 x 11 976) | within | 411 551 | 500 000 |
| 79^3 (from 2000^3) | within | 493 039 | 500 000 |
| 56^3 rev_3D | unchanged | 175 616 | under |
| 3^3 trimer | unchanged | 27 | under |

**`a576661`.** The checked pass restarts binning after
`put_back_into_SimulVolume()`, and the restart was `molItr = 0` inside a
`++molItr` loop, so it resumed at molecule 1. `clear_member_lists()` had just
emptied every list, so molecule 0 ended the step absent from the grid along
with every pair it belongs to. The counter is now signed and the restarts set
it to -1.

The message that marks this path, "Attempting to fit back into box", appears
**zero times across all 18 cases**, which is why it survived and why the fix
changes nothing measurable.

### 21.4 Verdict

All three are bitwise identical across all 18 cases, so the tip carries
sections 19 and 20's results unchanged. The rebuilt tip hashes to
`fe381b55964209a7`, matching the binary these checks were run against.

## 22. Flattening `Molecule` to structure-of-arrays: built, measured, reverted

A review of this branch objected that

> base data structure classes such as `Molecule` [are] a pointer-heavy
> array-of-structures with many nested `std::vector`s and an `unordered_map`,
> unsuitable for efficient device access without conversion to flat
> structure-of-arrays storage.

Two conditions were put on any such conversion: **bitwise identity**, and
**similar or higher efficiency**. Both were tested by building a flattening
rather than by arguing about it, together with a control that measures the same
axis from the opposite side, and a size sweep because the answer turned out to
depend on size.

Summary of what was found:

| | verdict |
| --- | --- |
| Bitwise identity, host-side flattening | **achievable** -- the build here is byte-identical on all 13 cases |
| Bitwise identity, device | **already measured unachievable**, and not for layout reasons (22.2) |
| Layout headroom exists | **yes, and it grows with size** -- inflating the stride costs 1.21x at 40,000 molecules (22.6) |
| The flattening built here captures it | **no** -- 0.90x to 0.98x everywhere (22.5, 22.6), for the reason in 22.7 |
| Device offload at NERDSS's current sizes | **three to five orders of magnitude below break-even** (22.8) |

Measured 2026-08-27, same host and settings as section 21: Apple M5, 10
physical cores, L1D 128 KB and L2 16 MB on the performance cores, Apple clang
21.0.0, `-O3 -std=c++0x`, GSL 2.8, `make serial`, seed 20260810, at `cc1d954`
plus the working tree's whitespace-only reformat of `EXEs/nerdss.cpp`. The
reference build hashes to `c1db1694d0724b0c`. Everything below was reverted
afterwards and the tree rebuilt to that same hash, so none of it is still in the
source.

### 22.1 What the class actually holds

| | bytes | cache lines |
| --- | ---: | ---: |
| `Molecule` | 656 | 10.25 |
| `Molecule::Iface` | 72 | 1.125 |
| `Complex` | 312 | 4.875 |

`Molecule` owns 22 heap-allocating members: `interfaceList`, `tmpICoords` and
20 `std::vector`s of `int`, `double` or `std::array<int,3>`. At 24 bytes per
header that is 528 of the 656 bytes -- **80% of the object is pointer triples**.

The fields the pairwise search reads before deciding a candidate pair is worth
evaluating are spread over four of those ten lines:

| field | offset | line |
| --- | ---: | ---: |
| `myComIndex`, `molTypeIndex` | 0, 4 | 0 |
| `comCoord` | 32 | 0 |
| `interfaceList` | 56 | 0-1 |
| `isImplicitLipid` | 91 | 1 |
| `freelist` | 160 | 2 |
| `bndlist` | 184 | 2-3 |
| `bndpartner` | 208 | 3 |

So the rejection cascade in `check_bimolecular_reactions()` touches about four
cache lines per molecule, eight per candidate pair, to read 40 bytes of
information. That part of the objection is accurate, and 22.6 shows it is not
free.

**The `unordered_map` is not accurate.** `Molecule::mapIdToIndex` and
`Complex::mapIdToIndex` are `static` members -- one map per process, not one per
object -- and the only translation unit that references them is
`src/mpi/id_index.cpp`. `src/mpi` is in the Makefile's `INCLUDE_FOLDERS` only
for `make mpi`; the serial and OpenMP builds never compile it. The maps
contribute zero bytes per molecule and are touched on no hot path measured in
this document.

### 22.2 What bitwise identity constrains, and what it does not

A layout change moves bytes; it does not change arithmetic. If a flattened
container holds bit-exact copies of the same values, and the same code performs
the same operations in the same order on them, the output is identical. That is
why 22.4 works.

Two things are constrained, and neither is a property of the container.

**Append order is load-bearing.** `determine_if_reaction_occurs()` walks
`mol.crossbase` in list order and draws one `rand_gsl64()` **per entry**,
returning on the first draw that lands under that entry's probability. The
number of random numbers consumed therefore depends on *where in the list* the
accepted entry sits. Any scheme that changes the order in which the sweep
appends to `crossbase`, `probvec`, `mycrossint` and `crossrxn` changes the
random stream from that point onward and ends bitwise comparison. That is
exactly why `openmp-production-lane` had to reconstruct the serial append order
with conflict waves instead of simply parallelising the loop. It applies to
array-of-structures storage just as much as to flat storage.

**Device arithmetic is not host arithmetic.** This repository's own CUDA proof
of concept (`gpu_poc/`, commit `303e620`) flattened precisely the independent
numerical prefix of the 3D path and compared it against the same
double-precision CPU functions. It reports zero discrete mismatches and a
**maximum relative error of 9.55e-15** for reaction checking, 8.88e-16 for
propagation. Nonzero. The probability path calls libm `erfc` and `exp` and
`Faddeeva::erfcx`, whose last-ulp results are implementation-defined, and a
device math library is a different implementation. **Flattening for device
execution and bitwise identity are already measured to be incompatible here,
and the obstacle is the math library, not the container.** A device lane has to
be validated statistically, the way sections 5 and 10.2 validate the
RNG-changing group.

### 22.3 The two builds

**The flattening.** `MolHotView`, 32 bytes, two to a cache line:

```cpp
struct MolHotView {
    double x, y, z;          // comCoord, bit-exact copy
    int molTypeIndex;
    unsigned int flags;      // implicitLipid | hasBonds | hasFree | comOnSurface
};
```

`comOnSurface` folds in `complexList[myComIndex].OnSurface`, so the gates stop
following the molecule-to-complex indirection as well. The mirror is rebuilt
once per timestep by an O(N) pass placed beside `shellIndex.rebin()`, which
carries the same precondition: both must describe the molecule set the search is
about to walk. That placement is sound because the search itself never moves a
molecule, changes a bond list, or moves a complex on or off the surface --
`check_implicit_reactions()` and `check_compartment_reaction()` compute
probabilities and nothing else, and every association runs later in the
timestep.

Five gates were rewritten to read the mirror -- the implicit-lipid test, the
`rxnPartners` type scan, the `excludeVolumeBound` test, the already-bound test
and the `rMaxLimit` distance cutoff -- plus the four `typeMask` tests in the
cell walk in `EXEs/nerdss.cpp`. Everything past the gates falls through to the
existing code unchanged. The cascade's per-molecule footprint drops from about
four cache lines to one.

**The control.** A second build appended `char[1024]` to the tail of `Molecule`.
No hot field's offset moves; only the array stride, from 656 to **1,680 bytes**.
If the search is limited by molecule-data cache traffic, this build pays for it,
and the size of what it pays is the headroom any layout change is competing for.

### 22.4 Bitwise

`run_suite.sh` over `cases.tsv`, then `compare_suites.sh` against the reference
build: **13 of 13 cases byte-identical, for both builds.** So a host-side
flattening satisfies the first condition, and this is the demonstration.

### 22.5 Timing at the suite's sizes

`interleaved_timing.sh` over `debugremoval_cases.tsv`, three builds interleaved
per case, 7 repetitions, CPU time (user+sys) rather than wall clock, medians:

| case | nItr | base (s) | pad (s) | mirror (s) | pad/base | mirror/base | base sd |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `clathrin` | 150,000 | 3.920 | 3.880 | 3.800 | 1.010 | 1.032 | 0.314 |
| `hexamer` | 25,000 | 7.800 | 7.960 | 8.070 | 0.980 | 0.967 | 0.345 |
| `implicit_lipid` | 100,000 | 3.750 | 3.750 | 3.820 | 1.000 | 0.982 | 0.092 |
| `mem_localization` | 2,000 | 0.790 | 0.840 | 0.850 | 0.940 | 0.929 | 0.366 |
| `michaelis_menten` | 150,000 | 4.770 | 4.790 | 4.810 | 0.996 | 0.992 | 0.095 |
| `rev_3D` | 20,000 | 6.890 | 7.250 | 7.140 | 0.950 | 0.965 | 0.143 |
| **total** | | **27.920** | **28.470** | **28.490** | **0.981** | **0.980** | |

Multiplying the stride by 2.56 costs 1.9%; dividing the hot footprint by about
20 costs 2.0%. Same magnitude, same sign, and every per-case difference inside
one or two standard deviations. Read on its own this says the layout axis does
not move the program.

**That reading would have been wrong.** These cases hold 100 to 3,955
molecules, 66 KB to 2.6 MB of `Molecule`, against a 128 KB L1D and a 16 MB L2.
The suite cannot see an effect that only appears when the array stops fitting.
`cases.tsv` was built to span the sizes the branch's other changes cared about,
and it is not a size sweep.

### 22.6 Timing above the suite: the headroom is real and it grows

Four larger models. `scale_10k` and `scale_40k` are `rev_3D` with the copy
numbers multiplied and the box side scaled by the cube root of the same factor,
so the density, the reactions and the timestep are unchanged and only the size
moves. `enzyme` is `sample_inputs/enzyme/parms_clat_enzyme.inp`, the largest
model in `sample_inputs` that actually runs. Seven repetitions for the first
three, five for `scale_40k`; medians of CPU time:

| case | molecules | `Molecule` bytes | base (s) | pad/base | mirror/base | base sd | pad sd |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `scale_2k` | 2,000 | 1.3 MB | 4.630 | 0.953 | 0.971 | 0.217 | 0.269 |
| `enzyme` | 6,410 | 4.2 MB | 5.160 | 0.942 | 0.974 | 0.286 | 0.396 |
| `scale_10k` | 10,000 | 6.6 MB | 10.000 | 0.864 | 0.899 | 1.175 | 2.585 |
| `scale_40k` | 40,000 | 26.2 MB | 29.540 | **0.827** | 0.979 | 1.293 | 1.436 |

**The padding penalty rises monotonically with the working set, from within
noise at 2,000 molecules to 1.21x at 40,000.** At `scale_40k` the gap is 6.2 s
against a 1.29 s baseline standard deviation -- more than four standard
deviations, and the one measurement in this section that is unambiguous on its
own. Molecule-array memory traffic *is* on the critical path once the array is
large, which is the substance of the review's complaint.

The mirror does not capture it. It is 0.90x to 0.98x at every size, never
faster, and its deficit does not shrink as the headroom grows.

### 22.7 Why the mirror fails to collect headroom that exists

The mirror shadows the class; it does not replace it. Two consequences.

It only covers the gates. A pair that survives them still walks the full
656-byte `Molecule` for `freelist`, `interfaceList` and `bndpartner`, and the
312-byte `Complex` for `D`, `Dr` and `radius`. The bytes the padding control
penalises are mostly still being touched.

And the refresh is O(N) while the search is O(pairs). Rebuilding the mirror
writes 32 bytes per molecule per timestep whether or not that molecule is ever
examined. Measured pair counts (22.8) put `scale_40k` at 19,841 candidate pairs
per timestep -- about 39,700 molecule visits -- against 40,000 refresh writes.
**The refresh touches more molecules than the search it is accelerating.** At
`scale_10k`, 1,132 pairs per timestep against 10,000 refresh writes, it is worse
by an order of magnitude, and `scale_10k` is indeed the worst case in the table
at 0.899x.

That is a verdict on this design, not on flattening. A genuine conversion --
where the flat arrays *are* the storage, so there is nothing to refresh and the
whole object is narrow rather than just its gate fields -- carries neither cost
and would be competing for the headroom the padding control measures. It is also
a far larger change: `interfaceList` alone has 835 uses across `src`, `include`
and `EXEs`, and each one would need the bitwise verification 22.4 gives here
cheaply. Two smaller variants sit between the two and are not tested here:
refreshing the mirror lazily per cell rather than per timestep, and simply
reordering `Molecule`'s declarations so the gate fields share one line, which
costs nothing at runtime and needs no refresh at all.

### 22.8 The device question is a problem-size question

A third build counted calls to `check_bimolecular_reactions()` -- one per
candidate pair that survives the cell list and the type mask:

| case | molecules | nItr | pair checks | per timestep |
| --- | ---: | ---: | ---: | ---: |
| `clathrin` | 100 | 150,000 | 9,365,681 | **62** |
| `rev_3D` | 2,000 | 20,000 | 1,914,860 | **96** |
| `hexamer` | 1,000 | 25,000 | 29,993,419 | **1,200** |
| `scale_10k` | 10,000 | 5,000 | 5,661,633 | **1,132** |
| `enzyme` | 6,410 | 2,000 | 11,363,595 | **5,682** |
| `scale_40k` | 40,000 | 2,000 | 39,681,693 | **19,841** |

The GPU proof of concept measured its end-to-end break-even -- kernel plus
transfers, for the arithmetic prefix alone, with no packing or scatter cost --
at about **1,000,000 pairs** per launch: 1.02x at 1,048,576, 0.92x at 262,144,
0.97x at 65,536, 0.24x at 1,024. NERDSS produces 62 to 19,841 pairs per
timestep across everything measurable here, and the pairs cannot be batched
across timesteps because each step's positions depend on the previous step's
outcome. `clathrin`'s **entire 150,000-iteration run** produces 9.4 million
pairs, nine break-even batches in total.

The dependency structure gives the same answer independently:
`openmp-production-lane` measures about 41 pairs per conflict wave, against the
parallel-region costs in section 13.2 and the break-even widths in 13.3.

Flattening changes none of this. It is a prerequisite for device execution, not
a reason for it.

### 22.9 What flattening the ragged lists could buy, at most

The 20 per-molecule vectors could be replaced by CSR-style offset and value
arrays, and that can be made order-preserving and therefore bitwise identical.
The ceiling on what it would win is already measured.

`crossbase`, `mycrossint`, `crossrxn` and `probvec` are emptied with `.clear()`
at `EXEs/nerdss.cpp:1808-1811`, and `clear_reweight_vecs()` swaps the six
`curr*` lists into the six `prev*` lists and clears the `curr*` side. `clear()`
retains capacity, so after the first few timesteps the sweep performs **no
allocations at all**, and section 16.2 measured the allocator at **0.94%** of
leaf samples once the `std::string` temporary was removed. That is the whole
budget a perfect CSR conversion is competing for -- separately from the cache
headroom in 22.6, which is about the 656-byte stride, not about the heap blocks.

Against it: rows grow by data-dependent amounts during the sweep, so a flat
layout needs either a per-row capacity cap -- a new failure mode when a dense
model exceeds it -- or a counting pass, which doubles the sweep. And rows are
emptied individually mid-timestep (`associate_box.cpp:1033` clears exactly the
two reacting molecules while their partners' rows still name them), so the
representation must support single-row truncation too.

### 22.10 How large the models actually are, which is the crux

The size at which this question changes answer is between 6,410 and 10,000
molecules on this host, and `cases.tsv` tops out at 3,955. What the repository
ships is a wider range than the suite covers:

| model | explicit molecules | status |
| --- | ---: | --- |
| `clathrin_coat/flat_clat-ap2-pip2.dir/parms.inp` | 500,400 | does not reach the timestep loop, see below |
| `clathrin_pioneer/parmsMacroRate.inp` | 10,310 | `pip2.mol` absent from the directory |
| `enzyme/parms_clat_enzyme.inp` | 6,410 | runs; used above |
| `clathrin_coat/.../parms_clath_ap_pip.inp` | 5,400 | not measured |
| `mem_localization` | 3,955 | largest in `cases.tsv` |
| `gagsphere` | 2,500 | used in section 13.4 |

The 500,400-molecule model is explicit, not implicit: `pip2.mol` sets `isLipid`
but not `isImplicitLipid`, and the run's own output confirms it with
`NUMBER OF MOLECULES IN GEN COORDS: 500400`. That is 328 MB of `Molecule` and
156 MB of `Complex` before any interface storage. It did not get past setup in
ten minutes on this host, and the reason is not the layout: the line after that
message is `fixOverlappingMolecules()`, whose fallback branch at
`generate_coordinates.cpp:354` is an all-pairs O(N^2) double loop repeated up to
`MAX_ITERATIONS = 50`. At N = 500,400 that is up to 1.25e13 pair tests before
the first timestep.

So the honest position on the review is: it is right that the layout costs
something, it is right that the cost grows with size, and the sizes at which it
starts to matter are ones this repository ships models for and this branch's
suite does not measure. It is wrong about the `unordered_map`, and the
conclusion it draws -- convert to flat storage for device access -- runs into
22.2 and 22.8 rather than into any property of the containers. The first thing
worth doing is not a conversion; it is a benchmark at 10,000 molecules and
above, because as 22.5 and 22.6 together show, a suite that stops at 3,955
returns "no effect" for a change that has a 1.21x effect at 40,000.

## 23. Narrowing `Molecule`, stage 1: the reweighting vectors

Section 22 established that `Molecule`'s 656-byte stride costs real time at
scale -- a control build that inflated it to 1,680 bytes measured 0.827x at
40,000 molecules -- and that the `MolHotView` proof of concept failed to collect
that headroom because it shadowed the class instead of narrowing it. This is the
first stage of narrowing it for real, on branch `molecule-layout`. The plan for
the remaining stages is in
[`docs/molecule_layout_plan.md`](../../docs/molecule_layout_plan.md).

Measured 2026-08-27, Apple M5, Apple clang 21.0.0, `-O3 -std=c++0x`, GSL 2.8,
seed 20260810, parent `cc1d954`, baseline binary `c1db1694d0724b0c`.

### 23.1 What changed

Twelve `std::vector`s -- `prevlist`, `currlist`, `prevmyface`, `currmyface`,
`prevpface`, `currpface`, `prevnorm`, `currprevnorm`, `ps_prev`, `currps_prev`,
`prevsep`, `currprevsep` -- were two groups of six **parallel arrays**, not
twelve independent lists. Every read indexes all six at one position:
`determine_3D_bimolecular_reaction_probability()` scans `prevlist` for a
(partner, myFace, partnerFace) match and then reads `prevsep`, `ps_prev` and
`prevnorm` at that same subscript. Every write pushes six values in lockstep.

They are now one `Molecule::ReweightEntry` and two vectors of it. 288 bytes of
vector header become 48, and the lookup follows one pointer instead of six.

**`sizeof(Molecule)`: 656 -> 416 bytes**, 10.25 cache lines to 6.5.

### 23.2 Bitwise

| table | result |
| --- | --- |
| `cases.tsv` | 13 of 13 byte identical |
| `coverage_cases.tsv` | 5 of 5 byte identical |
| restart read path | 15 files byte identical |
| `make mpi` | builds clean |

`coverage_cases.tsv` was run because `cases.tsv` does not enter every edited
path: `gagsphere` is the only model reaching the `excludeVolumeBound` half of
`check_bimolecular_reactions()` and `compartment` the only one reaching the
transmission path.

The restart **read** path needed its own check. The suite hashes `RESTARTS/`, so
it proves restart files are *written* identically and says nothing about reading
them -- and the on-disk format still carries six separate length-prefixed
arrays, which the reader now has to reassemble. Both builds were continued from
the same `rev_3D/RESTARTS/restart10000.dat`; every file they then produced
matched.

One edited file is covered by no test at all.
`determine_1D_bimolecular_reaction_probability()` is called only from the
`com1.onFiber && com2.onFiber` arm of `check_bimolecular_reactions()`, and no
input under `sample_inputs` sets `onFiber`. Section 15's rule applies: a bitwise
suite that never reaches a function reports "identical" in the same words it
uses for code it does reach. It was missing from `known_uncovered.tsv`, which
listed the fiber *sweep* but not the probability function behind the same gate;
it and its two exclusive callees are now listed there, and its change is argued
by reading rather than by the suite.

### 23.3 A measurement change, because this host is never idle

The first timing pass reported 1.196x and 1.200x against baselines that were
three times their own stage-0 values, with three unrelated `nerdss_mpi`
processes at 99% CPU. Section 11.6 records this failure mode and says to check
`uptime` first -- which tells you when to distrust a number, not how to take
one, on a desktop that is rarely idle.

`/usr/bin/time -l` on Apple silicon reports **retired instructions** and
**elapsed cycles**. Four repetitions of each build on `scale_10k` under that
same load:

| | instructions | spread | cycles | spread | wall |
| --- | ---: | ---: | ---: | ---: | ---: |
| base | 89.83 G | 0.1% | 104.9 G | 1.5% | inflated ~3x |
| stage 1 | 87.68 G | 0.05% | 85.5 G | 0.6% | inflated ~3x |

The counters are nearly load-proof where wall and CPU time are not. They also
separate the two things a change can do, which is the distinction this whole
line of work turns on: **instructions falling means less work; cycles falling at
an unchanged instruction count means fewer stalls.**
`interleaved_timing.sh` now records CPU time, wall time, instructions, cycles
and peak RSS, and prints instructions per cycle.

### 23.4 Result

> **Superseded, and by how much matters.** The table below was measured on a
> host carrying three unrelated `nerdss_mpi` processes at 99% CPU plus a Python
> at 98.6%. Section 24 re-measures all three builds together under far lighter
> load and gets 1.058x to 1.172x where this says 1.119x to 1.389x. **Contention
> inflated stage 1's apparent gain by up to 30 points.** 23.7 warned that it
> might; it did. The numbers here are kept because the warning is only worth
> anything next to what it was warning about -- use section 24's.

Interleaved, medians; `scale_2k`, `enzyme` and `scale_10k` at 5 repetitions,
`scale_40k` at 7:

| case | molecules | instructions | cycles | CPU time | RSS | IPC base -> stage 1 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `scale_2k` | 2,000 | 1.032x | 1.119x | 1.113x | 1.018x | 1.908 -> 2.069 |
| `enzyme` | 6,410 | 1.022x | 1.389x | 1.383x | 1.182x | 1.783 -> 2.422 |
| `scale_10k` | 10,000 | 1.024x | 1.246x | 1.253x | 1.032x | 0.846 -> 1.029 |
| `scale_40k` | 40,000 | 1.014x | 1.264x | 1.223x | 1.095x | 1.018 -> 1.270 |

The one thing here that section 24 confirms rather than overturns is the shape:
**instructions fall by 1.4% to 3.2% while cycles fall by much more.** Merging six
`push_back`s into one is the whole instruction saving, and it is far too small to
explain the cycle saving; IPC rises to make up the difference. That is a stall
reduction, which is what a narrower object is supposed to produce and what no
amount of doing-less-work could imitate.

### 23.5 RSS confirms the change landed where it was aimed

240 bytes per molecule should disappear. Measured against predicted:

| case | predicted | measured | ratio |
| --- | ---: | ---: | ---: |
| `scale_2k` | 0.48 MB | 0.56 MB | 1.16 |
| `enzyme` | 1.54 MB | 2.03 MB | 1.32 |
| `scale_10k` | 2.40 MB | 2.83 MB | 1.18 |
| `scale_40k` | 9.60 MB | 11.17 MB | 1.16 |

Consistently 16-32% *more* than the struct-size prediction, which is the second
half of the change showing up: twelve heap-owning members became two, so a
molecule that ever recorded an encounter also stops carrying up to ten separate
allocator blocks. RSS is barely affected by load, and section 24 reproduces
these figures within a few percent.

### 23.6 A size-dependence claim that did not survive re-measurement

This section originally read: "The size dependence is not the whole story",
and argued from the table in 23.4 that `enzyme` at 6,410 molecules gained more
(1.389x) than `scale_40k` at 40,000 (1.264x), so part of the win had to come
from the reweighting *lookup* -- six pointer chases becoming one -- which scales
with reweighting traffic rather than with molecule count.

**The premise was an artifact of the contaminated measurement.** Under the
lighter load of section 24 the ordering is monotone in size after all -- 1.058x
at 2,000, 1.059x at 6,410, 1.101x at 10,000, 1.172x at 40,000 -- which is the
same shape section 22.6's padding control produced. `enzyme` is not special; it
was simply the case whose measurement contention distorted most.

The lookup argument may still be true in part, but nothing here measures it, and
the evidence that was offered for it has evaporated. Recorded rather than
deleted because the reasoning error is the useful part: a mechanism was inferred
from a ranking, and the ranking was noise.

### 23.7 Caveat on the magnitudes, and how it turned out

Every number in 23.4 was taken on a host carrying one to two cores of unrelated
load. Interleaving means both builds met the same conditions and the counters
are robust to it, but memory contention plausibly *amplifies* a locality win:
with other cores pressing on a shared L2, the build with the larger working set
suffers more than it would alone. The load-proof floor is the instruction
reduction, about 2%.

That caveat was correct and the amplification was large. Section 24's re-measure,
with all three builds interleaved under lighter load, puts stage 1 at 1.058x to
1.172x rather than 1.119x to 1.389x. **The right lesson is not that the
counters failed -- instructions and RSS reproduced to within a few percent --
but that cycles are load-proof only in the sense of being measurable, not in the
sense of measuring the same thing.** Stalls are a property of the machine's
memory system at that moment, and that moment included four busy cores.

### 23.8 Verdict

Stage 1 meets every exit criterion its plan set: byte identical on both case
tables and across a restart, `sizeof(Molecule)` at 416, MPI still building, and
a measured gain rather than the flat result section 22 recorded for the proof of
concept. Stage 2 -- merging `crossbase`, `mycrossint`, `crossrxn` and `probvec`
into one vector, 416 -> 344 bytes -- is warranted on this evidence.

## 24. Narrowing `Molecule`, stage 2: the crossing lists

Stage 2 of the plan in
[`docs/molecule_layout_plan.md`](../../docs/molecule_layout_plan.md). Measured
2026-08-27, same host and build settings as section 23, with all three builds
interleaved in one pass so they meet identical conditions.

### 24.1 What changed

`probvec`, `crossbase`, `mycrossint` and `crossrxn` were four parallel vectors
indexed by one subscript everywhere they are read. They are now one
`Molecule::CrossEntry` and one vector of it: 96 bytes of vector header become 24,
and the traversals follow one pointer instead of four.

**`sizeof(Molecule)`: 416 -> 344 bytes.** With stage 1, 656 -> 344 -- 10.25 cache
lines to 5.38.

Two pieces of scaffolding went with it. `record_crossing_pair()` carried an
`alsoInitProbvec` flag whose only purpose was keeping the four lists in step;
with one entry there is nothing to keep in step. And
`determine_1D_bimolecular_reaction_probability()` held a seventh hand-written
copy of that helper -- the same two crossbase pushes, two mycrossint pushes, two
crossrxn pushes, two ncross bumps and two probvec seeds, in the same order --
which is now a call.

### 24.2 What had to be established first

The merge looked blocked. **Sixteen sites call `crossbase.clear()` alone**
mid-timestep, after an association, deliberately leaving the other three
populated for the rest of the timestep. One merged vector clears all four.

Three findings settle it.

*Nothing loops on the other three's sizes.* Every traversal in the program --
`determine_if_reaction_occurs()`, `perform_bimolecular_reactions()`,
`zero_partner_probvec()`, the four overlap sweeps, `Cluster` -- bounds itself by
`crossbase.size()` and reads the rest at that subscript. A search for
`probvec.size()`, `mycrossint.size()` and `crossrxn.size()` returns nothing but
a commented-out debug line. So after `crossbase.clear()` the leftovers are
unreachable: dead data until the end-of-timestep clear. Dropping them with the
same call changes nothing that can be read.

*Alignment holds while they are populated.* Checked rather than argued: an
instrumented build reported any molecule whose four lists disagreed in length at
the end of the pairwise sweep. **Zero reports across all 13 cases in
`cases.tsv` and all 5 in `coverage_cases.tsv`.**

*One site could have broken it, and is the one real semantic change.*
`determine_2D_implicitlipid_reaction_probability()` pushed `probvec`
unconditionally at entry, while the matching `crossbase`, `mycrossint` and
`crossrxn` pushes sit inside its `!isDissociated && rate > 0` branch. On any path
that skips that branch, the molecule ended the timestep with one more `probvec`
entry than `crossbase` entry -- and since every traversal reads them at a common
subscript, a later entry would have been paired with the *wrong* probability.
One entry carrying all four values cannot express that state, so the merge
removes the possibility.

The instrumented run says it never happens: either that branch is not reached by
any case here, or it is never taken early. So nothing measurable changes. But
this is a latent defect the merge closes rather than a mechanical rewrite, and
it is marked as such at the site.

### 24.3 Bitwise

Against the **pre-branch** baseline `c1db1694d0724b0c`, so this covers stages 1
and 2 together:

| table | result |
| --- | --- |
| `cases.tsv` | 13 of 13 byte identical |
| `coverage_cases.tsv` | 5 of 5 byte identical |
| restart read path | 15 files byte identical |
| `make mpi` | builds clean |

The MPI wire format is unchanged: `serialize_back()` and `deserialize_back()`
split and reassemble the same four arrays in the same order.

### 24.4 Result

Three builds interleaved, 5 repetitions, medians. Host load average 4.25 at the
start and 4.34 at the end -- two busy cores, against the four-plus that
contaminated section 23.4. Baseline standard deviations are 0.5% to 3%.

**Cycles**, the load-robust measure:

| case | molecules | base | stage 1 | stage 2 | stage 1/base | stage 2/base |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `scale_2k` | 2,000 | 14.29 G | 13.51 G | 13.17 G | 1.058x | 1.085x |
| `enzyme` | 6,410 | 18.31 G | 17.29 G | 16.98 G | 1.059x | 1.078x |
| `scale_10k` | 10,000 | 25.52 G | 23.18 G | 22.14 G | 1.101x | 1.153x |
| `scale_40k` | 40,000 | 85.52 G | 72.97 G | 64.52 G | 1.172x | **1.325x** |
| **total** | | **143.6 G** | **126.9 G** | **116.8 G** | **1.132x** | **1.230x** |

CPU time agrees closely: 1.134x and 1.226x aggregate, and 1.317x for
`scale_40k`.

**Instructions barely move**, which is the point:

| case | stage 1/base | stage 2/base |
| --- | ---: | ---: |
| `scale_2k` | 1.033x | 1.038x |
| `enzyme` | 1.022x | 1.025x |
| `scale_10k` | 1.023x | 1.028x |
| `scale_40k` | 1.013x | 1.016x |

A 2.2% instruction saving accompanies a 23.0% cycle saving. The difference is
stalls, and **instructions per cycle** shows it directly:

| case | base | stage 1 | stage 2 |
| --- | ---: | ---: | ---: |
| `scale_2k` | 3.734 | 3.825 | 3.904 |
| `enzyme` | 5.279 | 5.472 | 5.558 |
| `scale_10k` | 3.499 | 3.765 | 3.925 |
| `scale_40k` | 2.867 | 3.316 | 3.741 |

`scale_40k` gains 30% of IPC across the two stages while executing 1.6% fewer
instructions.

**RSS**, against a prediction of 312 bytes per molecule for the two stages:

| case | predicted | measured | ratio |
| --- | ---: | ---: | ---: |
| `scale_2k` | 0.62 MB | 0.77 MB | 1.24 |
| `enzyme` | 2.00 MB | 2.59 MB | 1.29 |
| `scale_10k` | 3.12 MB | 3.78 MB | 1.21 |
| `scale_40k` | 12.48 MB | 14.99 MB | 1.20 |

Again 20-29% beyond the struct arithmetic, from sixteen heap-owning members
becoming three.

### 24.5 The size dependence, restated

Under lighter load both stages are monotone in model size, which is the shape
section 22.6's padding control produced and which section 23.6 mistakenly
thought it had refuted:

| molecules | stage 1 | stage 2 |
| ---: | ---: | ---: |
| 2,000 | 1.058x | 1.085x |
| 6,410 | 1.059x | 1.078x |
| 10,000 | 1.101x | 1.153x |
| 40,000 | 1.172x | 1.325x |

This is what section 22 predicted and what the `MolHotView` proof of concept
failed to deliver. The control build that *widened* `Molecule` by 2.56x cost
0.827x at 40,000 molecules; narrowing it by 1.91x returns 1.325x there. The two
are the same effect measured in opposite directions, and the headroom the
padding identified is now being collected -- by narrowing the object rather than
shadowing it.

### 24.6 Verdict

Stage 2 meets its exit criteria: byte identical on both case tables and across a
restart, `sizeof(Molecule)` at 344, MPI building, and 1.230x aggregate against
the pre-branch baseline in cycles.

Stages 3 and 4 of the plan -- reordering the survivors so the rejection
cascade's fields share one cache line, and deciding about the association
scratch -- remain. Neither has been attempted.
