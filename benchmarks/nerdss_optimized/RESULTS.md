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
