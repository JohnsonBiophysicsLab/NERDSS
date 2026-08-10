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
