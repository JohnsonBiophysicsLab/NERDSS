# The `nerdss-optimized` branch

Branched from `master` at `260f6e2`. It implements the five profiling findings
filed as issues #8-#12 and keeps input files, output files and command-line
behavior unchanged, so existing models run without modification.

The changes fall into two groups, and they need different kinds of testing:

| Group | Issues | Effect on results |
| --- | --- | --- |
| Result-preserving | #8 (optimization part), #9, #11 | Bit-for-bit identical output |
| Stream-changing | #8 (correctness part), #10, #12 | Output differs by construction; validated statistically |

Issues #10 and #12 replace random samplers, so every subsequent random draw
changes and bitwise comparison stops being a meaningful test. That is why the
work is split into two builds during validation, and why the results below report
bitwise identity for one group and statistical agreement for the other.

## What changed

### Issue #8 - `matchList` in `find_which_reaction.cpp`

`find_which_reaction()` and `find_reaction_rate_state()` each built a
`std::vector<std::array<int, 2>>` of every matching rate state, then rescanned it
to pick the match with the most required ancillary interfaces. Both functions run
once per candidate interface pair per timestep, so that vector cost a heap
allocation on the hottest path in the program to hold something a running best
already carries.

Both now use one shared streaming selector, `best_matching_rate()` in
[`shared_reaction_functions.hpp`](../include/reactions/shared_reaction_functions.hpp).
The selection comparison is strict (`>`), which preserves the original
first-match-on-tie behavior: the old loop seeded its best index to 0 and only
replaced it on a strictly larger ancillary count. Removing the vector also
removes the pointer arithmetic that recovered rate indices from element
addresses.

The forward and symmetric paths are exactly result-preserving, including the
`totMatches` counter, which is still incremented on precisely the same four
return paths.

**Reverse bimolecular state-change branch.** The issue lists five defects here
and recommends handling them separately, because they cannot be preserved: each
one leaves the branch dead or undefined.

1. The gate read `conjBackRxnIndex > 0`. Zero is a valid index and `-1` is the
   sentinel, so a model whose conjugate back reaction sits at `backRxns[0]`
   never entered the branch at all.
2. The rate index was recovered as `&oneRate - &oneRxn.rateList[0]` where
   `oneRate` came from `backRxns[...].rateList`. Subtracting pointers into
   unrelated containers is undefined behavior.
3. The single-match case returned without assigning `rateIndex`, and every
   caller requires `rateIndex != -1`, so the reaction was discarded.
4. The multiple-match case returned without selecting a match at all.
5. The zero-match case read element 0 of an empty match list.

All five are fixed. Fixing them exposed a sixth defect, in the callers rather
than in `find_which_reaction()`: `perform_bimolecular_reactions()` selects the
state-change path by testing `forwardRxns[rxnIndex[0]].rxnType`, so `rxnIndex[0]`
is an index into `forwardRxns`, but `perform_bimolecular_state_change_box()`,
`..._sphere()` and the two implicit-lipid variants used that same value to index
`backRxns`. That only works when both lists happen to have one entry. Those four
files now resolve the conjugate back reaction through the forward reaction
(`backRxns[forwardRxns[rxnIndex].conjBackRxnIndex]`), which is an identity in
every model where the path worked before.

No input under `sample_inputs/` contains a bimolecular state change, so this path
has its own regression model at
[`benchmarks/nerdss_optimized/regression/bimol_state_change`](../benchmarks/nerdss_optimized/regression/bimol_state_change).
That model showed the path cannot be exercised end to end on this codebase: the
machinery around it fails on `master` for three further reasons, all outside the
scope of issues #8-#12 and all left alone here. See the regression model's README
and section 6 of `RESULTS.md`. The result-preserving build fails on that model at
exactly the same molecule and iteration as `master`, which confirms these
corrections are inert for anything that runs today.

### Issue #8 - `hasIntangibles.cpp`

`hasIntangibles()` scanned the whole interface list of a molecule for each
ancillary-interface requirement, so it cost
O(requirements x interfaces). `RxnIface::relIfaceIndex` is by definition the
position of the interface in `Molecule::interfaceList` - every site that
populates a molecule assigns `interfaceList[i].relIndex = i` - so the
requirement is now looked up directly, giving O(requirements).

The lookup is guarded twice so behavior is preserved even for a malformed model:
an out-of-range index simply fails to match, and if a molecule ever violates the
`relIndex == position` invariant the original linear scan runs instead.

### Issue #8 - redundant checks audited in `check_bimolecular_reactions.cpp`

Three redundancies, all result-preserving:

- `canInteract` scanned `molTemplateList[pro1MolType].rxnPartners` and only
  afterwards discarded the pair if `pro2` was an implicit lipid. The
  implicit-lipid test is a single flag read that rules the pair out on its own,
  so it now runs first and the partner scan is skipped entirely for those pairs.
- The already-bound test scanned both molecules' `bndpartner` lists before
  combining the two results with `&&`. The second scan now only runs if the
  first found a match.
- Each molecule's `myComIndex` was re-read from `moleculeList` about twenty
  times per candidate pair. Because the surrounding calls take `moleculeList` by
  non-const reference, the compiler cannot prove the value is stable and must
  reload it every time. It is now read once, as a commented-out attempt in the
  original code already suggested, along with references to the two complexes.
- One dead local (`molTypeIndex1` in the exclude-volume block) removed.

### Issue #9 - redundant `rotQuat.unit()` in `class_Molecule_Complex.cpp`

In `Complex::propagate()` the rotation quaternion is built from half-angle
products of the ZYX Euler rotation, which is a unit quaternion up to double
rounding. The `rotQuat.unit()` call that followed was also a no-op, because
`Quat::unit()` returns the normalized quaternion rather than normalizing in
place. Removed; the two remaining `rotQuat = rotQuat.unit()` sites in that file
were assignments, and are handled by issue #10 instead.

### Issue #10 - biased random quaternion generation

Random molecular orientations were drawn as

```cpp
Quat rotQuat{rand_gsl() * 2 - 1, rand_gsl() * 2 - 1, rand_gsl() * 2 - 1, rand_gsl() * 2 - 1};
rotQuat = rotQuat.unit();
```

which samples uniformly inside a 4-cube and then projects onto the unit
3-sphere. Directions towards the corners of the cube therefore receive more
probability mass than directions towards its face centers, and the resulting
rotations are not uniform over orientations.

`rand_unit_quat()` in [`class_Quat.cpp`](../src/classes/class_Quat.cpp) replaces
it with Shoemake's subgroup algorithm: exactly uniform on the 3-sphere, three
uniform variates, no rejection, and unit by construction because
`r1^2 + r2^2 = (1 - u1) + u1 = 1`. Both call sites in
`Molecule::create_random_coords()` use it, and neither needs the normalization
step any more.

This affects the initial orientation of every molecule, so it changes initial
conditions and therefore trajectories. It does not change the dynamics.

### Issue #11 - `Coord` class

[`class_Coord.hpp`](../include/classes/class_Coord.hpp) is now the definition
site for the arithmetic, which was the substantive change: `operator+`,
`operator-`, `operator+=`, `zero_crds()`, `isOutOfBox()`, `get_magnitude()`,
`round()` and `roundv()` all lived in `class_Coord.cpp`, so every use inside the
propagation, reflection and association loops was an out-of-line call that the
compiler could not inline across translation units. The build uses no LTO, so
none of them could be inlined at all.

Alongside that:

- read-only operands taken by `const` reference, including `operator-=`, which
  took `Coord&` and so could not accept a temporary;
- compound assignments return `Coord&`. The free
  `operator+=(Coord&, const std::array<double, 3>&)` did not even modify its
  operand - it returned a new `Coord` and discarded the result, so `c += arr`
  silently did nothing. It now assigns. There were no call sites, so nothing
  observable changes;
- `operator/=` took `double&`, which rejected literals; it now takes `double`;
- `const` and `noexcept` where applicable, so `get_magnitude()` and
  `isOutOfBox()` can be called on a `const Coord`;
- `magnitude_squared()` added, and `is_co_linear()` rewritten around it. Heron's
  formula needed three square roots for the side lengths and is badly
  conditioned for exactly the sliver triangles the test cares about. The cross
  product gives twice the triangle area directly, so comparing squared
  quantities applies the same 1E-8 area threshold with no square root and better
  accuracy near degeneracy;
- the `std::vector<double>` constructor validated before indexing. It used to
  initialize `x`, `y`, `z` from `vals[0..2]` in its member-initializer list and
  only then check the size, so a short vector was read out of range before the
  check could reject it;
- `serialize()`/`deserialize()` use `memcpy` instead of a typed store through
  `arrayRank + nArrayRank`, which carries no alignment guarantee. The bytes
  written are unchanged, so serialized buffers stay compatible;
- `operator==` combines its three rounded comparisons with `&&` rather than
  multiplying `bool`s, and `operator!=` is `!(a == b)`.

### Issue #12 - `GaussV()`

The Marsaglia polar method drew about 2.55 uniforms per call plus a logarithm, a
square root and a division, and branched unpredictably because it rejects
samples outside the unit disc. `GaussV()` is now
`gsl_ran_gaussian_ziggurat(r, 1.0)`, which answers from a table lookup roughly
98% of the time.

This changes both the values and the number of uniforms consumed per call, so
trajectories differ.

## Reproducing the measurements

```bash
make serial
# output hashes and per-build timings
./benchmarks/nerdss_optimized/run_suite.sh bin/nerdss <label> 3
./benchmarks/nerdss_optimized/compare_suites.sh <baseline-label> <label>
# timings with the builds interleaved, which is what the reported speedups use
./benchmarks/nerdss_optimized/interleaved_timing.sh 3 baseline=<bin> candidate=<bin>
```

- [`cases.tsv`](../benchmarks/nerdss_optimized/cases.tsv) lists the models and
  the iteration count used for each.
- [`known_broken.tsv`](../benchmarks/nerdss_optimized/known_broken.tsv) lists the
  sample inputs excluded because they already fail on `master`.
- [`rng_quality/run.sh`](../benchmarks/nerdss_optimized/rng_quality/run.sh)
  measures the samplers from issues #10 and #12 directly.
- [`statistical_check.sh`](../benchmarks/nerdss_optimized/statistical_check.sh)
  compares seed-averaged copy numbers between two builds.

Results are recorded in
[`benchmarks/nerdss_optimized/RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md).

## Results

Measured on an Apple M5 with Apple clang 21 at `-O3`, serial build, 13 models
from `sample_inputs/VALIDATE_SUITE`. Full tables, iteration counts and
methodology in
[`benchmarks/nerdss_optimized/RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md).

- **Bitwise:** all 13 cases byte-identical between `master` and the
  result-preserving subset (issues #8 optimization, #9, #11).
- **Speed:** 1.091x over the suite for the result-preserving subset alone, 1.168x
  with issues #10 and #12 included. Per case, 1.13x to 1.56x for models with
  multiple interfaces per reaction; parity for the two single-interface models,
  which have nothing for the reaction-matching work to save.
- **Issue #10:** the old orientation sampler gives chi^2/dof of 9,810 against
  uniform with individual bins off by up to 59%; the replacement gives 1.02 and
  0.70%.
- **Issue #12:** both samplers match the standard normal to four digits; the
  ziggurat is 3.19x faster in isolation.
- **Statistical equivalence:** across three models and six seeds each, the
  largest Welch `|z|` between `master` and the full branch is 0.61, so the
  sampler replacements did not move the simulated physics.

## Defects found but deliberately not changed

Each of these is a live or pre-existing defect discovered while doing the work
above. None is part of issues #8-#12, and every one needs its own model and its
own validation, so all are left untouched here.

- `find_which_state_change_reaction()` (unimolecular state change) assigns
  `rateIndex = backRxns[...].relRxnIndex` on its reverse path - a reaction index
  written into a rate index - and never assigns `rxnIndex` there at all. This is
  the unimolecular sibling of the issue #8 reverse-branch defects, but unlike
  that branch it is reachable from existing models.
- `set_rMaxLimit()` inspects only `ReactionType::bimolecular`. A model whose only
  bimolecular reaction is a state change therefore keeps `rMaxLimit == 0`, and
  `SimulVolume::Dimensions` divides the box length by it, producing a negative
  cell count and an endless `CELL PAIR MAX EXCEEDED` rescale that makes no
  progress. `VALIDATE_SUITE/unimolecular_reverse` hits the same fault through the
  same route, having no bimolecular reaction at all.
- Executing a bimolecular state change produces NaN coordinates within a few
  hundred iterations with `nan` association angles, and corrupts the copy
  counters (negative and above-population counts) with concrete angles.
- `VALIDATE_SUITE/create_destroy` segfaults during zeroth-order creation starting
  from zero molecules. `VALIDATE_SUITE/clock_model` segfaults immediately after
  setup, seed-dependently - it completes with seeds 1-5 and crashes with seed
  20260810.

The optimized build reproduces every one of these identically where the random
stream is unchanged.
