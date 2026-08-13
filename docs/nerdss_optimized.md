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

No input under `sample_inputs/` contains a *reversible* bimolecular state change,
so this path has its own regression model at
[`benchmarks/nerdss_optimized/regression/bimol_state_change`](../benchmarks/nerdss_optimized/regression/bimol_state_change).

> **Correction.** This section previously said no input contains a bimolecular
> state change at all. That is wrong: `sample_inputs/enzyme/parms_clat_enzyme.inp`
> line 116, `syn(pi) + pip2(head~U) -> syn(pi) + pip2(head~P)`, parses as one, as
> the build's own reaction dump confirms. It is irreversible, so
> `conjBackRxnIndex == -1`: the reverse branch is gated off by the corrected
> `conjBackRxnIndex != -1` test, and the `perform_bimolecular_state_change_*`
> change only applies when `isStateChangeBackRxn` is true, which that reaction can
> never set. So the conclusion that these corrections are inert for existing
> models still holds; only the justification given for it was wrong. `enzyme` is
> in neither `cases.tsv` nor `known_broken.tsv`, so the forward bimolecular
> state-change path it does exercise is untested by the suite - which is how the
> mistake survived. It is byte-identical over 20,000 iterations across the builds
> compared here.
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

> `class_Coord.hpp` has since been merged with `class_Vector.hpp` into
> [`class_Vec3D.hpp`](../include/classes/class_Vec3D.hpp); see
> [the `Vec3D` merge](#coord-and-vector-merged-into-vec3d) below. Everything in
> this section still holds, under the new names.

`class_Coord.hpp` is now the definition
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

### Reverse unimolecular state change

This one is not part of issues #8-#12. It was recorded as a known defect and left
alone when they were done, on the assumption that the path was reachable from
existing models and that touching it could move results. Establishing that it is
*not* reachable is what made it safe to fix, and that is the first result below.

`find_which_state_change_reaction()` has the same forward/reverse shape as
`find_which_reaction()`, and its reverse branch had the unimolecular versions of
the same defects:

1. `rateIndex` was computed correctly inside the rate loop and then overwritten
   after it with `backRxns[...].relRxnIndex` - a reaction index written into a
   rate index - while `rxnIndex` was never assigned on that path at all. Both
   callers require `rxnIndex != -1`, so every reverse unimolecular state change
   was discarded before it could fire; had one survived, it would have indexed a
   rate list with a reaction index. Both the `matches > 1` and `matches == 1`
   branches did this, and the two were otherwise identical, so they are now one
   branch.
2. `backRxns[oneRxn.conjBackRxnIndex]` was read without testing the `-1`
   sentinel. The branch is entered on a product-state match, which does not imply
   the reaction is reversible, so an irreversible state change whose product
   state matched read out of bounds. The gate is now `conjBackRxnIndex != -1`,
   the same correction made to the bimolecular gate.

**`rxnIndex` indexes `backRxns` here, unlike the bimolecular case.** The two
paths have genuinely different conventions, so this was established from the
callers rather than assumed. `check_for_unimolecular_reactions.cpp` and
`check_for_unimolstatechange_reactions.cpp` both read
`backRxns[rxnIndex].rateList[rateIndex]`, `backRxns[rxnIndex].productListNew[0]`
and `backRxns[rxnIndex].observeLabel` when `isStateChangeBackRxn` is set, and
nothing on this path tests `forwardRxns[rxnIndex].rxnType`. The bimolecular path
needs the opposite because `perform_bimolecular_reactions()` selects the
state-change case by testing `forwardRxns[rxnIndex[0]].rxnType`, which is why
`perform_bimolecular_state_change_*` was corrected to resolve the back reaction
through `forwardRxns[rxnIndex].conjBackRxnIndex` instead. `relRxnIndex` on the
`BackRxn` equals `conjBackRxnIndex` by construction, so the value the original
code computed was right; only the variable it was assigned to was wrong.

**A third defect, in the parser.** `class_MolTemplate.hpp` documents
`stateChangeRxns` as `(forward, back)` and `find_which_state_change_reaction()`
indexes `forwardRxns[rxnItr.first]` on that basis, but
`populate_reaction_lists.cpp` stored `(back, forward)` for the state on the
*product* side, in all six places it builds the list. With a single reversible
pair the forward and back indices are both 0 and the two orderings coincide,
which is why this never surfaced. With the index spaces offset it silently
selects an unrelated reaction. Fixed at the population site so the documented
invariant holds; the direction is recovered where it always was, by testing the
molecule's state against the reactant and the product. That also makes the two
branches identical, so each is now a single condition.

**A fourth defect blocked the model from running at all.**
`BackRxn::display()` read `rate.otherIfaceLists[1]` whenever the list was
non-empty. `otherIfaceLists` holds one entry per reactant, so the conjugate of a
unimolecular state change has one entry and element 1 is past the end. It
crashed four runs in five on the same seed, depending on heap layout.
`ForwardRxn::display()` avoids this by skipping the block for
`uniMolStateChange`; the back version now iterates the list it actually has, so
its output is unchanged wherever the old code was well defined. Display-only, and
the bitwise suite confirms it.

**No existing model is affected, and this was checked rather than argued.** No
input under `sample_inputs/` declares a unimolecular state change with `<->`:
`michaelis_menten`, `unimolecular_reverse` and `auto_phos` write theirs as
separate `->` reactions, which leaves `conjBackRxnIndex` at `-1` so the reverse
branch is never entered. All 13 suite cases stay byte-identical, including
`michaelis_menten`. Two models that use states but are not in the suite were
checked separately and are also byte-identical: `auto_phos` over 200,000
iterations and `enzyme` over 20,000.

Because nothing existing reaches the path, it has its own model at
[`sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible`](../sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible),
with a pass/fail check at
[`benchmarks/nerdss_optimized/regression/unimol_state_change_reversible`](../benchmarks/nerdss_optimized/regression/unimol_state_change_reversible).
It deliberately offsets the forward and back index spaces, which a single
reversible pair cannot do. Unlike the bimolecular model, this path does work end
to end once corrected: both state changes reach `kf/(kf+kb)` within statistical
error, against a factor-of-four different answer before. See section 8 of
`RESULTS.md`.

### `Coord` and `Vector` merged into `Vec3D`

`Coord` and `Vector` were the same three doubles twice over: `Vector` derived
from `Coord` and added a cached `magnitude`. They are now one type,
[`class_Vec3D.hpp`](../include/classes/class_Vec3D.hpp), and
`class_Coord.{hpp,cpp}` and `class_Vector.{hpp,cpp}` are gone.

Having two types meant the type of an expression depended on which header its
operands came from, and the two halves of the API had drifted apart. Both
`operator+(Vector, Coord)` overloads existed, one returning `Vector` and one
returning `Coord`, and which one a call site got was decided by whether its left
operand happened to be const. Four spellings covered two length operations.
`Vector::cross()` normalized its result, which a cross product does not do.
`vector_projection()` returned the rejection, which is the other half of the
decomposition. And reaching `dot`, `cross` or `normalize` from a `Coord`
required a conversion, so the association code is full of `Vector { a - b }`
round trips that existed only to change the type name.

The standardized API is `length()` / `length_squared()`, `dot()`, `cross()` (the
cross product) and `unit_cross()` (normalized, which is what `Vector::cross()`
did), `normalize()` / `normalized()`, `angle_between()`, `rejection_from()`,
`zero()`, and the full operator set - `+ - * /` with scalars on either side,
compound assignment, unary minus, rounded `==` / `!=`.

**The cached magnitude is gone, and that is the part that needed care.** It was
written by `calc_magnitude()` and maintained by nothing: no operator that
changed x, y or z updated it. Almost every read sat one or two lines below the
`calc_magnitude()` that served it, and `length()` gives those the same bits.
Four call sites depended on the cache holding something *other* than the current
length - an unmeasured zero, or a pre-rotation length on a vector since rotated
- and each of them decided a real branch, including one whole coordinate
transform that `check_bases.cpp` has always skipped without meaning to. Those
four now pass the length they mean as an argument to `angle_between()`, so the
dependency is in the call rather than in an object's history. None of them was
"corrected", because correcting any of them would change results; see section
14.2 of `RESULTS.md` for what each one was.

`Vector` was 32 bytes and `Vec3D` is 24; `Coord` was already 24, so the
coordinate arrays that dominate a real run are unchanged in size and only
transient vectors shrink. `normalize()` drops one of its two square roots - the
second only refreshed the cache.

All 13 suite cases stay byte-identical. Whole-simulation timing is unchanged
(0.999x, inside the run-to-run spread); the operations themselves are 5-11%
faster in isolation, which is a tenth off a fraction of a percent of runtime and
so does not surface at the suite level. Sections 14.3-14.5 of `RESULTS.md`.

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
- `benchmarks/vec3d_benchmark.cpp` measures the vector operations the `Vec3D`
  merge changed against the pre-merge `Coord`/`Vector`, and checks that they
  agree bit for bit. Build it with
  `g++ -O3 -std=c++11 -Iinclude $(gsl-config --cflags) src/classes/class_Vec3D.cpp benchmarks/vec3d_benchmark.cpp -o vec3d_benchmark`.
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
  result-preserving subset (issues #8 optimization, #9, #11), and again across
  the `Vec3D` merge.
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
- **`Vec3D` merge:** 13 of 13 byte-identical, runtime 0.999x over the suite -
  neutral, and within the run-to-run spread. In isolation `normalize()` is
  1.096x, the angle between two vectors 1.109x, and streaming an array of
  vectors 1.08x, from one square root removed and 32 bytes down to 24.

## Defects found but deliberately not changed

Each of these is a live or pre-existing defect discovered while doing the work
above. None is part of issues #8-#12, and every one needs its own model and its
own validation, so all are left untouched here.

- ~~`find_which_state_change_reaction()` (unimolecular state change) assigns
  `rateIndex = backRxns[...].relRxnIndex` on its reverse path.~~ **Fixed** - see
  [Reverse unimolecular state change](#reverse-unimolecular-state-change) above.
  The reason it was deferred, that the path is reachable from existing models,
  turned out not to hold: no input under `sample_inputs/` declares a unimolecular
  state change with `<->`, so the reverse branch is unreachable and all 13 suite
  cases stay byte-identical.
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

## Candidate corrections and optimizations

Filed issues that are not part of #8-#12 and have not been worked yet. Each entry
records what was established by reading the code, so the work can be scoped
without rediscovering it. Nothing here is implemented.

### Suite coverage gap - `sample_inputs/enzyme` (validation)

`enzyme` is in neither `cases.tsv` nor `known_broken.tsv`, so it is neither
validated nor recorded as broken. It is the only sample input holding a
bimolecular state change (irreversible), so it is the only one exercising the
forward `perform_bimolecular_state_change_*` path. It runs to completion for
20,000 iterations and is byte-identical across the builds compared here, so the
cheap fix is to add it to `cases.tsv` with a calibrated `nItr`. Doing so would
have caught the incorrect claim corrected above at the time it was written.

### Issue #4 - false out-of-bounds abort (correction)

[Issue #4](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/4) reports a run
aborting with `Cannot fit complex 1 into simulation volume. Exiting...` for a
molecule at `[-63.8805, 81.449, -250]` in a `[166.9, 166.9, 500.0]` box, which is
inside the box. The reporter notes it is intermittent across replicates.

The message the reporter quotes is the one without a named dimension, so the
branch that fired is the `currBin` test at
[`class_SimulVolume.cpp:416`](../src/classes/class_SimulVolume.cpp), not one of
the three coordinate tests. That is the whole problem, and it does not need the
model to see:

- `SimulVolume::update_memberMolLists()` tests the three dimensions first and the
  bin index last, as an `if / else if` chain. So the `currBin` branch is reached
  **only when all three coordinate tests have already passed**, that is, only for
  a molecule the same function considers inside the box.
- All four branches then call the same corrector,
  `Complex::put_back_into_SimulVolume()`, whose translation is derived purely
  from the three coordinate tests: `transVec.x` is nonzero only if
  `x > waterBox.x/2` or `x < -waterBox.x/2`, and likewise for y and z. Those are
  the conditions the caller just ruled out.
- So on the `currBin` path the translation is exactly `(0, 0, 0)`. The complex
  does not move, `update_properties()` recomputes the same bin, the caller
  restarts its scan, and this repeats until the corrector's own counter reaches
  1000 and calls `exit(1)`. The abort is guaranteed, not probabilistic; what is
  intermittent is only whether a molecule ever lands in that state.

`z = -250` is exactly `-waterBox.z/2`, the lower face, and the reporter's model
pins lipids there (`isLipid = true`, `D = [1.0, 1.0, 0.0]`), which is consistent
with a boundary-coincidence case rather than a genuine escape.

Two further discrepancies in the same code are worth folding into the fix, since
they concern the same predicate:

- The coordinate tests allow `1E-6` of slack on the lower face
  (`comCoord.z + 1E-6 < -(waterBox.z / 2)`) but none on the upper face
  (`comCoord.z > (waterBox.z / 2)`), and the binning adds its own `1E-6` only in
  z. The corrector's `else if (zDiff < -waterBox.z)` is strict, so a coordinate
  sitting exactly on a face is treated as inside by one predicate and outside by
  another depending on which one asks.
- The bin guard is `currBin > numSubCells.tot`, so `currBin == numSubCells.tot`
  passes the guard and is then used to index `subCellList`, one past the end. The
  `xItr`/`yItr`/`zItr` clamps just above only repair the exact values `-1` and
  `numSubCells.<dim>`, not larger overshoots.

Scoping note: a fix has to make the four out-of-bounds predicates and the
corrector agree on a single convention, and needs a model that puts a molecule
exactly on a face. Because it changes what happens on the boundary, it is a
result-changing correction and needs the statistical comparison rather than the
bitwise one.

### Issue #7 - `Quat` const-correctness and batch rotation (optimization)

[Issue #7](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/7) asks for
`const`/`noexcept`/`[[nodiscard]]` on the read-only `Quat` operations, separation
of copy-returning from in-place forms, conventional compound operators, explicit
zero-quaternion handling, and - the one change with a measurable payoff - hoisting
the inverse out of loops that rotate many vectors with one quaternion, optionally
behind a prepared-rotation type.

Two things to weigh before starting:

- The API cleanups touch every call site of `Quat` and are mostly mechanical, but
  `CMakeLists.txt` sets `CMAKE_CXX_STANDARD 11` and the Makefile passes
  `-std=c++0x`, so `[[nodiscard]]` needs a standard bump or has to be dropped from
  the scope.
- Issue #9, already done, is a symptom of exactly what issue #7 asks for. The
  `rotQuat.unit()` call it removed was a no-op because `Quat::unit()` returns the
  normalized quaternion instead of normalizing in place - the confusion between
  copy-returning and in-place forms that issue #7's `normalized()` / `normalize()`
  split is meant to prevent. So the API change has a demonstrated defect behind
  it, not just style.
- The inverse-hoist has a concrete target. `Quat::rotate()` builds
  `this->inverse()` on every call
  ([`class_Quat.cpp:57`](../src/classes/class_Quat.cpp)), and it is called once
  per interface inside loops that reuse a single `rotQuat`: `rotate.cpp:18`,
  `class_Molecule_Complex.cpp:409`, `:540` and `:1087`. So the redundant work is
  proportional to interfaces per complex, which is also the axis along which
  issue #8 produced its largest wins.
- It should still be measured before being assumed worthwhile. The per-case
  speedups this branch reports come from reaction matching, not from rotation, and
  the profiling behind issues #8-#12 did not isolate the rotation paths.

Both parts are expected to be bit-for-bit result-preserving, so the 13-case
bitwise suite is the right check, which makes this the cheaper of the two to
validate.
