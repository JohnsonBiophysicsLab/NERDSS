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

Later work on the branch falls into the same two groups for the same reason. The
[pairwise cell grid](#the-pairwise-cell-grid---how-it-is-sized-and-how-it-is-walked)
changes are six result-preserving commits and three that resize the grid; a
different grid orders the candidate pairs differently, which draws the random
numbers differently, so those three are validated statistically as well.

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

### The pairwise cell grid - how it is sized, and how it is walked

The serial bimolecular search bins every molecule into a uniform Cartesian grid
laid over `Membrane::waterBox` and walks a 13-neighbour half-stencil, so each
pair is offered once. The grid is built once at parse time by
`SimulVolume::create_simulation_volume()` and never rebuilt; every timestep,
`update_memberMolLists()` re-bins all molecules by centre of mass, and the
pairwise loop in [`EXEs/nerdss.cpp`](../EXEs/nerdss.cpp) hands each same-cell
and stencil-neighbour pair to `check_bimolecular_reactions()`.

An audit of how that grid gets its dimensions found three correctness bugs -
one silent one in the sizing itself, one in the MPI variant of the same cap, and
one the branch had introduced in its own bookkeeping - and five ways the search
did more work than it needed to. Nine commits follow, each verified against its
own immediate parent, with the measurements in section 19 of
[`RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md).

#### Motivation

The grid has exactly one physical parameter; everything else that decides its
shape is a cap, a floor or a rescale. The physical one is
`params.rMaxLimit`: the largest, over every bimolecular reaction, of
`3 sqrt(6 Dtot dt) + bindRadius` plus both interface-to-COM arms. Sub-volumes
of that edge, with a stencil reaching one sub-volume in each direction, are
exactly what makes the search complete - a pair further apart than `rMaxLimit`
cannot react, so it need not be offered.

Four rules then overrode that, and none of them consulted the interaction range:

| rule | where | value |
| --- | --- | --- |
| z floor | `Dimensions()` | `z = max(4, floor(Lz / rMaxLimit))` |
| z floor, again | `check_dimensions()` | `z = max(4, floor(Lz / rMaxLimit) / 2)` |
| per-dimension cap | `check_dimensions()` | 30 (serial); total 27000 (MPI) |
| sub-volume budget | `check_dimensions()` | `max(4000, 0.5 * N^2)`, `minCells` 64, `scale` from 2 |

Measured per geometry, the cost of those rules shows up as *over-inclusion* -
candidate pairs the stencil offers per pair actually within `rMaxLimit`. The
floor for any cell list with sub-volumes exactly `rMaxLimit` wide is
`27 / (4 pi / 3)`, about 6.4x. Cases well above it are paying for a guard rail:

| case | grid | edge / rMaxLimit | over-inclusion | what set it |
| --- | --- | --- | --- | --- |
| clathrin | 14^3 | 1.05 | 5.0x | nothing; correct |
| mem_localization | 19x19x30 | 1.01 | 2.9x | nothing; correct |
| rev_3D | 30^3 | 1.88 | 11.8x | the 30 cap |
| rev_3Dto2D | 30^3 | 1.97 | 11.0x | the 30 cap |
| rev_2D | 30x30x4 | 5.75 in-plane, 0.43 in z | 52.3x | the 30 cap and the z floor |
| compartment | 17x17x8 | 2.07, 2.07, 4.41 | 97.6x | the `maxPairs` budget |

Two further findings framed the work. The pairwise loop walked *all* of
`subCellList` however few molecules the system held - and occupancy is low
almost everywhere: 97 of 2744 sub-volumes for clathrin, about 1800 of 27000 for
rev_3D, 448 of 8000 for sphere. And the search was blind to what the two
molecules were: `mem_localization` offered about 165000 pairs per step of which
699, 0.42%, got past the molecule-type test inside
`check_bimolecular_reactions()`, the rest being lipid-lipid pairs no reaction
names.

#### The changes, and whether each is a fix or an optimization

| commit | change | kind |
| --- | --- | --- |
| `ebac7db` | Register created molecules in `occupiedSubCells` | **fix** (branch-introduced) |
| `459bd83` | Size sub-volumes from the interaction range in every dimension | **fix** (pre-existing) |
| `cf646a9` | Visit only the occupied sub-volumes in the pairwise search | optimization |
| `a81c7a6` | Drop the `max(4000, N^2/2)` sub-volume budget | optimization |
| `77daa3a` | Bound sub-volumes by a memory budget, not 30 per dimension | optimization |
| `1a704f1` | Reject out-of-range pairs before the interface search | optimization |
| `51e5d33` | Skip candidate pairs whose molecule types cannot react | optimization |
| `808272b` | Track occupied sub-volumes with a bit per sub-volume | optimization (reimplements `cf646a9`) |
| `c25b0a0` | Stop the MPI sub-volume cap from raising a dimension | **fix** (pre-existing, MPI only) |

**`ebac7db` - duplicate cell membership after any creation reaction.** This one
is a defect the branch introduced, at `0784dd8`. `clear_member_lists()` empties
only the sub-volumes listed in `occupiedSubCells`, but
`create_molecule_and_complex_from_rxn()` binned a new molecule by pushing
straight into `subCellList[currBin].memberMolList`. Creation runs before
`update_memberMolLists()`, so a molecule created into a sub-volume the previous
re-binning pass had left empty was invisible to the next clear: its entry
survived and the re-binning pass added the same molecule a second time. That
sub-volume then stayed off the registry, was never cleared again, and its member
list only grew. `VALIDATE_SUITE/create_destroy` reports untracked non-empty
sub-volumes and duplicate members from `simItr 1`; with the fix it reports none.
`master` is unaffected, because it clears every sub-volume unconditionally.
`add_member()` is now the only way to grow a member list.

**`459bd83` - sub-volumes thinner than the interaction range.** Pre-existing;
the arithmetic is byte-identical on `master`. `Dimensions()` gave z a
`max(4, ...)` floor, so a box of any thickness got at least four sub-volumes
along z. Below `waterBox.z = 4 * rMaxLimit` that made them thinner than the
interaction range, and the plus-or-minus-one stencil stopped spanning it:
molecules two sub-volumes apart in z could be within `rMaxLimit` and were never
offered to `check_bimolecular_reactions()` at all. `check_dimensions()` opened
with a repair for exactly this, which re-applied the same floor and so could not
fix it. Counting every pair within `rMaxLimit` by brute force, on the rev_3D
model with the box thickness varied:

| box z | cells z | edge z | pairs in range | reached | missed |
| --- | --- | --- | --- | --- | --- |
| 30 nm | 4 | 7.50 nm | 29123 | 26951 | 2172 (7.5%) |
| 50 nm | 4 | 12.50 nm | 20475 | 20391 | 84 (0.41%) |
| 100 nm | 5 | 20.00 nm | 11350 | 11350 | 0 |
| 940 nm | 30 | 31.33 nm | 1328 | 1328 | 0 |

Zero missed whenever the edge is at or above `rMaxLimit`, non-zero the moment
the floor forces it below. All three dimensions now come from
`floor(L / rMaxLimit)`, which is the finest grid whose sub-volumes still cover
the range; a box shorter than the range gets one sub-volume, which has no
neighbours and so cannot miss any. `create_simulation_volume()` asserts the
resulting edges, so a later change to the arithmetic cannot bring the missed
reactions back silently. It also rejects `rMaxLimit <= 0` with a message rather
than dividing by it - see the defects list below.

**`cf646a9` and `808272b` - walking only the occupied sub-volumes.** The loop
cost time proportional to the sub-volume count regardless of the molecule count.
Timed on its own, a bare pass over every sub-volume was 2 to 6% of total
runtime. `occupiedSubCells` already listed the non-empty sub-volumes for
`clear_member_lists()`; in ascending order it is exactly the non-empty
subsequence of the old walk, so the pairs come out in the same order and the
random stream is untouched.

The first version got that order by sorting the list, which `add_member()` fills
in molecule order. That sort cost more than the walk it replaced - on rev_3D at
nItr 60000, median of five interleaved repetitions on the same
27000-sub-volume grid, the full walk took 9.127 s and the sorted version 9.335 s,
2.3% slower. `808272b` replaced it with one bit per sub-volume: registration is
a single OR, and the ascending list is one pass over the mask, a word per 64
sub-volumes, plus one `push_back` per occupied sub-volume. Setting a bit is also
idempotent, which is what the sorted version needed a `std::unique` for.

**`a81c7a6` - the `max(4000, N^2/2)` budget.** The budget existed so that
walking the sub-volumes could not cost more than testing all `N^2/2` pairs
outright, and `cf646a9` removed that reason. What it did cost was resolution, on
whichever axes its loop happened to pick: it divided x and y by `scale` but z by
`2 * scale`, so the grid it left behind was not cubic, which the 13-neighbour
stencil assumes. On `compartment`, 100 molecules set the budget at 5000
sub-volumes, the 27000 of the capped grid exceeded it, and the loop settled on
17 x 17 x 8 - sub-volumes of 58.8 x 58.8 x 125 nm against a 28.35 nm interaction
range, an over-inclusion of 98x, the worst of any sample, on the case whose own
molecule count triggered the guard meant to make it faster. `N` also came from
`Molecule::numberOfMolecules` read once at setup, and the grid is never rebuilt,
so a system that grew kept whatever resolution its starting count had bought.

**`77daa3a` - the 30-per-dimension cap.** It took no account of the interaction
range. Candidate pairs grow with the cube of the sub-volume edge, and the
measured over-inclusion tracked it: rev_3D wants 56 sub-volumes per side, got
30, and doubled from about 5x to 11.8x; rev_2D wants 172 per side and reached
52x. Memory is the only thing that still justifies a cap, so the cap is now
stated in sub-volumes and applied to the total: one costs about 145 bytes, so
the 500000 budget is roughly 72 MB of grid, and the largest sample asks for
59^3 = 205379. Over budget, all three axes shrink by the same cube-root factor,
which keeps the sub-volumes cubic and can only lower a count, so it cannot push
an edge back below `rMaxLimit`.

**`1a704f1` - rejecting out-of-range pairs early.** `check_bimolecular_reactions()`
opened the full `freelist x rxnPartners x freelist` search, and
`find_which_reaction()` for every combination, on any pair the cell list handed
it. The cell list is deliberately over-inclusive, so most of what it offers is
out of range: 78% on clathrin, 99% on compartment. Testing the centre-of-mass
separation first costs three subtractions and three multiplies, and rests on the
bound the cell list already assumes,

```
|COM1 - COM2| <= |iface1 - iface2| + r1 + r2 < Rmax + r1 + r2 <= rMaxLimit
```

for any pair close enough to react. Two exemptions, both measured rather than
assumed: the volume-exclusion path is not gated because it carries its own radii
and its 2D branch reaches to `RMax * 10`, and pairs with both complexes on the
surface are not gated because `rMaxLimit` is not in fact a bound for them - see
below.

**`51e5d33` - skipping pairs whose types cannot react.** A bit table settles the
question before the call: bit `t2` of `interactionMask[t1]` is set iff a
`(t1, t2)` pair can do anything inside the function, meaning `t2` appears in
`t1`'s `rxnPartners` or either template carries `excludeVolumeBound`. The second
term is deliberately loose, because the volume-exclusion path turns on an
interface being bound, which is runtime state. Each sub-volume also carries the
set of types it holds, so a whole neighbouring sub-volume is skipped when
nothing in it can pair with the molecule being tested - which is where most of
the saving is, since with 95% of `mem_localization`'s molecules being lipids
most sub-volumes hold nothing else. The table exists only for models with at
most 64 molecule types; past that it is left empty and every pair is offered, as
before.

**`c25b0a0` - the MPI cap could raise a dimension.** The MPI branch of
`check_dimensions()` reacts to a grid over 27000 sub-volumes by assigning y and
z the same per-axis figure, which can raise them: a 200 x 200 x 1 grid comes out
200 x 11 x 11. Raising a count shortens that edge, which is what makes the
stencil stop covering the interaction range - so what used to be silently missed
pairs would now hit the assertion added in `459bd83`. `std::min` keeps the cap
doing only what a cap should. The change is inside `#ifdef mpi_` and a clean
serial rebuild gives the same SHA-256, so the serial results carry over
unchanged; `make mpi` compiles, and this build has no MPI validation beyond
that.

#### Why five cases are not bitwise identical

Six of the nine commits are exactly result-preserving: all 18 benchmark cases,
13 in `cases.tsv` and 5 in `coverage_cases.tsv`, are byte-identical across
them. Three commits change the grid, and a different grid means a different
order of candidate pairs. The order matters because the search fills each
molecule's `crossbase` and `probvec` in the order it offers the pairs, and
`determine_if_reaction_occurs()` then walks that list drawing one
`rand_gsl64()` per entry and returning on the first draw under the entry's
probability. Reorder the list and both the random numbers consumed and
the partner selected change. Bitwise comparison stops being a meaningful test at
that point, exactly as it does for issues #10 and #12 above, and statistical
agreement takes over.

Five cases move, over seven case-and-commit pairs, and they are precisely the
ones whose grid one of those rules was distorting - no others. rev_2D and
compartment each appear twice because two different rules were acting on them:

| case | commit | grid before | grid after | why it moved |
| --- | --- | --- | --- | --- |
| trimer | `459bd83` | 3x3x4, 29.60 nm in z | 3^3, 39.47 nm | z floor forced sub-volumes below the 37.58 nm range |
| rev_2D | `459bd83` | 30x30x4, 2.50 nm in z | 30x30x1, 10 nm | same floor, on a 10 nm slab |
| compartment | `a81c7a6` | 17x17x8 | 30^3 | the only grid the `maxPairs` budget touched |
| rev_3D | `77daa3a` | 30^3 | 56^3 | the 30 cap |
| rev_3Dto2D | `77daa3a` | 30^3 | 59^3 | the 30 cap |
| rev_2D | `77daa3a` | 30x30x1 | 172x172x1 | the 30 cap |
| compartment | `77daa3a` | 30^3 | 35^3 | the 30 cap |

Each lands on the finest grid that still covers its interaction range - edges of
16.79 against 16.70, 16.95 against 16.94, 5.81 against 5.79, 28.57 against
28.35, 39.47 against 37.58 - which is the intended end state, not an accident of
the new arithmetic.

That the physics does not move is expected for a different reason in each case.

- **`77daa3a`, the finer grids.** A finer grid drops only pairs that were
  over-included and would have failed the distance test anyway. Nothing that
  could react stops being offered, so only the order changes. Over 12 seeds,
  plateau-averaged copy numbers agree within seed scatter: largest Welch `|z|`
  1.24 on rev_3Dto2D, 0.58 on rev_3D and rev_2D.
- **`a81c7a6`, compartment.** Same argument, and the case is additionally
  degenerate: it cannot be run past about 8000 iterations (see
  `known_broken.tsv`) and at 4000 its observables barely move. Over 16 seeds the
  largest `|z|` is 1.00.
- **`459bd83`, rev_2D.** Its trajectory does not change at all. Every file under
  `DATA/` and `PDB/` is byte-identical; the restart files differ only in each
  molecule's recorded `mySubVolIndex`, all 1600 of them by the same
  2700 = 3 x 900, the z-row offset. Every molecule sits on the membrane, so both
  grids put all of them in one z-row behind the same in-plane stencil - the
  sub-volume indices are renumbered and nothing else.
- **`459bd83`, trimer.** This is the one case where the physics *should* move,
  because the old grid was missing reactions. It moves in the right direction
  and by too little to resolve; see below.

#### Three things the work turned up

**1. The missed-pair bug is not a thin-box bug.** The z floor was written as a
guard for flat geometries, and the failure mode reads that way - a 30 nm slab
loses 7.5% of its in-range pairs. But the condition that triggers it is
`floor(Lz / rMaxLimit) < 4`, which a perfectly cubic box meets whenever the
interaction range is a large fraction of the box. `trimer` is a cubic 118.41 nm
box against a 37.58 nm range: `floor(118.41 / 37.58) = 3`, so the floor forced
four z sub-volumes of 29.60 nm - thinner than the range, in a cubic box.
Reconstructing both grids over 101 trajectory frames, 300 molecules and 454619
pairs within `rMaxLimit`, the old 3x3x4 grid could not reach 1808 of them,
0.398% +- 0.009%; the new 3x3x3 grid reaches all of them.

That is small enough that the equilibrium comparison does not resolve it. Over
28 seeds at nItr 30000, trimer's total bound pairs move from 46.194 +- 1.250 to
48.025 +- 1.278, `z = +1.02`; the direction is right, all three bound species up
and all six free species down, but the seed-to-seed scatter is 2.7% against a
0.4% mechanism. The +3.96% figure is scatter, not signal, and the pair count is
what carries the claim. Any model whose box is under four interaction ranges on
a side was affected, which is not a rare shape.

**2. `rMaxLimit` does not bound 2D reactions.** `1a704f1` rests on `rMaxLimit`
bounding how far apart two molecules can be and still react - which is the same
thing the cell list assumes when it sizes sub-volumes by it. The first version
of that commit broke rev_2D, so the assumption was measured directly, counting
every reacting pair:

| model | reacting pairs | beyond rMaxLimit | worst |
| --- | --- | --- | --- |
| rev_2D | 20220 | 1550 (7.7%) | 1.109x, a 2D pair |
| rev_3D | 1338617 | 0 | |
| rev_3Dto2D | 728096 | 0 | |
| clathrin | 13362366 | 0 | |
| mem_localization | 1255726 | 0 | |

`set_rMaxLimit()` estimates a 2D reaction's reach as `3 sqrt(6 Dtot dt)`, while
`determine_2D_bimolecular_reaction_probability()` uses `3.5 sqrt(4 Dtot dt)` over
a `Dtot` that `add_2D_rotational_diffusion()` and `discretize_2D_Dtot()` have
both revised. rev_2D is the one model whose `rMaxLimit` is itself set by a 2D
reaction; everywhere else a 3D reaction sets it and covers the 2D one.

`1a704f1` therefore exempts pairs with both complexes on the surface, which is
enough to make that commit result-preserving. The underlying gap is wider and
pre-existing: with sub-volumes 5.81 nm wide and a reach to 6.43 nm, rev_2D's
cell list already cannot guarantee it offers every reacting pair. Closing it
means changing `set_rMaxLimit()`, hence `rMaxLimit`, hence the grid, so it is
left as its own change - see the defects list below.

**3. Two cases get nothing back, for structural reasons.** rev_2D is 2% slower
and sphere is flat.

rev_2D is pure 2D, so `1a704f1`'s filter is exempted for almost all of its
pairs, and its two molecule types do react, so `51e5d33`'s mask never fires. Its
runtime is dominated by the 2D reaction-probability tables rather than by the
search, so cutting candidate pairs does not help it much. What it does get is
the finer grid, whose cost it pays without the benefit. Both of the reasons it
misses out point at the same root cause as finding 2 above.

sphere is flat because none of these nine commits touches what is wrong with it.
The grid is laid over the `(2R)^3` bounding cube, and for surface-bound molecules
only the shell sub-volumes can ever be occupied - a fraction `pi h / 2R` of the
box, which on the `sphere` sample is 1257 reachable sub-volumes of 8000, of
which 448 are occupied. 94.4% of the grid is swept for nothing before `cf646a9`
and skipped for nothing after it, but the resolution that the empty 94.4% costs
is not recovered either way. `get_distance()` also measures surface pairs by
geodesic arc length while the grid bins on Cartesian coordinates; arc is never
shorter than chord, so the grid stays conservative and no pairs are lost, but it
over-includes on small spheres. A surface decomposition for `OnSurface`
complexes - latitude-longitude bands over the shell, with the Cartesian grid
kept for the interior - is
[the next section](#the-spherical-membrane---a-latitude-longitude-index).

#### What it costs and what it buys

`interleaved_timing.sh`, median of 3 repetitions, builds interleaved per case so
both meet the same machine conditions, on an otherwise idle Apple M5:

| case | before | after | |
| --- | --- | --- | --- |
| mem_localization | 2.457 s | 0.784 s | **3.13x** |
| rev_3Dto2D | 5.967 s | 3.437 s | 1.74x |
| michaelis_menten | 4.365 s | 2.597 s | 1.68x |
| clathrin | 3.362 s | 2.082 s | 1.62x |
| closed_homoTrimer | 6.364 s | 4.409 s | 1.44x |
| homoTrimer | 6.389 s | 4.441 s | 1.44x |
| hetTrimer | 6.152 s | 4.488 s | 1.37x |
| trimer | 5.757 s | 4.206 s | 1.37x |
| implicit_lipid | 2.631 s | 2.108 s | 1.25x |
| hexamer | 5.112 s | 4.183 s | 1.22x |
| rev_3D | 3.572 s | 3.386 s | 1.06x |
| sphere | 5.957 s | 5.865 s | 1.02x |
| rev_2D | 41.332 s | 42.105 s | 0.98x |
| **total** | **99.417 s** | **84.091 s** | **1.18x** |

`mem_localization` is the case `51e5d33` was aimed at, and it is the largest
single win. `rev_3D` gains little because its two effects nearly cancel: it gets
2.3x more sub-volumes from `77daa3a`, which cuts its candidate pairs, and pays
for the larger grid in `refresh_occupied_cells()` and in cache footprint.

#### How the measurements were made

The suite numbers come from the committed scripts listed under
[Reproducing the measurements](#reproducing-the-measurements). Three of the
figures above do not, and were taken with temporary instrumentation that was
reverted before each commit; redoing them means re-adding it.

- **Missed pairs against box thickness, and the funnel counts per geometry.**
  Counters in `EXEs/nerdss.cpp` sampling every 50 steps: occupied sub-volumes,
  stencil pair count, and the same pairs filtered to those within `rMaxLimit`;
  plus, behind an environment flag, an `O(N^2)` sweep of every pair within
  `rMaxLimit` regardless of sub-volume. The difference between the stencil count
  and the brute-force count is the missed pairs. The box thickness was varied by
  editing `WaterBox` in the rev_3D model.
- **trimer's 0.398%.** Reconstructed in Python from
  `DATA/trajectory.xyz`, which writes each molecule's COM followed by its
  interfaces, so every third line is a COM. Both grids' indices are recomputed
  from the COMs with the same binning arithmetic as `update_memberMolLists()`,
  and a pair counts as unreachable when its index differs by two or more on any
  axis. This needs no instrumented build, only a run with `trajWrite` turned
  down.
- **Reacting pairs beyond `rMaxLimit`.** A global holding the current pair's COM
  separation, set at the top of `check_bimolecular_reactions()` and tested in
  `get_distance()` on the branch that records a crossing.

Because the first and third need an instrumented binary, the numbers they
produced are recorded in section 19 of `RESULTS.md` rather than being
re-derivable from the tree.

#### Building: `make serial` and `make mpi` share `obj/`

Both targets compile into `obj/`, and `make mpi` adds `-Dmpi_`. Building one
after the other therefore leaves object files compiled for the wrong target,
which `make` will not rebuild because they are newer than their sources. That
matters here specifically: `check_dimensions()` has an `#ifdef mpi_` branch, so
a serial binary linked against an MPI-flavoured `class_SimulVolume.o` sizes its
grid by the MPI rules. Run `make clean` between the two targets. `make clean`
also removes `bin/`, so anything kept there for comparison has to live
elsewhere.

### The spherical membrane - a latitude-longitude index

`ShellIndex` covers surface-to-surface pairs on a spherical system with bands
one cutoff wide in colatitude, each divided into cells one cutoff wide along
its own narrowest parallel. Molecules stay in the Cartesian grid as well, so
surface-to-interior pairs are unaffected; the Cartesian pass skips a pair only
when both of its molecules are in this index, and the second pass offers it
there. Commit `336322d`.

#### What was measured before building it

Section 19's work had already removed most of what the sphere was paying. What
was left, on the two spherical samples:

| sample | pairs offered / step | of which surface-surface |
| --- | --- | --- |
| `sphere`, R=100 | 0 | 0 |
| `gagsphere`, R=70 | 1915 | 70 (3.6%) |

`sphere` offers nothing: its only reaction binds A to an implicit lipid, so the
molecule-type mask discards every pair before the search reaches it. On the
R=70 sphere, 96.4% of what is offered has at least one molecule off the
surface, which a shell index does not touch. And a sampling profile of that
model puts 11655 of 11662 samples inside
`determine_2D_bimolecular_reaction_probability` -> `create_survMatrix` /
`create_pirMatrix` -> `integrator` -> Bessel `j0/y0/j1/y1`. That is 2D
reaction-table construction, it happens only for pairs already within `Rmax`,
and no neighbour structure can reduce it.

So the index was expected to be neutral on every current sample before a line
of it was written, and it is. What it is for is spheres large enough that the
`(2R)^3` bounding cube exceeds the sub-volume budget and the interior starts
taking resolution from the shell - at R around 2000 nm against a 30 nm range,
`(2R/h)^3` passes the 500 000 budget from `77daa3a`.

#### Why the stencil is complete

The sizing is the whole correctness argument, and it comes from the haversine
identity: for two points at colatitudes `tA`, `tB` with longitude difference
`dp` and central angle `g`,

```
sin^2(g/2) = sin^2((tA - tB)/2) + sin(tA) sin(tB) sin^2(dp/2)
```

Both terms on the right are non-negative, so each is separately bounded by the
left. A pair within `g` therefore has `|tA - tB| <= g`, and
`sin(|dp|/2) <= sin(g/2) / sqrt(sin tA sin tB)`. Bands `g` wide in colatitude
put such a pair at most one band apart; cells wide enough to cover the second
bound put it at most one cell apart within a band. That is what makes a
forward-only, plus-or-minus-one stencil complete, and it is why the longitude
count has to be computed per band rather than globally - the bound blows up
near the poles, where the count collapses to one and the band becomes a single
cap.

[`benchmarks/shell_index_test.cpp`](../benchmarks/shell_index_test.cpp) checks
this against an `O(N^2)` sweep over seven radius-and-cutoff combinations, from
a 1000 nm sphere against a 30 nm cutoff to a 30 nm sphere against 19.4 nm:
every pair within the cutoff offered, none offered twice, in all seven.

#### Three ways it stays out of the way

The index is inactive for any non-spherical system, so nothing outside a sphere
changes at all.

It is also inactive when no two explicit molecules can react with each other.
That case is not hypothetical and not free: the first version activated for any
sphere, and `sphere` came out **9% slower while staying bitwise identical** -
448 surface molecules rebinned every step, an `acos` and an `atan2` apiece, to
fill an index whose every pair the type mask then discarded. With the gate,
`sphere` is 0.997x.

And a molecule below `0.9 * sphereR` is left to the Cartesian grid rather than
binned here, because a cell only covers a cutoff of arc for molecules far
enough out. Surface complexes are pinned by the reflectors at `sphereR - RS3D`,
which measures 0.971 and 0.995 of `sphereR` on the two spherical samples, so
the floor admits all of them. Nothing about correctness rests on the fraction:
a molecule below it is handled by the Cartesian grid exactly as it was before
this index existed, so the constant trades one search for another and cannot
lose a pair. That is the property the hyperparameters section 19 removed did
not have.

#### Verdict

17 of the 18 benchmark cases are bitwise identical and time at 0.99x to 1.00x.
`cluster_gagsphere` is the only case where the index activates and the only one
that moves; over 12 seeds its copy numbers agree within seed scatter, largest
`|z|` 0.29, including all three gag-gag surface reactions.

Its own timing is neutral. A single seed reads 0.842x, but that model runs
anywhere from 32 s to 72 s depending on seed - its cost follows how many
distinct 2D tables the trajectory happens to need - and across five seeds the
medians are 44.89 s before and 45.11 s after. A 2.2x spread between seeds
against a 0.5% difference in medians.

### Re-audit - cumulative drift, and three residual defects

Sections above verified each commit against its own parent. Two questions that
does not answer: does the sequence as a whole drift anywhere it was not meant
to, and is anything left in the decomposition worth fixing. Section 21 of
`RESULTS.md` has the tables.

#### Cumulative drift

Running the suite at the tip against the manifests from the build before any of
this work, six of the eighteen cases differ, and they are exactly the union of
the cases each stream-changing commit was expected to move: `trimer` and
`rev_2D` from the z floor, `compartment` from the maxPairs budget, `rev_3D`,
`rev_3Dto2D`, `rev_2D` and `compartment` from the 30 cap, and
`cluster_gagsphere` from the shell index. The other twelve are byte-identical
from the original baseline to the tip. Nothing moved that was not meant to.

One timing in that run needs discounting rather than explaining: `homoTrimer`
recorded 1698.693 s against a 6.060 s baseline while still producing
byte-identical output. Its grid is 10^3 SubBoxes of 32.1 nm against a 29.7 nm
range, which is correct and small, and the same binary, input and seed re-runs
at 3.97, 3.95 and 3.91 s wall against 3.86, 3.88 and 3.85 s user. It is the
host, which the preamble to `RESULTS.md` already records doing this to the same
case. Bitwise verdicts are load-independent and the interleaved timings are
load-tolerant; a single-shot suite timing on this host is not evidence.

#### Three residual defects

All three are bitwise identical across all 18 cases, so the results above carry
over unchanged.

**`3d253cc` - the shell cutoff angle was the small-angle limit of the right
bound.** `ShellIndex` sized its cells from `rMaxLimit / radiusFloor`, which is
the small-angle limit of `2*asin(rMaxLimit / 2*radiusFloor)` and smaller than
it, since `asin(u) >= u`. rMaxLimit is a *chord* bound -- it is built so that
two molecules close enough to react have `|COM1 - COM2| <= rMaxLimit` -- and a
chord subtends its widest angle at the smallest radius, so the radius floor is
the worst case and the exact form is the one the guarantee needs. The old form
was correct only by the margin between where surface complexes actually sit,
0.971 and 0.995 of `sphereR`, and the 0.9 floor. Making it exact costs 0.4% in
cell width on the R=70 sample and leaves the band and cell counts unchanged.

The test now checks the property the simulation relies on, rather than the
index's internal consistency: a pair within `rMaxLimit` of each other by
straight-line distance must be offered, with the points placed at the radius
floor because that is the worst case. It keeps the angular check as a separate
count, so a change that breaks one and not the other is distinguishable.

**`391a1a3` - the sub-volume budget could fail to bind.** It scaled all three
axes by one cube-root factor and clamped each at a minimum of one SubBox, and
those two rules disagree: the factor assumes every axis absorbs its share, and
an axis already at one cannot. A 1 x 1 x 1197604 grid came out of the single
pass at 1 x 1 x 895000, still nearly twice the 500 000 cap. Iterating fixes it
-- each pass either lowers the product or leaves it alone, and the loop stops
on the latter. No sample input comes near the cap, so nothing in the suite
moves; a cap that does not cap is the same class of defect as the guards this
work removed, which is why it is worth the loop anyway.

**`a576661` - the re-bin restart skipped molecule 0.** The checked pass calls
`put_back_into_SimulVolume()` for a molecule found outside the box, empties
every member list and restarts the binning. The restart was `molItr = 0` inside
a `++molItr` loop, so it resumed at molecule 1 -- and since
`clear_member_lists()` had just thrown away molecule 0's membership, nothing
put it back. For that step molecule 0 was absent from the grid, and so was
every pair it belongs to. The counter is now signed and the restarts set it to
-1.

The message that marks this path, "Attempting to fit back into box", appears
zero times across all 18 cases. That is why it survived, and why fixing it
changes nothing measurable. The `MpiContext` overload carries the same four
restarts and is left alone, as with the other MPI-only paths here.

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
- **Pairwise cell grid:** one silent correctness bug fixed - sub-volumes could
  be thinner than the interaction range, losing up to 7.5% of the pairs within
  it in a thin box and 0.40% in a cubic one - plus four sizing rules replaced.
  1.18x over the suite, from 3.13x on `mem_localization` down to 0.98x on
  `rev_2D`; 13 of the 18 benchmark cases byte-identical across the whole
  sequence, and the five that move are exactly the ones whose grid the old rules
  distorted. Largest Welch `|z|` across the stream-changing commits is 1.24.
  Section 19 of `RESULTS.md`.
- **Spherical shell index:** a latitude-longitude index for surface-to-surface
  pairs on a sphere, proved complete against a brute-force sweep over seven
  radius-and-cutoff combinations. 17 of 18 cases byte-identical, the one that
  moves agrees statistically (largest `|z|` 0.29), and the whole thing is
  neutral on every current sample by construction - the measurements that say
  why were taken before it was built. Section 20 of `RESULTS.md`.
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
- ~~`set_rMaxLimit()` inspects only `ReactionType::bimolecular`. A model whose
  only bimolecular reaction is a state change therefore keeps `rMaxLimit == 0`,
  and `SimulVolume::Dimensions` divides the box length by it, producing a
  negative cell count and an endless `CELL PAIR MAX EXCEEDED` rescale that makes
  no progress.~~ **Diagnosed, not repaired** - `create_simulation_volume()` now
  rejects `rMaxLimit <= 0` with a message naming the cause instead of looping
  (`459bd83`), because an assertion on the sub-volume edges is worth little
  while the divisor can be zero. The underlying gap is unchanged: such a model
  still cannot run, it just says so. `VALIDATE_SUITE/unimolecular_reverse` hits
  it through the same route, having no bimolecular reaction at all. The obvious
  repair is to guard on the value rather than on `forwardRxns` being empty,
  which is the condition the existing 40 nm fallback tests.
- `set_rMaxLimit()` also under-estimates the reach of 2D reactions, so
  `rMaxLimit` is not the upper bound the cell list assumes it is. It estimates
  `3 sqrt(6 Dtot dt)` where
  `determine_2D_bimolecular_reaction_probability()` uses
  `3.5 sqrt(4 Dtot dt)` over a `Dtot` that `add_2D_rotational_diffusion()` and
  `discretize_2D_Dtot()` have both revised. Measured on
  `VALIDATE_SUITE/bimolecular_reversible/rev_2D`, 1550 of 20220 reacting pairs
  (7.7%) lie beyond `rMaxLimit`, by up to 1.109x; on four models whose
  `rMaxLimit` is set by a 3D reaction, none of 16.7 million does. Only models
  whose `rMaxLimit` is itself set by a 2D reaction are exposed, and for those
  the cell list can silently skip reacting pairs - rev_2D runs with 5.81 nm
  sub-volumes against a reach of about 6.43 nm. Left here because repairing it
  changes `rMaxLimit`, hence the grid, hence the random stream for every 2D
  model, which needs its own change and its own statistical validation.
  `1a704f1` works around it by exempting surface-surface pairs from its
  centre-of-mass filter; that exemption can go once the bound is real.
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
