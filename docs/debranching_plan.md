# De-branching the NERDSS main loop

NERDSS grew by branching. Each new capability - a spherical boundary, an
implicit-lipid membrane, a compartment, diffusion on a fiber - arrived as an
`if` beside the code that came before it, and, where the bodies diverged far
enough, as a second file with `_sphere` or `_implicitlipid` appended to the
name. The result works and is well tested, but the branch count now grows as the
product of the capabilities rather than their sum, and the same decision is
written down in many places, where it drifts.

This document is a plan for removing most of that branching without changing
what the program computes. It is written to be executed in order; each step is
small enough to validate on its own.

## The problem is not the number of `if`s

It is that three unrelated kinds of variation are all expressed the same way, so
they multiply instead of composing.

| Axis | Examples | Decided | Sites |
| --- | --- | --- | --- |
| Static model configuration | `isSphere` / `isBox`, `hasCompartment`, `implicitLipid`, `clusterOverlapCheck` | once, at parse time; never changes again | 65 + 22 + ~40 |
| Dynamic per-complex state | `OnSurface`, `onFiber`, `D.z < 1e-16` - i.e. 1D / 2D / 3D | every timestep | pervasive in `src/reactions` |
| Entity kind | `isImplicitLipid`, `insideCompartment`, `bindToSurface` | per molecule, per type | 251 |

Each axis wants a different remedy, and applying one remedy to all three is why
the combinatorics feel unmanageable. A static configuration flag should not be a
runtime branch at all. A genuinely dynamic per-complex property should be a
value stored on the complex, not a predicate re-derived at each use. An entity
kind that the loop must skip should be filtered out of the loop's input, not
tested inside its body.

### Two symptoms already visible in the code

[`class_Membrane.hpp:134-135`](../include/classes/class_Membrane.hpp) declares
both `isBox` and `isSphere` as independent booleans. That is a two-state choice
encoded in a way that permits two nonsensical states, and nothing in the type
system prevents them.

More seriously,
[`reflect_dispatch.cpp`](../src/boundary_conditions/reflect_dispatch.cpp)
documents the drift in its own header comment: of the six reflection entry
points, *"the two that take a compartment flag consult it, the other four look
only at `membraneObject.isSphere()`."* That divergence was invisible while the six
lived in six files. It is the characteristic failure of branching-by-copy: the
decision is duplicated, so the copies fall out of step, and no single place can
be corrected.

A third instance: `check_bimolecular_reactions.cpp` re-derives the
sphere-versus-box geodesic separation inline in its volume-exclusion path rather
than calling the one already written in
[`get_distance.cpp`](../src/reactions/get_distance.cpp).

## Axis 1 - static geometry becomes a constraint set

The obvious refactor is a `Geometry` base class with a virtual `reflect()`. It
does not work directly, because the box and sphere routines are not the same
algorithm with a different distance function.
[`reflect_traj_complex_rad_rot_box.cpp`](../src/boundary_conditions/reflect_traj_complex_rad_rot_box.cpp)
loops over three axes tracking a wall extent on each;
[`reflect_traj_complex_rad_rot_sphere.cpp`](../src/boundary_conditions/reflect_traj_complex_rad_rot_sphere.cpp)
finds a single farthest radial point. Different shapes of computation.

One level up, they are the same. Both routines are:

> For each boundary constraint, find the complex's worst violation. Reflect by
> twice the overshoot along the outward normal. If that correction pushes the
> complex out through an opposing constraint, resample.

Under that description a box is six planar constraints, a sphere is one radial
constraint violated when `|r| > R`, and a compartment is one radial constraint
with the sign reversed.

```cpp
struct Constraint {                          // signed gap: > 0 inside, < 0 violated
    virtual double gap(const Vec3D& p) const = 0;
    virtual Vec3D  outwardNormal(const Vec3D& p) const = 0;
};
struct PlaneConstraint  : Constraint { int axis; double wall; double sign; };
struct RadialConstraint : Constraint { double R; double sign; };  // +1 inside, -1 outside

class Boundary {
    std::vector<std::unique_ptr<Constraint>> constraints;  // built once, from parsed input
    double surfaceOffset(const Complex&, double RS3Dinput) const;
};
```

The payoff is that **`hasCompartment` stops being a branch**. A compartment is
one more entry in `constraints`. The six reflection functions with inconsistent
compartment handling collapse into a single loop that cannot forget the
compartment, because there is no longer anything to forget. The same holds for
`check_if_spans`, the `sweep_separation_*` family, and the association-time span
checks.

### What survived contact with the code

Most of the above did not, and the part that did paid off more than expected.
Recorded here rather than quietly rewritten, because the reasoning is the useful
part.

**The scan was already extracted.**
[`complex_extent.hpp`](../include/boundary_conditions/complex_extent.hpp)
had already collected the doubly-nested walk over member COMs and interfaces
that all sixteen routines shared. What was left to unify was the *decision*, not
the traversal.

**The box is not a constraint list.** A box reflection treats its three axes
independently: each axis is corrected on its own, and the routine re-samples
only if a correction pushes the complex out through the opposing wall *of that
axis*. Expressing that as a loop over six half-space constraints changes what
"the opposing constraint" means and loses the per-axis pairing. The box
reflectors are left as they are.

**Sphere and compartment are one reflector.** This is where the win actually
was. The three `*_compartment` files were their `*_sphere` twins with a handful
of signs flipped, and every one of those flips is the same flip:

| | sphere (contains) | compartment (excludes) |
| --- | --- | --- |
| boundary radius | `radius - RS3D` | `compartmentR + RS3D` |
| bounding test | `L + r > R` | `L + r < R` |
| point score | `+` | `-` |
| reflection tail | *identical* | *identical* |

All three are `sign * x` for `sign = +1` inside and `-1` outside, so a single
`RadialSide` parameter covers them - including the reflecting-surface offset,
which always shrinks the region the complex may occupy and therefore *adds* to
an excluding radius and *subtracts* from a containing one. The reflection tail
needs no sign at all: `lamda = -2 (targR - R) / targR` is already negative
outside a containing boundary and positive inside an excluding one.

So `hasCompartment` does stop being a separate code path - three files, and the
duplication in them, are gone - even though it remains a branch in the
dispatcher. That is the honest version of the claim above.

**Two scoring conventions, deliberately kept apart.** The radial routines rank
candidate points in one of two ways: by how far past the boundary a point has
got (`|p| - R`, seeded at zero), or by signed radius (`±|p|`, seeded at `±R`).
They agree mathematically and not in floating point - `|p| - R` can round two
distinct radii to the same double and hand a tie to whichever point was scanned
first. Each routine keeps the convention it was written with. Merging them would
be a silent, untestable change in which point wins a tie.

**Two compartment omissions preserved, not fixed.** The compartment reflector
never skipped surface-bound complexes and never re-checked the span after
reflecting ("assume the complex is not huge enough" - unchecked). Both are now
explicit arguments at the call site rather than a missing branch, so the
asymmetry is visible; correcting either changes results and needs its own
commit and its own argument.

The cost is one indirect call per complex per timestep at the reflection level -
not per interface and not per candidate pair. Set against the `sqrt` and
`GaussV()` calls already on that path, it is not measurable.

Geometry does reach two genuinely hot inner loops: the pair-distance calls in
`get_distance` and the exclusion path of `check_bimolecular_reactions`. Those
should stay non-virtual. Give `Boundary` an inlineable `pairDistance()` that
switches on a small enum, or template that one kernel. Do not template the
simulation loop as a whole: the code bloat is not repaid, and the branch is
perfectly predicted anyway, since it takes the same direction on every one of
the billion iterations of a run.

## Axis 2 - dimensionality becomes a value

The 1D / 2D / 3D decision is currently re-derived from `onFiber`, `OnSurface`
and `D.z` at each use site, and each site then rewrites the `Dtot` average by
hand. The same three-way ladder appears twice inside one file -
`check_bimolecular_reactions.cpp:180` for the binding path and again near
`:271` for volume exclusion - with bodies that are similar but not identical.

Give the rule and the sum one home each:

```cpp
enum class Dim : int { Fiber1D = 1, Surface2D = 2, Bulk3D = 3 };

Dim    pair_dim(const Complex&, const Complex&);
double weighted_D_sum(const Vec3D& D1, const Vec3D& D2, Dim);
```

**The pair rule is not `min(dim(c1), dim(c2))`.** An earlier draft of this
document said it was, and that is wrong: `onFiber` and `OnSurface` are set
independently in `Complex::update_properties()`, from `isPromoter` and
`isLipid`, so a fiber complex meeting a membrane complex agrees on neither and
must fall through to 3D. `min` would route that pair to 1D. Only a pair that
*agrees* drops a dimension, and the fiber case is tested first:

```cpp
if (c1.onFiber   && c2.onFiber)   return Dim::Fiber1D;
if (c1.OnSurface && c2.OnSurface) return Dim::Surface2D;
return Dim::Bulk3D;
```

This is a correctness change as much as a brevity one. The `Dtot` averaging rule
is presently written out four times, and the commented-out
`std::abs(D.z) < 1E-16` predecessors still sitting beside the live `OnSurface`
tests are evidence that it has already drifted once.

### The floating-point constraint is FMA contraction, not reassociation

`weighted_D_sum` must reproduce each original expression **as a single
statement**, not as a loop or a running accumulator over the three axes. The
loop form is algebraically identical and numerically is not: at `-O3` the
compiler may fuse `w*x + w*y + w*z` written as one expression differently from
the same terms accumulated across statements. Measured over 4 million random
triples, a loop version disagreed with the original expression on 3.2% of the 3D
cases, in the last ulp. That is enough - `Dtot` feeds `sqrt()` and the 2D
probability tables, so one ulp changes which random draws a run consumes and the
trajectory diverges from there.

Because no input file in the tree exercises the 1D fiber path, this cannot be
validated by running models alone; it needs a direct comparison of the helper
against the expressions it replaced.

## Axis 3 - implicit lipid and compartment become partners

This is the largest single reduction and the least obvious. Today the implicit
lipid is a molecule that every loop must remember to skip: `if
(mol.isImplicitLipid) continue;` appears roughly fifteen times in
[`nerdss.cpp`](../EXEs/nerdss.cpp) alone, alongside parallel
`associate_implicitlipid_*` and `check_implicit_reactions` paths.

Note that `check_implicit_reactions` and `check_compartment_reaction` are
already invoked back to back, at `nerdss.cpp:1035-1057`, with near-identical
signatures. They are one concept: binding to a continuous surface field rather
than to a discrete partner molecule.

```cpp
struct SurfaceField {          // implicit-lipid membrane; compartment inner / outer
    virtual void checkBinding(int molIndex, ..., ReactionContext&) = 0;
    virtual void checkDissociation(int molIndex, ...) = 0;
};
std::vector<std::unique_ptr<SurfaceField>> fields;   // built once, from parsed input
```

The main loop becomes:

```cpp
for (auto& field : fields) field->checkBinding(targMolIndex, ...);
```

An explicit-lipid box system has zero fields; implicit lipid adds one; a
compartment adds another. The `params.implicitLipid == true` guards disappear,
because an empty vector is its own guard.

### Filter the input, do not test in the body

The fastest way to remove a branch from a hot loop is to arrange that the loop
never visits the objects the branch would reject. NERDSS already does this -
`occupiedSubCells` and the SubBox `typeMask` in the main sweep are exactly that
idea. Extend it. Build, once per step for the dynamic sets and once per run for
the static ones:

```cpp
std::vector<int> activeMols;         // not isEmpty, not isImplicitLipid
std::vector<int> surfaceBinders;     // molTemplate.bindToSurface
std::vector<int> compartmentBinders;
```

`for (int m : activeMols)` then replaces `for (all molecules) { if (isEmpty ||
isImplicitLipid) continue; }` at every one of those sites. The branch is gone
outright rather than merely made cheaper, and the loop touches less memory.
For the innermost loops this is strictly better than any form of polymorphism.

## Order of work

NERDSS is a validated scientific code, so the sequence is chosen to keep every
step independently checkable. `run_code_tests/` provides golden trajectories:
pin them with a fixed RNG seed first and make **bit-identical output** the
acceptance criterion for every step below that claims to be result-preserving.
This mirrors the split already used on the `nerdss-optimized` branch, where
result-preserving and stream-changing commits are validated by different means
(see [`nerdss_optimized.md`](nerdss_optimized.md)).

| # | Step | Risk | Result-preserving |
| --- | --- | --- | --- |
| 1 | `isBox` / `isSphere` -> `enum class BoundaryShape` | very low | yes, bitwise |
| 2 | `pair_dim()` + a single `weighted_D_sum()` | low | yes, if the expressions are kept verbatim |
| 3 | Merge the sphere/compartment reflectors behind `RadialSide` | medium | yes, bitwise |
| 4 | `SurfaceField` for implicit lipid and compartment | medium | yes |
| 5 | Participant lists (`activeMols` and friends) | very low | yes |
| 6 | Merge `associate_box` / `associate_sphere` | high | to be determined |

**Step 1** is mechanical and eliminates the impossible states. Do it alone, so
the diff stays reviewable.

**Step 2** removes the duplicated ladders. It is bitwise-identical only if each
expression is transcribed exactly and kept as one statement - see the
FMA note above. Note also that the two volume-exclusion blocks in
`check_bimolecular_reactions.cpp` are not the same ladder: the `pro2` block has
no fiber branch, so a pair on a fiber is measured there as a 3D reaction. That
asymmetry looks like copy-paste drift rather than intent, but correcting it is a
behaviour change and belongs in its own commit.

**Step 3** turned out to be result-preserving after all - see *What survived
contact with the code*. The compartment inconsistency is now visible as
arguments (`skipOnSurface`, `recheckSpan`) instead of as a missing branch, which
is the prerequisite for deciding it, but the decision itself is deferred to its
own commit. Note also that the header comment in `reflect_dispatch.cpp` miscounts
it: three of the six routines consult the compartment, not two.

**Step 6 is last on purpose.** `associate_box.cpp` (1073 lines) and
`associate_sphere.cpp` (598) are the least mergeable pair in the codebase -
genuinely different geometry throughout, not a shared body with a swapped
distance. Attempting them first would produce one large, risky diff. After steps
1 through 5 the shared primitives exist to merge *into*, and the remaining
difference is small enough to see.

Steps 1, 2 and 5 carry almost no risk and deliver most of the readability. Steps
3 and 4 are where the architecture actually changes.

## Open items

Carried forward deliberately, with the reasoning, so they are decisions rather
than omissions.

### The compartment span re-check - a decision for a domain expert

The compartment reflector does not re-check the span after reflecting. Its
comment reads "assume the complex is not huge enough", and that assumption is
unchecked. It is now the explicit `recheckSpan` argument of
`reflect_traj_complex_radial` rather than a branch that is simply absent, so it
can be decided; this document does not decide it, because it is a physics and
numerics question rather than a structural one.

The gap is real but bounded. The compartment reflection runs *after* the outer
box or sphere reflection and pushes the complex radially outward, i.e. toward
the outer wall, and nothing re-checks. What bounds it is that
[`nerdss.cpp`](../EXEs/nerdss.cpp) already refuses to start unless
`waterBox/2 - compartmentR > rMaxLimit`, so there is at least `rMaxLimit` of
clearance, and the reflection displacement is at most about twice the
penetration depth, which is itself bounded by the trajectory step. Escaping the
box therefore needs an unusually large step or an unusually large complex.

Fixing it is not simply "call the span check afterwards": that check tests the
*outer* boundary, so re-running it can push the complex back into the
compartment, and a correct fix needs an iterate-until-consistent loop with a
convergence argument. The spherical reflector already has such a loop, capped at
100 attempts, and it calls `exit(1)` when it fails to converge - which is the
shape of the problem, and the reason this wants a domain judgement about which
constraint should win.

By contrast the compartment reflector's *other* apparent omission - that it does
not skip surface-bound complexes - was checked and is correct, not an
oversight. `OnSurface` means bound to the implicit-lipid membrane, which is a
different surface from the compartment; there is no reason to exempt those
complexes from compartment reflection. It is `skipOnSurface = false` at the call
site and stays that way.

### Two pre-existing crashes found while looking for validation models

Neither is caused by the de-branching work; both were found because this work
needed a compartment model to validate against.

1. **Fixed.** The overlap loop propagated destroyed molecules, so an emptied
   complex reached `create_complex_propagation_vectors()`, which reads
   `memberList[0]` unconditionally. `sample_inputs/compartment/` segfaulted a
   few iterations in.
2. **Open.** Restarting the compartment model segfaults in
   `initialize_paramters_for_implicitlipid_and_compartment_model()`, on the
   unmodified code as well. Not investigated further; it is on the compartment
   restart path, which no test covers.

### Coverage gaps that make some of this unvalidatable by running models

Worth knowing before trusting a bitwise A/B:

* **No fiber model exists.** Nothing under `sample_inputs/` or
  `run_code_tests/` sets `isPromoter`, so `Dim::Fiber1D` is unreachable in every
  shipped model and the 1D paths cannot be exercised. `weighted_D_sum` is
  covered instead by a direct comparison against the expressions it replaced
  (`run_code_tests/weighted_D_sum_check.cpp`).
* **One compartment model**, and it needed a crash fix before it would run.
* **No model reaches the compartment restart path** at all, per crash 2 above.

## What this does not address

MPI. [`nerdss_mpi.cpp`](../EXEs/nerdss_mpi.cpp) is a second 1058-line main with
its own copy of much of this structure. Every step above should be applied to
the shared code it calls rather than to either main, so that the two mains
converge rather than diverge further; but unifying the two mains themselves is a
separate piece of work with its own risks, and is out of scope here.
