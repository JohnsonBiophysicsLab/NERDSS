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
only at `membraneObject.isSphere`."* That divergence was invisible while the six
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

Store the answer on the object instead:

```cpp
// on Complex, assigned once per step where OnSurface / onFiber are already maintained
enum class Dim : uint8_t { Fiber1D = 1, Surface2D = 2, Bulk3D = 3 };
Dim dim;
```

The pair rule and the diffusion sum then each have one home, and the ladder
becomes a table:

```cpp
const Dim d = std::min(com1.dim, com2.dim);              // the pair's dimensionality
const double Dtot = weighted_D_sum(com1.D, com2.D, d);   // one rule, one place

kProbabilityKernel[int(d)](rxnIndex, rateIndex, biMolData, ctx);
```

`ctx` carries the extras only the 2D kernel needs - `tableIDs`, `DDTableIndex`,
`normMatrices`, `survMatrices`, `pirMatrices` - so all three kernels share a
signature.

This is a correctness change as much as a brevity one. The `Dtot` averaging rule
is presently written out four or five times, and the commented-out
`std::abs(D.z) < 1E-16` predecessors still sitting beside the live `OnSurface`
tests are evidence that it has already drifted once.

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
| 2 | `Complex::dim` + a single `weighted_D_sum` | low | yes, if FP order is kept |
| 3 | `Boundary` constraint set replacing the six dispatchers | medium | **no** - see below |
| 4 | `SurfaceField` for implicit lipid and compartment | medium | yes |
| 5 | Participant lists (`activeMols` and friends) | very low | yes |
| 6 | Merge `associate_box` / `associate_sphere` | high | to be determined |

**Step 1** is mechanical and eliminates the impossible states. Do it alone, so
the diff stays reviewable.

**Step 2** removes the duplicated ladders. It is bitwise-identical only if the
arithmetic order is preserved: `1/2*(a+b) + 1/2*(c+d)` is not `(a+b+c+d)/2` in
floating point, and the existing expressions must be transcribed exactly.

**Step 3** is where the compartment inconsistency quoted above actually gets
fixed, which means a deliberate behavior change. Before writing the constraint
set, decide which of the six reflection routines was right about the compartment
and record the decision; the four that ignore it are either a latent bug or an
intentional optimisation, and the code does not currently say which.

**Step 6 is last on purpose.** `associate_box.cpp` (1073 lines) and
`associate_sphere.cpp` (598) are the least mergeable pair in the codebase -
genuinely different geometry throughout, not a shared body with a swapped
distance. Attempting them first would produce one large, risky diff. After steps
1 through 5 the shared primitives exist to merge *into*, and the remaining
difference is small enough to see.

Steps 1, 2 and 5 carry almost no risk and deliver most of the readability. Steps
3 and 4 are where the architecture actually changes.

## What this does not address

MPI. [`nerdss_mpi.cpp`](../EXEs/nerdss_mpi.cpp) is a second 1058-line main with
its own copy of much of this structure. Every step above should be applied to
the shared code it calls rather than to either main, so that the two mains
converge rather than diverge further; but unifying the two mains themselves is a
separate piece of work with its own risks, and is out of scope here.
