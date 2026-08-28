# Two defects in `Molecule`'s bond bookkeeping

Found while investigating whether `bndlist`, `bndpartner` and `bndRxnList` are
index-parallel arrays that could be merged the way section 23 and 24 of
[`RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md) merged the reweighting
and crossing lists. They are not, and the reasons are two separate defects.

Both **predate this branch**: `break_interaction.cpp`,
`break_interaction_implicitlipid.cpp`, `correct_structutre.cpp` and
`src/mpi/delete_disappeared.cpp` are untouched by it, as are the `bnd*` pushes
in `associate_box.cpp`. Neither is caused by the layout work; both were found by
it.

Measured 2026-08-28, Apple M5, Apple clang 21.0.0, `-O3 -std=c++0x`, seed
20260810, over all 13 cases in `cases.tsv` and all 5 in `coverage_cases.tsv`.

## Summary

| | defect 1 | defect 2 |
| --- | --- | --- |
| what | `bndRxnList` is never erased and has no live reader | `bndlist` and `bndpartner` are maintained by three mutually inconsistent erase paths |
| status | **actively wrong in 13 of 18 cases** | **latent: never fires in any tested model** |
| observable today | no | no |
| why it does not bite | its only reader is dead code | no tested model double-bonds one pair of molecules |
| what would make it bite | re-enabling `correct_structure()` | any model where two molecules bind through two interfaces |

Neither changes any output today. Both are traps for the next person.

## The three lists

```cpp
std::vector<int> bndlist;     // relative indices of this molecule's bound interfaces
std::vector<int> bndpartner;  // the molecule each of those is bound to
std::vector<int> bndRxnList;  // the reaction that formed each bond
```

`bndpartner`'s declaration carries the comment *"Make this have the same
numbering !!"*. That reads like a statement of invariant. It is a wish.

The authoritative record of what is bound to what is not these lists at all: it
is `interfaceList[i].isBound` together with
`interfaceList[i].interaction.partnerIndex`. The three vectors are a denormalised
index over that, maintained by hand at every association and dissociation site.

## Method

Rather than reason about it, a build was instrumented to compare all three lists
against the interfaces themselves at the end of every timestep, for every live
molecule, and to count:

- `bndlist.size() != bndpartner.size()` — the two disagreeing with each other;
- either disagreeing with the set of interfaces actually flagged `isBound`;
- the same partner appearing twice in `bndpartner`, which is the precondition
  for defect 2;
- `bndRxnList.size()` against the true bond count, and the high-water mark of
  each.

## Defect 1 — `bndRxnList` is write-only, unbounded, and its reader is dead

### What the code does

**Pushed in exactly one place.** `associate_box.cpp:1000-1001`, for both
reactants. Not by `associate_sphere.cpp`, not by
`associate_ImplicitLipid_box.cpp`, not by `associate_ImplicitLipid_sphere.cpp`.
So three of the four association paths form a bond and record nothing.

**Never erased.** The code that would erase it on dissociation is commented out
at `break_interaction.cpp:110-128`.

**Read in exactly one place, which is dead.**
`correct_structutre.cpp:41` reads `bndRxnList[0]`. Its only two call sites, at
`break_interaction.cpp:209` and `:216`, are commented out. `correct_structure()`
is unreachable.

### Measured

`rxnVsTruth` counts molecule-timesteps where `bndRxnList.size()` differs from
the number of interfaces actually bound. `maxRxnList` is the largest
`bndRxnList` ever seen; `maxBnd` the largest number of simultaneous bonds.

| case | `rxnVsTruth` | `maxBnd` | `maxRxnList` |
| --- | ---: | ---: | ---: |
| `sphere` | 16,427,893 | 1 | **0** |
| `implicit_lipid` | 15,220,974 | 1 | **0** |
| `hetTrimer` | 3,089,810 | 2 | 6 |
| `homoTrimer` | 1,816,669 | 2 | **24** |
| `closed_homoTrimer` | 1,816,669 | 2 | **24** |
| `rev_3D` | 831,904 | 1 | 3 |
| `trimer` | 185,226 | 2 | 7 |
| `unimol_state_change` | 127,128 | 1 | 2 |
| `clathrin` | 106,372 | 3 | 3 |
| `michaelis_menten` | 98,562 | 1 | 1 |
| `cluster_gagsphere` | 52,429 | 2 | **0** |
| `rev_3Dto2D` | 12,976 | 1 | 2 |
| `rev_2D` | 620 | 1 | 2 |
| `hexamer`, `mem_localization`, `gagsphere`, `cluster_mem_loc`, `compartment` | 0 | 0-2 | 0-2 |

Three regimes, and the mechanism explains all three exactly:

1. **Grows without bound.** `homoTrimer` reaches **24 entries for a molecule
   that is never bound more than twice**. Every association appends; no
   dissociation removes. The list is a log of every bond the molecule has ever
   formed, not a description of the bonds it has.
2. **Stays empty while bonds exist.** `sphere` and `implicit_lipid` bind through
   `associate_sphere` and `associate_ImplicitLipid_*`, which never push, so
   `bndRxnList` is empty for 15-16 million molecule-timesteps in which the
   molecule is bound.
3. **Happens to track.** `hexamer` and `gagsphere` bind through `associate_box`
   and, over the iteration counts used, dissociate rarely enough that the list
   is never left stale.

### Why it does not bite today, and what would make it

`correct_structure()` is unreachable, so nothing reads the list. If the
commented-out calls at `break_interaction.cpp:209` are restored — and they are
marked *"This is commented out when merging"*, i.e. disabled by accident during
a merge rather than retired on purpose — then `bndRxnList[0]` becomes live, and:

- on `sphere` or any implicit-lipid model, `bndRxnList` is **empty**, so
  `bndRxnList[0]` indexes past the end of a zero-length vector. That is
  undefined behaviour, and `operator[]` will not diagnose it;
- on a dissociating box model, `[0]` is whichever reaction bound the molecule
  **first**, quite possibly a bond broken long ago. The value is used as
  `forwardRxns[...].bindRadius`, so a structure correction would be applied at
  the wrong radius.

### A third defect in the same function

`correct_structure()` cannot work even if it were called. It takes
`const std::vector<Molecule>& moleculeList`, then does:

```cpp
Molecule movingMol{moleculeList[proteinIndex]};   // a copy
movingMol.comCoord = dispVec + movingMol.comCoord;
for (unsigned int j{0}; j < movingMol.interfaceList.size(); ++j)
    movingMol.interfaceList[j].coord = dispVec + movingMol.interfaceList[j].coord;
```

It computes the displacement correctly and applies it to a **local copy** that
is then discarded. Nothing is written back, and the `const` parameter means
nothing could be. Re-enabling the call sites would therefore restore a function
that reads a stale or out-of-bounds reaction index in order to compute a
correction it throws away.

### Suggested resolution

In increasing order of ambition:

1. **Delete `bndRxnList` and `correct_structure()`.** Nothing reads either. This
   removes 24 bytes from `Molecule` and a trap, and is the only option that is
   verifiable today — the bitwise suite proves a no-op because nothing reads
   them.
2. **Or repair the list and leave the reader dead**, by restoring the
   index-parallel erase described under defect 2 and adding the missing pushes
   to the other three association paths.
3. **Or repair the list, fix `correct_structure()` to take a non-const
   `moleculeList` and write back, and re-enable the calls.** This changes
   trajectories and needs a model that exercises it — the call sites are gated
   on `complexList[...].onFiber`, which no input under `sample_inputs` sets, so
   it also needs a fiber model before it can be tested at all.

Option 1 is what the evidence supports. Options 2 and 3 are only worth it if
somebody knows what `correct_structure()` was for.

## Defect 2 — `bndlist` and `bndpartner` have three inconsistent erase paths

### What the code does

They are pushed together, always, at four sites: `associate_box.cpp:996-999`,
`associate_sphere.cpp:525-528`, `associate_ImplicitLipid_box.cpp:362-363` and
`:444-445`, `associate_ImplicitLipid_sphere.cpp:280-281`. On the push side the
invariant holds.

They are erased by three different mechanisms.

**(a) Explicit dissociation — `break_interaction.cpp`.** The two lists are
erased hundreds of lines apart, on different branches, by different keys, both
using the erase-remove idiom:

```cpp
// line 20
reactMol1.bndpartner.erase(remove(begin, end, reactMol2.index), end);   // by PARTNER
// line 280
reactMol1.bndlist.erase(remove(begin, end, relIface1), end);           // by INTERFACE
```

`std::remove` drops **every** matching element, not one. An interface can be
bound to at most one partner, so `relIface1` occurs at most once and `bndlist`
loses exactly one entry. A molecule can be bound to the same partner through
**several** interfaces, so `reactMol2.index` can occur several times and
`bndpartner` loses **all** of them.

Dissociating one bond between a doubly-bonded pair therefore removes one
`bndlist` entry and two `bndpartner` entries. The lists diverge in length, and
the molecule stops listing a partner it is still bound to.

The same function has a third asymmetry: `break_interaction.cpp:177`, the
cancel-dissociation path, pushes `bndpartner` alone without a matching
`bndlist` push.

**(b) Implicit-lipid dissociation — `break_interaction_implicitlipid.cpp:46-47`.**
Different semantics again: `find_if` plus `erase` removes exactly one from each,
which is the correct count. But `erase(end())` is undefined behaviour, and
neither call checks that `find_if` found anything.

**(c) MPI cleanup — `src/mpi/delete_disappeared.cpp:133-135`.** This one
*assumes* the invariant:

```cpp
for (int i = 0; i < partner.bndpartner.size(); i++) {
  if (partner.bndpartner[i] == targMolIndex) {
    partner.bndpartner.erase(partner.bndpartner.begin() + i);
    partner.bndlist.erase(partner.bndlist.begin() + i);
  }
}
```

It erases both at the same index, which is only correct if they are parallel —
the property (a) does not preserve. It also erases while iterating without
decrementing `i`, so two consecutive matching entries leave the second in place.

### Measured

Across all 18 cases: `sizeMismatch=0`, `listVsTruth=0`, `partnerVsTruth=0`,
**`dupPartner=0`**.

`dupPartner` is the one that matters. It counts molecules holding the same
partner index twice in `bndpartner` — the precondition for (a)'s asymmetry. It
is zero everywhere, including `gagsphere` and `cluster_gagsphere`, the two models
with ring closure and therefore the likeliest to form a second bond between an
already-bound pair.

**So the two lists agree with each other and with the interfaces in every tested
model, and this defect never fires.** The maximum number of simultaneous bonds
on any molecule reaches 3 (`clathrin`), but never twice to the same partner.

### Why it is still worth fixing

The invariant is not enforced anywhere; it holds by accident of which models are
run. A model in which two molecules bind through two interface pairs — which the
reaction language permits, and which `bindRadSameCom` and `irrevRingClosure`
exist to support — breaks it on the first dissociation of such a pair. The
consequences are not subtle:

- `check_bimolecular_reactions()` decides whether a pair is already bound by
  scanning `bndpartner`. A molecule that has dropped a partner it is still bound
  to would be offered that binding again, and could form a second bond on an
  interface that already has one.
- `bndlist.size() > 0` gates the `excludeVolumeBound` path.
- `delete_disappeared.cpp` would erase mismatched index pairs, corrupting the
  partner list of a third molecule.

### The fix already exists, commented out

`break_interaction.cpp:110-128` contains exactly the right approach, dated
`2023-01-04`, with the comment *"For the sake of safty, all indeces are found
according to bndlist (unique elements)"*:

```cpp
auto it = std::find_if(reactMol1.bndlist.begin(), reactMol1.bndlist.end(),
                       [&](const size_t& iface) { return iface == relIface1; });
int eraseIndex1 = std::distance(reactMol1.bndlist.begin(), it);
reactMol1.bndlist.erase(it);
reactMol1.bndpartner.erase(reactMol1.bndpartner.begin() + eraseIndex1);
reactMol1.bndRxnList.erase(reactMol1.bndRxnList.begin() + eraseIndex1);
```

Locate by the unique key (the interface), take its position, erase all three at
that position. That is the correct algorithm and it maintains all three lists
including `bndRxnList`, which would also resolve defect 1.

Someone found this problem, wrote the fix, and it was disabled — the block is
introduced by *"This is commented out when merging"*.

**The disabled code has a bug in its own second half.** For `reactMol2` the order
is inverted:

```cpp
auto it = std::find_if(reactMol2.bndlist.begin(), reactMol2.bndlist.end(), ...);
reactMol2.bndlist.erase(it);                                   // erases first
int eraseIndex2 = std::distance(reactMol2.bndlist.begin(), it); // then uses `it`
```

`std::distance` is called on an iterator that `erase` has already invalidated.
The `reactMol1` half above it gets the order right. So the block cannot simply
be uncommented; the second half needs the two statements swapped, and both
halves need a guard for `find_if` returning `end()`.

### Suggested resolution

1. Restore the `2023-01-04` erase in `break_interaction.cpp`, with the
   `reactMol2` statement order corrected and an `end()` guard on both halves.
2. Give `break_interaction_implicitlipid.cpp:46-47` the same treatment; it needs
   the `end()` guard regardless of anything else.
3. Fix the erase-while-iterating in `delete_disappeared.cpp` (decrement `i`, or
   iterate backwards), or rewrite it in terms of the same locate-then-erase
   helper.
4. Factor all three into one `break_bond(mol, relIface, partnerIndex)` helper so
   there is one erase path rather than three. The three sites already disagree;
   a shared helper is the only way that stays fixed.
5. Add a model with a double bond between one pair of molecules to
   `coverage_cases.tsv`. Without one, none of the above is testable — the suite
   returns "identical" for a correct fix and for a wrong one alike, which is the
   rule section 15 of `RESULTS.md` records.

Step 5 is the prerequisite for trusting steps 1-4, and is the reason this
document recommends investigation before repair.

## Reproducing

The instrumented build is not committed; it was a temporary probe on
`EXEs/nerdss.cpp`, in the end-of-timestep bookkeeping loop, guarded by
`NERDSS_CHECK_BONDS`. It walked `moleculeList` and, per molecule, built the true
bound set from `interfaceList[i].isBound` and
`interfaceList[i].interaction.partnerIndex`, then compared `bndlist`,
`bndpartner` and `bndRxnList` against it and against each other, accumulating the
five counters tabulated above and printing a `BOND_SUMMARY` line on the final
iteration.

```bash
make clean && make serial CXXFLAGS='-std=c++0x -DNERDSS_CHECK_BONDS'
CASES_FILE=cases.tsv          ./benchmarks/nerdss_optimized/run_suite.sh bin/nerdss bond_main 1
CASES_FILE=coverage_cases.tsv ./benchmarks/nerdss_optimized/run_suite.sh bin/nerdss bond_cov  1
grep -rh BOND_SUMMARY benchmarks/nerdss_optimized/results/bond_*/runs/*/rep1/stderr.log
```
