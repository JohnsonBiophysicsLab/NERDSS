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
| status | **actively wrong in 13 of 18 cases** | **actively wrong on `homoTrimer` and `closed_homoTrimer`** |
| observable in serial output today | no | no |
| why it does not bite | its only reader is dead code | nothing in the serial path pairs the two lists by index |
| what would make it bite | re-enabling `correct_structure()` | the two MPI sites that do pair them by index |

Neither changes any **serial** output today: defect 1's reader is dead code, and
nothing in the serial path pairs `bndlist` with `bndpartner` by index. Both are
wrong in the data structures themselves, both are reachable by MPI code that
does rely on them, and both are traps for the next person.

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

- `interfaceList[bndlist[k]].interaction.partnerIndex != bndpartner[k]` — the
  pairing itself, which is the invariant the MPI code relies on and the only one
  of these that a permutation violates;
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

**Correction to the first version of this document.** It reported this defect as
latent -- "never fires in any tested model" -- on the strength of a probe that
compared list *sizes*, set membership, and duplicate partners, all of which came
back clean across 18 cases. That probe could not see a **permutation**: two lists
of the same length, holding the same elements, paired up wrongly. Re-measured
against the actual invariant, the defect is live.

The invariant the MPI sites rely on is that for every `k`,
`interfaceList[bndlist[k]].interaction.partnerIndex == bndpartner[k]`. Checked
every timestep for every live molecule:

| case | broken slots | molecule-timesteps affected |
| --- | ---: | ---: |
| `homoTrimer` | **88,490** | 44,245 |
| `closed_homoTrimer` | **88,490** | 44,245 |
| the other 11 in `cases.tsv` | 0 | 0 |
| all 5 in `coverage_cases.tsv` | 0 | 0 |

A representative instance, and it persists once it happens:

```
BOND_PAIR_BROKEN itr=1150 mol=725 slot=0 bndlist=0 bndpartner=890 truePartner=472 nBonds=2
BOND_PAIR_BROKEN itr=1150 mol=725 slot=1 bndlist=1 bndpartner=472 truePartner=890 nBonds=2
BOND_PAIR_BROKEN itr=1151 mol=725 slot=0 bndlist=0 bndpartner=890 truePartner=472 nBonds=2
```

Molecule 725's interface 0 is bound to 472 and interface 1 to 890; `bndpartner`
says the opposite. The two entries are cleanly transposed.

The size-based counters remain zero everywhere -- `sizeMismatch=0`,
`listVsTruth=0`, `partnerVsTruth=0`, `dupPartner=0`, including both `gagsphere`
models with ring closure. So mechanism (a) below, the `std::remove` asymmetry,
genuinely does not fire: no tested model binds one pair of molecules through two
interfaces. What fires is the cancel-dissociation path.

### Cause: erase-then-append in the cancel path

`break_interaction.cpp:20` removes the partner from wherever it sits, because
`determine_parent_complex_IL` downstream reads `bndpartner` and must not see the
bond being broken:

```cpp
reactMol1.bndpartner.erase(remove(begin, end, reactMol2.index), end);
```

If the dissociation is then cancelled, `break_interaction.cpp:177` puts it back
-- at the **end**:

```cpp
reactMol1.bndpartner.push_back(reactMol2.index);
```

`bndlist` is never touched on this path. A molecule with two bonds whose *first*
bond has a cancelled dissociation therefore ends up with `bndpartner` permuted
against `bndlist`, exactly the transposition above. `homoTrimer` and
`closed_homoTrimer` are the two models here that combine multiply-bonded
molecules with enough loop closure to exercise the cancellation.

### Why no serial output is wrong today

Nothing in the serial path reads the two lists at a common subscript.
`check_bimolecular_reactions()` scans `bndpartner` for membership, which is
order-independent; `bndlist.size()` gates `excludeVolumeBound`; and `bndlist[k]`
is read as an interface index in `check_dissociation*` and
`check_bimolecular_reactions` without pairing it to `bndpartner`.

The two MPI sites named above *do* pair them by index, so a run that reaches them
with `homoTrimer`-like topology gets a partner list entry rewritten for the wrong
interface, or a third molecule's lists erased at mismatched positions.

It is also a trap for any new code that assumes what the comment promises.

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


## Resolution

Both defects are fixed on branch `bond-bookkeeping-fix`. Filed as
[#18](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/18) and
[#19](https://github.com/JohnsonBiophysicsLab/NERDSS/issues/19).

### Defect 1 — deleted

`bndRxnList` and `correct_structure()` are removed, along with the two pushes in
`associate_box.cpp`, the two clears in `class_Molecule_Complex.cpp`, and the
commented-out call sites. Nothing read either, so the bitwise suite proves the
deletion is a no-op. **`sizeof(Molecule)`: 328 -> 304 bytes.**

This is option 1 of the three the investigation listed. Options 2 and 3 keep a
list nothing reads, or revive a function that cannot write back through its own
`const` parameter; neither is worth doing before somebody establishes what
`correct_structure()` was for. The deletion does not foreclose either: the git
history holds both, and the reaction that formed a bond is recoverable from
`interfaceList[i].interaction.conjBackRxn`.

### Defect 2 — one erase path instead of three

A shared `erase_bond(mol, relIface)` in `reaction_bookkeeping.cpp` locates the
bond through `bndlist` -- the only key that names it uniquely -- and removes the
entry from both vectors at that position, or does nothing if the interface is
not bound. `find_bond_slot()` exposes the lookup for the one caller that needs
the slot without erasing.

The three sites now go through it:

- **`break_interaction.cpp`.** The early `bndpartner` erase has to stay early,
  because `determine_parent_complex_IL` downstream must not see the bond, so the
  slot is located once at the top and reused. The cancel path
  **`insert`s at that slot instead of `push_back`ing**, which is the actual
  cause of the measured breakage. The committed path erases `bndlist` at the
  same slot. The `std::remove` calls, which dropped every matching entry rather
  than the one being broken, are gone.
- **`break_interaction_implicitlipid.cpp`.** Two unguarded `find_if` + `erase`
  calls selecting by different keys become one `erase_bond()`. The
  `erase(end())` undefined behaviour is gone with them.
- **`src/mpi/delete_disappeared.cpp`.** Now walks backwards, so erasing does not
  skip the entry after a match, and bounds the `bndlist` access.

### Verified

| check | result |
| --- | --- |
| pairing probe, `cases.tsv` | **88,490 broken slots -> 0**, all 13 cases clean |
| `cases.tsv`, non-restart output | 13 of 13 byte identical |
| `cases.tsv`, restart files | differ on `homoTrimer` and `closed_homoTrimer` only |
| `coverage_cases.tsv` | 5 of 5 byte identical |
| `make mpi` | builds clean |

**The physics does not move.** Every `DATA/` trajectory, observable, histogram
and species file is byte identical across all 18 cases. The only files that
change are `RESTARTS/*.dat` and `DATA/restart.dat`, on exactly the two models
where the pairing was measured broken -- and they change because
`write_restart.cpp:508-513` serialises `bndlist` and `bndpartner`, so they are
recording the data that was wrong. A seed-averaged comparison over six
independent seeds finds the copy numbers bit-identical, `|z| = 0.00`.

That the set of changed files is exactly `{restart}` and the set of changed
cases is exactly `{homoTrimer, closed_homoTrimer}` -- the two the probe flagged
-- is the strongest evidence available that the fix changes what it was meant to
and nothing else.

### One consequence worth stating

A restart file written by an older build carries the permuted lists. The fix
corrects how they are maintained, not what a stale file contains; a run resumed
from such a file inherits the permutation until the affected bonds break. Since
no serial path reads the two lists pairwise, this is only a concern for the MPI
sites named above.
