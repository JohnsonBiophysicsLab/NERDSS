# Narrowing `Molecule`: implementation plan

Branch `molecule-layout`, from `nerdss-optimized` at `cc1d954`.

This turns section 22 of
[`RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md) from a proof of
concept into a change that can ship. Section 22 established four things this
plan is built on:

1. A host-side layout change can be **bitwise identical** -- the `MolHotView`
   proof of concept was byte-identical on all 13 cases.
2. The layout **does** cost time, and the cost grows with the model: a control
   build that inflated `Molecule`'s stride from 656 to 1,680 bytes measured
   0.953x at 2,000 molecules, 0.942x at 6,410, 0.864x at 10,000 and **0.827x at
   40,000**.
3. The proof of concept did not collect that headroom, because it *shadowed*
   the class instead of narrowing it and paid an O(N)-per-timestep refresh
   against an O(pairs) search.
4. The existing suite cannot see any of this -- `cases.tsv` tops out at 3,955
   molecules, and the effect appears between 6,410 and 10,000.

## The target is a hot/cold split, not per-field SoA

The review asked for "flat structure-of-arrays storage". Taken literally -- one
array per field -- that would make this program **slower**, and it is worth
being explicit about why before writing any code.

The pairwise search reaches molecules by *gather*: it walks
`SubBox::memberMolList`, which holds unsorted molecule indices, and for each one
reads five or six fields of that single molecule. Per-field SoA puts those six
fields on six different cache lines. Today they sit on four. The access pattern
is wrong for SoA and right for a compact struct.

What the padding control actually measured is **stride**, not field grouping.
Inflating the object hurt; the lever is therefore to make the object small. So:

- keep `Molecule` an array of structs;
- merge groups of parallel `std::vector`s that are always indexed together into
  one vector of a small struct -- fewer headers in the object, and one pointer
  chase instead of six at the use site;
- leave the gather-hot scalars where they are and, last, reorder them so they
  land on one line.

This is the same transformation SoA is usually reaching for -- fewer bytes
touched per useful byte -- applied in the direction this access pattern rewards.

## Inventory

`Molecule` is 656 bytes. 528 of them are the 24-byte headers of 22 heap-owning
members.

| group | members | bytes | when read |
| --- | --- | ---: | --- |
| gather-hot | `myComIndex`, `molTypeIndex`, `index`, `comCoord`, flags, `interfaceList`, `freelist`, `bndlist`, `bndpartner` | ~148 | every candidate pair |
| crossings | `crossbase`, `mycrossint`, `crossrxn`, `probvec` | 96 | only pairs that survive the gates |
| reweighting | `prevlist`, `currlist`, `prevmyface`, `currmyface`, `prevpface`, `currpface`, `prevnorm`, `currprevnorm`, `ps_prev`, `currps_prev`, `prevsep`, `currprevsep` | 288 | only inside the reweighting lookups |
| association scratch | `tmpComCoord`, `tmpICoords` | 48 | only during association |
| misc / MPI | `mass`, `transmissionProb`, `id`, `complexId`, five bools | ~45 | rarely |

The reweighting group is **44% of the object** and is never touched by the
rejection cascade. It is the first target.

## The transformation, concretely

The twelve reweighting vectors are two groups of six parallel arrays. Every
read indexes all six at the same position -- `get_prevNorm()`-style lookups scan
`prevlist` for a partner match and then read `prevnorm[i]`, `ps_prev[i]`,
`prevsep[i]` at that same `i`. They are an array of structs written as six
arrays.

```cpp
struct ReweightEntry {
    double norm;      // prevnorm  / currprevnorm
    double survProb;  // ps_prev   / currps_prev
    double sep;       // prevsep   / currprevsep
    int partner;      // prevlist  / currlist
    int myFace;       // prevmyface / currmyface
    int partnerFace;  // prevpface / currpface
};                    // 40 bytes

std::vector<ReweightEntry> prevReweight;  // 24
std::vector<ReweightEntry> currReweight;  // 24
```

288 bytes of headers become 48. Six allocations per group become one, and the
lookup does one pointer chase instead of six.

The same shape applies to the crossing lists, which are also four parallel
arrays indexed together:

```cpp
struct CrossEntry {
    double prob;                 // probvec
    std::array<int, 3> rxn;      // crossrxn
    int partner;                 // crossbase
    int myIface;                 // mycrossint
};                               // 32 bytes

std::vector<CrossEntry> crossings;  // 24
```

96 bytes become 24.

Together: **656 -> 344 bytes, a 1.9x narrowing** -- the opposite direction from
the control build's 2.56x widening, which cost 21% at 40,000 molecules.

## Stages

Each stage is independently buildable, independently bitwise-verifiable and
independently revertible. No stage begins until the previous one is byte
identical on `cases.tsv` and measured on the size sweep.

### Stage 0 -- a benchmark that can see the effect

Nothing later can be judged without this, and its absence is why section 22.5
nearly reached the wrong conclusion.

- `benchmarks/nerdss_optimized/scale_cases.tsv`: constant-density `rev_3D` at
  2,000 / 10,000 / 40,000 molecules, plus `enzyme` at 6,410, the largest model
  in `sample_inputs` that runs.
- `benchmarks/nerdss_optimized/make_scale_inputs.sh`: generates the scaled
  inputs, so the table is reproducible rather than checked-in output. Copy
  numbers scale by *k*, the box side by the cube root of *k*, so density,
  reactions and timestep are untouched.
- `interleaved_timing.sh`: honour `CASES_FILE` and `ROOT_DIR`. It currently
  hardcodes `$SCRIPT_DIR/cases.tsv`, which is why section 11.6 recorded a
  measurement being lost to it.
- Measure in CPU time (user+sys), not wall clock, per section 11.6.

**Exit criterion:** the sweep reproduces section 22.6's padding control --
roughly 0.95x at 2,000 rising to 0.83x at 40,000 -- so we know the harness can
detect the effect we are about to chase.

**Result: met, with one caveat.** Baseline against the 1,680-byte padded build,
5 repetitions, medians of CPU time:

| case | pad/base, committed harness | pad/base, section 22.6 |
| --- | ---: | ---: |
| `scale_2k` | 0.971 | 0.953 |
| `enzyme` | 1.077 | 0.942 |
| `scale_10k` | 0.887 | 0.864 |
| `scale_40k` | **0.834** | 0.827 |

`scale_40k` reproduces tightly and `scale_10k` closely, so the harness detects
the effect where the effect is large. `enzyme` does not reproduce and should not
be read: at roughly five seconds a run its spread is 0.498 s on the baseline and
0.810 s on the padded build, which swallows the difference. Either raise its
`nItr` or treat it as a correctness case rather than a timing one.

### Stage 1 -- merge the reweighting vectors

`Molecule` 656 -> 416 bytes. 123 member references across 16 files.

Three constraints found before starting:

- **The restart file format writes the six `prev*` lists as six separate
  length-prefixed sequences** (`write_restart.cpp:532-554`,
  `read_restart.cpp:1051-1091`). The on-disk format must not change. The writer
  loops the merged vector six times; the reader fills it field by field.
- `clear_reweight_vecs()` swaps `curr*` into `prev*` and clears `curr*`. That
  becomes one swap and one clear instead of six of each -- and must stay a swap,
  because `.clear()` retaining capacity is what keeps the sweep allocation-free
  in steady state (section 16.2, allocator at 0.94%).
- **Correction to this plan as first written.** It said the reweighting fields
  are absent from the MPI serialisation. That is true only of the six `prev*`
  fields, whose calls are commented out. `serialize_back()` and
  `deserialize_back()` **do** carry all six `curr*` fields. They now assemble the
  same six arrays from `currReweight` and read them back the same way, so the
  wire format is byte-identical to the pre-merge build. The temporaries this
  costs are free in context: `serialize_primitive_vector()` takes its vector by
  value, so all six were already being copied there.
- `src/mpi/id_index.cpp` rewrites `prevlist` between molecule IDs and indices in
  both directions. Those two loops now walk `prevReweight` and touch `.partner`.

**Exit criteria:** `sizeof(Molecule) == 416`; 13 of 13 cases byte identical;
restart round-trip byte identical; `make mpi` builds; size sweep recorded.
`coverage_cases.tsv` is run too, because `cases.tsv` does not enter every path
stage 1 edits -- `gagsphere` is the only model that reaches the
`excludeVolumeBound` half of `check_bimolecular_reactions()`, and `compartment`
the only one that reaches the transmission path.

**Result: all correctness criteria met.** `sizeof(Molecule)` is 416, down from
656 -- 10.25 cache lines to 6.5.

| check | result |
| --- | --- |
| `cases.tsv` | 13 of 13 byte identical |
| `coverage_cases.tsv` | 5 of 5 byte identical |
| restart read path | 15 files byte identical |
| `make mpi` | builds clean |

The restart *read* path needed its own check, because the suite only proves
restart files are *written* identically: both builds were continued from the
same `rev_3D/RESTARTS/restart10000.dat` and every file they then produced
matched.

**One of the six edited sites is covered by no test.** Section 15's rule applies:
a bitwise suite that never reaches a function reports "identical" for any change
to it, in the same words it uses for code it does reach.
`determine_1D_bimolecular_reaction_probability()` is called only from the
`com1.onFiber && com2.onFiber` arm of `check_bimolecular_reactions()` at line
193, and no input under `sample_inputs` sets `onFiber` -- the same gate that
already keeps `sweep_separation_complex_rot_fiber` dark. It was missing from
`known_uncovered.tsv`, which listed the fiber sweep but not the probability
function behind the same gate; it and its two exclusive callees, `passocF_1D`
and `pirr_pfree_ratio_psF_1D`, are now listed there.

So that file's change is argued by reading, not by the suite. It is a
transliteration: the loop bound, the three-way match, the `>= RMax` test, both
`pirr_pfree_ratio_psF_1D()` calls and the `break` are unchanged except for where
the six values come from, and the six-argument `emplace_back` is identical in
form and order to the five other sites, three of which the suite does exercise.

**Timing: first pass contaminated, do not quote it.** Seven repetitions,
interleaved, medians of CPU time:

| case | base (s) | base sd | stage1 (s) | stage1 sd | stage1/base |
| --- | ---: | ---: | ---: | ---: | ---: |
| `scale_2k` | 24.370 | 8.405 | 27.350 | 5.496 | 0.891 |
| `enzyme` | 57.600 | 19.174 | 52.250 | 17.849 | 1.102 |
| `scale_10k` | 30.880 | 0.401 | 25.810 | 0.511 | 1.196 |
| `scale_40k` | 69.950 | 3.906 | 58.280 | 3.831 | 1.200 |

The absolute times give it away: `scale_10k`'s baseline is 30.88 s here against
9.12 s for the same binary on the same input in the stage 0 pass, and `enzyme`'s
is 57.6 s against 5.44 s. `ps` during the run showed **three `nerdss_mpi`
processes at 99% CPU from an unrelated session, plus a Python at 98.6%** --
four of ten cores taken throughout. This is the failure mode RESULTS.md section
11.6 records, and the remedy it gives is to check the load before trusting a
number.

Interleaving means both builds met the same conditions, so the ratios are not
meaningless -- `scale_10k`'s 5.07 s gap against 0.4-0.5 s standard deviations is
far outside noise, and `scale_40k` agrees at 1.200x. But contention inflates
everything roughly threefold and plausibly *amplifies* the effect being measured:
with four other cores hammering the shared L2, the build with the larger working
set suffers more.

**Waiting for an idle host was the wrong fix.** This machine is a desktop that
is rarely idle, and section 11.6's advice -- check `uptime` before trusting a
number -- tells you when to distrust a measurement, not how to take one.
`/usr/bin/time -l` on Apple silicon reports **retired instructions** and
**elapsed cycles**, and both are nearly load-proof. Measured against three
unrelated processes at 99% CPU, four repetitions of the same binary on the same
input:

| | instructions | spread | cycles | spread |
| --- | ---: | ---: | ---: | ---: |
| base | 89.83 G | 0.1% | 104.9 G | 1.5% |
| stage1 | 87.68 G | 0.05% | 85.5 G | 0.6% |

Wall time over the same eight runs was inflated roughly threefold against the
stage 0 pass; the counters were not. They also separate the two things a change
can do, which is exactly the distinction this work turns on: **instructions
falling means less work was done, cycles falling at an unchanged instruction
count means fewer stalls.** `interleaved_timing.sh` now records both, plus peak
RSS, and prints instructions per cycle.

Both are in the harness, so this applies to every stage from here on.

**Result.** The first pass of these numbers was taken against three unrelated
`nerdss_mpi` processes at 99% CPU and reported 1.119x to 1.389x. Re-measured
with all three builds interleaved under lighter load, stage 1 is **1.058x to
1.172x** in cycles, monotone in model size. Contention had inflated it by up to
30 points -- the caveat recorded at the time was right, and the corrected
figures are in section 24.4 of
[`RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md).

Instructions fall about 2% while cycles fall much more, so the win is stalls
rather than work; RSS falls 16-32% beyond the 240 bytes per molecule the struct
shrank, from twelve heap-owning members becoming two.

**Stage 1 passes.**

### Stage 2 -- merge the crossing lists: done

`Molecule` 416 -> 344 bytes; with stage 1, 656 -> 344.

Three things had to be established before the merge was safe, because sixteen
sites clear `crossbase` alone mid-timestep and one merged vector clears all four:

* **Nothing loops on the other three's sizes.** Every traversal bounds itself by
  `crossbase.size()`, so after that clear the leftovers are unreachable.
* **Alignment holds while populated.** An instrumented build reported any
  molecule whose four lists disagreed in length: zero across all 18 cases.
* **One site could have broken it.**
  `determine_2D_implicitlipid_reaction_probability()` seeded `probvec`
  unconditionally while pushing the other three only inside a branch, so a
  skipped branch left a stray entry that would pair later entries with the wrong
  probability. The merge removes the possibility. Inert on every tested case,
  but a latent defect closed rather than a mechanical rewrite, and marked as such
  at the site.

Byte identical on both tables and across a restart, against the **pre-branch**
baseline. **1.230x aggregate in cycles, 1.325x at 40,000 molecules.**

Two constraints from the plan as first written held up: `serialize_back()` and
`deserialize_back()` do carry all four lists, so both sides split and reassemble
the same four arrays and the wire format is unchanged; and
`determine_if_reaction_occurs()` compares reaction triples with `==` on
`std::array<int,3>`, so `CrossEntry::rxn` stays an `std::array`. 190 member
references were rewritten across 29 files.

**Stage 2 passes. Stages 3 and 4 remain unattempted.**

### Stage 3 -- reorder the survivors: done

`Molecule` 344 -> 328 bytes; across the three stages, **656 -> 328, exactly
half**.

Only the order of the member declarations changed. The rejection cascade's
fields -- `comCoord`, `myComIndex`, `molTypeIndex`, `isImplicitLipid`,
`freelist`, `bndlist`, `bndpartner` -- sat on four cache lines and now sit on
two. The 16-byte size drop is a side effect: the twelve `bool`s cluster into one
8-byte run instead of being scattered between wider fields.

A pure reorder is safe here because `Molecule` has user-provided constructors
(never aggregate-initialised), `operator==` names its fields through `std::tie`,
and `serialize()`/`deserialize()` keep their own fixed order so the MPI wire
format does not move. The one hazard is that members initialise in *declaration*
order regardless of a constructor's init list; that list was reordered to match
and the build is clean under `-Wreorder`.

Byte identical on both tables and across a restart. **1.032x over stage 2**,
1.046x at 40,000 molecules -- small, but free at runtime, and reproduced within
half a point by a second, noisier pass.

**Stage 3 passes.**

**A range, not a number.** Stage 2 measured 1.230x aggregate against baseline on
one pass and 1.389x on another; both were interleaved with tight spreads.
Dividing cycles by CPU seconds shows the first ran at 3.98 GHz and the second at
2.22 GHz -- performance cores versus efficiency cores, 128 KB/16 MB of cache
versus 64 KB/6 MB. The narrower the memory hierarchy, the more a narrower object
is worth. Cumulatively the three stages are worth roughly **1.23x on a
performance core and 1.43x on an efficiency core**, rising to 1.33x and 1.67x on
the 40,000-molecule case, with instruction counts unchanged to within 2%
throughout. See section 25.4 of
[`RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md).

### Stage 4 -- measured, and declined

Two blocks were candidates. Only one survived inspection, and it did not survive
measurement.

**The bond lists cannot be merged.** `bndlist`, `bndpartner` and `bndRxnList`
look parallel -- `bndpartner`'s comment even says *"Make this have the same
numbering !!"* -- but `break_interaction.cpp` erases `bndpartner` by partner
index through `std::remove` (which drops every match) and `bndlist` by interface
index, on different branches, and line 177 re-pushes `bndpartner` alone. The
comment is an intention, not an invariant.

Separately, `bndRxnList` is close to vestigial and looks like a latent defect:
pushed only by `associate_box`, read only at `[0]`, and its erase-on-dissociation
commented out, so on a box model it grows and never shrinks. A correctness
question, not a layout one; left alone and recorded.

**The association scratch was measured before being written.** `tmpComCoord` and
`tmpICoords` are 48 bytes across 290 references in 41 mostly-geometry files.
Instead of paying that to find out, section 22.4's control was run in reverse:
48 inert bytes appended to the current 328-byte struct, taking the stride to 376
and moving no field. Widening by exactly what stage 4 would remove bounds what
removing it could return.

**The ceiling is 1.019x aggregate** -- 1.4-1.6% on the two tightest cases,
instructions flat at 1.000x confirming the probe measured stride alone. And it
is a ceiling: any side container adds an indirection to all 290 accesses.

Stage 3 removed 16 bytes for a measured 1.032x, because those bytes were fields
the cascade reads. Stage 4 would remove three times as many for half the gain,
because they already sit past `interfaceList` at offset 144 -- the cascade reads
two cache lines today and would read two afterwards.

**Closed as measured-and-declined rather than untried.** See section 26 of
[`RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md). If `Molecule` ever
needs to be narrower, this is where the next 48 bytes are and what they are
worth.

## What is explicitly not being done

- **No per-field SoA**, for the gather reason above.
- **No device port.** Section 22.2 measured bitwise identity to be unreachable
  on a device regardless of layout -- the CUDA proof of concept's arithmetic
  prefix differed from the CPU by up to 9.55e-15 -- and section 22.8 measured
  the available batch at 62 to 19,841 pairs per timestep against a break-even
  near 1,000,000. Narrowing the object is a prerequisite for that work, not a
  step into it.
- **No CSR conversion of the ragged rows.** Section 22.9 bounds it at the 0.94%
  the allocator still costs, against a per-row capacity cap or a counting pass.
- **No change to append order anywhere.** `determine_if_reaction_occurs()` draws
  one `rand_gsl64()` per `crossbase` entry in list order, so append order is the
  random stream. Every stage preserves it exactly; that is what makes the
  bitwise check meaningful.

## Verification, every stage

```bash
make serial
./benchmarks/nerdss_optimized/run_suite.sh bin/nerdss <label> 1
./benchmarks/nerdss_optimized/compare_suites.sh <baseline> <label>
CASES_FILE=scale_cases.tsv ./benchmarks/nerdss_optimized/interleaved_timing.sh 7 \
    base=<baseline-bin> cand=bin/nerdss
```

A stage that is not byte identical on all 13 cases does not proceed. A stage
that is byte identical but shows no gain on the sweep is recorded and reverted,
the way section 11 recorded `find_which_reaction()`.
