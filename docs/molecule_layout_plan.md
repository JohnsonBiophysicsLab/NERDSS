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

**Result, with the counters.** Interleaved medians; `scale_40k` at 7
repetitions, the rest at 5:

| case | molecules | instructions | cycles | CPU time | RSS | IPC base -> stage 1 |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| `scale_2k` | 2,000 | 1.032x | 1.119x | 1.113x | 1.018x | 1.908 -> 2.069 |
| `enzyme` | 6,410 | 1.022x | **1.389x** | 1.383x | 1.182x | 1.783 -> 2.422 |
| `scale_10k` | 10,000 | 1.024x | 1.246x | 1.253x | 1.032x | 0.846 -> 1.029 |
| `scale_40k` | 40,000 | 1.014x | 1.264x | 1.223x | 1.095x | 1.018 -> 1.270 |

Instructions fall 1.4-3.2%, cycles fall 11-28%, and IPC rises 8-36% to make up
the gap. Merging six `push_back`s into one is the entire instruction saving and
is far too small to explain the cycles: this is a stall reduction. RSS falls by
16-32% *more* than the 240-bytes-per-molecule the struct shrank, which is twelve
heap-owning members becoming two.

Full tables, the RSS-against-prediction check, and the caveat about measuring on
a loaded host are in section 23 of
[`RESULTS.md`](../benchmarks/nerdss_optimized/RESULTS.md).

**Stage 1 passes. Proceed to stage 2.**

### Stage 2 -- merge the crossing lists

`Molecule` 416 -> 344 bytes. 190 member references.

Larger and riskier than stage 1, and it must wait for stage 1's measurement,
because if stage 1 shows no gain at 40,000 molecules the premise is wrong and
this stage should not be attempted.

Constraints:

- `serialize_back()` / `deserialize_back()` **do** carry all four lists, so the
  MPI wire format is affected. It has no on-disk compatibility requirement, but
  both sides must change together.
- `determine_if_reaction_occurs()` compares `crossrxn` entries with `==` on
  `std::array<int,3>`; keeping the member an `std::array` preserves that.
- Rows are truncated individually mid-timestep (`associate_box.cpp:1033` and
  eight similar sites clear exactly the two reacting molecules). One vector
  still supports that; four separate `.clear()` calls become one.

### Stage 3 -- reorder the survivors

Once the object is 344 bytes, move the gather-hot fields to the front so the
rejection cascade reads one line instead of four. Free at runtime, no refresh,
no new state. Measured last because it is only worth doing on the narrowed
object.

### Stage 4 -- decide about the rest

`tmpComCoord` and `tmpICoords` (48 bytes) are association-only scratch; the MPI
fields are 16 bytes. Both are candidates, neither is obviously worth its diff.
Decide from stage 3's numbers, not now.

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
