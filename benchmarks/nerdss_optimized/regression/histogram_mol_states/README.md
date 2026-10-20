# `histogramMolStates`: interface states in the complex histogram

`DATA/histogram_complexes_time.dat` counts the complexes present at every
`timeWrite` by composition.  Until version 1.3 the composition was the number of
molecules of each *type*, so a phosphorylated `pip2(head~P)` and an
unphosphorylated `pip2(head~U)` were one species, `pip2`, and a complex of two
P lipids and one U lipid was the same line as one of one P and two U.  For a
model whose chemistry is the state change (`sample_inputs/enzyme`,
`sample_inputs/auto_phos`) that loses the quantity of interest.

`histogramMolStates = true` in the parameters block makes the histogram tell
members of one type apart by the states of their interfaces, named as the
species of `copy_numbers_time.dat` are:

```
Time (s): 0.001
2	A(ser~P): 1. Phos: 1.
2	A(ser~P): 1. A(ser~U): 1.
3	A(ser~P): 1.
7	Phos: 1.
99	A(ser~U): 1.
```

where the default writes `2  A: 1. Phos: 1.`, `2  A: 2.`, `7  Phos: 1.` and
`102  A: 1.` for the same frame.  Only interfaces with more than one state are
named, `A(tag~U,ser~P)` when a molecule has two; an implicit lipid is always the
plain `name: n.` entry, since its links to a complex carry no state.  The switch
changes this file only.  It is off by default so that existing runs, and the
tools that parse the file, see what they always did.

The switch is recorded in the JSON restart file when set, so a restart keeps it
without an add file.  The `.dat` format cannot carry it (its layout is fixed), so
a restart from a `.dat` file needs `histogramMolStates = true` in the add file.

Both histogram writers, the serial one and `nerdss_mpi`'s per-rank one, build
the line through `complex_composition_string()` in
[`src/io/print_complex_hist.cpp`](../../../../src/io/print_complex_hist.cpp).

## The check

```bash
make serial
./check.sh ../../../../bin/nerdss [work-dir]
```

Three state-change models are each run three times with one seed: as is
(`plain`), with the switch written out as `false` (`off`) and with it `true`
(`on`):

| case | model | why |
| --- | --- | --- |
| `auto_phos` | `sample_inputs/auto_phos` | A(ser~U/~P) binds through the interface that changes state, so bound members are left out of the copy-number columns |
| `unimol_state_change` | `VALIDATE_SUITE/unimol_state_change_reversible` | A carries two state interfaces, which exercises the `A(tag~U,ser~P)` naming |
| `enzyme` | `sample_inputs/enzyme` | the motivating model: pip2(head~U/~P) on a membrane, phosphorylated by syn and bound through head by ap2 |

For each case:

1. `off` equals `plain` on every output file, byte for byte: the explicit
   default changes nothing.
2. `on` equals `plain` on every output file but the histogram, and its restart
   files differ from `plain`'s only by the parameter: the switch changes the
   histogram and nothing else.
3. `check_histogram.py` reads the two histograms and `copy_numbers_time.dat`:
   collapsing the state-aware histogram onto types (dropping the `(...)` part
   of every name and adding the counts of members that then share a name)
   must give exactly the plain histogram, frame by frame; every type with a
   multi-state interface must always be written with its states; and the
   histogram's totals must agree with the copy-number columns.  A column
   `Type(iface~S)` counts molecules whose interface is in state S *and
   unbound*; the histogram does not record bonds, so the case table declares
   which state interfaces bind (`--binding-site A.ser`) and members of such a
   type are then only counted on monomer lines, which is exact when the
   interface is the type's only binding interface.
4. A restart from the `on` run's checkpoint (with `rng_state` copied in, as
   `restart_continuation/check.sh` does) must carry states in its histogram and
   reproduce the `on` run's frames at every shared time.

The suite's `coverage_cases.tsv` also runs `auto_phos` with the switch on
(`hist_mol_states`), so the state-aware path stays in the bitwise comparison
between builds.
