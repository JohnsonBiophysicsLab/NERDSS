# Startup overlaps in the states molecules start in

A molecules block can start a type's copies in different interface states:

```
Q : 150 (q~U), 150 (q~P)
```

`initialize_states()` applies those states, by index among the type's
molecules, and it ran only after `generate_coordinates()` in both startup paths
(`EXEs/nerdss.cpp`, `parse_input_for_a_new_simulation()`). The overlap pass
therefore saw every molecule in its interfaces' first states: a pair that can
bind only in a later state was never tested and could start inside sigma.

`generate_coordinates()` now applies the starting states
(`initialize_molecule_states()`) before the overlap passes. A re-placed molecule
gets its states back too: `create_random_coords()` puts every interface in its
first state, and the next pass has to test the molecule as it will start.
`initialize_states()` still runs later, for the implicit lipid's counts, and
sets the same states again.

The model: 300 P, all `p~P`, and 300 Q, the first 150 `q~U` and the last 150
`q~P`, in a 40 nm box; one interface each, 0.5 nm from the center. Only
`P(p~P) + Q(q~P)` binds, at sigma 2 nm. `U` is each interface's first state.

```bash
./check.py <path-to>/nerdss
```

runs seeds 1-5 (`SEEDS` overrides), a few seconds. Each start must have no P
within sigma of a phosphorylated Q, and over all seeds there must be at least
6 per seed of P within sigma of an unphosphorylated Q, which binds nothing and
which chance puts there 15-21 times per seed.

Before the fix (on PR #91's c1f7ec32, and with only the per-pair fix):

```
seed 1: 300 P, 300 Q, re-placements 0, P-Q(q~P) pairs < sigma 20, P-Q(q~U) pairs < sigma 21
  FAIL  seed 1: no P within sigma of a phosphorylated Q
...
```

With the states applied but not kept across a re-placement, 1-5 pairs per seed
are left (every seed fails). With the fix:

```
seed 1: 300 P, 300 Q, re-placements 23, P-Q(q~P) pairs < sigma 0, P-Q(q~U) pairs < sigma 21
...
PASS
```
