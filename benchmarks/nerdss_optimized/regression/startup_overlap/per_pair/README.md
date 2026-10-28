# Startup overlaps by interface pair, not by molecule type

`checkAndResolveOverlap()` re-places a randomly started molecule when one of its
interfaces is too close to a partner's. It used to decide "too close" by
molecule type: for each interface of the first molecule it took the first
reaction between the two types, then held that interface to that reaction's
sigma against *every* interface of the second molecule, whichever interface or
state the reaction named. An interface that reacts with nothing excluded
partners as if it were the binding site.

It now uses `interfacesOverlap()`, the test molecule creation
(`moleculeOverlaps()`) and add files (`moleculeOverlapsForRestart()`) already
used: a pair overlaps when `isReactant()` matches both interfaces (type,
interface, state, bound) to the two reactants of one bimolecular reaction, and
they are closer than that reaction's sigma.

The model:

* 300 P and 300 Q in a 60 nm box, nothing else; `nItr = 1`, read from
  `DATA/initial_crds.xyz`.
* P has one interface `p`, 0.5 nm from its center.
* Q has `q`, 0.5 nm from its center, which binds `p` at sigma 2 nm, and `r`,
  3 nm out on the other side, which reacts with nothing.

```bash
./check.py <path-to>/nerdss
```

runs seeds 1-5 (`SEEDS` overrides), about 10 s in all. Each start must have no
`p`-`q` pair within sigma, and over all seeds there must be at least 4 per seed
of `p` within sigma of an `r`, which chance puts there about 10-20 times per
seed.

Before the fix (on PR #91's c1f7ec32):

```
seed 1: 300 P, 300 Q, re-placements 27, p-q pairs < sigma 0, p-r pairs < sigma 0
...
  FAIL  p-r pairs within sigma over 5 seeds: 0 (>= 20; r reacts with nothing and must not exclude p)
```

With it:

```
seed 1: 300 P, 300 Q, re-placements 16, p-q pairs < sigma 0, p-r pairs < sigma 9
...
seed 5: 300 P, 300 Q, re-placements 11, p-q pairs < sigma 0, p-r pairs < sigma 20
PASS
```
