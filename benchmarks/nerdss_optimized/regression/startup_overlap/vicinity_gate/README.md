# Startup overlaps beyond the 2 nm vicinity margin

This model exists for `are_molecules_in_vicinity()`, the gate
`checkAndResolveOverlap()` puts in front of every pair it compares while
`generate_coordinates()` removes overlaps from a new simulation's random start.
An overlap is an interface pair of two types that share a bimolecular reaction,
closer than that reaction's sigma; the second molecule is then drawn again.

The gate used to pass a pair only when the centers were closer than
r1 + r2 + 2 nm, where r is a template's largest interface distance from its
center. Two interfaces can be closer than sigma while the centers are anywhere
up to r1 + r2 + sigma apart, so for sigma > 2 nm every overlap in between was
never looked at, and the log still said "No overlapping proteins found". The
gate is now r1 + r2 + the largest bindRadius of a bimolecular reaction between
the two types, the bound the triangle inequality gives for any orientation.

The model makes those overlaps common and easy to count:

* 20 P and 20 Q in a 30 nm box, nothing else.
* One interface per molecule, 0.5 nm from the center, so r1 + r2 + 2 = 3 nm.
* `P(p) + Q(q)` at `sigma = 5`, so an overlapping pair's centers can be up to
  6 nm apart.
* `nItr = 1`: only the start is under test, and it is read from
  `DATA/initial_crds.xyz`, written before the first step.

Run it with:

```bash
./check.sh <path-to>/nerdss
```

`check.sh` starts the model at seeds 1-10 (`SEEDS` overrides) and counts P-Q
interface pairs closer than sigma in each start; any pair fails the test. It
takes about a second.

## What it is built to show

On `nerdss-optimized-1.4` at 4017ece8, before the fix:

```
seed 1: 20 P, 20 Q, 3 P-Q interface pairs within 5.0 nm (3 with centers >= 3.0 nm apart), closest 3.328 nm
seed 2: 20 P, 20 Q, 5 P-Q interface pairs within 5.0 nm (5 with centers >= 3.0 nm apart), closest 3.404 nm
...
seed 7: 20 P, 20 Q, 10 P-Q interface pairs within 5.0 nm (10 with centers >= 3.0 nm apart), closest 3.149 nm
...
FAIL: a P and a Q start within sigma of each other
```

Every pair left behind has its centers at least 3 nm apart, the old margin:
the gate let none of them through, and the pairs it did let through were
resolved. With the fix:

```
seed 1: 20 P, 20 Q, 0 P-Q interface pairs within 5.0 nm (0 with centers >= 3.0 nm apart), closest 5.417 nm
...
PASS
```
