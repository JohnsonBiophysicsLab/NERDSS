# Creating molecules into a crowded box

This model exists for `moleculeOverlaps()`, the test that decides where
`create_molecule_and_complex_from_rxn()` may leave a molecule a creation
reaction has made. No case table runs a model where that test has anything to
reject: `create_destroy` has no bimolecular reaction, so nothing a created
molecule can overlap; `clock_model` creates into 1612 nm of mostly empty box,
where a partner lands within 8 nm of a place about once in 34000; and
`compartment` reaches the test three times in a whole run.

So the model is built to make the test bite, and to make its verdict readable
off the final coordinates:

* 1000 A and, from `NULL -> B(b)`, about 36000 B in a 100 nm box. 1000
  exclusion spheres of 5 nm radius cover 39% of that volume, so about two
  places in five drawn at random have to be rejected.
* `A(a) + B(b)` at `sigma = 5`. Both templates put their one interface on the
  center, so an interface pair and a center pair are the same distance apart and
  a plain distance count over the final coordinates is the whole verdict.
* `D = [0,0,0]` for both. Nothing moves after it is placed, so the final
  coordinates *are* the places the test accepted, and no association can consume
  an overlapping pair before it is counted. It also puts `rMaxLimit` at exactly
  `sigma`: `set_rMaxLimit()` adds `3 * sqrt(6 * Dtot * dt)` and both interfaces'
  distances from their centers, and all three are zero here.

`rMaxLimit = 5` over a 100 nm box gives a 20 x 20 x 20 grid of SubBoxes 5 nm on
a side, which is the point. A 5 nm exclusion sphere is 523.6 nm^3 and a SubBox
is 125 nm^3, so a test that walks only the SubBox a place falls in can see at
most 23% of the sphere it has to keep clear -- at most, because the sphere is
centred wherever in that SubBox the place happened to land, and only the part
of it inside the SubBox is looked at.

Run it with:

```bash
./check.sh <path-to>/nerdss
```

`check.sh` runs the model and counts A-B pairs closer than `sigma` in
`DATA/final_coords.xyz`. It takes about 20 s on a build that places the
molecules properly and several minutes on one that does not, because the pairs
such a build leaves behind are offered to the reaction machinery on every step
for the rest of the run.

## What it is built to show

On `nerdss-optimized` with the 27-SubBox walk, at the seed `check.sh` defaults
to:

```
1000 A, 36234 B, 0 A-B pairs within 5.000 nm, closest 5.000032 nm
PASS
```

Before it -- the same build with the walk covering one SubBox instead of 27:

```
1000 A, 36263 B, 13512 A-B pairs within 5.000 nm, closest 0.507452 nm
FAIL: a created B was left within 5.0 nm of an A
```

17866 pairs is what placing every B at random would give: 0.49267 A lie within
5 nm of a point drawn uniformly in the box, integrated over those 1000 A by
Monte Carlo, and slightly below the 0.5236 a boundless box would give because a
ball centred near a face reaches outside, where there are no A. So the
one-SubBox walk removed 24.4% of the pairs -- against the 23.9% of a 523.6 nm^3
exclusion sphere that a 125 nm^3 SubBox can hold, which is the ceiling on what
it could ever remove -- and the 27-SubBox walk removes all of them.

The B counts differ between the two because rejecting a place consumes random
numbers, which moves the Poisson draw that decides how many molecules each step
creates. The run took 19 s against 9 min 47 s.
