# Bond site separation

A bond is written at the binding radius `sigma`, and a complex moves as a rigid
body, so the two interface sites of a bonded pair stay `sigma` apart for as long
as the bond lives. This test asserts that, and exists because the serial build
used to violate it.

## What it runs

`parms.inp` and `A.mol` are a 75-molecule actin model built from PDB 6BNO by
`ionerdss.build_system_from_pdb` (interface cutoff 1.0, overlap separation limit
3.0). Two reversible binding reactions give each molecule a forward and a
backward site on two different interfaces, so the model grows filaments that
readily branch and meet again — which is what it takes to get two associations
into one timestep with overlapping complexes.

```bash
./run_test.sh [path-to-nerdss] [work-dir]
```

`run_test.sh` runs 40000 iterations at five seeds and passes each run to
`check_bond_site_separation.py`. It takes a few seconds per seed. The work
directory is a temporary one that is removed on success and kept on failure.

`check_bond_site_separation.py` can be pointed at any run that was made with
`bondedComplexWrite` set:

```bash
python3 check_bond_site_separation.py path/to/run --tolerance 1.15
```

It reconstructs each bonded site's position from `DATA/COMPLEXES/*.json` as
`com + R(q) . template_site`, using the `.mol` files next to the run, and reports
the largest separation it finds in units of the bond's `sigma`.

## The tolerance

1.15 sigma, not 1.0, because a loop closure is legitimately allowed to bond at a
separation up to `bindRadSameCom * sigma`, and `bindRadSameCom` defaults to 1.1.
A reaction between two molecules that are already in one complex bonds them in
place without moving anything, so whatever separation they are at becomes the
bond length, and that bond then persists. On this model the largest legitimate
separation observed is 1.063 sigma.

The margin above 1.1 has a cost: the defect also produces bonds just over the
limit, at 1.110 and 1.115 sigma on this model (see
[`docs/bond_site_separation_defect.md`](../../docs/bond_site_separation_defect.md)),
and those pass this check. `../loop_closure_separation_check.cpp`, run by
`make checks`, pins the criterion itself at its exact boundary, with no
dependence on a trajectory.

## The defect this guards against

`check_bimolecular_reactions()` admits a pair in two different complexes on a
diffusion-based `Rmax`, which is far looser than `sigma`, because associating
them would first move the two complexes together. That crossing is recorded once
per timestep and consumed later in the same timestep. If an earlier association
in that timestep merged the two complexes the pair belongs to, the crossing
arrives at `associate()` as a loop closure — and the loop-closure branch bonds in
place, moving nothing, so it wrote a bond between sites that were still many
sigma apart. The result was a single reported complex made of filaments that
never touched: on this model, 15% of runs at 100000 iterations and 23% at
1000000, with separations up to 34 sigma.

`associate()` now re-applies the same-complex criterion,
`R1 < bindRadSameCom * bindRadius`, against the complexes as they are at the
moment of association.

## Seeds

The five seeds in `run_test.sh` are those in 1..60 whose unfixed run first shows
a far-apart bond before iteration 40000, measured on `nerdss-optimized` at
8025fd7:

| seed | first frame with the bad bond | worst separation |
| ---: | ---: | ---: |
| 32 | 16800 | 7.5 sigma |
| 21 | 23700 | 16.2 sigma |
| 46 | 28200 | 19.5 sigma |
| 54 | 33400 | 6.4 sigma |
| 15 | 36200 | 34.4 sigma |

All five fail on the unfixed build and pass on the fixed one. The defect needs a
coincidence, so most seeds never reach it, and the seeds only keep reaching it
while their trajectories stay the same: a later change that moves the random
stream can quietly turn this into a test that passes for the wrong reason. That
is what `loop_closure_separation_check.cpp` is for. `SEEDS` and `TOLERANCE` in
the environment override the defaults.
