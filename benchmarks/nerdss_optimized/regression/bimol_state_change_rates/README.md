# Reversible bimolecular state change: equilibrium vs the two rates

`K(k) + A(a~X) <-> K(k) + A(a~Y)` with `ka` (X -> Y) = 10 and `kb` (Y -> X) = 50
nm^3/us. K is unchanged by either direction, so at equilibrium Y/X = ka/kb = 0.2,
a fraction of 1/6 in Y. `check.sh` runs the model and checks that fraction, and
that the copy counters add up.

```bash
./check.sh <path-to>/nerdss
```

No input under `sample_inputs/` has a bimolecular state change, so neither suite
table reaches this path.
[`../bimol_state_change`](../bimol_state_change) is the issue #8 reproducer for
the lookup. It cannot run past a few hundred steps: its `nan` assocAngles reach
`theta_rotation()` unguarded and send the coordinates to NaN. Here every interface
sits on its molecule's COM. Both molecules are then points, so a state change
applies no rotation and the run stays finite.

## Four defects stand between the model and the right answer

| # | where | effect | fixed by |
| --- | --- | --- | --- |
| 1 | `perform_bimolecular_state_change_{box,sphere}` | counts under the relative state index, so K(k) and K(kd) move and A(a~X), A(a~Y) never do | `9c95f918` (this branch) |
| 2 | `populate_reaction_lists()` | the facilitator never lists the product partner, so the back direction is tested only for pairs with A first | `c0355bc4` (reversible-lipid branch) |
| 3 | `determine_if_reaction_occurs()` | a chosen back state change is matched against the forward reactants, so it fails whenever K is `mol` | `73bf19ed` (reversible-lipid branch) |
| 4 | `determine_*_probability()`, `evaluate_binding_within_complex()` | the back direction runs at the forward rate | `9fe61395` (this branch) |

## Results

`check.sh`, seed 20260810 (clang -O3 builds on `nerdss-optimized-1.4` a998d6b4):

| build | fraction in Y | Y/X | verdict |
| --- | --- | --- | --- |
| a998d6b4 | 0.0000 (counters) | - | FAIL, counters do not add up |
| + 1 | 0.9823 | 55.6 | FAIL, back direction never fires |
| + 1, 4 | 0.9800 | 49.0 | FAIL, back direction never fires |
| + 1, 2, 3 | 0.4908 | 0.96 | FAIL, both directions at ka |
| + 1, 2, 3, 4 | 0.1780 | 0.217 | PASS |

Six seeds (11-16) per rate pair, each a mean over 2-10 ms:

| ka/kb | ka/(ka+kb) | + 1, 2, 3 | + 1, 2, 3, 4 |
| --- | --- | --- | --- |
| 10/50 | 0.1667 | 0.5043 +- 0.0044 | 0.1659 +- 0.0021 |
| 50/10 | 0.8333 | 0.4977 +- 0.0024 | 0.8338 +- 0.0030 |
| 10/10 | 0.5000 | 0.5043 +- 0.0044 | 0.5043 +- 0.0044 |

The per-run standard deviation at 10/50 is 0.005, so the 0.04 tolerance in
`check.sh` is 8 sigma. The ka = kb runs are byte-identical between the two
builds. Correcting each rate for diffusion on its own, with
k_eff = 1/(1/k + 1/kD) and kD = 4 pi sigma Dtot = 628 nm^3/us, would predict
0.1753 and 0.8247. The six-seed means are 4.4 and 3.0 sigma from those values,
and 0.4 and 0.2 sigma from ka/(ka+kb).
