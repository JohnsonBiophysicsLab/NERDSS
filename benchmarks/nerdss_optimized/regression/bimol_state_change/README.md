# Bimolecular state-change reproducer

This model exists for the reverse bimolecular state-change corrections in
issue #8. No input under `sample_inputs/` contains a bimolecular state change,
so nothing in the validation suite reaches that code.

Run it with:

```bash
mkdir -p /tmp/bimol && cp *.mol parms.inp /tmp/bimol && cd /tmp/bimol
<path-to>/nerdss -f parms.inp -s 20260810
```

## What it is built to show

The reaction order gives `forwardRxns[1].conjBackRxnIndex == 0`, the arrangement
master could not handle. The irreversible `A(a~X) -> A(a~Y)` at rate 0 occupies
`forwardRxns[0]` without creating a `BackRxn`, so the forward and back index
spaces diverge and the state change lands on `forwardRxns[1]` / `backRxns[0]`.
Master's gate (`conjBackRxnIndex > 0`) rejects index 0, so the reverse `Y -> X`
direction was never considered; and had it passed, the caller would have indexed
`backRxns[1]`, which does not exist.

## Status on nerdss-optimized-1.4

When this model was written it ran on no build, master or optimized. Three
defects stood in the way: the first is worked around, the other two are fixed.

1. `set_rMaxLimit()` only inspects `ReactionType::bimolecular`. A model whose
   only bimolecular reaction is a state change therefore keeps `rMaxLimit == 0`,
   and `SimulVolume::Dimensions` divides the box length by it, yielding a
   negative cell count and an endless `CELL PAIR MAX EXCEEDED` rescale. This is
   still true. The `K(kd) + K(kd)` association at the end of `parms.inp` is
   scaffolding that works around it; it is listed last so the state-change
   reaction keeps indices 1 and 0.
2. The state change turned coordinates into NaN (`ERROR: [Before overlap
   checking] nan in the coordinates of molecule 22 at iteration 139` on master).
   `apply_state_change_rotations()` ran `theta_rotation()` on the unspecified
   (`nan`) theta1 and theta2 of `assocAngles`; it now skips them, as
   `associate_box` does. Before that fix every K and the A it last touched
   froze at NaN after the first event, so the run ended at A~X 170 / A~Y 30 for
   every seed (30 being the K count), with 60 molecules holding `null`
   coordinates in `DATA/restart.json`.
3. With concrete `assocAngles`, the copy counters went wrong instead: `K(k)`
   went negative and `K(kd)` exceeded the 30 molecules present. The state
   change was counted under the relative state index (PR #84), and the
   reverse `Y -> X` direction hardly ran (PR #83) and ran at the forward rate
   (PR #84).

With all of them, the two states settle where they should. Equilibrium is
`A~Y / A~X = onRate3Dka / offRatekb = 1/5`, so A~Y averages 200/6 = 33.3. Mean
A~Y over the second half of the run, seeds 1-3: 31.0, 32.4, 33.3, and no
`null` coordinates in `DATA/restart.json`. Measured with the NaN fix on a tree
without PRs #83 and #84, the reverse direction almost never fired (198-200 of
200 A ended in Y) and `K(k)` went negative.
