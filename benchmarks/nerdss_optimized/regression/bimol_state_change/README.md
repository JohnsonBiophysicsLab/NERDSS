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

## Result: the path cannot be validated end to end on this codebase

The model does not run to completion on **any** build, master or optimized,
because the surrounding bimolecular state-change machinery has its own defects
that are outside the scope of issues #8-#12:

1. `set_rMaxLimit()` only inspects `ReactionType::bimolecular`. A model whose
   only bimolecular reaction is a state change therefore keeps `rMaxLimit == 0`,
   and `SimulVolume::Dimensions` divides the box length by it, yielding a
   negative cell count and an endless `CELL PAIR MAX EXCEEDED` rescale. The
   `K(kd) + K(kd)` association at the end of `parms.inp` is scaffolding that
   works around this; it is listed last so the state-change reaction keeps
   indices 1 and 0.
2. With that worked around, executing the state change produces NaN coordinates
   within a few hundred iterations:
   `ERROR: [Before overlap checking] nan in the coordinates of molecule 22 at
   iteration 139`. Setting the state-change rate to zero removes the NaN and the
   run completes, which localizes the fault to the state-change execution rather
   than to the scaffolding.
3. Supplying concrete `assocAngles` instead of `nan` avoids the NaN but corrupts
   the copy counters instead: `K(k)` goes negative (-54) and `K(kd)` exceeds the
   30 molecules present (114), while no reaction is recorded as having fired.

Measured on the three builds, same seed:

| build | outcome |
| --- | --- |
| master (`260f6e2`) | NaN, molecule 22, iteration 139 |
| optimized, result-preserving subset | NaN, molecule 22, iteration 139 (identical) |
| optimized, full branch | NaN, molecule 0, iteration 210 (same fault, different random stream) |

The result-preserving build failing at exactly the same molecule and iteration as
master is the useful signal here: the issue #8 changes are inert on this path,
which is consistent with the corrected branch being dead on master. Making the
bimolecular state change actually work needs the three defects above fixed first,
each with its own validation.
