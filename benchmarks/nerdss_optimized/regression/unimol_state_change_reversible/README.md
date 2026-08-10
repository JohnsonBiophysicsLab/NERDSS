# Reversible unimolecular state-change reproducer

This is the validation case for the reverse (product-matching) path of
`find_which_state_change_reaction()`. The model itself lives in
[`sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible`](../../../../sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible);
`check.sh` runs it against a given binary and compares both state changes to the
equilibrium fraction they must reach.

```bash
./check.sh <path-to>/nerdss
```

Exit status 0 means both reversible state changes settled on `kf/(kf+kb)`.

## Why a new model was needed

No input under `sample_inputs/` declares a unimolecular state change with `<->`.
`michaelis_menten`, `unimolecular_reverse`, `auto_phos` and `enzyme` all write two
separate `->` reactions, which leaves `conjBackRxnIndex == -1` and never reaches
the reverse path. So nothing that existed could exercise it, in either direction.

The model also arranges for the forward and back index spaces to *differ*, which
a single reversible pair cannot do (it would put the pair at forward 0 / back 0,
where confusing the two lists is invisible). An irreversible rate-0 reaction at
forward index 1 creates no `BackRxn`, so the two state changes sit at forward
2 and 3 but back 1 and 2.

## What the model measures

For `U <-> P` with rates `kf` and `kb`, the equilibrium fraction in P is
`kf/(kf+kb)`, here `1000/(1000+3000) = 0.25`. Relaxation time is `1/(kf+kb)` =
250 us and the run is 2000 us, so `check.sh` time-averages the second half.

A build whose reverse path never fires makes P absorbing, so P climbs toward 1.00
instead of settling at 0.25. That is a factor-of-four separation, not a subtle
statistical one.

## Results

Three builds, seed 20260810, time-averaged over the equilibrated half:

| build | A(ser~P) | B(thr~P) | expected | verdict |
| --- | --- | --- | --- | --- |
| `nerdss-optimized` at `dbb5cbf` (`2b23b2da`) | - | - | 0.25 | SIGSEGV during reaction display |
| plus `BackRxn::display()` fix only (`20c79d91`) | 0.765 | 0.753 | 0.25 | FAIL, P is absorbing |
| all four fixes (`1c044a13`) | 0.229 | 0.269 | 0.25 | PASS |

The middle row is the point of interest: it isolates the index-assignment defect
from the display crash. With the display crash fixed but the index assignments
untouched, the model runs and the reverse direction simply never happens - P
decays one-way from 200 U toward 200 P. The 0.765 is not an equilibrium, it is
where a one-way decay happens to be after 8 time constants.

Across seeds 20260810, 1, 2, 3 and 7 the corrected build gives P fractions between
0.224 and 0.269 against a predicted 0.25, with a binomial spread of 0.031 on 200
copies, so the agreement is not seed-specific.

## Three defects had to be fixed for this model to run

1. **`find_which_state_change_reaction()`, reverse path.** `rateIndex` was
   computed correctly in the loop and then overwritten with a *reaction* index,
   and `rxnIndex` was never assigned. Both callers require `rxnIndex != -1`, so
   every reverse unimolecular state change was discarded. This is the defect
   recorded but deliberately left alone under issues #8-#12.
2. **The same branch read `backRxns[conjBackRxnIndex]` without testing for -1**,
   so an irreversible state change whose product state matched read out of
   bounds.
3. **`stateChangeRxns` violated its own documented ordering.**
   `class_MolTemplate.hpp` declares the pair as `(forward, back)`, and
   `find_which_state_change_reaction()` indexes `forwardRxns[rxnItr.first]`
   accordingly, but `populate_reaction_lists.cpp` stored `(back, forward)` for
   the state on the product side. With one reversible pair the two happen to be
   equal, which is why this never showed up; with the index spaces offset it
   selects the wrong reaction. Fixed at the population site, restoring the
   documented invariant.

A fourth defect blocks the model before the simulation starts and is fixed here
too, since without it the model cannot run at all:

4. **`BackRxn::display()` read `rate.otherIfaceLists[1]`** whenever the list was
   non-empty. `otherIfaceLists` holds one entry per reactant, so a unimolecular
   back reaction has one entry and element 1 is out of bounds. It crashed 4 runs
   in 5 with the same seed, depending on heap layout.
   `ForwardRxn::display()` avoids it by skipping the block for `uniMolStateChange`
   entirely. Display-only, so it cannot affect results; and the 13-case bitwise
   suite confirms it does not.
