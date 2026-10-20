# B8 on the fiber model: master against nerdss-optimized-1.3

Measurement behind the B8.3 fix (`1a5d164a`, "Exempt fiber pairs from the
rMaxLimit pre-cut"), the grid sizing that followed from it (`22828ec3`), and the master-versus-branch comparison of B8 as a whole
(B8.2, the pro2 exclusion block's missing fiber branch, `07d2b8d5`; B8.3,
the pre-cut). Run on 2026-10-14. Builds: `master` at `c507ae30`; the 1.3
branch at `5f56ede1` (before), `1a5d164a` (fixed) and `22828ec3` (grid); Apple clang, `make
serial`, `-O3`.

## The defect (B8.3)

The centre-of-mass pre-cut in front of `check_bimolecular_reactions()` drops
any pair further apart than `rMaxLimit`, exempting only pairs with both
complexes on the membrane. `set_rMaxLimit()` estimates every reach as
3 sqrt(6 Dtot dt) with Dtot averaged over three axes, while the 1D reaction
and exclusion rules reach 4 sqrt(2 (D1.x + D2.x) dt), 4/3 as far when
D.y = D.z = 0. On this model `rMaxLimit` is 62 nm against a 1D reach of
82 nm (free site against free site) and 71.3 nm (bound site against free
site), so fiber pairs between the two radii, which `master` evaluates, were
never evaluated on the branch. The fix exempts fiber pairs the way surface
pairs are exempted.

## The model

`model/`: the two-species fiber model of the branch report's appendix
(`reference_2026-10-09/fiber_b8/model_pro2_order`): a 2000 nm fiber along x
with 50 `A` and 200 `B` sites, all `isPromoter = true`, `D = [100, 0, 0]`,
dt = 1 us; `A(l) + B(r) -> A(l!1).B(r!1)` irreversibly with sigma = 2 nm and
`excludeVolumeBound = true`. All 50 `A` bind within the first 30 steps and
the 150 free `B` keep sliding past the 50 bound sites. 300 000 steps, seeds
20260811 to 20260814, trajectory every 20 steps. `model_kin/parms.inp` is
the same model for 200 steps with copy numbers every step (association
kinetics, seeds 20260901 to 20260948).

## Builds and what each isolates

| Build | Source | Purpose |
| --- | --- | --- |
| `master` | `c507ae30` | baseline |
| `b13_before` | `5f56ede1` | branch with B8.2, before B8.3 |
| `b13_fix` | `1a5d164a` | branch with B8.2 and B8.3 (the shipped state) |
| `b13_fix_legacy_nob82` | `1a5d164a` + `legacy_stream.patch` + `revert_B8_2.patch` | master's random stream and grid, B8.2 off, B8.3 on: must equal `master` byte for byte |
| `b13_before_legacy` | `5f56ede1` + `legacy_stream.patch` | B8.2 alone under master's stream |
| `b13_fix_legacy` | `1a5d164a` + `legacy_stream.patch` | B8.2 + B8.3 under master's stream |
| `b13_fix_count` | `1a5d164a` + `scripts/count_B8_3.patch` | counts, for every fiber pair, whether it lies beyond `rMaxLimit` and whether the 1D rules admit it |
| `b13_grid` | `22828ec3` | `b13_fix` plus `set_rMaxLimit()` sized by the 1D reach for promoter pairs (section below) |
| `b13_grid_count` | `22828ec3` + `scripts/count_B8_3.patch` | the same, instrumented |

The two patches are in `reference_2026-10-09/` (`scripts/legacy_stream_patch.py`,
`fiber_b8/scripts/revert_B8_2.patch`).

## Results

**Isolation.** `b13_fix_legacy_nob82` reproduces `master`'s trajectory byte
for byte on 4 of 4 seeds over 300 000 steps (`runs/*/seed*.trajectory.md5`).
Those two builds differ only in B8.3, so the fix restores exactly the pairs
`master` evaluates. Without the fix (the reference's `optlegacy_nob8` runs)
the trajectories differ from the first frame.

**Realized overlaps** (A-site / B-site pairs closer than 0.9 sigma, counted
on every 20th step of a 300 000-step run; a bound pair sits at exactly sigma):

| Build | per run, 4 seeds | mean +- SEM |
| --- | --- | ---: |
| `master` | 237, 220, 215, 228 | 225.0 +- 4.8 |
| `b13_before` (B8.2 only) | 14, 16, 15, 17 | 15.5 +- 0.6 |
| `b13_fix` (B8.2 + B8.3) | 17, 9, 15, 9 | 12.5 +- 2.1 |
| `b13_grid` (B8.2 + B8.3 + grid by the 1D reach) | 11, 12, 9, 14 | 11.5 +- 1.0 |
| `b13_fix_legacy_nob82` (= `master`) | 237, 220, 215, 228 | 225.0 +- 4.8 |
| `b13_before_legacy` | 11, 10, 8, 13 | 10.5 +- 1.0 |
| `b13_fix_legacy` | 9, 11, 11, 9 | 10.0 +- 0.6 |

The 15 to 20x reduction is B8.2; B8.3 is within seed scatter on this
observable, as the branch report predicted. The legacy-stream rows
reproduce the reference's numbers exactly.

**What the pre-cut dropped** (`b13_fix_count`, 4 seeds, `runs/b13_fix_count/*.counters.txt`):

| | pair evaluations | of which beyond `rMaxLimit` |
| --- | ---: | ---: |
| fiber pairs offered by the cell list | 1 056 246 044 | 368 654 892 (34.9%) |
| admitted by the 1D reaction rule (free A vs free B, first ~30 steps) | 21 762 | 5 156 (23.7%) |
| admitted by the 1D exclusion rule (bound site vs free site) | 1 032 727 583 | 75 799 864 (7.3%) |

Before the fix, a quarter of the pairs the 1D reaction rule admits and 7% of
the pairs the exclusion rule admits were never reached. All of them sit
between 62 and 82 nm, 3 to 4 standard deviations of the pair's relative
step, which is why the realized effect is small.

**Association kinetics** (48 seeds, copy numbers every step): `master`,
`b13_before` and `b13_fix` agree within SEM at every step; at step 10 the
mean bond count is 38.9, 38.8 and 39.6 of 50.

**Run time** (100 000 steps, one run at a time, 2 repeats): `master` 29.4 s,
`b13_before` 24.9 s, `b13_fix` 25.1 s, `b13_grid` 27.0 s. Exempting fiber
pairs costs about 1% on this model, where every pair is a fiber pair; the
wider cells of `b13_grid` cost a further 8% (the cell list offers 33% more
fiber pairs per run, 352 million against 264 million).

**Standard suite** (`suite_*.tsv`, the 13 cases of `cases.tsv`): byte
identical between `b13_before` and `b13_fix`, and between `b13_fix` and
`b13_grid`. No case in the tree sets `isPromoter`. The timing columns of the
second file are not comparable: that suite ran alongside the fiber batch.

## Figures

- `figs/fig1_b8_master_vs_1.3.png`: overlaps per build, separation
  distribution next to the exclusion zone, association kinetics, run time.
- `figs/fig2_b8_3_precut.png`: where the admitted fiber pairs sit relative to
  `rMaxLimit`, and the share the pre-cut dropped.
- `figs/fig3_grid_1d_reach.png`: what the cell list offers with the grid
  sized by the 3D estimate and by the 1D reach.

## The cell grid sized by the 1D reach (`22828ec3`)

Panel A of figure 2 (and of figure 3) shows a second consequence of the
same 3D estimate: the cell grid is `rMaxLimit` wide with a +-1 stencil, so
a fiber pair further apart than one cell width (62.5 nm here) is offered to
the pair search with a probability falling linearly to zero at two widths,
while the 1D rules reach 82 nm. That holds on `master` and on the branch
alike, and the B8.3 exemption cannot help with it: a pair the cell list
never offers never reaches the pre-cut.

`set_rMaxLimit()` now takes, for a reaction whose two reactant templates are
both promoters, the larger of the 3D estimate and the 1D reach
4 sqrt(2 (D1.x + D2.x) dt) + sigma. On this model `rMaxLimit` becomes 82 nm
and the grid 24 cells of 83.3 nm instead of 32 of 62.5 nm. Instrumented
(`b13_grid_count`, 4 seeds):

| | `b13_fix_count` (62 nm) | `b13_grid_count` (82 nm) |
| --- | ---: | ---: |
| fiber pairs offered by the cell list, per run | 264 million | 352 million |
| admitted by the 1D reaction rule, beyond `rMaxLimit` | 5 156 of 21 762 (23.7%) | 0 of 21 962 |
| admitted by the 1D exclusion rule, beyond `rMaxLimit` | 75 799 864 of 1 032 727 583 (7.3%) | 0 of 1 037 576 221 |

Figure 3 shows the offered histogram flat out to the 1D reach with the new
grid, where before it had already lost a quarter of the pairs at 82 nm.
Realized overlaps are unchanged within scatter (11.5 +- 1.0 against
12.5 +- 2.1), as expected for pairs 3 to 4 standard deviations apart. The
fiber exemption in the pre-cut (`1a5d164a`) is kept: it still covers a
fiber pair whose templates are not both promoters, and it is now redundant
rather than wrong for the rest.

## Reproducing

```bash
export FIBER_SCRATCH=/some/scratch            # holds build/<name>/bin/nerdss and fiber/
cp -r model $FIBER_SCRATCH/fiber/model; cp -r model_kin $FIBER_SCRATCH/fiber/model_kin
scripts/run_one.sh b13_fix 20260811            # 300 000 steps, overlaps + md5 + wall time
scripts/run_kin.sh b13_fix 20260901            # 200 steps, copy numbers every step
scripts/run_time.sh b13_fix 1                  # 100 000 steps, wall time only
python3 scripts/aggregate.py $FIBER_SCRATCH/fiber results.json
python3 scripts/make_figures.py results.json figs   # needs matplotlib
```
