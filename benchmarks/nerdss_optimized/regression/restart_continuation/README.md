# Exact continuation from a checkpoint

A restart from a checkpoint, with the RNG state restored, must continue the run
that wrote the checkpoint bit for bit. `check.sh` tests exactly that:

```bash
./check.sh <path-to>/nerdss [work-dir]
```

For each case in `continuation_cases.tsv` it runs the model for 2N steps with
`checkPoint = N` and `RNGwrite = true`, restarts in a clean directory from
`RESTARTS/restart<N>.dat` with `RESTARTS/rng_state<N>` copied in as `rng_state`,
and requires

- the two final `DATA/restart.dat` files to be identical apart from line 2
  (`numItr`) and line 4 (`currSimTime`, which the restart reaches by a different
  but equivalent expression), and
- every record the restart wrote to `DATA/*_time.dat` to equal the uninterrupted
  run's record for the same time. That covers what the restart file does not
  hold: the copy numbers `init_counterCopyNums()` recounts at a restart, for
  one, are compared through the record the restart writes as it starts.

Exit status 0 means every case continued exactly. The whole table takes well
under a minute, and needs bash, awk and python3.

This is stronger than a read->write round trip of a restart file, which shows
that the reader consumes what the writer wrote but says nothing about state the
file never held.

## What the cases found

On 01bf678 every implicit-lipid model and every box model with a creation
reaction diverged. Four defects were responsible, each fixed in its own commit;
each case in the table is the first to expose one of them.

| Defect | Cases that expose it |
| --- | --- |
| The implicit lipid's single molecule holds `trajStatus = propagated` from its first step on, because the end-of-step reset skips it. `trajStatus` is not in the file, so a restart propagated it again, drawing random numbers the uninterrupted run did not. | `sphere`, `implicit_lipid`, `mem_localization_IL` |
| The 2D binding table for the implicit lipid is built entry by entry from the free-lipid count at the step each entry is first needed. A restart rebuilt it from the count at the restart. | `two_site` |
| `numberOfProteinEachState`, the count of lipid-binding interfaces at step 0 that new table entries use, was recounted from the molecules at the restart. | `two_site_depleting` |
| `read_restart()` left `waterBox.xLeft` and `xRight` at zero, and a molecule created in a box is placed at `x = xLeft + (xRight - xLeft) * rand`: after a restart every created molecule landed on `x = 0`. | `create_destroy`, `pucadyil` |

`rev_3D` and `compartment` continued exactly before any of this and are kept as
controls.

Measured with the case table's seeds, "diverged" meaning a differing final
`restart.dat`:

| Build | Controls | trajStatus cases | `two_site` | `two_site_depleting` | creation cases |
| --- | --- | --- | --- | --- | --- |
| 01bf678 | exact | diverged | diverged | diverged | diverged |
| + trajStatus | exact | exact | diverged | diverged | diverged |
| + 2D table | exact | exact | exact | diverged | diverged |
| + protein counts | exact | exact | exact | exact | diverged |
| + water box bounds | exact | exact | exact | exact | exact |

## Why two models of its own

No sample model exposes the two table defects within a short run. Only two
bind a second lipid in 2D in one: `gagsphere` builds its table after the
checkpoint, at the same step and from the same lipid count in both runs, and
`secretion_catalyzed_assembly/membrane/kd10`, which builds it before, never
draws a random number between the two runs' binding probabilities. Both
continue exactly once the trajStatus fix is in.

- `two_site`: C has two sites that bind the lipid, so a C bound by one binds by
  the other in 2D, and fast unbinding keeps it cycling. The table is consulted
  every step, and a wrong entry reaches the trajectory within the run.
- `two_site_depleting`: the same C with fewer lipids than binding sites, so an
  entry's block distance is set by the protein count, and with C destroyed over
  time, so the count at a checkpoint differs from the count at step 0. With seed
  777 the table is still empty at the checkpoint at step 10, so it is first
  built after the restart; a wrong count shows in the reweighting survival
  probabilities of the molecules that then attempt 2D binding.
