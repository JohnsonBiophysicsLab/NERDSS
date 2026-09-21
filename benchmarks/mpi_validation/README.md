# Serial vs MPI validation harness

`FINDINGS.md` is the write-up. This directory holds what produced it.

| file | role |
|---|---|
| `screen2.sh` | feasibility screen; refuses spherical inputs by reading the `.inp` |
| `bench2.sh` | runs the benchmark, keeps the pathway histograms |
| `cases.txt` | the selected cases: `name|input|nItr` |
| `analyze.py` | species means, Welch + Holm + Benjamini-Hochberg |
| `analyze_pathway.py` | assembly pathway: permutation test + per-composition tests |
| `check_integrity.py` | mass conservation and compositions the network cannot make |
| `standalone_cmp.sh` | merged vs upstream `origin/mpi`, same input and seeds |
| `stress.sh` | longer-run replicate of that comparison |
| `bench2_species.tsv`, `standalone_cmp.tsv` | the collected results |
| `headers/` | species column names per case |

The raw per-run histograms (~24 MB) are not committed; re-running `bench2.sh`
regenerates them.

Note: these scripts contain absolute paths to the scratchpad they were written
in (`SP=...`) and to prebuilt binaries in `bench_bins/`. Point `SP` at a working
directory and put `nerdss_serial` / `nerdss_mpi` in `$SP/bench_bins` before
reuse. They are recorded as the exact procedure that produced `FINDINGS.md`,
not as a turnkey tool.
