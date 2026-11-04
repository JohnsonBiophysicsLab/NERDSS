# Repository cleanup: redundant cache and generated files

Date 2026-11-04, branch `nerdss-optimized-1.4` at `9200537`.

This note records what was removed, why it is safe, and what else is
proposed but deliberately **not** removed yet because it needs an owner's call.
Removal from the tip does not shrink the clone: every file stays in history and
can be restored with `git checkout 9200537 -- <path>`. Shrinking the history
would need a rewrite (`git filter-repo`) of a shared branch, which is not
proposed.

## Tier 0 - removed in this change: `.cache/clangd/index/` (442 files, 1.05 MB)

**What it is.** The background index that the `clangd` language server writes
for an editor (VS Code, Neovim, CLion with clangd, ...). Each `*.idx` file is a
binary RIFF serialization of the symbols of one source file, named
`<source>.<16-hex hash of the source's absolute path>.idx`. clangd rebuilds it
on its own the first time a project is opened, so a checked-in copy is never
needed.

**How it got here.** It arrived in `22828ec` (2026-10-14), the bulk import
described in [`REPLACEMENT_REPORT.md`](../REPLACEMENT_REPORT.md): the import was
an `rsync --delete` of a local working directory, and the report states
"Did not skip other source hidden paths. The source includes
`.cache/clangd/...`, so those files were copied into the working tree." It was
never intended to be tracked; the `.gitignore` simply did not mention it.

**Why it is safe to remove** - each point was checked:

| Check | Result |
| --- | --- |
| Referenced by the build? | No. `Makefile`, `CMakeLists.txt` and `Dockerfile` never mention `.cache` or `clangd` (the one `CACHE` hit in `CMakeLists.txt` is the CMake cache-variable keyword). |
| Referenced by any script, notebook or doc? | No. `git grep` finds `.cache`/`clangd` only in `REPLACEMENT_REPORT.md` (which lists it as an accidental copy) and inside the `.idx` files themselves. |
| Is it current? | No. It indexes 239 distinct source names; **51 of them no longer exist** in `src/`, `include/` or `EXEs/` - e.g. `class_Coord.cpp`, `class_Vector.cpp`, `associate.cpp`, `reflect_traj_complex_rad_rot.cpp`, `math_functions.cpp` - files removed or renamed by the 1.2-1.4 work (Coord+Vector -> `Vec3D`, the box/sphere splits). |
| Is it one machine's index? | No. Every source has **two** `.idx` files with different path hashes, i.e. indexes from two checkouts at different absolute paths were mixed. Neither matches any other clone, so clangd ignores them and re-indexes anyway. |
| Does anything at run time read it? | No. NERDSS never opens files under `.cache/` (`git grep` over `src/`, `include/`, `EXEs/` finds no such path). Checked 2026-11-04: a clean `make serial -j4` in a fresh worktree of this commit (no `.cache/`) succeeds and the binary runs `rev_3D`; the 19-case `cases.tsv` + `coverage_cases.tsv` suite behind `benchmarks/nerdss_optimized/profile_suite_2026-11-04/` was run after the directory was deleted. |

**Guard against recurrence.** `.gitignore` now ignores `/.cache/` and
`.ccls-cache/`. A root `.clangd` *configuration* file, if anyone adds one, is
not ignored - that one is meant to be shared.

## Tier 1 - proposed: build products and editor leftovers (safe, not yet removed)

These are the same kind of thing as the clangd index - generated, or editor
backups - and nothing reads them. They are listed rather than removed because
some of them sit in `sample_inputs/`, which users browse, and the owner may
want to keep a particular one.

| Item | Files | Size | Evidence it is redundant |
| --- | --- | --- | --- |
| Compiled `nerdss` executables inside `sample_inputs/`: `nerdss_clathrin`, `nerdss_clock_model`, `nerdss_IL` (x2), `nerdss_EL`, `nerdss_SlowDsol`, `nerdss_michaelis_menten`, `nerdss_unimolecular_reverse`, and five named `nerdss` (`Gag/`, `create_destroy/`, `unimolecular_reverse/`, `clathrin_coat/puckered_clat.dir/`, `membrane_rev_localization/SmallBox/SlowDsol/`) | 13 | 263 MB | Old builds: eight are the same 32 MB debug ELF (one build ID), two other x86-64 ELFs, three x86-64 Mach-O. They do not match the current source, the Mach-O ones do not run on Linux or Apple silicon without Rosetta, and `make serial` rebuilds the program. Only the matching `OUTPUT` logs and `REPLACEMENT_REPORT.md` mention them. |
| `benchmarks/nerdss_optimized/rng_quality/rng_quality` | 1 | 35 KB | arm64 Mach-O build of `rng_quality.cpp`; `run.sh` next to it rebuilds it. |
| Editor backups `*~` in `sample_inputs/` | 9 | 30 KB | Already matched by `.gitignore` (`*~`) but tracked from before. Eight have a newer sibling without `~` (2-57 lines different). `compartment/parms.inp~` has **no** sibling - check whether `clath_compartment.inp` superseded it before deleting. |
| Jupyter `.ipynb_checkpoints/` in `sample_inputs/` | 9 | 27 KB | Already matched by `.gitignore`; autosave copies of files that exist next to them. |
| `REPLACEMENT_REPORT.md` at the root | 1 | 30 KB | Describes itself as "a local documentation artifact ... not part of the source directory". Consider moving it to `docs/` as history rather than deleting it. |

Procedure, if approved: one commit per row, `git rm` only; then
`make serial && benchmarks/nerdss_optimized/run_suite.sh bin/nerdss after` and
`compare_suites.sh` against a run made before the commit. None of these files is
an input to any case, so every case must stay byte-identical.

## Tier 2 - keep (looked at, not redundant)

| Item | Why it stays |
| --- | --- |
| `sample_inputs/VALIDATE_SUITE/*/OUTPUT` (9 files, 60 MB) | stdout logs of the 2022 reference runs. Generated, but they are the only record of the seed, start date and command behind the reference figures in those folders. If space matters, compress them (`xz` reduces such logs >10x) rather than delete. |
| `docs/html/` (831 files, 9.8 MB) | Doxygen output, but `docs/index.html` redirects to `html/index.html`, so it is the published API documentation. Regenerating it in CI instead would be a separate change. |
| `NERDSS_USER_GUIDE_MASTER.docx` / `.pdf` | Source of the user guide. A Markdown/HTML edition is proposed on the website branch; the `.docx` remains the archival copy until that edition is adopted. |
