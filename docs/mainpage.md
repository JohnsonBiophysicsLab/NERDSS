\mainpage
\tableofcontents

See \ref changelog for a list of changes.

---

\section quick Quickstart

\subsection Dependencies

- C++ compiler, one of the following
  - GCC >= 4.9.0 
  - Intel >= 16.0
  - Clang >= 3.3
- GSL >= 1.5
- CMake >= 3.0
- Make >= 4.2.1

\subsection build Building NERDSS

1. Create a build directory and move into it: `mkdir build; cd build`
2. Create the Makefile using CMake: `cmake ..`
3. Make the executable: `make`
  - Use the `-j $NUMCORES` flag for faster building, if you wish, where `$NUMCORES` is the number of cores your CPU has (probably 2 or 4).
4. The executable is now in the build directory, entitled `vector_rd_reweight_NOPBC`

\subsection run Running NERDSS

- See \ref input for information on input file syntax and format.
- Each molecule `X` in the `molecules` block is read from a mol file, by default `X.mol` in the working directory (see \ref molfiles for other places)
- The parameters file can be in any directory, as long as you pass the path to the binary as a command line argument.
- To run, simply type `./$BINARY -f $PARAMFILE`, where \$BINARY is the path to the binary and \$PARAMFILE is the path to the parameter file.

- Can use several command line flags to modify at runtime:

  | Flag | Description|
  |------|------------|
  | `-f <filename>`, `--parmfile <filename>`    | declares the parameter file (required unless restarting)                  |
  | `-s <n>`, `--seed <n>`                      | manually declare a seed for RNG, a whole number (optional)                |
  | `-r <filename>`, `--restart <filename>`     | declares the restart file (optional)                                      |
  | `-a <filename>`, `--add <filename>`         | declares an add file applied on restart (optional)                        |
  | `-c <filename>`, `--coordinate <filename>`  | declares a file of initial coordinates (optional)                         |
  | `-m <dir>`, `--mol-dir <dir>`               | read `X.mol` from `<dir>` instead of the working directory (optional)     |
  | `-M <dir>`, `--mol-dir-recursive <dir>`     | read `X.mol` from anywhere under `<dir>` (optional)                       |
  | `-v`, `--verbose`                           | more log output; `-vv` for the most (optional)                            |
  | `--debug-force-dissoc`                      | force dissociation to occur whenever possible (optional, for debug only)  |
  | `--debug-force-assoc`                       | force association to occur whenever possible (optional, for debug only)   |
  | `--print-system-info`                       | print the system information (optional)                                   |
  | `-h`, `--help`                              | list the flags and exit                                                   |

- Flags are read with `getopt_long`, so a long flag also takes its value after `=` (`--seed=12`) and short flags can be combined (`-vv`). An unknown flag, a flag missing its value, a seed that is not a whole number, or a stray argument stops the run with a message

\subsection molfiles Where mol files are read from

- For each molecule `X`, the first of these that applies gives its mol file:
  1. An entry for `X` in a `molFiles` block of the input file (below)
  2. `X.mol` in the `-m`/`--mol-dir` directory, or the one file named `X.mol` anywhere under the `-M`/`--mol-dir-recursive` directory
  3. `X.mol` in the working directory
- `-M`/`--mol-dir-recursive` stops with an error if `X.mol` is in more than one place under the directory; name the one to use in a `molFiles` block. It skips hidden directories (`.git`) and does not follow symbolic links to directories. The two flags cannot be combined
- A `molFiles` block maps molecule names to mol files, which may then have any name:

      start molFiles
          A = mols/clathrin_v2.mol      # relative to this input file
          B = /shared/models/AP2.mol
      end molFiles

  - A relative path is relative to the directory of the file the block is in, not the working directory, so a model folder can be run from anywhere with `-f`
  - One `name = path` per line; `#` starts a comment, so paths cannot contain `#`. Molecules without an entry are looked up as in steps 2 and 3
  - Put the block before `start reactions` (the top of the file is simplest). Every NERDSS version before this one ignores the block, so such a file still runs on them, provided `X.mol` is also in the working directory. Placed after the reactions block, it is read as a reaction by every version, this one included
  - An add file (`-a`) can have its own `molFiles` block for the molecules it adds; the flags apply to add files too. A restart without an add file reads no mol files

\subsection hist The complex histogram

- `DATA/histogram_complexes_time.dat` lists, at every `timeWrite`, how many complexes of each composition exist, one line per composition: `4	A: 2.` means four complexes made of two A molecules
- By default the composition counts molecule types only, so a phosphorylated `pip2(head~P)` and an unphosphorylated `pip2(head~U)` are one species, `pip2`. With `histogramMolStates = true` in the parameters block they are two, named as in `copy_numbers_time.dat`: `pip2(head~P): 2. pip2(head~U): 1.` and `pip2(head~P): 1. pip2(head~U): 2.` are then different complexes. The switch changes this file only; the trajectory and every other output are the same either way

\subsection restr Restarting a Simulation

- With each write of the restart file, the current state of the RNG will also be written as a binary file (called rng_state)
- To restart a simulation, make sure the trajectory, RNG state file, and restart file are in the same directory, then run `./$EXE -r $RESTART`, where `$EXE` is the executable and $RESTART is the restart file (`DATA/restart.json` by default; checkpoints go to `RESTARTS/restart<step>.json`)
- Restart files are JSON. A `.dat` restart file written by an earlier build is still read: the format is recognised from the file's contents, not its name
- `legacyRestartFormat = true` in the parameters block makes a run write the old `.dat` format instead (`DATA/restart.dat`, `RESTARTS/restart<step>.dat`). A restart keeps the format of the file it started from; an add file can set the parameter either way
- As of right now, some limitations
  - The restart will fail if the trajectory and rng_state files are not in the directory.
  - The restart blindly appends to the trajectory, i.e. it does not yet check to make sure the last written iteration matches the restart file.

\subsection examples Examples

Example simulations can be found in \ref examplesPage.

---

\section help Troubleshooting

\subsection genDebug General Debugging

- To do general debugging (e.g. with GDB/LLDB), cmake with the following command: `cmake .. -DCMAKE_BUILD_TYPE=Debug`, in order to activate the debug flags.
  - `-g`: generate debugging tables for GDB/LLDB
  - `-fsanitize=address`: Uses a sanitizer to give information about memory access issues, such as out-of-range element access
  - `-fno-omit-frame-pointer`: stores the stack frame pointers for functions (needed for `-fsanitize`)

\subsection buildHelp Build Errors

### Issue: FindGSL.cmake is not in CMAKE_MODULE_PATH

- CMake can have trouble finding GSL on machines on which it is installed in non-standard directories (e.g. supercomputers/clusters using `module` environments)
- Solution: Set the `GSL_ROOT_DIR` environment path to the root directory containing GSL.
  - For example, on MARCC type `export GSL_ROOT_DIR=/software/apps/gsl/2.5`

### Issue: Make has no recipe for all

- Full error: 

      make[1]: *** No rule to make target '/some/path'.  Stop.
      Makefile:83: recipe for target 'all' failed
      make: *** [all] Error 2

- Solution: Explicitly state the target executable: `make $EXEC`, where `$EXEC` is `nerdss` or `vector_rd_reweight_NOPBC` (depending on the version)

---

\section cite Citations

If you use NERDSS, please cite the software publication:
  - Varga, M.V.; ...; Johnson, M.E. *In preparation*.

\subsection FPR Algorithm Papers
  - Johnson, M. E. [Modeling the Self-Assembly of Protein Complexes through a Rigid-Body Rotational Reaction--Diffusion Algorithm.](https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.8b08339) *J. of Phys. Chem. B* 2018, **122**, 11771–11783.
  - Johnson, M.E., and G. Hummer. [Free propagator reweighting integrator for single-particle dynamics in reaction-diffusion models of heterogeneous protein-protein interactions systems.](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.4.031037) *Phys. Rev. X* **4**, 031037 (2014)
  - Yogurtcu, O.N., and M.E. Johnson. [Theory of bi-molecular association dynamics in 2D for accurate model and experimental parameterization of binding rates.](https://www.ncbi.nlm.nih.gov/pubmed/26328828) *J. Chem. Phys.* **143**, 084117 (2015)
  - Yogurtcu, O.N., and M.E. Johnson. [Cytoplasmic proteins can exploit membrane localization to trigger functional assembly.](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006031) *PLoS Comput. Biology* **14**, e1006031 (2018)
