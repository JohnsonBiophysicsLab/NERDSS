# Deterministic OpenMP production lane

## Outcome

This change adds an actual OpenMP execution lane to the production NERDSS
serial executable. It parallelizes independent 3D bimolecular candidate checks
and the per-molecule transform inside sufficiently large complex propagation,
while preserving the mutation order required by the existing reaction data
structures.

The fixed-seed production validation is stronger than a tolerance check. On the
measured system, all scientific output from the new serial and OpenMP builds is
byte-for-byte identical at 1, 2, 4, and 8 OpenMP threads. The only raw-file
difference is the expected wall-clock `CREATED` timestamp in the first line of
`PDB/0.pdb`; the PDB body and every trajectory, restart, observable, and system
file are identical.

The measured whole-program median speedup at eight physical-core workers is
1.074x for the supplied 10,000-molecule production workload. One to four
workers regress slightly because deterministic conflict scheduling and barriers
cost more than the available parallel pair work. This is a deliberately honest
production result: the 3.4-5.7x isolated-kernel POC scaling does not translate
directly to an application in which candidate discovery, conflict construction,
reaction execution, small-complex propagation, and output remain serial.

## Scope of the implementation

The production path now covers these requested functions and related POC lanes:

- `measure_separations_to_identify_possible_reactions()` retains the original
  cell and neighbor traversal order, materializes that ordered stream as
  candidate-pair records, and submits the batch to the OpenMP lane.
- `check_bimolecular_reactions()` remains the authoritative stateful evaluator.
  The new `check_bimolecular_reaction_pairs()` batch entry point calls it in
  deterministic conflict waves.
- `determine_3D_bimolecular_reaction_probability()` uses an extracted pure
  `calculate_3D_bimolecular_parameters()` prefix for rotational diffusion and
  reaction-radius calculation. The serial and OpenMP paths therefore execute
  the same arithmetic implementation.
- `get_distance()` uses an extracted pure `calculate_bimolecular_distance()`
  geometry calculation, then performs the existing ordered vector appends and
  complex-counter updates.
- `Complex::propagate()` uses static OpenMP work sharing for the molecule-member
  loops in its spherical, translation-only, and quaternion rotation branches.
- `EXEs/nerdss.cpp` now calls the shared separation function for ordinary
  production systems. Compartment systems retain the original inline traversal
  because they require `transmissionRxns`, which the shared API does not own.

## Why the outer loop cannot be a naive `omp parallel for`

The original candidate check is not a pure map. A successful check appends to
per-molecule vectors including `crossbase`, `mycrossint`, `crossrxn`, and
`probvec`; increments per-complex `ncross`; and, for specialized lanes, can
touch 2D lookup tables or global counters. Two candidate pairs that share a
molecule or complex would race, and even a lock-based implementation could
change append order. Changed append order can later change association-event
selection and the random-number stream.

The production implementation therefore separates discovery from evaluation:

1. Enumerate pairs in exactly the original cell/member/neighbor order.
2. Assign each pair to the earliest conflict wave that is later than the last
   wave used by either molecule or either complex.
3. Add a global ordered dependency for cases that are not independent 3D pairs:
   same-complex, 1D, 2D, implicit-lipid, ghosted/invalid, and exclude-volume
   cases.
4. Evaluate one wave with `omp for schedule(static)`.
5. Use the implicit barrier before starting the next wave.

Within a wave, no two pairs share a molecule or complex. Consequently, each
worker owns every vector and counter it mutates. Across waves, all mutations to
a particular molecule or complex occur in the same order as the serial
enumeration. Specialized cases retain their global serial ordering.

The OpenMP build falls back to the original serial traversal for a timestep if
an active template binds to an implicit surface. This preserves the original
interleaving of implicit-lipid checks with explicit pair checks. Compartment
models also retain the serial production path.

## Pure arithmetic extraction

The POC demonstrated that distance, diffusion-radius, and quaternion arithmetic
are data-parallel and bitwise stable when each record keeps its original
floating-point operation order. Production now makes that property explicit:

- `calculate_bimolecular_distance()` only reads molecule, complex, reaction,
  and membrane geometry and returns separation, distance, and cutoff status.
- `get_distance()` consumes that result and performs the legacy mutations.
- `calculate_3D_bimolecular_parameters()` computes the two rotational diffusion
  contributions and `Rmax` in the original statement order.
- `determine_3D_bimolecular_reaction_probability()` consumes those values and
  continues through the existing probability/history logic.

There are no OpenMP reductions and no reassociation of floating-point
expressions. A candidate is evaluated from start to finish by one worker.

## Propagation lane

`Complex::propagate()` constructs the motion and quaternion exactly once, as it
did before. Its member transform is then a static index-based OpenMP loop. Each
iteration writes one distinct molecule and its distinct interface list. The
implicit loop barrier completes all molecule transforms before
`update_properties()` reads the complex.

The default threshold is 32 complex members. Small complexes stay serial to
avoid repeatedly creating worker teams for only a handful of coordinates. The
threshold can be changed for workload tuning or validation:

```bash
NERDSS_OMP_MIN_COMPLEX_MEMBERS=64 OMP_NUM_THREADS=8 bin/nerdss_omp \
  -f parm.inp -s 20260806
```

The benchmark includes a forced-path smoke test with the threshold set to one,
so the OpenMP propagation branch is exercised even though the homotrimer
performance case contains only small complexes.

## Build integration

The Makefile has independent serial, OpenMP, and MPI object directories. This
is necessary because make cannot otherwise distinguish an object compiled with
`NERDSS_USE_OPENMP` from one compiled without it.

```bash
make serial   # obj_serial -> bin/nerdss
make omp      # obj_omp, -fopenmp, NERDSS_USE_OPENMP -> bin/nerdss_omp
make mpi      # obj_mpi -> bin/nerdss_mpi
```

Invoke these as separate make commands. The benchmark runner does so. `make
clean` removes all four historical/current object directories and `bin`.

The OpenMP linker flag is present in both compilation and linking, and the
OpenMP preprocessor definition is applied while compiling the executable source
as well as library sources.

## Runtime controls

- `OMP_NUM_THREADS`: worker count.
- `OMP_DYNAMIC=FALSE`: prevents the runtime from silently changing team size.
- `OMP_WAIT_POLICY=active`: used by the benchmark to reduce wake-up latency.
- `OMP_PROC_BIND=spread`, `OMP_PLACES=cores`: spreads workers across the eight
  core places exposed by WSL.
- `NERDSS_OMP_MIN_PAIRS`: minimum candidate batch size before conflict waves
  and a worker team are used; default 64.
- `NERDSS_OMP_MIN_COMPLEX_MEMBERS`: minimum complex size for parallel member
  propagation; default 32.

Both NERDSS thresholds are read once per process. Values below one are clamped
to one.

## Correctness methodology

The benchmark uses `benchmarks/openmp_production/parm.inp`, a fixed-density 3D
homotrimer case with 10,000 molecules, 200 iterations, and RNG seed 20260806.
For each executable/core-count run, the harness creates an isolated directory
and hashes every file under `DATA`, `PDB`, and `RESTARTS`.

Two manifests are compared:

- the raw manifest hashes every byte;
- the simulation-state manifest hashes every non-PDB file unchanged and hashes
  each PDB after omitting only line 1, which contains the wall-clock creation
  timestamp.

There are 23 output files in each run. The serial/OpenMP comparison covered all
four worker counts and three repetitions per count. Every simulation-state
manifest matched. Raw manifests differed only by the two diff entries (old and
new hash) for `PDB/0.pdb`; inspection showed that its sole content difference
was the timestamp in the title line.

Three additional checks were performed:

1. The pre-change production executable from commit `44efa1c` was built in an
   isolated worktree and run with the same input and seed. Its `DATA` and
   `RESTARTS` trees exactly matched the new serial build, and its normalized PDB
   hash matched. This verifies that helper extraction and moving the ordinary
   traversal into the shared function did not alter serial simulation state.
2. OpenMP propagation was forced with
   `NERDSS_OMP_MIN_COMPLEX_MEMBERS=1` and
   `parm_propagation_smoke.inp`. Serial, OpenMP-1, and OpenMP-8 outputs matched.
3. Both executables were built from separate object trees and completed without
   compiler or linker errors; `git diff --check` also passed.

This establishes exact equality for the tested compiler, flags, platform,
input, and seed. Bitwise reproducibility should be revalidated after compiler,
math-library, optimization-flag, scheduling, or arithmetic changes.

## Performance methodology and result

Measured environment:

- Intel Core i5-13450HX host;
- WSL reports 8 cores, 16 logical CPUs, and 2 threads per core;
- Linux 6.18.33.2-microsoft-standard-WSL2;
- GCC 15.2.0;
- `-O3 -std=c++0x`, plus `-fopenmp` for the OpenMP build;
- three fresh process runs at each worker count;
- `/usr/bin/time` wall-clock measurement around the complete executable.

Median speedup is relative to the 10.56-second serial median:

| Variant | Workers | Trials (s) | Median (s) | Speedup | Efficiency |
|---|---:|---|---:|---:|---:|
| serial | 1 | 10.66, 10.40, 10.56 | 10.56 | 1.000x | 100.0% |
| OpenMP | 1 | 11.06, 11.29, 11.59 | 11.29 | 0.935x | 93.5% |
| OpenMP | 2 | 11.44, 11.20, 11.65 | 11.44 | 0.923x | 46.2% |
| OpenMP | 4 | 10.76, 11.46, 12.02 | 11.46 | 0.921x | 23.0% |
| OpenMP | 8 | 19.87, 9.83, 8.76 | 9.83 | 1.074x | 13.4% |

The 19.87-second eight-worker observation is retained in the raw data. It is a
clear WSL/hybrid-host scheduling outlier relative to the other eight-worker
runs, but removing it would conceal the variability a user can encounter.

The scaling shape shows that correctness-preserving production scheduling is
not free. At low worker counts, candidate materialization, wave construction,
barriers, and serial application work dominate. Eight workers finally recover
that overhead and reduce median total time by about 6.9%. The best observed
eight-worker run is 1.205x faster than the serial median, but the median is the
reported result.

This workload does not exercise parallel propagation at its default threshold,
because a homotrimer has at most three members. Systems with large polymers or
assemblies may benefit independently from the propagation lane; they need a
representative benchmark before changing the threshold.

## Reproduction

From the repository root:

```bash
BENCH_REPS=3 BENCH_THREADS="1 2 4 8" \
  benchmarks/openmp_production/run_benchmark.sh
```

The runner prints its timestamped result directory. It contains raw logs,
timings, environment metadata, SHA-256 manifests, and per-run diffs. Generated
result directories are ignored by git; the measured raw table is retained in
`benchmark_results_i5_13450hx.txt`.

To force the propagation correctness lane:

```bash
NERDSS_OMP_MIN_COMPLEX_MEMBERS=1 \
BENCH_INPUT=parm_propagation_smoke.inp \
BENCH_REPS=1 BENCH_THREADS="1 8" \
  benchmarks/openmp_production/run_benchmark.sh
```

## Known limits and next optimization target

- Parallel evaluation is intentionally limited to conflict-free production
  work. It does not attempt concurrent mutation of 2D tables, implicit-lipid
  state, exclude-volume state, or compartment transmission state.
- Candidate discovery and greedy wave construction are serial and allocate a
  batch every timestep. Reusing capacity and replacing `vector<vector<Pair>>`
  with flat wave offsets are likely low-risk overhead reductions.
- Each wave has a barrier. A dependency-aware task graph could reduce barrier
  cost, but reproducible append order would require careful validation.
- Parallel propagation is inside one large complex, not across many small
  complexes. A higher-level complex batch could help assembly-poor workloads if
  its surrounding state updates can be proven independent.
- The production result should guide deployment: use the serial binary for
  small or highly conflicting systems, and benchmark representative large 3D
  systems before selecting an OpenMP core count.
