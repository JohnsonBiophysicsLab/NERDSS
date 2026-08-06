# OpenMP reaction-checking and propagation proof of concept

## Executive summary

This proof of concept replaces the CUDA execution backend in `gpu_poc` with
OpenMP for the same flattened reaction-checking and propagation calculations.
It keeps the serial and parallel arithmetic in the same C++ executable, uses
the same input arrays, and validates every output after each timed run.

On the measured WSL host, eight physical-core workers provide approximately
4.7-5.0x reaction-checking speedup and 3.2-3.6x propagation speedup for the
largest batch in repeated performance sweeps. The dedicated bitwise-validation
sweep measured 5.70x and 3.40x, respectively. This is materially better than
the transfer-inclusive CUDA POC because no host/device packing or transfers are
required.

For every tested size and thread count, the OpenMP results were **exactly
bitwise identical** to serial for every numeric output. They were not merely
within a floating-point tolerance.

Production NERDSS code is not modified by this POC. The current production
reaction traversal performs shared, stateful mutations, so integrating OpenMP
safely requires thread-local or index-aligned results followed by a
deterministic merge.

## Repository changes

All OpenMP-specific files are in `omp_poc/`:

| File | Purpose |
|---|---|
| `reaction_propagation_omp_bench.cpp` | Serial/OpenMP implementations, deterministic data generation, timing, relative-error checks, bitwise checks, and correctness gates |
| `Makefile` | GCC/OpenMP build with `-O3 -march=native -std=c++17 -fopenmp` |
| `run_benchmark.sh` | Reproducible CPU metadata capture, core affinity, build, execution, and raw-result capture |
| `benchmark_results_i5_13450hx.txt` | Representative performance sweep used for the original scaling analysis |
| `benchmark_results_bitwise_i5_13450hx.txt` | Full sweep with explicit numeric-value counts and bitwise mismatch counts |
| `.gitignore` | Excludes the generated binary and timestamped result directory |

The POC is additive. Existing NERDSS source and build targets are unchanged.

## Relationship to the production hot path

The production timestep performs the following relevant work:

1. `measure_separations_to_identify_possible_reactions()` enumerates molecule
   pairs in the current and neighboring simulation-volume cells.
2. `check_bimolecular_reactions()` traverses template partner data and free
   interfaces.
3. `determine_3D_bimolecular_reaction_probability()` computes the rotational
   diffusion contribution and reaction radius.
4. `get_distance()` computes interface distance and separation and contributes
   to the cutoff decision.
5. The surrounding reaction code mutates `crossbase`, `mycrossint`, `crossrxn`,
   `probvec`, probability-reweighting history, and per-complex `ncross`.
6. Reaction execution and overlap handling eventually call
   `Complex::propagate()`, which constructs a quaternion and transforms the COM
   and interfaces of every molecule in the complex.

The POC covers the independent numerical part of steps 3-4 and the
quaternion/coordinate calculations in step 6. It deliberately excludes the
irregular pair search and the shared mutations in steps 1-2 and 5.

## OpenMP implementation

### Reaction checking

Each flattened candidate-pair record contains only the coordinate, diffusion,
rotational-diffusion, magnitude, and binding-radius values required by the
reaction equations. Each iteration writes one distinct output record and does
not read or update another iteration's output.

The parallel form is therefore a static work-sharing loop:

```cpp
#pragma omp parallel for schedule(static) num_threads(threads)
for (std::int64_t k = 0; k < count; ++k) {
  output[k] = checkReaction(input[k], dt);
}
```

Static scheduling is appropriate because every record performs approximately
the same work. It also makes item assignment repeatable and avoids dynamic
scheduler overhead.

### Propagation

Propagation has two ordered phases:

1. construct one normalized quaternion per complex;
2. transform every point using the quaternion belonging to its complex.

Both phases are individually data parallel, but phase 2 depends on completion
of phase 1. The implementation uses one parallel region with two `omp for`
loops. The implicit barrier at the end of the first loop establishes the
dependency without creating a second worker team:

```cpp
#pragma omp parallel num_threads(threads)
{
#pragma omp for schedule(static)
  for (...) {
    quaternions[k] = makeQuat(motions[k]);
  }

#pragma omp for schedule(static)
  for (...) {
    output[k] = propagate(input[k], motions[complex], quaternions[complex]);
  }
}
```

Combining the loops in one region matters for smaller batches because OpenMP
team creation and synchronization can otherwise be comparable to the useful
work.

### Thread placement

The runner uses:

```text
OMP_BENCH_THREADS=1,2,4,8
OMP_WAIT_POLICY=active
OMP_PROC_BIND=spread
OMP_PLACES=cores
```

WSL reports eight cores and 16 logical CPUs. Spreading workers across core
places avoids placing the main benchmark on sibling threads first. An
exploratory 16-logical-thread sweep frequently regressed relative to eight
workers because barrier and scheduling costs outweighed SMT gains.

## Numerical and bitwise equivalence

### What is checked

The benchmark performs two independent correctness checks after the timed
serial and OpenMP regions:

1. maximum relative error for every floating-point field;
2. an exact comparison of each `double`'s 64-bit IEEE-754 representation.

The bit pattern is copied into `std::uint64_t` with `std::memcpy`, avoiding
strict-aliasing undefined behavior. Fields are compared individually rather
than comparing whole structures, because structure padding bytes are not part
of the numerical result and need not be initialized.

Reaction output checks cover:

- distance;
- separation;
- effective diffusion;
- the `withinRmax` flag.

Propagation checks cover:

- all four quaternion components;
- all three output coordinates;
- the complex ID.

The program exits with status 3 if relative error exceeds `1e-12`, if any
numeric bit pattern differs, or if any flag/ID differs.

### Result

The complete validation sweep covered seven problem sizes and four core
counts. In total it compared:

- 156,577,792 floating-point values bit for bit;
- 44,736,512 reaction flags or complex IDs.

Results:

```text
maximum relative error:       0
numeric bit mismatches:       0 / 156,577,792
flag or ID mismatches:        0 / 44,736,512
```

Thus the tested OpenMP block is exactly bitwise identical to the serial block
under this compiler, build configuration, machine, and input set.

### Why exact equality is expected here

OpenMP changes which worker executes an item, but it does not change the order
of floating-point operations *within* an item. There are no floating-point
reductions, no shared accumulators, and no cross-item dependencies. Serial and
OpenMP call the same inline functions compiled once with the same flags.

This result should not be generalized to different compilers, compiler flags,
math libraries, CPU instruction sets, or a future implementation that adds a
parallel reduction. Full production integration should keep the same
per-candidate arithmetic and deterministic merge order if reproducibility is a
requirement.

## Benchmark methodology

### Measured environment

- CPU: 13th Gen Intel Core i5-13450HX
- WSL topology: 8 cores, 16 logical CPUs, 2 threads per exposed core
- Compiler: GCC 15.2.0
- Flags: `-O3 -march=native -std=c++17 -fopenmp -DNDEBUG`
- Precision: double
- Scheduling: static
- Propagation layout: four points per complex

### Timing protocol

- Deterministic pseudo-random inputs are regenerated with fixed seeds.
- Sizes range from 1,024 to 4,194,304 records.
- Each row contains seven trials; the median time is reported.
- Each trial repeats the operation 300, 80, 15, or 4 times depending on size.
- Serial timing uses a dedicated ordinary C++ loop.
- OpenMP timing includes entry to and exit from the parallel region.
- Correctness and bitwise scans occur after timing and are not included in the
  reported execution time.
- Three complete performance sweeps were used to assess run-to-run variation,
  followed by a complete sweep with explicit bitwise counters.

The host CPU is hybrid, while WSL exposes a synthetic symmetric topology. Turbo
frequency and P/E-core placement can produce apparently superlinear two- or
four-worker results. The large-batch multi-run range is more reliable than an
individual small-batch efficiency value.

## Core scaling

The following values are from the full bitwise-validation sweep at 4,194,304
items. Speedup is relative to the dedicated serial loop.

| Stage | Workers | Serial ms | OpenMP ms | Speedup | Efficiency |
|---|---:|---:|---:|---:|---:|
| reaction | 1 | 118.82 | 98.88 | 1.20x | 120% |
| reaction | 2 | 118.82 | 38.70 | 3.07x | 154% |
| reaction | 4 | 118.82 | 28.36 | 4.19x | 105% |
| reaction | 8 | 118.82 | 20.85 | 5.70x | 71% |
| propagation | 1 | 60.78 | 60.15 | 1.01x | 101% |
| propagation | 2 | 60.78 | 27.19 | 2.24x | 112% |
| propagation | 4 | 60.78 | 21.49 | 2.83x | 71% |
| propagation | 8 | 60.78 | 17.90 | 3.40x | 42% |

The superlinear rows reflect the hybrid-core/turbo caveat above, not an
algorithmic claim. Across the three preceding complete sweeps, eight-worker
reaction speedup at this size was 4.66-4.99x (median 4.93x), while propagation
was 3.17-3.58x (median 3.54x). These ranges are the recommended planning values.

At smaller sizes, parallel overhead is visible and the best worker count can be
below eight. The raw files contain all sizes and thread counts.

## Comparison with the CUDA POC

The prior transfer-inclusive CUDA POC on the same host measured:

| Stage, 4,194,304 items | CUDA end-to-end | Representative 8-worker OpenMP |
|---|---:|---:|
| reaction checking | 76.10 ms | 20.05 ms |
| propagation | 42.57 ms | 17.70 ms |

This is a directional comparison between separate benchmark runs, not a
single-process head-to-head measurement. Nevertheless, the margin is large and
consistent with the implementation difference: OpenMP operates directly on
host-resident data, while the naive CUDA path pays packing and transfer costs.

## Production integration design

### Reaction checking

A pragma must not be placed directly around the existing stateful traversal.
The production code updates shared vectors, per-complex counters, and
probability history. A safe design is:

1. retain current same/neighbor-cell candidate enumeration;
2. flatten immutable candidate inputs into an index-stable array;
3. evaluate the independent numerical prefix with a static OpenMP loop;
4. write one result per candidate or use per-thread reaction-record buffers;
5. merge accepted records in original candidate order;
6. perform `crossbase`, `mycrossint`, `crossrxn`, `probvec`, reweighting, and
   `ncross` mutations in that deterministic merge.

The merge preserves ordering and avoids locks in the arithmetic-heavy loop.
Profiling should verify that flattening plus merging does not consume the saved
time at realistic candidate counts.

### Propagation

Propagation can be parallelized across complexes when each molecule belongs to
exactly one propagated complex and no concurrent phase accesses its mutable
coordinates. Construct one motion/quaternion record per complex, then either:

- parallelize the outer complex loop, keeping each complex's member loop
  serial; or
- flatten all points and use the two-phase implementation demonstrated here.

The first option is less invasive and likely preferable when complexes contain
enough members to amortize scheduling. The second provides better load balance
for highly uneven complex sizes.

### Determinism requirements

To retain the POC's exact serial equivalence in production:

- keep floating-point reductions out of the parallel region;
- give each iteration exclusive output storage;
- preserve random-number assignment independently of worker scheduling;
- merge reaction records in the original serial order;
- avoid `fast-math` changes between serial and OpenMP builds;
- add regression tests that compare output bits for representative timesteps.

## Build and reproduce

```bash
cd /home/yying/NERDSS/omp_poc
make clean all
./run_benchmark.sh

# Optional explicit sizes:
./run_benchmark.sh 65536 1048576 4194304

# Optional logical-thread experiment:
OMP_BENCH_THREADS=1,2,4,8,16 \
OMP_WAIT_POLICY=active OMP_PROC_BIND=spread OMP_PLACES=threads \
./reaction_propagation_omp_bench 1048576
```

The CSV columns `numeric_values_checked`, `numeric_bit_mismatches`,
`discrete_values_checked`, and `id_or_flag_mismatches` make the exact-equality
result auditable for every row.

## Limitations

This proof of concept does not include:

- simulation-volume cell enumeration;
- partner/template and free-interface traversal;
- 1D, 2D, or spherical reaction paths;
- random association sampling;
- probability-history and reweighting mutations;
- association/dissociation side effects;
- overlap resampling and boundary reflection;
- MPI communication;
- production flattening and deterministic-merge cost;
- load imbalance from real complex-size and candidate distributions.

The benchmark establishes that the independent arithmetic is OpenMP-friendly
and exactly reproducible. It does not by itself establish whole-application
speedup.

## Recommendation

Use OpenMP as the next production prototype before pursuing GPU integration.
Target the independent numerical prefix and per-complex propagation, retain a
deterministic serial merge for state mutations, and profile with real timestep
distributions. The measured CPU speedup, exact bitwise result, and lack of data
transfer make this the lower-risk acceleration path.
