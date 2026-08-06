# OpenMP reaction-checking and propagation proof of concept

This benchmark replaces the CUDA execution backend in `gpu_poc` with OpenMP
for the same flattened numerical work. It does not modify production NERDSS
simulation behavior.

## Exact blocks tested

The CPU serial and OpenMP paths call the same double-precision functions for:

- the independent 3D reaction prefix: rotational diffusion contribution,
  effective diffusion, `Rmax`, distance, separation, and cutoff;
- nonspherical propagation: one quaternion per complex followed by the COM and
  interface-coordinate transform.

Reaction checking uses `#pragma omp parallel for schedule(static)`. Propagation
uses one parallel region containing two static `omp for` loops; the implicit
barrier between them ensures all complex quaternions exist before points are
transformed.

The production pair enumeration, object-graph traversal, probability-history
updates, and reaction-vector mutations remain outside this POC because a naive
pragma around those operations would introduce races.

## Build and reproduce

```bash
cd /home/yying/NERDSS/omp_poc
make clean all
./run_benchmark.sh

# Optional explicit problem sizes:
./run_benchmark.sh 65536 1048576
```

The runner binds OpenMP workers across the physical cores exposed by WSL:

```text
OMP_BENCH_THREADS=1,2,4,8
OMP_WAIT_POLICY=active
OMP_PROC_BIND=spread
OMP_PLACES=cores
```

Set `OMP_BENCH_THREADS` manually when running the binary to test another list.

## Method and measured environment

- CPU: 13th Gen Intel Core i5-13450HX
- WSL topology: 8 cores, 16 logical CPUs, 2 threads per exposed core
- Compiler: GCC 15.2.0 with `-O3 -march=native -fopenmp`
- Sizes: 1,024 through 4,194,304 items
- Seven timing trials per row; the reported value is the median
- Adaptive repetitions: 300, 80, 15, or 4 per trial by problem size
- Deterministic inputs and static scheduling
- Correctness gate: zero ID/flag mismatches and relative error <= `1e-12`
- Three complete benchmark sweeps were run to check run-to-run variation

The tracked raw output is `benchmark_results_i5_13450hx.txt`.

## Representative scaling

The table uses the final full sweep. Speedup is versus the dedicated serial
loop, not the one-thread OpenMP loop.

| Stage | Items | Serial ms | 2 cores | 4 cores | 8 cores |
|---|---:|---:|---:|---:|---:|
| reaction | 1,024 | 0.0154 | 1.81x | 2.09x | 2.71x |
| reaction | 65,536 | 1.5043 | 2.43x | 4.53x | 4.93x |
| reaction | 1,048,576 | 22.3292 | 2.33x | 3.62x | 4.48x |
| reaction | 4,194,304 | 98.8404 | 2.54x | 3.58x | 4.93x |
| propagation | 1,024 | 0.0142 | 2.16x | 3.54x | 2.76x |
| propagation | 65,536 | 0.9353 | 2.43x | 2.90x | 4.67x |
| propagation | 1,048,576 | 16.9798 | 2.59x | 3.28x | 3.88x |
| propagation | 4,194,304 | 62.7034 | 2.34x | 2.75x | 3.54x |

For 4,194,304 items across all three full sweeps, eight-core reaction speedup
was 4.66-4.99x (median 4.93x), and propagation speedup was 3.17-3.58x
(median 3.54x). Every row had zero mismatches and zero numerical error because
each item executes the same instruction sequence without reductions.

Some two- and four-core rows appear superlinear. The i5-13450HX is a hybrid CPU
whose P/E-core distinctions are hidden by WSL's synthetic topology, so worker
placement and turbo frequency can make the serial baseline slower than a
multi-worker run. Treat the large-batch three-run ranges as more reliable than
individual small-batch efficiency values.

An exploratory 16-logical-thread sweep frequently regressed relative to eight
workers, especially for small and medium batches. Barrier and scheduling costs
on the eight exposed cores outweigh SMT gains for these kernels.

## Comparison with the CUDA POC

On the same host, the prior transfer-inclusive CUDA POC took 76.10 ms for
4,194,304 reaction items and 42.57 ms for 4,194,304 propagation points. The
representative eight-core OpenMP times were 20.05 ms and 17.70 ms,
respectively. The comparison is directional rather than a single-process
head-to-head run, but it confirms that avoiding packing and PCIe transfers is
decisive for this block.

## Recommendation for production

Prefer OpenMP over the naive GPU offload for this hot path. Do not simply place
a pragma around the current production pair traversal: `crossbase`,
`mycrossint`, `crossrxn`, `probvec`, reweighting history, and `ncross` are
mutated during checking.

A safe integration should:

1. enumerate candidate pairs as today;
2. evaluate independent numerical probability inputs with a static OpenMP loop;
3. write index-aligned results or thread-local reaction records;
4. merge stateful mutations deterministically after the parallel region;
5. parallelize propagation across complexes only when molecule ownership is
   disjoint for the timestep.

This keeps the simple CPU implementation while preserving determinism and
avoiding the GPU POC's transfer break-even requirement.
