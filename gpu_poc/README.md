# GPU reaction-checking and propagation proof of concept

This standalone CUDA benchmark tests whether the data-parallel numerical work in
NERDSS's reaction-checking and propagation hot path can benefit from GPU
execution. It does not change production simulation behavior.

## Production hot path and dependencies

The timestep loop in `EXEs/nerdss_mpi.cpp` calls:

1. `measure_separations_to_identify_possible_reactions()`, which enumerates
   molecule pairs in the same and neighboring `SimulVolume` cells.
2. `check_bimolecular_reactions()`, which traverses molecule/template partner
   and free-interface vectors.
3. For 3D pairs,
   `determine_3D_bimolecular_reaction_probability()` calculates rotational
   diffusion contributions and `Rmax`, then `get_distance()` calculates
   interface separation.
4. `get_distance()` and the remaining probability code append to per-molecule
   `crossbase`, `mycrossint`, `crossrxn`, `probvec`, and reweighting
   vectors, and increment per-complex `ncross`.
5. `perform_bimolecular_reactions()` creates motion vectors. Then
   `check_overlap()` eventually calls `Complex::propagate()`, which constructs
   a quaternion per complex and transforms every member COM and interface.

Steps 1, 2, and the vector mutations in step 4 are irregular, stateful
object-graph operations. The POC flattens only the independent numerical prefix
of the 3D path (the production rotational-diffusion, `Rmax`, distance,
separation, and cutoff equations) and the nonspherical quaternion transform from
`Complex::propagate()`. Both CPU and GPU execute the same double-precision
functions.

## Build and reproduce

Requirements are a CUDA-capable NVIDIA GPU, `nvidia-smi`, a C++ compiler, and
`nvcc`.

```bash
cd /home/yying/NERDSS/gpu_poc
make clean all
./run_benchmark.sh
# Optional explicit sizes:
./run_benchmark.sh 1024 65536 1048576
```

The Makefile detects compute capability and emits a native `sm_XX` cubin. This
is essential on the measured host: driver 581.95 advertises CUDA 13.0, while the
installed nvcc is CUDA 13.3. Default PTX compilation fails at launch with
"the provided PTX was compiled with an unsupported toolchain"; every launch is
therefore error-checked.

## Method and measured environment

- GPU: NVIDIA GeForce RTX 4050 Laptop GPU, compute capability 8.9, 6 GiB
- Driver: 581.95 (WSL reports 580.112 / CUDA 13.0)
- Toolkit: CUDA 13.3, nvcc 13.3.73
- Host compiler: GCC 15.2.0; optimization: `-O3`
- Sizes: 1,024 through 4,194,304 items
- One untimed warmup; 200/100/30/10 repetitions depending on size
- Deterministic inputs; propagation uses four points per complex
- Resident timing measures kernels with inputs already on-device
- End-to-end timing includes H2D inputs, kernel execution, and D2H results
- Correctness gate: zero discrete mismatches and maximum relative error <=1e-12

The full output is in `benchmark_results_rtx4050.txt`.

| Stage | Items | CPU ms | GPU resident ms | GPU end-to-end ms | Resident speedup | End-to-end speedup |
|---|---:|---:|---:|---:|---:|---:|
| reaction | 1,024 | 0.0289 | 0.0108 | 0.1213 | 2.69x | 0.24x |
| reaction | 65,536 | 1.2887 | 0.1022 | 1.3271 | 12.60x | 0.97x |
| reaction | 1,048,576 | 19.6600 | 1.4735 | 19.2934 | 13.34x | 1.02x |
| reaction | 4,194,304 | 79.8094 | 5.8736 | 76.1024 | 13.59x | 1.05x |
| propagation | 1,024 | 0.0135 | 0.0230 | 0.0852 | 0.59x | 0.16x |
| propagation | 65,536 | 0.7485 | 0.0578 | 0.7507 | 12.96x | 1.00x |
| propagation | 262,144 | 3.3471 | 0.1896 | 2.4915 | 17.66x | 1.34x |
| propagation | 4,194,304 | 50.9950 | 4.4885 | 42.5696 | 11.36x | 1.20x |

All runs had zero flag/ID mismatches. Maximum relative error was
`9.55e-15` for reaction checking and `8.88e-16` for propagation.

## Interpretation and recommendation

Resident data is strongly favorable: large batches are about 11-14x faster.
Transfers erase nearly all of the reaction-checking benefit, with a measured
end-to-end break-even around one million pairs and only a 1.05x gain at 4.19
million. Propagation reaches reliable end-to-end benefit around 262,000 points,
but production currently invokes it one complex at a time, far below that
threshold.

Do not integrate per-pair or per-complex GPU launches. A useful production GPU
design would first flatten a whole timestep's candidate pairs and propagation
points, keep coordinates and parameters resident across stages/timesteps, and
scatter only compact accepted-candidate/state updates. Proceed only if
production profiles show batches near these break-even sizes after packing
costs. Otherwise CPU data-layout, SIMD, and thread-level parallelism are the
lower-risk optimization.

This POC excludes cell enumeration, partner/interface searches, 1D/2D and
spherical paths, random sampling, `passocF`/reweighting history, association
side effects, overlap resampling, boundary reflection, MPI communication, and
the packing/scatter cost required by a full integration.
