# NERDSS GPU Acceleration Guide

## Overview

NERDSS supports optional CUDA GPU acceleration to parallelize computationally
intensive parts of the simulation. When built with GPU support, the following
phases are offloaded to the GPU:

1. **Propagation vector generation** -- Random translational and rotational
   Brownian displacements are computed in parallel for all complexes on the GPU
   using cuRAND, replacing the per-complex sequential CPU computation.

2. **Pairwise distance calculations** -- The cell-based neighbor search for
   identifying candidate bimolecular reaction pairs is parallelized on the GPU.
   Molecule pairs within `rMaxLimit` are identified in bulk and transferred back
   to the host for the sequential reaction-decision logic.

The reaction decision logic, association/dissociation mechanics, and I/O remain
on the CPU, as they involve complex branching and sequential state mutations.

## Requirements

- **CUDA Toolkit >= 10.0** (nvcc compiler, CUDA runtime, cuRAND library)
- **A CUDA-capable NVIDIA GPU** (compute capability >= 6.0, i.e. Pascal or newer)
- **GSL** (same requirement as the CPU build)
- **CMake >= 3.10** (for the CMake build) or **GNU Make** (for the Makefile build)

## Building with CMake

```bash
mkdir build && cd build
cmake .. -DNERDSS_CUDA=ON -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

This produces two executables:
- `nerdss` -- CPU-only (always built)
- `nerdss_cuda` -- GPU-accelerated (built when `-DNERDSS_CUDA=ON`)

### Specifying GPU architectures

By default, the build targets compute capabilities 60, 70, 75, 80, and 86.
Override with:

```bash
cmake .. -DNERDSS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="80;86;89;90"
```

## Building with Make

```bash
make cuda
```

This produces `bin/nerdss_cuda`. The Makefile automatically detects `nvcc` and
targets common GPU architectures. Edit the `CUDA_ARCH` variable in `Makefile`
to customize.

## Running

The GPU-accelerated binary is a drop-in replacement for the CPU binary:

```bash
./bin/nerdss_cuda -f input.inp
```

At startup, the program prints GPU device information and whether GPU
acceleration is active. If no CUDA device is detected at runtime, the program
falls back to CPU-only execution automatically.

## Architecture

### Data Layout

The GPU uses a Structure-of-Arrays (SoA) layout for coalesced memory access,
converting from the host's Array-of-Structures (AoS) layout (`std::vector<Complex>`,
`std::vector<Molecule>`) at transfer boundaries. Key structures:

- `gpu::ComplexSoA` -- per-complex diffusion constants, coordinates, propagation vectors
- `gpu::MoleculeSoA` -- per-molecule coordinates, interface coordinates, type indices
- `gpu::CellListGPU` -- flattened cell-list for GPU neighbor search

### Kernel Overview

| Kernel | Purpose | Parallelism |
|--------|---------|-------------|
| `propagationKernel` | Generate Brownian displacements | 1 thread per complex |
| `distanceKernel` | Pairwise distance search | 1 thread per molecule |
| `translateKernel` | Apply translational motion | 1 thread per molecule |
| `initRandStatesKernel` | Initialize cuRAND states | 1 thread per complex |

### Host-Device Transfer

The `gpu::GPUManager` class manages all device memory and transfers:

1. **Upload** complex/molecule data from host vectors to device SoA buffers
2. **Launch** CUDA kernels
3. **Download** results (propagation vectors, candidate pairs) back to host

Buffers are over-allocated (2x current size) to accommodate system growth
without frequent reallocation.

## Performance Notes

- GPU acceleration is most beneficial for large systems (>10,000 molecules)
  where the propagation and distance kernels dominate wall time.
- For small systems, CPU-only may be faster due to host-device transfer overhead.
- The GPU is used between the bimolecular reaction search and the reaction
  decision phases of each timestep.

## File Structure

```
include/gpu/
  gpu_structs.hpp    -- GPU SoA data structures and kernel parameter types
  gpu_manager.hpp    -- GPUManager class interface

src/gpu/
  gpu_kernels.cu     -- CUDA kernel implementations
  gpu_manager.cu     -- GPUManager class implementation

docs/
  GPU_BUILD.md       -- This file
```
