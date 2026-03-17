/*! \file gpu_structs.hpp
 * \brief GPU-friendly Structure-of-Arrays (SoA) data layouts for CUDA kernels.
 *
 * These flat, contiguous arrays replace the host-side AoS (std::vector<Complex>,
 * std::vector<Molecule>) for efficient GPU memory access and coalesced reads/writes.
 *
 * Created for NERDSS GPU parallelization.
 */
#pragma once

#include <cstddef>

#ifdef NERDSS_CUDA
#include <cuda_runtime.h>
#include <curand_kernel.h>
#endif

namespace gpu {

/// Flat SoA representation of Complex data needed on the GPU.
/// Each array is indexed by complex index [0..numComplexes).
struct ComplexSoA {
    // Center-of-mass coordinates
    double* comX = nullptr;
    double* comY = nullptr;
    double* comZ = nullptr;

    // Translational diffusion constants
    double* Dx = nullptr;
    double* Dy = nullptr;
    double* Dz = nullptr;

    // Rotational diffusion constants
    double* Drx = nullptr;
    double* Dry = nullptr;
    double* Drz = nullptr;

    // Output: translational propagation vectors
    double* trajTransX = nullptr;
    double* trajTransY = nullptr;
    double* trajTransZ = nullptr;

    // Output: rotational propagation vectors
    double* trajRotX = nullptr;
    double* trajRotY = nullptr;
    double* trajRotZ = nullptr;

    // Flags
    int* isEmpty = nullptr;
    int* onSurface = nullptr;
    int* onFiber = nullptr;

    // Number of active complexes in this batch
    int numComplexes = 0;
};

/// Flat SoA representation of Molecule coordinate data needed on the GPU.
/// Each array is indexed by molecule index [0..numMolecules).
struct MoleculeSoA {
    // Center-of-mass coordinates
    double* comX = nullptr;
    double* comY = nullptr;
    double* comZ = nullptr;

    // Interface coordinates (flattened: mol i's j-th interface at offset[i] + j)
    double* ifaceX = nullptr;
    double* ifaceY = nullptr;
    double* ifaceZ = nullptr;

    // Per-molecule offset into the flattened interface arrays
    int* ifaceOffset = nullptr;
    // Per-molecule count of interfaces
    int* ifaceCount = nullptr;

    // Molecule type index
    int* molTypeIndex = nullptr;

    // Which complex this molecule belongs to
    int* myComIndex = nullptr;

    // Flags
    int* isEmpty = nullptr;
    int* isImplicitLipid = nullptr;

    int numMolecules = 0;
    int totalIfaces = 0; // total number of interfaces across all molecules
};

/// Cell-list structure for GPU neighbor search.
/// Uses a sorted-by-cell approach for efficient parallel iteration.
struct CellListGPU {
    // For each cell: start index in sortedMolIndices
    int* cellStart = nullptr;
    // For each cell: count of molecules in that cell
    int* cellCount = nullptr;
    // Molecule indices sorted by cell membership
    int* sortedMolIndices = nullptr;
    // For each cell: neighbor cell indices (flattened, max 26 neighbors + self)
    int* neighborOffset = nullptr;
    int* neighborCount = nullptr;
    int* neighborList = nullptr;

    int numCells = 0;
    int numMolecules = 0;
    int totalNeighborEntries = 0;
};

/// Candidate reaction pair produced by the GPU distance kernel.
/// Transferred back to host for the sequential reaction-decision logic.
struct CandidatePair {
    int mol1Index;
    int mol2Index;
    double dist2; // squared distance between interfaces or COMs
    int iface1;   // interface index on mol1
    int iface2;   // interface index on mol2
};

/// Simulation parameters needed on the GPU (small, constant memory).
struct GPUSimParams {
    double timeStep;
    double rMaxLimit;
    double boxX, boxY, boxZ; // water box dimensions (half-extents)
    int isSphere;
    double sphereR;
};

} // namespace gpu
