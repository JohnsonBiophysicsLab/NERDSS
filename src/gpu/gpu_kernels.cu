/*! \file gpu_kernels.cu
 * \brief CUDA kernels for NERDSS GPU-accelerated simulation.
 *
 * Contains kernels for:
 *   1. Propagation vector generation (random translational + rotational motion)
 *   2. Pairwise distance calculation for neighbor-cell reaction search
 *   3. cuRAND state initialization
 */

#include "gpu/gpu_structs.hpp"

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cmath>
#include <cstdio>

namespace gpu {

// ---------------------------------------------------------------------------
// Constant memory for simulation parameters (fast broadcast to all threads)
// ---------------------------------------------------------------------------
__constant__ GPUSimParams c_params;

// ---------------------------------------------------------------------------
// Helper: Box-Muller Gaussian from cuRAND uniform
// (curand_normal_double is also available but this keeps it explicit)
// ---------------------------------------------------------------------------
__device__ inline double gaussianRand(curandState* state)
{
    return curand_normal_double(state);
}

// ---------------------------------------------------------------------------
// Kernel: Initialize cuRAND states
// ---------------------------------------------------------------------------
__global__ void initRandStatesKernel(curandState* states, int n, unsigned long seed)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    // Each thread gets a unique sequence number; offset = 0
    curand_init(seed, idx, 0, &states[idx]);
}

// ---------------------------------------------------------------------------
// Kernel: Generate propagation vectors for all complexes
//
// Each thread handles one complex. Generates 3 translational + 3 rotational
// Gaussian-distributed random displacements scaled by diffusion constants.
// This is the GPU equivalent of create_complex_propagation_vectors().
// ---------------------------------------------------------------------------
__global__ void propagationKernel(
    ComplexSoA complexes,
    curandState* randStates,
    int numComplexes)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numComplexes) return;

    // Skip empty or invalid complexes
    if (complexes.isEmpty[idx]) return;

    curandState localState = randStates[idx];
    double dt = c_params.timeStep;

    double dx = complexes.Dx[idx];
    double dy = complexes.Dy[idx];
    double dz = complexes.Dz[idx];
    double drx = complexes.Drx[idx];
    double dry = complexes.Dry[idx];
    double drz = complexes.Drz[idx];

    // On-surface complexes in sphere mode: handled by host (rare path).
    // For box systems and interior-of-sphere complexes, generate standard
    // Brownian displacements.
    // On-surface complexes (membrane-bound): D.z is typically ~0, producing
    // near-zero z displacement -- consistent with the CPU code path.

    // Translational displacement: dx_i = sqrt(2 * dt * D_i) * N(0,1)
    complexes.trajTransX[idx] = sqrt(2.0 * dt * dx) * curand_normal_double(&localState);
    complexes.trajTransY[idx] = sqrt(2.0 * dt * dy) * curand_normal_double(&localState);
    complexes.trajTransZ[idx] = sqrt(2.0 * dt * dz) * curand_normal_double(&localState);

    // Rotational displacement: dtheta_i = sqrt(2 * dt * Dr_i) * N(0,1)
    complexes.trajRotX[idx] = sqrt(2.0 * dt * drx) * curand_normal_double(&localState);
    complexes.trajRotY[idx] = sqrt(2.0 * dt * dry) * curand_normal_double(&localState);
    complexes.trajRotZ[idx] = sqrt(2.0 * dt * drz) * curand_normal_double(&localState);

    // Write back cuRAND state
    randStates[idx] = localState;
}

// ---------------------------------------------------------------------------
// Kernel: Pairwise distance calculation using cell lists
//
// Each thread handles one molecule. For that molecule, it iterates over
// molecules in the same cell and neighbor cells, computing squared distances.
// Pairs within rMaxLimit are written to the candidate-pair output buffer
// using atomicAdd on a global counter.
//
// This replaces the nested cell-iteration loop in the CPU code that calls
// check_bimolecular_reactions() for each pair.
// ---------------------------------------------------------------------------
__global__ void distanceKernel(
    MoleculeSoA molecules,
    CellListGPU cellList,
    CandidatePair* candidatePairs,
    int* numCandidates,
    int maxCandidates,
    int numMolecules)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numMolecules) return;

    // Get the molecule index from the sorted list
    int mol1Idx = cellList.sortedMolIndices[idx];

    // Skip empty or implicit-lipid molecules
    if (molecules.isEmpty[mol1Idx] || molecules.isImplicitLipid[mol1Idx])
        return;

    double m1x = molecules.comX[mol1Idx];
    double m1y = molecules.comY[mol1Idx];
    double m1z = molecules.comZ[mol1Idx];

    double rMax2 = c_params.rMaxLimit * c_params.rMaxLimit;

    // Determine which cell this molecule belongs to
    // (We iterate the sorted list, so we can binary-search for cell boundaries,
    //  but a simpler approach: each molecule knows its cell via the sorted order.)
    // Find cell index for this position in the sorted array
    int cellIdx = -1;
    for (int c = 0; c < cellList.numCells; ++c) {
        if (idx >= cellList.cellStart[c] && idx < cellList.cellStart[c] + cellList.cellCount[c]) {
            cellIdx = c;
            break;
        }
    }
    if (cellIdx < 0) return;

    // Check molecules in the same cell (only those with higher sorted index to avoid double-counting)
    int cellBegin = cellList.cellStart[cellIdx];
    int cellEnd = cellBegin + cellList.cellCount[cellIdx];

    for (int j = idx + 1; j < cellEnd; ++j) {
        int mol2Idx = cellList.sortedMolIndices[j];
        if (molecules.isEmpty[mol2Idx] || molecules.isImplicitLipid[mol2Idx])
            continue;

        double dx = molecules.comX[mol2Idx] - m1x;
        double dy = molecules.comY[mol2Idx] - m1y;
        double dz = molecules.comZ[mol2Idx] - m1z;
        double dist2 = dx * dx + dy * dy + dz * dz;

        if (dist2 <= rMax2) {
            int slot = atomicAdd(numCandidates, 1);
            if (slot < maxCandidates) {
                candidatePairs[slot].mol1Index = mol1Idx;
                candidatePairs[slot].mol2Index = mol2Idx;
                candidatePairs[slot].dist2 = dist2;
                candidatePairs[slot].iface1 = -1; // to be resolved on host
                candidatePairs[slot].iface2 = -1;
            }
        }
    }

    // Check molecules in neighbor cells
    int nStart = cellList.neighborOffset[cellIdx];
    int nCount = cellList.neighborCount[cellIdx];

    for (int ni = 0; ni < nCount; ++ni) {
        int neighCell = cellList.neighborList[nStart + ni];
        int neighBegin = cellList.cellStart[neighCell];
        int neighEnd = neighBegin + cellList.cellCount[neighCell];

        for (int j = neighBegin; j < neighEnd; ++j) {
            int mol2Idx = cellList.sortedMolIndices[j];
            if (molecules.isEmpty[mol2Idx] || molecules.isImplicitLipid[mol2Idx])
                continue;

            double dx = molecules.comX[mol2Idx] - m1x;
            double dy = molecules.comY[mol2Idx] - m1y;
            double dz = molecules.comZ[mol2Idx] - m1z;
            double dist2 = dx * dx + dy * dy + dz * dz;

            if (dist2 <= rMax2) {
                int slot = atomicAdd(numCandidates, 1);
                if (slot < maxCandidates) {
                    candidatePairs[slot].mol1Index = mol1Idx;
                    candidatePairs[slot].mol2Index = mol2Idx;
                    candidatePairs[slot].dist2 = dist2;
                    candidatePairs[slot].iface1 = -1;
                    candidatePairs[slot].iface2 = -1;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Kernel: Position update (propagate complexes)
//
// Each thread updates one complex's COM and all its member molecules'
// coordinates with the computed translational displacement.
// Rotation is more complex (Euler matrix) and handled on host for now.
// This kernel applies the simple translation component.
// ---------------------------------------------------------------------------
__global__ void translateKernel(
    MoleculeSoA molecules,
    ComplexSoA complexes,
    int numMolecules)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numMolecules) return;

    if (molecules.isEmpty[idx] || molecules.isImplicitLipid[idx])
        return;

    int comIdx = molecules.myComIndex[idx];
    if (comIdx < 0 || complexes.isEmpty[comIdx])
        return;

    double tx = complexes.trajTransX[comIdx];
    double ty = complexes.trajTransY[comIdx];
    double tz = complexes.trajTransZ[comIdx];

    // Update molecule COM
    molecules.comX[idx] += tx;
    molecules.comY[idx] += ty;
    molecules.comZ[idx] += tz;

    // Update interface coordinates
    int ifStart = molecules.ifaceOffset[idx];
    int ifCount = molecules.ifaceCount[idx];
    for (int i = 0; i < ifCount; ++i) {
        molecules.ifaceX[ifStart + i] += tx;
        molecules.ifaceY[ifStart + i] += ty;
        molecules.ifaceZ[ifStart + i] += tz;
    }
}

// ---------------------------------------------------------------------------
// Host-callable wrapper functions
// ---------------------------------------------------------------------------

void launchInitRandStates(void* states, int n, unsigned long seed)
{
    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / blockSize;
    initRandStatesKernel<<<numBlocks, blockSize>>>(
        reinterpret_cast<curandState*>(states), n, seed);
    cudaDeviceSynchronize();
}

void launchPropagation(ComplexSoA& complexes, void* randStates, int numComplexes)
{
    int blockSize = 256;
    int numBlocks = (numComplexes + blockSize - 1) / blockSize;
    propagationKernel<<<numBlocks, blockSize>>>(
        complexes,
        reinterpret_cast<curandState*>(randStates),
        numComplexes);
    cudaDeviceSynchronize();
}

void launchDistance(
    MoleculeSoA& molecules,
    CellListGPU& cellList,
    CandidatePair* d_candidatePairs,
    int* d_numCandidates,
    int maxCandidates,
    int numMolecules)
{
    // Reset candidate counter
    cudaMemset(d_numCandidates, 0, sizeof(int));

    int blockSize = 256;
    int numBlocks = (numMolecules + blockSize - 1) / blockSize;
    distanceKernel<<<numBlocks, blockSize>>>(
        molecules, cellList, d_candidatePairs, d_numCandidates,
        maxCandidates, numMolecules);
    cudaDeviceSynchronize();
}

void launchTranslate(MoleculeSoA& molecules, ComplexSoA& complexes, int numMolecules)
{
    int blockSize = 256;
    int numBlocks = (numMolecules + blockSize - 1) / blockSize;
    translateKernel<<<numBlocks, blockSize>>>(molecules, complexes, numMolecules);
    cudaDeviceSynchronize();
}

void uploadSimParams(const GPUSimParams& params)
{
    cudaMemcpyToSymbol(c_params, &params, sizeof(GPUSimParams));
}

} // namespace gpu
