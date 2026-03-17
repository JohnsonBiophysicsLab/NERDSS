/*! \file gpu_manager.cu
 * \brief Implementation of the GPUManager class for NERDSS CUDA acceleration.
 *
 * Handles device memory lifecycle, host<->device transfers, and kernel dispatch.
 */

#include "gpu/gpu_manager.hpp"
#include "gpu/gpu_structs.hpp"

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <iostream>

// Forward declarations from gpu_kernels.cu
namespace gpu {
void launchInitRandStates(void* states, int n, unsigned long seed);
void launchPropagation(ComplexSoA& complexes, void* randStates, int numComplexes);
void launchDistance(MoleculeSoA& molecules, CellListGPU& cellList,
                   CandidatePair* d_candidatePairs, int* d_numCandidates,
                   int maxCandidates, int numMolecules);
void launchTranslate(MoleculeSoA& molecules, ComplexSoA& complexes, int numMolecules);
void uploadSimParams(const GPUSimParams& params);
}

// Macro for CUDA error checking
#define CUDA_CHECK(call)                                                        \
    do {                                                                        \
        cudaError_t err = call;                                                 \
        if (err != cudaSuccess) {                                               \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,    \
                    cudaGetErrorString(err));                                    \
        }                                                                       \
    } while (0)

namespace gpu {

// ---------------------------------------------------------------------------
// Construction / Destruction
// ---------------------------------------------------------------------------

GPUManager::GPUManager() {}

GPUManager::~GPUManager()
{
    if (!initialized) return;

    freeComplexBuffers();
    freeMoleculeBuffers();
    freeCellListBuffers();
    freeCandidateBuffer();

    if (d_randStates) {
        cudaFree(d_randStates);
        d_randStates = nullptr;
    }
}

// ---------------------------------------------------------------------------
// Static utility methods
// ---------------------------------------------------------------------------

bool GPUManager::isGPUAvailable()
{
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    return (err == cudaSuccess && deviceCount > 0);
}

void GPUManager::printDeviceInfo()
{
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        std::cout << "No CUDA-capable GPU detected.\n";
        return;
    }

    for (int i = 0; i < deviceCount; ++i) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        std::cout << "GPU " << i << ": " << prop.name << "\n"
                  << "  Compute capability: " << prop.major << "." << prop.minor << "\n"
                  << "  Global memory: " << (prop.totalGlobalMem / (1024 * 1024)) << " MB\n"
                  << "  SM count: " << prop.multiProcessorCount << "\n"
                  << "  Max threads per block: " << prop.maxThreadsPerBlock << "\n"
                  << "  Warp size: " << prop.warpSize << "\n";
    }
}

// ---------------------------------------------------------------------------
// Initialization
// ---------------------------------------------------------------------------

void GPUManager::initialize(int maxComplexes, int maxMolecules, int maxIfaces,
                            const Parameters& params, const Membrane& membrane,
                            unsigned long seed)
{
    if (!isGPUAvailable()) {
        std::cerr << "GPU: No CUDA device available. GPU acceleration disabled.\n";
        return;
    }

    printDeviceInfo();

    // Allocate device buffers
    allocateComplexBuffers(maxComplexes);
    allocateMoleculeBuffers(maxMolecules, maxIfaces);

    // Candidate pair buffer: estimate max pairs as 10x number of molecules
    maxCandidatePairs = std::max(maxMolecules * 10, 100000);
    allocateCandidateBuffer(maxCandidatePairs);

    // Initialize cuRAND states
    initRandStates(maxComplexes, seed);

    // Upload simulation parameters to constant memory
    GPUSimParams gpuParams;
    gpuParams.timeStep = params.timeStep;
    gpuParams.rMaxLimit = params.rMaxLimit;
    gpuParams.boxX = membrane.waterBox.x / 2.0;
    gpuParams.boxY = membrane.waterBox.y / 2.0;
    gpuParams.boxZ = membrane.waterBox.z / 2.0;
    gpuParams.isSphere = membrane.isSphere ? 1 : 0;
    gpuParams.sphereR = membrane.sphereR;
    d_params = gpuParams;
    uploadSimParams(gpuParams);

    initialized = true;
    std::cout << "GPU: Initialized with capacity for "
              << maxComplexes << " complexes, "
              << maxMolecules << " molecules, "
              << maxIfaces << " interfaces.\n";
}

// ---------------------------------------------------------------------------
// Buffer allocation helpers
// ---------------------------------------------------------------------------

void GPUManager::allocateComplexBuffers(int n)
{
    freeComplexBuffers();
    allocatedComplexes = n;

    CUDA_CHECK(cudaMalloc(&d_complexes.comX, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.comY, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.comZ, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.Dx, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.Dy, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.Dz, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.Drx, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.Dry, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.Drz, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.trajTransX, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.trajTransY, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.trajTransZ, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.trajRotX, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.trajRotY, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.trajRotZ, n * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_complexes.isEmpty, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_complexes.onSurface, n * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_complexes.onFiber, n * sizeof(int)));
}

void GPUManager::freeComplexBuffers()
{
    if (allocatedComplexes == 0) return;
    cudaFree(d_complexes.comX);     cudaFree(d_complexes.comY);     cudaFree(d_complexes.comZ);
    cudaFree(d_complexes.Dx);       cudaFree(d_complexes.Dy);       cudaFree(d_complexes.Dz);
    cudaFree(d_complexes.Drx);      cudaFree(d_complexes.Dry);      cudaFree(d_complexes.Drz);
    cudaFree(d_complexes.trajTransX); cudaFree(d_complexes.trajTransY); cudaFree(d_complexes.trajTransZ);
    cudaFree(d_complexes.trajRotX); cudaFree(d_complexes.trajRotY); cudaFree(d_complexes.trajRotZ);
    cudaFree(d_complexes.isEmpty);  cudaFree(d_complexes.onSurface); cudaFree(d_complexes.onFiber);
    allocatedComplexes = 0;
}

void GPUManager::allocateMoleculeBuffers(int nMol, int nIface)
{
    freeMoleculeBuffers();
    allocatedMolecules = nMol;
    allocatedIfaces = nIface;

    CUDA_CHECK(cudaMalloc(&d_molecules.comX, nMol * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_molecules.comY, nMol * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_molecules.comZ, nMol * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_molecules.ifaceX, nIface * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_molecules.ifaceY, nIface * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_molecules.ifaceZ, nIface * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_molecules.ifaceOffset, nMol * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_molecules.ifaceCount, nMol * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_molecules.molTypeIndex, nMol * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_molecules.myComIndex, nMol * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_molecules.isEmpty, nMol * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_molecules.isImplicitLipid, nMol * sizeof(int)));
}

void GPUManager::freeMoleculeBuffers()
{
    if (allocatedMolecules == 0) return;
    cudaFree(d_molecules.comX);      cudaFree(d_molecules.comY);      cudaFree(d_molecules.comZ);
    cudaFree(d_molecules.ifaceX);    cudaFree(d_molecules.ifaceY);    cudaFree(d_molecules.ifaceZ);
    cudaFree(d_molecules.ifaceOffset); cudaFree(d_molecules.ifaceCount);
    cudaFree(d_molecules.molTypeIndex); cudaFree(d_molecules.myComIndex);
    cudaFree(d_molecules.isEmpty);   cudaFree(d_molecules.isImplicitLipid);
    allocatedMolecules = 0;
    allocatedIfaces = 0;
}

void GPUManager::allocateCellListBuffers(int nCells, int nNeighborEntries, int nMols)
{
    freeCellListBuffers();
    allocatedCells = nCells;
    allocatedNeighborEntries = nNeighborEntries;

    CUDA_CHECK(cudaMalloc(&d_cellList.cellStart, nCells * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellList.cellCount, nCells * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellList.sortedMolIndices, nMols * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellList.neighborOffset, nCells * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellList.neighborCount, nCells * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_cellList.neighborList, nNeighborEntries * sizeof(int)));
}

void GPUManager::freeCellListBuffers()
{
    if (allocatedCells == 0) return;
    cudaFree(d_cellList.cellStart);  cudaFree(d_cellList.cellCount);
    cudaFree(d_cellList.sortedMolIndices);
    cudaFree(d_cellList.neighborOffset); cudaFree(d_cellList.neighborCount);
    cudaFree(d_cellList.neighborList);
    allocatedCells = 0;
}

void GPUManager::allocateCandidateBuffer(int n)
{
    freeCandidateBuffer();
    maxCandidatePairs = n;
    CUDA_CHECK(cudaMalloc(&d_candidatePairs, n * sizeof(CandidatePair)));
    CUDA_CHECK(cudaMalloc(&d_numCandidates, sizeof(int)));
}

void GPUManager::freeCandidateBuffer()
{
    if (d_candidatePairs) { cudaFree(d_candidatePairs); d_candidatePairs = nullptr; }
    if (d_numCandidates) { cudaFree(d_numCandidates); d_numCandidates = nullptr; }
}

void GPUManager::initRandStates(int n, unsigned long seed)
{
    if (d_randStates) cudaFree(d_randStates);
    CUDA_CHECK(cudaMalloc(&d_randStates, n * sizeof(curandState)));
    launchInitRandStates(d_randStates, n, seed);
}

// ---------------------------------------------------------------------------
// Data upload: Host -> Device
// ---------------------------------------------------------------------------

void GPUManager::uploadComplexData(const std::vector<Complex>& complexList)
{
    int n = static_cast<int>(complexList.size());
    if (n > allocatedComplexes) {
        allocateComplexBuffers(n * 2); // over-allocate for growth
        initRandStates(n * 2, 42); // re-init with more states
    }
    d_complexes.numComplexes = n;

    // Pack host arrays
    std::vector<double> hComX(n), hComY(n), hComZ(n);
    std::vector<double> hDx(n), hDy(n), hDz(n);
    std::vector<double> hDrx(n), hDry(n), hDrz(n);
    std::vector<int> hEmpty(n), hOnSurface(n), hOnFiber(n);

    for (int i = 0; i < n; ++i) {
        hComX[i] = complexList[i].comCoord.x;
        hComY[i] = complexList[i].comCoord.y;
        hComZ[i] = complexList[i].comCoord.z;
        hDx[i] = complexList[i].D.x;
        hDy[i] = complexList[i].D.y;
        hDz[i] = complexList[i].D.z;
        hDrx[i] = complexList[i].Dr.x;
        hDry[i] = complexList[i].Dr.y;
        hDrz[i] = complexList[i].Dr.z;
        hEmpty[i] = complexList[i].isEmpty ? 1 : 0;
        hOnSurface[i] = complexList[i].OnSurface ? 1 : 0;
        hOnFiber[i] = complexList[i].onFiber ? 1 : 0;
    }

    CUDA_CHECK(cudaMemcpy(d_complexes.comX, hComX.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.comY, hComY.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.comZ, hComZ.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.Dx, hDx.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.Dy, hDy.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.Dz, hDz.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.Drx, hDrx.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.Dry, hDry.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.Drz, hDrz.data(), n * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.isEmpty, hEmpty.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.onSurface, hOnSurface.data(), n * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_complexes.onFiber, hOnFiber.data(), n * sizeof(int), cudaMemcpyHostToDevice));
}

void GPUManager::uploadMoleculeData(const std::vector<Molecule>& moleculeList,
                                     const std::vector<Complex>& complexList)
{
    int nMol = static_cast<int>(moleculeList.size());

    // Count total interfaces
    int totalIfaces = 0;
    for (int i = 0; i < nMol; ++i) {
        totalIfaces += static_cast<int>(moleculeList[i].interfaceList.size());
    }

    if (nMol > allocatedMolecules || totalIfaces > allocatedIfaces) {
        allocateMoleculeBuffers(nMol * 2, totalIfaces * 2);
    }

    d_molecules.numMolecules = nMol;
    d_molecules.totalIfaces = totalIfaces;

    // Pack host arrays
    std::vector<double> hComX(nMol), hComY(nMol), hComZ(nMol);
    std::vector<double> hIfaceX(totalIfaces), hIfaceY(totalIfaces), hIfaceZ(totalIfaces);
    std::vector<int> hIfaceOffset(nMol), hIfaceCount(nMol);
    std::vector<int> hMolType(nMol), hMyComIdx(nMol);
    std::vector<int> hEmpty(nMol), hIsIL(nMol);

    int ifOffset = 0;
    for (int i = 0; i < nMol; ++i) {
        const auto& mol = moleculeList[i];
        hComX[i] = mol.comCoord.x;
        hComY[i] = mol.comCoord.y;
        hComZ[i] = mol.comCoord.z;
        hMolType[i] = mol.molTypeIndex;
        hMyComIdx[i] = mol.myComIndex;
        hEmpty[i] = mol.isEmpty ? 1 : 0;
        hIsIL[i] = mol.isImplicitLipid ? 1 : 0;
        hIfaceOffset[i] = ifOffset;
        hIfaceCount[i] = static_cast<int>(mol.interfaceList.size());

        for (int j = 0; j < hIfaceCount[i]; ++j) {
            hIfaceX[ifOffset + j] = mol.interfaceList[j].coord.x;
            hIfaceY[ifOffset + j] = mol.interfaceList[j].coord.y;
            hIfaceZ[ifOffset + j] = mol.interfaceList[j].coord.z;
        }
        ifOffset += hIfaceCount[i];
    }

    CUDA_CHECK(cudaMemcpy(d_molecules.comX, hComX.data(), nMol * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.comY, hComY.data(), nMol * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.comZ, hComZ.data(), nMol * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.ifaceX, hIfaceX.data(), totalIfaces * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.ifaceY, hIfaceY.data(), totalIfaces * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.ifaceZ, hIfaceZ.data(), totalIfaces * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.ifaceOffset, hIfaceOffset.data(), nMol * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.ifaceCount, hIfaceCount.data(), nMol * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.molTypeIndex, hMolType.data(), nMol * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.myComIndex, hMyComIdx.data(), nMol * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.isEmpty, hEmpty.data(), nMol * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_molecules.isImplicitLipid, hIsIL.data(), nMol * sizeof(int), cudaMemcpyHostToDevice));
}

void GPUManager::uploadCellList(const SimulVolume& simulVolume,
                                 const std::vector<Molecule>& moleculeList)
{
    int nCells = static_cast<int>(simulVolume.subCellList.size());
    int nMol = static_cast<int>(moleculeList.size());

    // Compute total neighbor entries
    int totalNeighborEntries = 0;
    for (int c = 0; c < nCells; ++c) {
        totalNeighborEntries += static_cast<int>(simulVolume.subCellList[c].neighborList.size());
    }

    if (nCells > allocatedCells || totalNeighborEntries > allocatedNeighborEntries) {
        allocateCellListBuffers(nCells * 2, totalNeighborEntries * 2, nMol * 2);
    }

    d_cellList.numCells = nCells;
    d_cellList.numMolecules = nMol;
    d_cellList.totalNeighborEntries = totalNeighborEntries;

    // Build sorted molecule list and cell start/count arrays
    std::vector<int> hCellStart(nCells);
    std::vector<int> hCellCount(nCells);
    std::vector<int> hSortedMols;
    hSortedMols.reserve(nMol);

    int runningStart = 0;
    for (int c = 0; c < nCells; ++c) {
        hCellStart[c] = runningStart;
        hCellCount[c] = static_cast<int>(simulVolume.subCellList[c].memberMolList.size());
        for (int m : simulVolume.subCellList[c].memberMolList) {
            hSortedMols.push_back(m);
        }
        runningStart += hCellCount[c];
    }

    // Build flattened neighbor list
    std::vector<int> hNeighborOffset(nCells);
    std::vector<int> hNeighborCount(nCells);
    std::vector<int> hNeighborList;
    hNeighborList.reserve(totalNeighborEntries);

    int nOffset = 0;
    for (int c = 0; c < nCells; ++c) {
        hNeighborOffset[c] = nOffset;
        hNeighborCount[c] = static_cast<int>(simulVolume.subCellList[c].neighborList.size());
        for (int n : simulVolume.subCellList[c].neighborList) {
            hNeighborList.push_back(n);
        }
        nOffset += hNeighborCount[c];
    }

    int sortedSize = static_cast<int>(hSortedMols.size());

    CUDA_CHECK(cudaMemcpy(d_cellList.cellStart, hCellStart.data(), nCells * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cellList.cellCount, hCellCount.data(), nCells * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cellList.sortedMolIndices, hSortedMols.data(), sortedSize * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cellList.neighborOffset, hNeighborOffset.data(), nCells * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cellList.neighborCount, hNeighborCount.data(), nCells * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_cellList.neighborList, hNeighborList.data(), totalNeighborEntries * sizeof(int), cudaMemcpyHostToDevice));
}

// ---------------------------------------------------------------------------
// Kernel launch wrappers
// ---------------------------------------------------------------------------

void GPUManager::launchPropagationKernel()
{
    if (!initialized) return;
    launchPropagation(d_complexes, d_randStates, d_complexes.numComplexes);
}

void GPUManager::downloadPropagationResults(std::vector<Complex>& complexList)
{
    if (!initialized) return;

    int n = d_complexes.numComplexes;
    std::vector<double> hTrajTX(n), hTrajTY(n), hTrajTZ(n);
    std::vector<double> hTrajRX(n), hTrajRY(n), hTrajRZ(n);

    CUDA_CHECK(cudaMemcpy(hTrajTX.data(), d_complexes.trajTransX, n * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTrajTY.data(), d_complexes.trajTransY, n * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTrajTZ.data(), d_complexes.trajTransZ, n * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTrajRX.data(), d_complexes.trajRotX, n * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTrajRY.data(), d_complexes.trajRotY, n * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(hTrajRZ.data(), d_complexes.trajRotZ, n * sizeof(double), cudaMemcpyDeviceToHost));

    for (int i = 0; i < n && i < static_cast<int>(complexList.size()); ++i) {
        if (complexList[i].isEmpty) continue;
        complexList[i].trajTrans.x = hTrajTX[i];
        complexList[i].trajTrans.y = hTrajTY[i];
        complexList[i].trajTrans.z = hTrajTZ[i];
        complexList[i].trajRot.x = hTrajRX[i];
        complexList[i].trajRot.y = hTrajRY[i];
        complexList[i].trajRot.z = hTrajRZ[i];
    }
}

void GPUManager::launchDistanceKernel(std::vector<CandidatePair>& candidatePairs)
{
    if (!initialized) return;

    launchDistance(d_molecules, d_cellList, d_candidatePairs, d_numCandidates,
                  maxCandidatePairs, d_cellList.numMolecules);

    // Download result count
    int numFound = 0;
    CUDA_CHECK(cudaMemcpy(&numFound, d_numCandidates, sizeof(int), cudaMemcpyDeviceToHost));
    numFound = std::min(numFound, maxCandidatePairs);

    // Download candidate pairs
    candidatePairs.resize(numFound);
    if (numFound > 0) {
        CUDA_CHECK(cudaMemcpy(candidatePairs.data(), d_candidatePairs,
                              numFound * sizeof(CandidatePair), cudaMemcpyDeviceToHost));
    }
}

// ---------------------------------------------------------------------------
// Dynamic resizing
// ---------------------------------------------------------------------------

void GPUManager::resizeIfNeeded(int numComplexes, int numMolecules, int totalIfaces)
{
    if (numComplexes > allocatedComplexes) {
        int newSize = numComplexes * 2;
        allocateComplexBuffers(newSize);
        initRandStates(newSize, 42);
    }
    if (numMolecules > allocatedMolecules || totalIfaces > allocatedIfaces) {
        allocateMoleculeBuffers(numMolecules * 2, totalIfaces * 2);
    }
}

} // namespace gpu
