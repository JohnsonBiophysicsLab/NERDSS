/*! \file gpu_manager.hpp
 * \brief Host-side GPU memory manager for NERDSS CUDA acceleration.
 *
 * GPUManager handles:
 *   - Device memory allocation / deallocation
 *   - Host-to-device and device-to-host data transfers
 *   - Kernel launch orchestration
 *   - cuRAND state management
 *
 * The interface is designed so that the main simulation loop can call
 * high-level methods without touching CUDA directly.
 */
#pragma once

#include "gpu/gpu_structs.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Membrane.hpp"
#include "classes/class_SimulVolume.hpp"

#include <vector>

namespace gpu {

class GPUManager {
public:
    GPUManager();
    ~GPUManager();

    /// Initialize GPU device, allocate memory, set up cuRAND states.
    /// Call once after system setup, before the simulation loop.
    /// @param maxComplexes  upper bound on number of complexes
    /// @param maxMolecules  upper bound on number of molecules
    /// @param maxIfaces     upper bound on total number of interfaces
    /// @param params        simulation parameters
    /// @param membrane      membrane/boundary info
    /// @param seed          RNG seed for cuRAND
    void initialize(int maxComplexes, int maxMolecules, int maxIfaces,
                    const Parameters& params, const Membrane& membrane,
                    unsigned long seed);

    /// Upload complex data from host vectors to device SoA arrays.
    void uploadComplexData(const std::vector<Complex>& complexList);

    /// Upload molecule coordinate data from host vectors to device SoA arrays.
    void uploadMoleculeData(const std::vector<Molecule>& moleculeList,
                            const std::vector<Complex>& complexList);

    /// Upload cell-list data from SimulVolume to GPU for neighbor search.
    void uploadCellList(const SimulVolume& simulVolume,
                        const std::vector<Molecule>& moleculeList);

    /// Launch the propagation-vector kernel on GPU.
    /// Generates random translational + rotational displacements for every
    /// active complex in parallel.
    void launchPropagationKernel();

    /// Download computed propagation vectors back to host Complex objects.
    void downloadPropagationResults(std::vector<Complex>& complexList);

    /// Launch the pairwise-distance kernel on GPU.
    /// For each cell, computes distances between all molecule pairs in
    /// the cell and its neighbors, filtering by rMaxLimit.
    /// Returns candidate pairs on the host.
    void launchDistanceKernel(std::vector<CandidatePair>& candidatePairs);

    /// Resize GPU buffers if the system grew beyond current capacity.
    void resizeIfNeeded(int numComplexes, int numMolecules, int totalIfaces);

    /// Query whether a CUDA-capable GPU is available.
    static bool isGPUAvailable();

    /// Print GPU device info to stdout.
    static void printDeviceInfo();

    /// Get the number of active complexes currently on the GPU.
    int getNumComplexes() const { return d_complexes.numComplexes; }

    /// Get the number of active molecules currently on the GPU.
    int getNumMolecules() const { return d_molecules.numMolecules; }

private:
    // Device-side SoA data
    ComplexSoA d_complexes;
    MoleculeSoA d_molecules;
    CellListGPU d_cellList;
    GPUSimParams d_params;

    // cuRAND states (one per thread, enough for max complexes)
    void* d_randStates = nullptr; // curandState*

    // Candidate pair output buffer on device
    CandidatePair* d_candidatePairs = nullptr;
    int* d_numCandidates = nullptr;  // atomic counter on device
    int maxCandidatePairs = 0;

    // Current allocated capacities
    int allocatedComplexes = 0;
    int allocatedMolecules = 0;
    int allocatedIfaces = 0;
    int allocatedCells = 0;
    int allocatedNeighborEntries = 0;

    bool initialized = false;

    // Internal helpers
    void allocateComplexBuffers(int n);
    void allocateMoleculeBuffers(int nMol, int nIface);
    void allocateCellListBuffers(int nCells, int nNeighborEntries, int nMols);
    void allocateCandidateBuffer(int n);
    void freeComplexBuffers();
    void freeMoleculeBuffers();
    void freeCellListBuffers();
    void freeCandidateBuffer();
    void initRandStates(int n, unsigned long seed);
};

} // namespace gpu
