#include "debug/debug.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

bool is_ghosted(Molecule &mol, MpiContext &mpiContext,
                SimulVolume &simulVolume) {
  bool ghosted = false;
  int xBin = get_x_bin(mpiContext, mol);
  if (mpiContext.rank) {  // if not first rank
    if (xBin <= 0) {
      ghosted = true;
    }
  }
  if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
    if (xBin >= simulVolume.numSubCells.x - 1) {
      ghosted = true;
    }
  }
  return ghosted;
}

bool is_owned_by_processor(Molecule &mol, MpiContext &mpiContext,
                           SimulVolume &simulVolume) {
  bool ghosted = false;
  int xBin = get_x_bin(mpiContext, mol);
  if (mpiContext.rank) {  // if not first rank
    if (xBin <= 0) {
      ghosted = true;
    }
  }
  if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
    if (xBin >= simulVolume.numSubCells.x - 1) {
      ghosted = true;
    }
  }
  return !ghosted;
}

// TODO: update membraneObject.numberOfFreeLipidsEachState[0] when a implicit
// lipid-bound mol move cross ranks