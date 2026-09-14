/*! \file write_all_species.cpp
 *
 * \brief
 *
 * ### Created on 2019-06-05 by Matthew Varga
 */
#include "io/io.hpp"
#include "tracing.hpp"
#include "debug/debug.hpp"
#include "mpi/mpi_function.hpp"

// The implicit lipid is a system-wide pseudo-molecule, not a spatially owned
// one: every rank holds a copy of it and its own share of the free-lipid count,
// and its comCoord is not a position any rank owns.  is_ghosted() is purely
// positional, so it was gating each rank's lipid report on a coordinate that
// drifts as the pseudo-molecule is propagated.  Measured at np=2 on
// mem_localization/SmallBox/FastDsol/IL: rank 1 held 1779 free lipids from the
// first step and reported 0 for as long as that coordinate happened to land in
// its ghost zone (comCoord.x = -78.2 -> xBin = -3), then began reporting once
// it drifted out (comCoord.x = 210.2 -> xBin = 8).  The merged trajectory was
// therefore short by one rank's entire share -- about half the membrane at
// np=2 -- on an interval whose length varied with the seed.
//
// Every other site that touches the implicit lipid already exempts it:
// prepare() forces isGhosted = false for it, derive_ghost_from_ownership()
// skips it, and update_memberMolLists() skips it.  This is the one place that
// did not.
static bool is_ghosted_for_output(Molecule& mol, MpiContext& mpiContext,
                                  SimulVolume& simulVolume) {
  if (mol.isImplicitLipid) return false;
  return is_ghosted(mol, mpiContext, simulVolume);
}

void write_all_species(double simTime, std::vector<Molecule>& moleculeList,
                       std::ofstream& speciesFile, copyCounters& counterArrays,
                       const Membrane& membraneObject, MpiContext& mpiContext,
                       SimulVolume& simulVolume) {
  // FOR MONOMER, EACH PROCESSOR ONLY OUTPUT THE NON-GHOST MOLECULES; FOR PAIR,
  // EACH PROCESSOR OUTPUT THE PAIR IF THE LARGER ID OF THE PAIR IS NON-GHOST

  int i, j;

  int index;
  int p1, p2;
  for (i = 0; i < counterArrays.copyNumSpecies.size(); i++)
    counterArrays.copyNumSpecies[i] = 0;

  for (p1 = 0; p1 < moleculeList.size(); p1++) {
    int numIfaces = moleculeList[p1].interfaceList.size();
    bool isGhosted =
        is_ghosted_for_output(moleculeList[p1], mpiContext, simulVolume);
    for (j = 0; j < numIfaces; j++) {
      // partnerIndex is -1 for an unbound interface, so the partner may only be
      // looked up inside the PAIR branch below.  Reading it before the test
      // indexed moleculeList[-1] on every free interface -- which is most of
      // them -- and is_ghosted() dereferences comCoord.x, 624 bytes before the
      // start of the vector.  Whether that address happened to be mapped
      // decided whether the run crashed, so a parallel run died in about a
      // third of the seeds tried, in whatever unrelated place the corruption
      // surfaced.
      const int partnerIndex =
          moleculeList[p1].interfaceList[j].interaction.partnerIndex;
      if (partnerIndex == -1) {
        // MONOMER
        if (isGhosted == true) continue;
      } else {
        // PAIR
        bool isPatnerGhosted = is_ghosted_for_output(
            moleculeList[partnerIndex], mpiContext, simulVolume);
        if (moleculeList[p1].id >
            moleculeList[p1].interfaceList[j].interaction.partnerId) {
          if (isGhosted == true) continue;
        } else {
          if (isPatnerGhosted == true) continue;
        }
      }
      // find out which state each interface on the molecule is in, and
      // increment copyNumSpecies array.
      index = moleculeList[p1].interfaceList[j].index;
      if (moleculeList[p1].isImplicitLipid == false) {
        counterArrays.copyNumSpecies[index]++;
      } else {
        // For implicit lipid, set copy numbers based on read in template.
        // int molTypeIndex = moleculeList[p1].molTypeIndex;
        for (int tmpStateIndex = 0; tmpStateIndex < membraneObject.nStates;
             tmpStateIndex++) {
          counterArrays.copyNumSpecies[index + tmpStateIndex] =
              membraneObject.numberOfFreeLipidsEachState
                  [tmpStateIndex];  // IL mol must be the first place
        }
      }
    }  // end all interfaces
  }    // end all current molecules

  speciesFile << simTime;
  for (auto i = 0; i < counterArrays.copyNumSpecies.size(); i++) {
    if (counterArrays.singleDouble[i] == 2 &&
        counterArrays.implicitDouble[i] == false) {
      // product state, contains two proteins, so will be double counted above.
      speciesFile << ',' << counterArrays.copyNumSpecies[i] * 0.5;
    } else {
      speciesFile << ',' << counterArrays.copyNumSpecies[i];
    }
  }
  speciesFile << std::endl;
}
