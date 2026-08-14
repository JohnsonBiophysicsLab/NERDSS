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
    bool isGhosted = is_ghosted(moleculeList[p1], mpiContext, simulVolume);
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
        bool isPatnerGhosted =
            is_ghosted(moleculeList[partnerIndex], mpiContext, simulVolume);
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
