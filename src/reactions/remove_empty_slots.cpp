#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "boundary_conditions/reflect_functions.hpp"
#include "debug/debug.hpp"
#include "error/error.hpp"
#include "io/io.hpp"
#include "macro.hpp"
#include "math/constants.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#ifdef mpi_
#include "mpi.h"
#endif
#include "mpi/mpi_function.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "split.cpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void remove_empty_slots(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, MpiContext& mpiContext) {
  //------------------------------------------------------------------------------------
  // Remove empty complexes:
  if (DEBUG) {
    debug_firstEmptyIndex(mpiContext, "Before 1.10");
  }
  if (DEBUG) {
    DEBUG_FIND_MOL("5: before 1.10");
  }
  if (DEBUG) {
    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                     "// before 1.10");
  }
  // Put the last non-empty complex in the list to the first empty slot:
  sort(Complex::emptyComList.begin(), Complex::emptyComList.end());
  int lastNonEmptyIndex = complexList.size() - 1;

  // newComIndex[old] is where the complex that sat at `old` ends up, or -1 if it
  // is gone.  Announcing the move through the moved complex's own memberList,
  // which is what this loop used to do, reaches only the molecules that complex
  // claims; a molecule whose myComIndex names a complex that does not list it
  // back keeps the old index, and the pop_back() below puts that index past the
  // end of the list.
  const int oldComCount{static_cast<int>(complexList.size())};
  std::vector<int> newComIndex(oldComCount);
  for (int i{0}; i < oldComCount; ++i) newComIndex[i] = i;
  for (auto& emptied : Complex::emptyComList)
    if (emptied >= 0 && emptied < oldComCount) newComIndex[emptied] = -1;

  for (auto& firstEmptyIndex : Complex::emptyComList) {
    // Find the last non-empty complex:
    while (complexList[lastNonEmptyIndex].isEmpty) lastNonEmptyIndex--;
    if (lastNonEmptyIndex <= firstEmptyIndex) {
      break;
    }
    // Move last non-empty complex to the first empty position:
    // cout << "  Moving a complex from " << lastNonEmptyIndex << " to " <<
    // firstEmptyIndex << "... "
    //         << "complexList[firstEmptyIndex].isEmpty=" <<
    //         complexList[firstEmptyIndex].isEmpty << endl;
    complexList[firstEmptyIndex] = complexList[lastNonEmptyIndex];

    // Update complex index to match new position:
    complexList[firstEmptyIndex].index = firstEmptyIndex;

    newComIndex[lastNonEmptyIndex] = firstEmptyIndex;
    lastNonEmptyIndex--;
  }

  // Remove last elements in complexList
  for (int k = 0; k < Complex::emptyComList.size(); k++) complexList.pop_back();

  // Empty Complex::emptyComList:
  Complex::emptyComList.clear();

  // myComIndex is the only place a complex index is stored outside the complex
  // itself.  A molecule whose complex is simply gone gets -1, the value the rest
  // of the code already reads as "no complex": it is wrong, but it is the only
  // answer that is in range, and it cannot be written through.
  for (auto& mol : moleculeList) {
    if (mol.isEmpty) continue;
    mol.myComIndex = (mol.myComIndex < 0 || mol.myComIndex >= oldComCount)
                         ? -1
                         : newComIndex[mol.myComIndex];
  }

  if (DEBUG) {
    debug_firstEmptyIndex(mpiContext, "Before 3.10");
  }
  if (DEBUG) {
    DEBUG_FIND_MOL("6: before 3.10");
  }
  if (DEBUG) {
    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                     "// before 3.10");
  }

  // if(DEBUG) {debug_bndpartner_interface(mpiContext, "35: (before remove
  // empty molecules)");}

  //------------------------------------------------------------------------------------
  // Remove empty molecules:
  sort(Molecule::emptyMolList.begin(), Molecule::emptyMolList.end());

  // newMolIndex[old] is where the molecule that sat at `old` ends up, or -1 if
  // it is gone.  Moving molecules invalidates every stored molecule index, and
  // repairing only the bound partners reachable through the moved molecule's
  // own bndpartner -- which is what this loop used to do -- leaves the rest
  // pointing at a slot the pop_back() below takes away: a complex memberList,
  // a crossing list, or an interface whose partnerIndex the bndpartner list
  // does not mirror.  Record the move here and apply it everywhere afterwards.
  const int oldMolCount{static_cast<int>(moleculeList.size())};
  std::vector<int> newMolIndex(oldMolCount);
  for (int i{0}; i < oldMolCount; ++i) newMolIndex[i] = i;
  for (auto& emptied : Molecule::emptyMolList)
    if (emptied >= 0 && emptied < oldMolCount) newMolIndex[emptied] = -1;

  // Always copy last occupied element from moleculeList to the first empty
  // index:
  lastNonEmptyIndex = moleculeList.size() - 1;
  // While there is a molecule to delete:
  for (auto& firstEmptyIndex : Molecule::emptyMolList) {
    // cout << "firstEmptyIndex=" << firstEmptyIndex << endl;
    if (DEBUG && (firstEmptyIndex >= moleculeList.size()))
      error("10: firstEmptyIndex(" + to_string(firstEmptyIndex) +
            ") >= moleculeList.size() (" + to_string(moleculeList.size()) +
            ")");
    // Find the last non-empty molecule:
    while (moleculeList[lastNonEmptyIndex].myComIndex == -1)
      lastNonEmptyIndex--;
    if (lastNonEmptyIndex <= firstEmptyIndex) break;
    // Move last non-empty molecule to the first empty position:
    moleculeList[firstEmptyIndex] = moleculeList[lastNonEmptyIndex];
    // Update molecule index to match new position:
    moleculeList[firstEmptyIndex].index = firstEmptyIndex;
    newMolIndex[lastNonEmptyIndex] = firstEmptyIndex;
    lastNonEmptyIndex--;
  }

  // Remove last elements in moleculeList
  for (int k = 0; k < Molecule::emptyMolList.size(); k++)
    moleculeList.pop_back();
  Molecule::emptyMolList.clear();

  // Apply the map to every stored molecule index.  Destinations are slots that
  // were empty, so an index nothing moved into maps to -1 and the reference is
  // dropped rather than left naming whoever took the slot.
  auto mapped = [&newMolIndex, oldMolCount](int idx) {
    return (idx < 0 || idx >= oldMolCount) ? -1 : newMolIndex[idx];
  };
  for (auto& com : complexList) {
    if (com.isEmpty) continue;
    for (auto& memIdx : com.memberList) memIdx = mapped(memIdx);
    com.memberList.erase(
        std::remove(com.memberList.begin(), com.memberList.end(), -1),
        com.memberList.end());
  }
  for (auto& subBox : simulVolume.subCellList) {
    for (auto& memIdx : subBox.memberMolList) memIdx = mapped(memIdx);
    subBox.memberMolList.erase(std::remove(subBox.memberMolList.begin(),
                                           subBox.memberMolList.end(), -1),
                               subBox.memberMolList.end());
  }
  for (auto& mol : moleculeList) {
    if (mol.isEmpty) continue;
    // -1 already means "this interface has no partner visible on this rank",
    // so a partner that is gone lands on the sentinel the rest of the code
    // already understands.
    for (auto& iface : mol.interfaceList)
      iface.interaction.partnerIndex = mapped(iface.interaction.partnerIndex);
    // bndpartner and bndlist are parallel, so an entry leaves both.
    for (int i{static_cast<int>(mol.bndpartner.size()) - 1}; i >= 0; --i) {
      const int remapped{mapped(mol.bndpartner[i])};
      if (remapped == -1) {
        mol.bndpartner.erase(mol.bndpartner.begin() + i);
        if (i < static_cast<int>(mol.bndlist.size()))
          mol.bndlist.erase(mol.bndlist.begin() + i);
      } else {
        mol.bndpartner[i] = remapped;
      }
    }
    // A dropped crossing is uncounted from its complex, the way
    // remove_partner_crossings() uncounts the ones it removes.
    for (int i{static_cast<int>(mol.crossings.size()) - 1}; i >= 0; --i) {
      const int remapped{mapped(mol.crossings[i].partner)};
      if (remapped == -1) {
        mol.crossings.erase(mol.crossings.begin() + i);
        if (mol.myComIndex >= 0 &&
            mol.myComIndex < static_cast<int>(complexList.size()))
          --complexList[mol.myComIndex].ncross;
      } else {
        mol.crossings[i].partner = remapped;
      }
    }
  }
}