#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"
#include <algorithm>
#include <cstdint>

namespace {
/*! \brief Bit `molTypeIndex`, or every bit for a type the mask cannot number.
 *
 * SimulVolume::add_member()'s convention, for the same reason: a type past 63
 * matches everything, which costs a test that cannot be skipped, never a
 * partner that should have been looked at.
 */
uint64_t type_bit(int molTypeIndex)
{
    return (molTypeIndex >= 0 && molTypeIndex < 64) ? (uint64_t(1) << molTypeIndex) : ~uint64_t(0);
}
} // namespace

bool moleculeOverlaps(const Parameters& params, SimulVolume& simulVolume, Molecule& createdMol,
    std::vector<Molecule>& moleculeList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // TRACE();
    // get which box the Molecule belongs to
    int xItr { int((createdMol.comCoord.x + membraneObject.waterBox.x / 2) / simulVolume.subCellSize.x) };
    int yItr { int((createdMol.comCoord.y + membraneObject.waterBox.y / 2) / simulVolume.subCellSize.y) };
    int zItr { int(-(createdMol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) / simulVolume.subCellSize.z) };
    int currBin
        = xItr + (yItr * simulVolume.numSubCells.x) + (zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y);

    // if (molTemplateList[createdMol.molTypeIndex].D.z == 0
    if ((molTemplateList[createdMol.molTypeIndex].isImplicitLipid || 
        molTemplateList[createdMol.molTypeIndex].isLipid)
        && std::abs(createdMol.comCoord.z) - std::abs((membraneObject.waterBox.z / 2)) > 1E-6) {
        //std::cerr << "Molecule " << createdMol.index << " of type " << molTemplateList[createdMol.molTypeIndex].molName
        //          << " is off the membrane. Writing coordinates and exiting.\n";
        return true;
    }

    // Now make sure the Molecule is still inside the box in all dimensions
    if (createdMol.comCoord.z > (membraneObject.waterBox.z / 2) || createdMol.comCoord.z + 1E-6 < -(membraneObject.waterBox.z / 2)) {
        //std::cout << "Molecule " << createdMol.index
        //          << " is outside simulation volume in the z-dimension, with center of mass coordinates ["
        //          << createdMol.comCoord << "]. Attempting to fit back into box.\n";
        return true;
    } else if (createdMol.comCoord.y > (membraneObject.waterBox.y / 2)
        || createdMol.comCoord.y + 1E-6 < -(membraneObject.waterBox.y / 2)) {
        //std::cout << "Molecule " << createdMol.index
        //          << " is outside simulation volume in the y-dimension, with center of mass coordinates ["
        //          << createdMol.comCoord << "]. Attempting to fit back into box.\n";
        return true;
    } else if (createdMol.comCoord.x > (membraneObject.waterBox.x / 2)
        || createdMol.comCoord.x + 1E-6 < -(membraneObject.waterBox.x / 2)) {
        //std::cout << "Molecule " << createdMol.index
        //          << " is outside simulation volume in the x-dimension, with center of mass coordinates ["
        //          << createdMol.comCoord << "]. Attempting to fit back into box.\n";
        return true;
    } else if (currBin > (simulVolume.numSubCells.tot) || currBin < 0) {
        //std::cout << "Molecule " << createdMol.index
        //          << " is outside simulation volume with center of mass coordinates [" << createdMol.comCoord
        //          << "]. Attempting to fit back into box.\n";
        return true;
    } else {
        // A coordinate exactly on the far face of the box puts its index one
        // SubBox past the end -- (x + L/2) / subCellSize.x is then numSubCells.x
        // exactly -- which currBin above absorbs into the next row rather than
        // rejecting.  Brought into range the way update_memberMolLists() brings
        // the same indices in, so that the walk below has three real SubBoxes on
        // each axis and add_member() cannot write past subCellList.
        xItr = std::min(std::max(xItr, 0), simulVolume.numSubCells.x - 1);
        yItr = std::min(std::max(yItr, 0), simulVolume.numSubCells.y - 1);
        zItr = std::min(std::max(zItr, 0), simulVolume.numSubCells.z - 1);
        currBin = xItr + (yItr * simulVolume.numSubCells.x)
            + (zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y);

        // Inside the box: the new molecule overlaps when one of its interfaces
        // and one of another molecule's could react with each other and sit
        // closer than that reaction's bindRadius, the criterion
        // generate_coordinates() uses for a new simulation.
        //
        // Only a reaction between two interfaces can overlap; a unimolecular
        // state change has a single reactant, and reactantListNew[1] is one past
        // the end of its list.  Of those, only the ones this molecule is itself
        // a reactant in can bring it within a bindRadius of anything, so the
        // reach and the partner types are taken from those alone: a created
        // molecule with no reactive interface -- create_destroy's A, every
        // molecule in a model with no bimolecular reaction -- then skips the
        // walk outright.
        double maxBindRadius { 0.0 };
        uint64_t partnerTypeMask { 0 };
        for (const auto& oneRxn : forwardRxns) {
            if (oneRxn.reactantListNew.size() != 2)
                continue;
            bool isReached { false };
            for (const auto& createdIface : createdMol.interfaceList) {
                if (isReactant(createdIface, createdMol, oneRxn.reactantListNew[0])) {
                    partnerTypeMask |= type_bit(oneRxn.reactantListNew[1].molTypeIndex);
                    isReached = true;
                }
                if (isReactant(createdIface, createdMol, oneRxn.reactantListNew[1])) {
                    partnerTypeMask |= type_bit(oneRxn.reactantListNew[0].molTypeIndex);
                    isReached = true;
                }
            }
            if (isReached)
                maxBindRadius = std::max(maxBindRadius, oneRxn.bindRadius);
        }
        const double createdRadius { molTemplateList[createdMol.molTypeIndex].radius };

        // The SubBox the molecule landed in and the 26 around it.  One SubBox is
        // not enough: a SubBox edge is params.rMaxLimit or a little more, and
        // rMaxLimit is at least bindRadius + both interfaces' distances from
        // their centers for every bimolecular reaction, so a partner within a
        // bindRadius of this place is usually in a neighbouring SubBox and never
        // in the same one by more than the fraction of the exclusion sphere that
        // happens to fit.  The other way round, the 3 x 3 x 3 block holds every
        // point within one SubBox edge of this one in each axis, so it holds
        // every partner that can be within rMaxLimit and therefore every partner
        // that can overlap.
        //
        // SubVolume::neighborList is no use here: it lists 13 neighbours, the
        // ones forward and up, which is what a search that visits every SubBox
        // in turn needs to see each pair once -- this one visits one SubBox and
        // needs both halves.
        //
        // Walking members directly, rather than each member's whole complex,
        // reaches the same molecules: update_memberMolLists() bins every
        // molecule that is neither emptied nor an implicit lipid by its own
        // center, so a complex's members are in the block exactly when they are
        // in range, however far the complex reaches.  The complex walk also
        // visited a complex once per member of it the SubBox held.
        //
        // A molecule no two-reactant reaction names cannot be within a
        // bindRadius of anything, and leaves partnerTypeMask empty, so the walk
        // does not run at all: create_destroy's A, and every molecule in a model
        // that declares no bimolecular reaction.
        for (int zOffset { -1 }; partnerTypeMask != 0 && zOffset <= 1; ++zOffset) {
            const int zBin { zItr + zOffset };
            if (zBin < 0 || zBin >= simulVolume.numSubCells.z)
                continue;
            for (int yOffset { -1 }; yOffset <= 1; ++yOffset) {
                const int yBin { yItr + yOffset };
                if (yBin < 0 || yBin >= simulVolume.numSubCells.y)
                    continue;
                for (int xOffset { -1 }; xOffset <= 1; ++xOffset) {
                    const int xBin { xItr + xOffset };
                    if (xBin < 0 || xBin >= simulVolume.numSubCells.x)
                        continue;

                    const int neighBin { xBin + (yBin * simulVolume.numSubCells.x)
                        + (zBin * simulVolume.numSubCells.x * simulVolume.numSubCells.y) };
                    for (int memMol : simulVolume.subCellList[neighBin].memberMolList) {
                        if (memMol < 0 || memMol >= int(moleculeList.size()))
                            continue;
                        const Molecule& partMol = moleculeList[memMol];
                        // The molecule being placed is not yet in any member
                        // list, but it can hold a slot an emptied one left
                        // behind, and would then overlap itself at every place
                        // it was ever drawn.  An emptied slot has no interfaces
                        // and an implicit lipid no position to overlap; the time
                        // step's overlap checks skip both too.
                        if (partMol.index == createdMol.index || partMol.isEmpty || partMol.isImplicitLipid)
                            continue;
                        if ((partnerTypeMask & type_bit(partMol.molTypeIndex)) == 0)
                            continue;

                        // Each interface lies within its template's radius of
                        // its molecule's center, so from this far apart no pair
                        // of them can be within any bindRadius.
                        Vec3D comVec { createdMol.comCoord - partMol.comCoord };
                        if (comVec.length()
                            >= createdRadius + molTemplateList[partMol.molTypeIndex].radius + maxBindRadius)
                            continue;

                        for (const auto& oneRxn : forwardRxns) {
                            if (oneRxn.reactantListNew.size() != 2)
                                continue;
                            for (const auto& createdIface : createdMol.interfaceList) {
                                for (const auto& partIface : partMol.interfaceList) {
                                    bool canReact { (isReactant(createdIface, createdMol, oneRxn.reactantListNew[0])
                                                        && isReactant(partIface, partMol, oneRxn.reactantListNew[1]))
                                        || (isReactant(createdIface, createdMol, oneRxn.reactantListNew[1])
                                            && isReactant(partIface, partMol, oneRxn.reactantListNew[0])) };
                                    if (canReact
                                        && Vec3D { createdIface.coord - partIface.coord }.length() < oneRxn.bindRadius)
                                        return true;
                                }
                            }
                        }
                    }
                }
            }
        }

        createdMol.mySubVolIndex = currBin;
        // Through add_member(), not subCellList: creation runs before
        // update_memberMolLists(), so a Molecule binned here into a SubBox that
        // the last re-binning pass left empty would otherwise be invisible to
        // clear_member_lists() and end up in memberMolList twice.
        simulVolume.add_member(currBin, createdMol.index, createdMol.molTypeIndex);
        return false;
    }
}

void create_molecule_and_complex_from_rxn(int parentMolIndex, int& newMolIndex, int& newComIndex, bool createInVicinity,
    MolTemplate& createdMolTemp, Parameters& params, const CreateDestructRxn& currRxn, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject)
{
    newMolIndex = 0;
    newComIndex = 0;
    if (Molecule::emptyMolList.size() > 0) {
        // Check to make sure the object pointed to by the last element in the each list is actually empty
        // If it isn't delete it from the list and move on until you find one that is
        try {
            while (!moleculeList[Molecule::emptyMolList.back()].isEmpty)
                Molecule::emptyMolList.pop_back();

            // if there's an available empty Molecule spot, make the new Molecule in place
            newMolIndex = Molecule::emptyMolList.back();
            Molecule::emptyMolList.pop_back(); // remove the empty Molecule spot index from the list
        } catch (std::out_of_range& e) {
            newMolIndex = moleculeList.size();
            moleculeList.emplace_back(); // create new empty Molecule spot
        }
    } else {
        newMolIndex = moleculeList.size();
        moleculeList.emplace_back(); // create new empty Molecule spot
    }

    // if (Complex::emptyComList.size() > 0) {
    //     // Check to make sure the object pointed to by the last element in the each list is actually empty
    //     // If it isn't delete it from the list and move on until you find one that is
    //     try {
    //         while (!complexList[Complex::emptyComList.back()].isEmpty)
    //             Complex::emptyComList.pop_back();
    //         // if there's an available empty Complex spot, make the new Complex in place
    //         newComIndex = Complex::emptyComList.back();
    //         Complex::emptyComList.pop_back(); // remove the empty Complex spot index from the list
    //     } catch (std::out_of_range) {
    //         newComIndex = complexList.size();
    //         complexList.emplace_back(); // create new empty Complex spot
    //     }

    // } else {
    //     newComIndex = complexList.size();
    //     complexList.emplace_back(); // create new empty Complex spot
    // }
    newComIndex = complexList.size();
    complexList.emplace_back();  // create new empty Complex spot

    // Now create the new species, and move it for as long as it overlaps a
    // molecule already in place.  Only the coordinates are drawn again: the
    // initialize functions also count the molecule, in numberOfMolecules,
    // numTotalUnits, maxID and numEachMolType, so re-running one per attempt
    // counted every molecule once per place it was tried.
    if (createInVicinity) {
        moleculeList[newMolIndex] = initialize_molecule_after_uni_reaction(
            newMolIndex, moleculeList[parentMolIndex], params, createdMolTemp, currRxn);
        while (moleculeOverlaps(params, simulVolume, moleculeList[newMolIndex], moleculeList, forwardRxns,
            molTemplateList, membraneObject)) {
            draw_coords_after_uni_reaction(
                moleculeList[newMolIndex], moleculeList[parentMolIndex], createdMolTemp, currRxn);
        }
    } else {
        moleculeList[newMolIndex]
            = initialize_molecule_after_zeroth_reaction(newMolIndex, params, createdMolTemp, currRxn, membraneObject);
        while (moleculeOverlaps(params, simulVolume, moleculeList[newMolIndex], moleculeList, forwardRxns,
            molTemplateList, membraneObject)) {
            draw_coords_after_zeroth_reaction(moleculeList[newMolIndex], createdMolTemp, currRxn, membraneObject);
        }
    }

    moleculeList[newMolIndex].myComIndex = newComIndex;
    moleculeList[newMolIndex].trajStatus = TrajStatus::propagated;
    moleculeList[newMolIndex].isGhosted = false;
    moleculeList[newMolIndex].isDissociated = true;
    complexList[newComIndex] = Complex { newComIndex, moleculeList.at(newMolIndex), createdMolTemp };
    complexList[newComIndex].trajStatus = TrajStatus::propagated;
    moleculeList[newMolIndex].complexId = complexList[newComIndex].id;
    ++Complex::numberOfComplexes;

    // add to monomerList if canDestroy = true
    {
        Molecule& oneMol { moleculeList[newMolIndex] };
        MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
        if (oneTemp.canDestroy) {
            oneTemp.monomerList.emplace_back(oneMol.index);
        }
    }
}

void MPI_create_molecule_and_complex_on_rank(
    Molecule& mol, int& newMolIndex, int& newComIndex,
    MolTemplate& createdMolTemp, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject) {
  newMolIndex = 0;
  newComIndex = 0;
  if (Molecule::emptyMolList.size() > 0) {
    // Check to make sure the object pointed to by the last element in the each
    // list is actually empty If it isn't delete it from the list and move on
    // until you find one that is
    // TODO: OPTIMIZE AT THE END: no need to try-catch, but rather loop while
    // there are elements, and after check wheter empty
    try {
      while (!moleculeList[Molecule::emptyMolList.back()].isEmpty)
        Molecule::emptyMolList.pop_back();

      // if there's an available empty Molecule spot, make the new Molecule in
      // place
      newMolIndex = Molecule::emptyMolList.back();
      Molecule::emptyMolList
          .pop_back();  // remove the empty Molecule spot index from the list
    } catch (std::out_of_range& e) {
      newMolIndex = moleculeList.size();
      moleculeList.emplace_back();  // create new empty Molecule spot
    }
  } else {
    newMolIndex = moleculeList.size();
    moleculeList.emplace_back();  // create new empty Molecule spot
  }

  if (Complex::emptyComList.size() > 0) {
    // Check to make sure the object pointed to by the last element in the each
    // list is actually empty If it isn't delete it from the list and move on
    // until you find one that is
    // TODO: OPTIMIZE AT THE END: no need to try-catch, but rather loop while
    // there are elements, and after check wheter empty
    try {
      while (!complexList[Complex::emptyComList.back()].isEmpty)
        Complex::emptyComList.pop_back();
      // if there's an available empty Complex spot, make the new Complex in
      // place
      newComIndex = Complex::emptyComList.back();
      Complex::emptyComList
          .pop_back();  // remove the empty Complex spot index from the list
    } catch (std::out_of_range) {
      newComIndex = complexList.size();
      complexList.emplace_back();  // create new empty Complex spot
    }

  } else {
    newComIndex = complexList.size();
    complexList.emplace_back();  // create new empty Complex spot
  }

  moleculeList[newMolIndex].myComIndex = newComIndex;
  moleculeList[newMolIndex].trajStatus = TrajStatus::propagated;

  moleculeList[newMolIndex].isGhosted = false;
  // moleculeList[newMolIndex].id = Molecule::maxID++;

  complexList[newComIndex] =
      Complex{newComIndex, moleculeList.at(newMolIndex), createdMolTemp};
  complexList[newComIndex].trajStatus = TrajStatus::propagated;

  moleculeList[newMolIndex].complexId = complexList[newComIndex].id;
  ++Complex::numberOfComplexes;
}
