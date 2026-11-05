/*! \file state_change_bookkeeping.cpp
 * \brief What a state change does once its geometry is settled.
 *
 * The box and sphere versions of the bimolecular state change ended with the
 * same fifty lines, and so did the box and sphere versions of the
 * implicit-lipid state change: commit the new state, make the temporary
 * coordinates real, reflect, zero the pair's probabilities, move the
 * observable and species counts, clear the crossings.  Each pair now shares
 * one function.  The statements are the copies' own, in their order, so the
 * outputs do not change; box against sphere stays a difference of the callers.
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"

void commit_bimolecular_state_change(int stateChangeIface, int relStateIndex, const RxnIface& newState,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    bool isStateChangeBackRxn, int rxnIndex, int backRxnIndex, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList,
    copyCounters& counterArrays, const Membrane& membraneObject)
{
    // change the state
    // copyNumSpecies is indexed by absolute interface (state) index.  stateIndex
    // is the state's position within its own interface, so it used to move the
    // counts of whichever species happen to have absolute index 0 and 1 (in
    // regression/bimol_state_change, K(k) and K(kd)) and never the changed ones.
    int oldStateIndex { stateChangeMol.interfaceList[stateChangeIface].index };
    stateChangeMol.interfaceList[stateChangeIface].change_state(
        relStateIndex, newState.absIfaceIndex, newState.requiresState);
    stateChangeMol.refresh_iface_state_bits();

    // write temporary to real coords and clear temporary coordinates
    for (auto memMol : facilitatorCom.memberList) {
        moleculeList[memMol].commit_tmp_association_coords();
        moleculeList[memMol].trajStatus = TrajStatus::propagated;
    }
    facilitatorCom.update_properties(moleculeList, molTemplateList);

    for (auto memMol : stateChangeCom.memberList) {
        moleculeList[memMol].commit_tmp_association_coords();
        moleculeList[memMol].trajStatus = TrajStatus::propagated;
    }
    stateChangeCom.update_properties(moleculeList, molTemplateList);

    // Enforce boundary conditions
    reflect_complex_rad_rot(membraneObject, facilitatorCom, moleculeList, 0.0, false);

    zero_partner_probvec(stateChangeMol, moleculeList);

    // Set probability of this protein to zero in all reactions so it doesn't try to
    // react again but the partners still will avoid overlapping.
    zero_partner_probvec(facilitatorMol, moleculeList);

    // update observables, if applicable
    // TODO: Temporarily, if backRxn, iterate down, if forwardRxn, iterate up
    update_state_change_observable(isStateChangeBackRxn, rxnIndex, backRxnIndex, forwardRxns, backRxns, observablesList);

    /*Update species copy numbers*/
    counterArrays.copyNumSpecies[oldStateIndex]--; // decrement ifaceIndex1
    counterArrays.copyNumSpecies[stateChangeMol.interfaceList[stateChangeIface].index]++; // increment product state

    // Update the crossed molecule lists so that the current molecules won't avoid anything, but others will.
    stateChangeCom.ncross = -1;
    facilitatorCom.ncross = -1;
    stateChangeMol.crossings.clear();
    facilitatorMol.crossings.clear();
}

void commit_implicitlipid_state_change(int stateChangeIface, int relStateIndex, const ForwardRxn& currRxn,
    double RS3D, Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom,
    Complex& facilitatorCom, bool isStateChangeBackRxn, int rxnIndex, int backRxnIndex,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    std::vector<Molecule>& moleculeList, std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, Membrane& membraneObject)
{
    // change the state. for IL, update species copy numbers
    int oldStateIndex { stateChangeMol.interfaceList[stateChangeIface].stateIndex };
    //    stateChangeMol.interfaceList[stateChangeIface].change_state(
    //        relStateIndex, newState.absIfaceIndex, newState.requiresState);

    // A back state change runs currRxn from its products to its reactants, so
    // it consumes the forward products and releases the forward reactants.
    // Both used to follow the forward reaction whichever way it ran.
    const std::vector<RxnIface>& consumedSide { isStateChangeBackRxn ? currRxn.productListNew
                                                                     : currRxn.reactantListNew };
    const std::vector<RxnIface>& releasedSide { isStateChangeBackRxn ? currRxn.reactantListNew
                                                                     : currRxn.productListNew };
    counterArrays.copyNumSpecies[consumedSide[0].absIfaceIndex] -= 1;
    counterArrays.copyNumSpecies[consumedSide[1].absIfaceIndex] -= 1;
    counterArrays.copyNumSpecies[releasedSide[0].absIfaceIndex] += 1;
    counterArrays.copyNumSpecies[releasedSide[1].absIfaceIndex] += 1;

    //Update free number of lipids of IL for each state
    // The state list here is the changing molecule's own interface, not the
    // implicit lipid's, which is why implicit_lipid_state_list() is not used.
    const auto& implicitLipidStateList = molTemplateList[stateChangeMol.molTypeIndex].interfaceList[stateChangeIface].stateList;

    // These two lookups keep the previous relStateIndex when the state is absent,
    // as they always have; only the first lookup in this function resets it.
    int consumed { implicit_lipid_state_index(
        implicitLipidStateList, implicit_lipid_iface(consumedSide, molTemplateList)) };
    if (consumed != -1)
        relStateIndex = consumed;
    membraneObject.numberOfFreeLipidsEachState[relStateIndex] -= 1;

    int released { implicit_lipid_state_index(
        implicitLipidStateList, implicit_lipid_iface(releasedSide, molTemplateList)) };
    if (released != -1)
        relStateIndex = released;
    membraneObject.numberOfFreeLipidsEachState[relStateIndex] += 1;

    // write temporary to real coords and clear temporary coordinates
    for (auto memMol : facilitatorCom.memberList) {
        moleculeList[memMol].commit_tmp_association_coords();
        moleculeList[memMol].trajStatus = TrajStatus::propagated;
    }
    facilitatorCom.update_properties(moleculeList, molTemplateList);

    //for (auto memMol : stateChangeCom.memberList) {
    //    moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
    //    for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
    //        moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
    //    moleculeList[memMol].clear_tmp_association_coords();
    //    moleculeList[memMol].trajStatus = TrajStatus::propagated;
    //}
    //stateChangeCom.update_properties(moleculeList, molTemplateList);

    // clear implicit lipid properties
    moleculeList[stateChangeMol.index].clear_tmp_association_coords();
    stateChangeCom.memberList.clear();
    stateChangeCom.memberList.push_back(stateChangeMol.index);
    stateChangeCom.update_properties(moleculeList, molTemplateList);

    // Enforce boundary conditions
    reflect_complex_rad_rot(membraneObject, facilitatorCom, moleculeList, RS3D, false);

    //for (unsigned crossItr { 0 }; crossItr < stateChangeMol.crossings.size(); ++crossItr) {
    //    int skipMol { stateChangeMol.crossings[crossItr].partner };
    //    for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossings.size(); ++crossItr2) {
    //        if (moleculeList[skipMol].crossings[crossItr2].partner == stateChangeMol.index)
    //            moleculeList[skipMol].crossings[crossItr2].prob = 0;
    //    }
    //}

    // Set probability of this protein to zero in all reactions so it doesn't try to
    // react again but the partners still will avoid overlapping.
    zero_partner_probvec(facilitatorMol, moleculeList);

    // update observables, if applicable
    // TODO: Temporarily, if backRxn, iterate down, if forwardRxn, iterate up
    update_state_change_observable(isStateChangeBackRxn, rxnIndex, backRxnIndex, forwardRxns, backRxns, observablesList);

    /*Update species copy numbers*/
    //counterArrays.copyNumSpecies[oldStateIndex]--; // decrement ifaceIndex1
    //counterArrays.copyNumSpecies[stateChangeMol.interfaceList[stateChangeIface].stateIndex]++; // increment product state

    // Update the crossed molecule lists so that the current molecules won't avoid anything, but others will.
    stateChangeCom.ncross = -1;
    facilitatorCom.ncross = -1;
    stateChangeMol.crossings.clear();
    facilitatorMol.crossings.clear();
}
