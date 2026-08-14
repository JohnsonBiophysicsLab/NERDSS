/*! \file implicit_lipid_state_lookup.cpp
 * \brief Finding which implicit-lipid state a reaction consumes or releases.
 *
 * Whenever an implicit-lipid reaction fires, the free-lipid count for the state
 * it used has to be adjusted, and getting to that count takes three steps that
 * six call sites each wrote out in full:
 *
 *   1. of the reaction's two interfaces, decide which is the implicit lipid;
 *   2. fetch the implicit lipid's state list;
 *   3. turn the interface's absolute index into a position in that list.
 *
 * Step 2 is not always the same list - `perform_implicitlipid_state_change_*`
 * looks at the *changing* molecule's interface rather than the lipid's - so it
 * stays a separate function that those callers can skip.
 */
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

const RxnIface& implicit_lipid_iface(
    const std::vector<RxnIface>& ifacePair, const std::vector<MolTemplate>& molTemplateList)
{
    if (molTemplateList[ifacePair[1].molTypeIndex].isImplicitLipid == true)
        return ifacePair[1];
    return ifacePair[0];
}

int implicit_lipid_state_index(const std::vector<Interface::State>& stateList, const RxnIface& iface)
{
    for (auto& state : stateList) {
        if (state.index == iface.absIfaceIndex)
            return static_cast<int>(&state - &stateList[0]);
    }
    return -1;
}

const std::vector<Interface::State>& implicit_lipid_state_list(const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    return molTemplateList[moleculeList[membraneObject.implicitlipidIndex].molTypeIndex].interfaceList[0].stateList;
}
