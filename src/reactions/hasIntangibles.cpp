#include "reactions/shared_reaction_functions.hpp"

namespace {

/*!
 * \brief Tests one ancillary-interface requirement against a Molecule.
 *
 * RxnIface::relIfaceIndex is documented as an index into MolTemplate's (and
 * therefore Molecule's) interfaceList, and every site that populates a
 * Molecule sets interfaceList[i].relIndex = i.  So the requirement can be
 * looked up directly instead of scanning the whole interface list for it,
 * which turns the cost of hasIntangibles() from
 * O(requirements x interfaces) into O(requirements).
 *
 * The direct lookup is guarded two ways so a malformed model behaves exactly
 * as it did before: an out-of-range relIfaceIndex simply fails to match, and
 * if a Molecule ever violates the relIndex == position invariant the original
 * linear scan is used instead.
 */
bool hasAncillaryIface(const Molecule& reactMol, const RxnIface& anccIface)
{
    if (anccIface.relIfaceIndex >= 0
        && static_cast<std::size_t>(anccIface.relIfaceIndex) < reactMol.interfaceList.size()) {
        const Molecule::Iface& oneIface = reactMol.interfaceList[anccIface.relIfaceIndex];
        if (oneIface.relIndex == anccIface.relIfaceIndex) {
            // relIndex is unique within a Molecule, so this is the only
            // interface the scan below could ever have matched.
            return anccIface.molTypeIndex == oneIface.molTypeIndex
                && anccIface.requiresInteraction == oneIface.isBound
                && anccIface.requiresState == oneIface.stateIden;
        }
    }

    for (const auto& oneIface : reactMol.interfaceList) {
        if (anccIface.molTypeIndex == oneIface.molTypeIndex && anccIface.relIfaceIndex == oneIface.relIndex
            && anccIface.requiresInteraction == oneIface.isBound && anccIface.requiresState == oneIface.stateIden) {
            return true;
        }
    }
    return false;
}

} // namespace

bool hasIntangibles(int reactantIndex, const Molecule& reactMol, const RxnBase::RateState& currRxnState)
{
    for (const auto& anccIface : currRxnState.otherIfaceLists[reactantIndex]) {
        // if the ancillary interface was found in either of the two reactants, iterate the number of matches
        // TODO: Does it matter if there are matches in both reactants?
        if (!hasAncillaryIface(reactMol, anccIface))
            return false;
    }
    return true;
}

bool hasIntangibles(int reactIndex1, int reactIndex2, const Molecule& reactMol1, const Molecule& reactMol2,
    const RxnBase::RateState& currRxnState)
{
    if (currRxnState.otherIfaceLists[0].empty() && currRxnState.otherIfaceLists[1].empty()) {
        return true;
    }

    for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex1]) {
        if (!hasAncillaryIface(reactMol1, anccIface))
            return false;
    }
    for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex2]) {
        if (!hasAncillaryIface(reactMol2, anccIface))
            return false;
    }

    return true;
}
