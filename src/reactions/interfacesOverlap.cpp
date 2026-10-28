#include "reactions/shared_reaction_functions.hpp"

bool interfacesOverlap(const Molecule& mol1, const Molecule& mol2, const std::vector<ForwardRxn>& forwardRxns,
    double maxBindRadius)
{
    // Distance first: most interface pairs that get this far are still farther
    // apart than any bindRadius between the two molecules, and the reaction
    // scan below is the expensive part.  The distance is the one every caller
    // compared before this was shared, length() against bindRadius, so the
    // answer is the same to the last bit.
    for (const auto& iface1 : mol1.interfaceList) {
        for (const auto& iface2 : mol2.interfaceList) {
            const double distance { Vec3D { iface1.coord - iface2.coord }.length() };
            if (distance >= maxBindRadius)
                continue;

            // Only a reaction between two interfaces can overlap; a
            // unimolecular state change has a single reactant.
            for (const auto& oneRxn : forwardRxns) {
                if (oneRxn.reactantListNew.size() != 2 || distance >= oneRxn.bindRadius)
                    continue;
                if ((isReactant(iface1, mol1, oneRxn.reactantListNew[0])
                        && isReactant(iface2, mol2, oneRxn.reactantListNew[1]))
                    || (isReactant(iface1, mol1, oneRxn.reactantListNew[1])
                        && isReactant(iface2, mol2, oneRxn.reactantListNew[0])))
                    return true;
            }
        }
    }
    return false;
}
