#include "system_setup/system_setup.hpp"

#include <algorithm>

std::vector<std::vector<double>> max_bind_radius_between_types(const std::vector<MolTemplate>& molTemplates,
    const std::vector<ForwardRxn>& forwardRxns)
{
    // Only a reaction between two interfaces can overlap; a unimolecular
    // state change has a single reactant.
    std::vector<std::vector<double>> maxBindRadius(molTemplates.size(), std::vector<double>(molTemplates.size(), 0.0));
    for (const auto& oneRxn : forwardRxns) {
        if (oneRxn.reactantListNew.size() != 2)
            continue;
        int type0 { oneRxn.reactantListNew[0].molTypeIndex };
        int type1 { oneRxn.reactantListNew[1].molTypeIndex };
        maxBindRadius[type0][type1] = std::max(maxBindRadius[type0][type1], oneRxn.bindRadius);
        maxBindRadius[type1][type0] = maxBindRadius[type0][type1];
    }
    return maxBindRadius;
}

bool are_molecules_in_vicinity(const Molecule& mol1, const Molecule& mol2, const std::vector<MolTemplate>& molTemplates,
    double maxBindRadius)
{
    Vec3D tmpVec { mol2.comCoord - mol1.comCoord };
    double distance = tmpVec.length();

    // Each interface lies within its template's radius of the molecule's
    // center, so two interfaces can be closer than a bindRadius only when the
    // centers are closer than both radii plus that bindRadius.  A fixed margin
    // of 2 nm here used to hide every overlap of a pair whose reaction has a
    // larger sigma.
    double reach { molTemplates[mol1.molTypeIndex].radius + molTemplates[mol2.molTypeIndex].radius + maxBindRadius };

    return distance < reach;
}
