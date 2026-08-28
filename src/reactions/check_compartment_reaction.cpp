#include "classes/class_Membrane.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <algorithm>

// pro2Index is a lipid's index.

void check_compartment_reaction(int pro1Index, int pro2Index, int simItr,
    const Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<TransmissionRxn>& transmissionRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs)
{

    /*
    * Calculate distance between molecule and compartment.
    * Establish whether entering or exiting.  (Determined by molecule type.)
    */


    int pro1MolType = moleculeList[pro1Index].molTypeIndex;

    if (molTemplateList[pro1MolType].crossesCompartment == false) {
        return;
    }

    // Find distance between molecule's COM and compartment origin.
    bool isEntering = false;
    double distanceToOrigin = moleculeList[pro1Index].comCoord.length();
    if (distanceToOrigin > membraneObject.compartmentR) {
        isEntering = true;
    }

    // Entering and exiting were two ~30-line branches that differed in exactly
    // two expressions: which way round the signed distance to the compartment
    // surface is taken, and which probability routine is called.  Everything
    // between - the reaction lookup, the whole Dtot construction, Rmax, and the
    // multi-molecule boundary-condition escape - was duplicated character for
    // character.
    //
    // The distance is written as a ternary rather than as `sign * (a - b)`
    // because the two are not identical when the molecule sits exactly on the
    // surface: `-1.0 * (a - a)` is -0.0 where `a - a` is +0.0, and that value is
    // handed to the probability routines.
    const int rxnIndex { molTemplateList[pro1MolType].transmissionRxnIndex };

    // NOTE: this is the single-complex average `(x + y + z) / 3`, not the
    // pairwise weighted_D_sum() form.  They associate differently and so do not
    // agree to the last bit; do not substitute one for the other.
    double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x +
                               complexList[moleculeList[pro1Index].myComIndex].D.y +
                               complexList[moleculeList[pro1Index].myComIndex].D.z);
    double cf = cos(sqrt(4.0 * complexList[moleculeList[pro1Index].myComIndex].Dr.z * params.timeStep));
    double Dr1;
    int relIface1 { transmissionRxns[rxnIndex].reactantListNew[0].relIfaceIndex };
    double relIfaceDistance = moleculeList[pro1Index].interfaceList[relIface1].coord.length();

    Vec3D ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
        - complexList[moleculeList[pro1Index].myComIndex].comCoord };
    double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y + ifaceVec.z * ifaceVec.z };
    Dr1 = 2.0 * magMol1 * (1.0 - cf);
    Dtot += Dr1 / (6.0 * params.timeStep);
    Dtot += membraneObject.droplet.D;

    // this definition of Rmax may not need bindRadius depending on the binding model.
    const double distToCompartment { isEntering
            ? relIfaceDistance - membraneObject.compartmentR
            : membraneObject.compartmentR - relIfaceDistance };
    const double Rmax { 3.0 * sqrt(6.0 * Dtot * params.timeStep) + transmissionRxns[rxnIndex].bindRadius };

    // Kept as `< Rmax` with the else, rather than inverted to an early return:
    // `!(a < b)` and `a >= b` part company if either operand is ever NaN.
    if (distToCompartment < Rmax) {
        // Checking if molecule is part of a complex and if so applying boundary condition of compartment
        if (complexList[moleculeList[pro1Index].myComIndex].memberList.size() > 1) {
            moleculeList[pro1Index].transmissionProb = 0;
            moleculeList[pro1Index].enforceCompartmentBC = true;
            return;
        }

        moleculeList[pro1Index].transmissionProb = 0;
        if (isEntering) {
            determine_entering_compartment_probability(
                distToCompartment, transmissionRxns, rxnIndex, pro1Index, moleculeList, Dtot, params, membraneObject);
        } else {
            // TODO: This is to be defined
            determine_exiting_compartment_probability(
                distToCompartment, transmissionRxns, rxnIndex, pro1Index, moleculeList, Dtot, params, membraneObject);
        }
    } else {
        moleculeList[pro1Index].transmissionProb = -1;
    }
}
