/*! \file determine_compartment_probability.cpp
 * \brief Probability that a molecule crosses the compartment boundary this step.
 *
 * Entering and exiting were two files with the same forty lines in them,
 * differing in one call: `prob_entering_compartment` versus
 * `prob_exiting_compartment`.  Everything that decides whether the crossing is
 * even considered - the molecule must not have just dissociated, the reaction
 * must have a positive rate - and everything that builds the `paramsIL` handed
 * to it was written twice.
 */
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "tracing.hpp"

namespace {

//! \brief Which side of the compartment boundary the molecule is trying to cross.
enum class CrossingDirection { entering, exiting };

/*! \brief Stores the crossing probability on the molecule, or leaves it alone.
 *
 * `pro1Index` is the molecule trying to cross; the implicit lipid of the
 * droplet is stored in `membraneObject`.
 */
void determine_compartment_probability(CrossingDirection direction, double distToCompartment,
    const std::vector<TransmissionRxn>& transmissionRxns, int rxnIndex, int pro1Index,
    std::vector<Molecule>& moleculeList, double Dtot, const Parameters& parameters, Membrane& membraneObject)
{
    // TRACE();
    /*3D reaction*/

    // This movestat check is if you allow just dissociated proteins to avoid overlap
    if (moleculeList[pro1Index].isDissociated == true)
        return;
    if (!(transmissionRxns[rxnIndex].rateList[0].rate > 0))
        return;

    // declare intrinsic binding rate of 3D->2D case.
    double ktemp { 2.0 * transmissionRxns[rxnIndex].rateList[0].rate };

    paramsIL params3D {};
    params3D.R2D = 0.0;
    params3D.sigma = transmissionRxns[rxnIndex].bindRadius;
    params3D.Dtot = Dtot;
    params3D.ka = ktemp;
    params3D.dt = parameters.timeStep;
    params3D.compartmentR = membraneObject.compartmentR;
    params3D.compartSiteRho = membraneObject.droplet.rho; // This is to be defined

    double rxnProb { (direction == CrossingDirection::entering)
            ? prob_entering_compartment(distToCompartment, params3D)
            : prob_exiting_compartment(distToCompartment, params3D) };

    if (rxnProb > 1.000001) {
        std::cerr << "Error: prob of reaction is: " << rxnProb << " > 1. Avoid this using a smaller time step." << std::endl;
        // exit(1);
    }

    moleculeList[pro1Index].transmissionProb = rxnProb;
}

} // namespace

void determine_entering_compartment_probability(double distToCompartment, const std::vector<TransmissionRxn>& transmissionRxns,
    int rxnIndex, int pro1Index, std::vector<Molecule>& moleculeList, double Dtot, const Parameters& parameters,
    Membrane& membraneObject)
{
    determine_compartment_probability(CrossingDirection::entering, distToCompartment, transmissionRxns, rxnIndex,
        pro1Index, moleculeList, Dtot, parameters, membraneObject);
}

void determine_exiting_compartment_probability(double distToCompartment, const std::vector<TransmissionRxn>& transmissionRxns,
    int rxnIndex, int pro1Index, std::vector<Molecule>& moleculeList, double Dtot, const Parameters& parameters,
    Membrane& membraneObject)
{
    determine_compartment_probability(CrossingDirection::exiting, distToCompartment, transmissionRxns, rxnIndex,
        pro1Index, moleculeList, Dtot, parameters, membraneObject);
}
