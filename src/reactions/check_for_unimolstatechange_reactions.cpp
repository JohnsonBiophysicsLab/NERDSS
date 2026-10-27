#include "math/constants.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

namespace {

/*! \brief What one draw of a unimolecular state change decided.
 *
 * fired is false when no state change applies to the interface or the draw
 * did not come up; then nothing was counted and the other two fields are
 * unset.
 */
struct StateChangeDraw {
    bool fired { false };
    int relStateIndex { -1 }; //!< the product state's index in the interface's stateList
    const RxnIface* product { nullptr }; //!< the product state as the reaction names it
};

/*! \brief Offers one unbound interface with states its unimolecular state
 * change this timestep.
 *
 * Finds the reaction that takes currState away, draws it, and when it fires
 * looks the product up in stateList and counts the reaction's observable.
 * rand_gsl() is called exactly once when a reaction applies and not at all
 * when none does, which is what keeps this bit-identical to the two copies
 * of this sequence it replaces.
 *
 * The two scans in check_for_unimolstatechange_reactions() are the callers:
 * the implicit lipid once per free lipid in each state, an explicit molecule
 * once per interface with states.  What they do with a fired draw is not
 * shared and stays with them: the lipid moves one count between its
 * per-state tallies, the molecule changes the interface's state and marks
 * itself propagated.  Which copy-number species loses the count is theirs as
 * well - the lipid scan decrements the state's own index, the molecule scan
 * the state's absolute interface index - so this draws and reports, and
 * changes nothing itself.
 */
StateChangeDraw draw_unimolecular_state_change(int relIndex, const Molecule& mol,
    const Interface::State& currState, const std::vector<Interface::State>& stateList,
    const Parameters& params, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, std::map<std::string, int>& observablesList)
{
    StateChangeDraw draw {};

    int rxnIndex { -1 };
    int rateIndex { -1 };
    bool isStateChangeBackRxn { false };
    find_which_state_change_reaction(relIndex, rxnIndex, rateIndex, isStateChangeBackRxn,
        mol, currState, forwardRxns, backRxns);
    if (rxnIndex == -1 || rateIndex == -1)
        return draw;

    double rate { (!isStateChangeBackRxn) ? forwardRxns[rxnIndex].rateList[rateIndex].rate
                                          : backRxns[rxnIndex].rateList[rateIndex].rate };
    double prob { first_order_probability(rate, params.timeStep) };

    // make sure that the time step is resonable according to the prob of reaction
    if (prob > 1.000001) {
        std::cerr << "Error: prob of reaction > 1. Avoid this using a smaller time step." << std::endl;
        //exit(1);
    }

    if (!(prob > rand_gsl()))
        return draw;

    draw.fired = true;
    draw.product = (!isStateChangeBackRxn) ? &forwardRxns[rxnIndex].productListNew[0]
                                           : &backRxns[rxnIndex].productListNew[0];
    for (auto& state : stateList) {
        if (state.index == draw.product->absIfaceIndex) {
            draw.relStateIndex = static_cast<int>(&state - &stateList[0]);
            break;
        }
    }

    count_unimolecular_observable(isStateChangeBackRxn, rxnIndex, forwardRxns, backRxns, observablesList);
    return draw;
}

} // namespace

void check_for_unimolstatechange_reactions(unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    SimulVolume& simulVolume, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, Membrane& membraneObject)
{
    for (unsigned molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
        // only do checks if the Molecule exists
        if (moleculeList[molItr].isEmpty) continue;
        if (moleculeList[molItr].isGhosted) continue;
#ifdef mpi_
        // held back to a later pass in this step; see src/mpi/defer_edge_reactions.cpp
        if (moleculeList[molItr].firstOrderHeld) continue;
#endif

        //check for state changes
        if (params.hasUniMolStateChange == true) {
            if (molTemplateList[moleculeList[molItr].molTypeIndex].ifacesWithStates.empty() == false) {
                for (int relIndex : molTemplateList[moleculeList[molItr].molTypeIndex].ifacesWithStates) {
                    const auto& stateList = molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[relIndex].stateList;
                    if (moleculeList[molItr].isImplicitLipid == true) {
                        if (moleculeList[molItr].interfaceList[relIndex].isBound) {
                            // TODO
                        } else {
                            for (int tmpImplicitLipidStateIndex = 0; tmpImplicitLipidStateIndex < membraneObject.nStates; tmpImplicitLipidStateIndex++) {
                                int tmpOneStateLipid = membraneObject.numberOfFreeLipidsEachState[tmpImplicitLipidStateIndex];
                                for (int tmpItr = 0; tmpItr < tmpOneStateLipid; tmpItr++) {
                                    moleculeList[molItr].interfaceList[relIndex].index += tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].relIndex += tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].stateIndex += tmpImplicitLipidStateIndex;

                                    const StateChangeDraw draw { draw_unimolecular_state_change(relIndex, moleculeList[molItr],
                                        molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[0].stateList[tmpImplicitLipidStateIndex],
                                        stateList, params, forwardRxns, backRxns, observablesList) };
                                    if (draw.fired) {
                                        --counterArrays.copyNumSpecies[tmpImplicitLipidStateIndex];
                                        ++counterArrays.copyNumSpecies[draw.product->absIfaceIndex];
                                        --membraneObject.numberOfFreeLipidsEachState[tmpImplicitLipidStateIndex];
                                        ++membraneObject.numberOfFreeLipidsEachState[draw.relStateIndex];
                                    }
                                    moleculeList[molItr].interfaceList[relIndex].index -= tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].relIndex -= tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].stateIndex -= tmpImplicitLipidStateIndex;
                                }
                            }
                        }
                    } else {
                        if (moleculeList[molItr].trajStatus == TrajStatus::none) {
                            auto stateItr = std::find_if(stateList.begin(), stateList.end(), [&](const Interface::State& oneState) {
                                return oneState.iden == moleculeList[molItr].interfaceList[relIndex].stateIden;
                            });
                            if (stateItr == stateList.end()) {
                                std::cerr << "Error, state is invalid. Exiting...\n";
                                exit(1);
                            }

                            if (moleculeList[molItr].interfaceList[relIndex].isBound) {
                                // TODO
                            } else {
                                const StateChangeDraw draw { draw_unimolecular_state_change(relIndex, moleculeList[molItr],
                                    *stateItr, stateList, params, forwardRxns, backRxns, observablesList) };
                                if (draw.fired) {
                                    --counterArrays.copyNumSpecies[stateItr->index];
                                    ++counterArrays.copyNumSpecies[draw.product->absIfaceIndex];

                                    /*State change happens here.*/
                                    moleculeList[molItr].interfaceList[relIndex].change_state(
                                        draw.relStateIndex, draw.product->absIfaceIndex, draw.product->requiresState);
                                    moleculeList[molItr].refresh_iface_state_bits();

                                    moleculeList[molItr].trajStatus = TrajStatus::propagated;
                                    complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;
                                    break; // break out of interface loop
                                }
                            }
                        }
                    } // only try each pair dissociating once
                } // finished trying out stateChange reactions for this proteins
            }
        }
    }
}
