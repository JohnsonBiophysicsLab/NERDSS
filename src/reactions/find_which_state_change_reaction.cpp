#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

void find_which_state_change_reaction(int ifaceIndex, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Molecule& reactMol, const Interface::State& currState, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns)
{
    // TRACE();
    for (auto rxnItr : currState.stateChangeRxns) {
        const ForwardRxn& oneRxn = forwardRxns[rxnItr.first]; // rxnItr.first is the reaction's index in forwardRxns
        if (oneRxn.rxnType == ReactionType::uniMolStateChange) {
            if (reactMol.interfaceList[ifaceIndex].index == oneRxn.reactantListNew[0].absIfaceIndex) {
                int matches { 0 };
                for (const auto& oneRate : oneRxn.rateList) {
                    if (hasIntangibles(0, reactMol, oneRate)) {
                        ++matches;
                        rateIndex = static_cast<int>(&oneRate - &oneRxn.rateList[0]);
                    }
                }

                if (matches == 0) {
                    // TODO: either an error or the molecule doesn't match all the requirements. do what?
                } else {
                    // With several matches the last one found is used, which is
                    // what the loop above already left in rateIndex.
                    rxnIndex = oneRxn.relRxnIndex;
                }
            } else if (oneRxn.conjBackRxnIndex != -1
                && reactMol.interfaceList[ifaceIndex].index == oneRxn.productListNew[0].absIfaceIndex) {
                // REVERSE UNIMOLECULAR STATE CHANGE
                //
                // The molecule carries this reaction's product state, so the
                // conjugate BackRxn is the one that applies.  Two defects were
                // corrected here; both left the branch unable to fire, so no
                // model that works today can regress:
                //   1. rateIndex, correctly computed in the loop below, was
                //      then overwritten with a reaction index, and rxnIndex was
                //      never assigned at all.  Callers require rxnIndex != -1,
                //      so every reverse state change was discarded, and had one
                //      survived it would have indexed the rate list with a
                //      reaction index;
                //   2. backRxns[conjBackRxnIndex] was read without testing the
                //      -1 sentinel, so an irreversible state change whose
                //      product state matched read out of bounds.
                const BackRxn& backRxn = backRxns[oneRxn.conjBackRxnIndex];
                int matches { 0 };
                for (const auto& oneRate : backRxn.rateList) {
                    if (hasIntangibles(0, reactMol, oneRate)) {
                        ++matches;
                        rateIndex = static_cast<int>(&oneRate - &backRxn.rateList[0]);
                    }
                }

                if (matches == 0) {
                    // TODO: either an error or the molecule doesn't match all the requirements. do what?
                } else {
                    // Both callers read backRxns[rxnIndex].rateList[rateIndex]
                    // when isStateChangeBackRxn is set, so rxnIndex indexes
                    // backRxns here.  relRxnIndex equals conjBackRxnIndex by
                    // construction; it is spelled this way to mirror the
                    // forward branch above.  Note this differs from the
                    // bimolecular convention, where perform_bimolecular_reactions()
                    // tests forwardRxns[rxnIndex].rxnType and so requires a
                    // forward index.
                    rxnIndex = backRxn.relRxnIndex;
                    isStateChangeBackRxn = true;
                }
            } else {
                continue;
            }
        }
    }
}
