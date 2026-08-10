#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

void find_which_reaction(int ifaceIndex1, int ifaceIndex2, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Interface::State& currState, const Molecule& reactMol1, const Molecule& reactMol2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    // Both flags are properties of the Molecules, not of the reaction, so they
    // are read once instead of once per candidate reaction.
    const bool mol1IsImplicitLipid { reactMol1.isImplicitLipid };
    const bool mol2IsImplicitLipid { reactMol2.isImplicitLipid };
    const int absIface1 { reactMol1.interfaceList[ifaceIndex1].index };
    const int absIface2 { reactMol2.interfaceList[ifaceIndex2].index };

    for (auto rxnItr : currState.myForwardRxns) {
        const ForwardRxn& oneRxn = forwardRxns[rxnItr];
        // see if we can find both of the reactants in the reaction's reactantList
        int reactIndex1 { -1 };
        int reactIndex2 { -1 };
        for (std::size_t reactItr { 0 }; reactItr < oneRxn.reactantListNew.size(); ++reactItr) {
            const RxnIface& oneReactant = oneRxn.reactantListNew[reactItr];
            if (mol1IsImplicitLipid && molTemplateList[oneReactant.molTypeIndex].isImplicitLipid) {
                if (reactIndex1 == -1) {
                    reactIndex1 = static_cast<int>(reactItr);
                    continue;
                }
            } else {
                if (absIface1 == oneReactant.absIfaceIndex) {
                    if (reactIndex1 == -1) {
                        reactIndex1 = static_cast<int>(reactItr);
                        continue;
                    }
                }
            }
            if (mol2IsImplicitLipid && molTemplateList[oneReactant.molTypeIndex].isImplicitLipid) {
                if (reactIndex2 == -1) {
                    reactIndex2 = static_cast<int>(reactItr);
                    continue;
                }
            } else {
                if (absIface2 == oneReactant.absIfaceIndex) {
                    if (reactIndex2 == -1) {
                        reactIndex2 = static_cast<int>(reactItr);
                        continue;
                    }
                }
            }
        }

        if (mol2IsImplicitLipid) {
            reactIndex2 = 1;
        }
        if (mol1IsImplicitLipid) {
            reactIndex1 = 1;
        }

        // REVERSE BIMOLECULAR STATE CHANGE
        //
        // The pair did not match this reaction's reactants, so test whether it
        // matches its products instead, meaning the conjugate BackRxn applies.
        //
        // Five defects were corrected here (see issue #8).  Every one of them
        // made this branch either dead or undefined, so no working model can
        // regress:
        //   1. the gate required conjBackRxnIndex > 0, but 0 is a valid index
        //      and -1 is the sentinel.  With a single reversible reaction pair
        //      -- the only arrangement in which the downstream index
        //      convention is self-consistent -- the branch never ran at all;
        //   2. the rate index was recovered by subtracting a pointer into the
        //      BackRxn rate list from a pointer into the ForwardRxn rate list,
        //      which is undefined behavior for unrelated containers;
        //   3. the single-match case never assigned rateIndex, so callers
        //      rejected the reaction (they require rateIndex != -1);
        //   4. the multiple-match case returned without selecting a match;
        //   5. the zero-match case read element 0 of an empty match list.
        if (oneRxn.conjBackRxnIndex != -1) {
            if (oneRxn.rxnType == ReactionType::biMolStateChange && (reactIndex1 == -1 || reactIndex2 == -1)) {
                const std::vector<RxnBase::RateState>& backRateList = backRxns[oneRxn.conjBackRxnIndex].rateList;
                bool noBackRateMatched { false };

                /*Why is this allowed here, using the products?*/
                for (std::size_t prodItr { 0 }; prodItr < oneRxn.productListNew.size(); ++prodItr) {
                    if (absIface1 == oneRxn.productListNew[prodItr].absIfaceIndex) {
                        reactIndex1 = static_cast<int>(prodItr);
                    }
                    if (absIface2 == oneRxn.productListNew[prodItr].absIfaceIndex) {
                        reactIndex2 = static_cast<int>(prodItr);
                    }
                    if (reactIndex1 != -1 && reactIndex2 != -1) {
                        const RateMatch match {
                            best_matching_rate(backRateList, reactIndex1, reactIndex2, reactMol1, reactMol2)
                        };

                        if (match.matches == 0) {
                            // No back rate state accepts this pair, so this
                            // reaction does not apply.  Leave rxnIndex and
                            // rateIndex untouched and try the next reaction
                            // rather than reading an empty match list.
                            noBackRateMatched = true;
                            break;
                        }

                        // rateIndex indexes backRxns[conjBackRxnIndex].rateList,
                        // which is what isStateChangeBackRxn tells the caller.
                        rateIndex = match.bestRateIndex;
                        rxnIndex = oneRxn.relRxnIndex;
                        isStateChangeBackRxn = true;
                        return;
                    }
                }

                if (noBackRateMatched)
                    continue;
            }
        } //This loop should not be attempted if there is no conjBackRxn

        // ORDINARY FORWARD REACTION
        if (reactIndex1 != -1 && reactIndex2 != -1) {
            RateMatch match { best_matching_rate(oneRxn.rateList, reactIndex1, reactIndex2, reactMol1, reactMol2) };

            if (match.matches == 0) {
                // if there are no matching rates and the reaction symmetric (i.e. interface 1 binding to interface 1,
                // swap the reactants and check their ancillary interfaces
                if (oneRxn.rxnType == ReactionType::bimolecular
                    && (oneRxn.intReactantList[0] == oneRxn.intReactantList[1])) {
                    match = best_matching_rate(oneRxn.rateList, reactIndex2, reactIndex1, reactMol1, reactMol2);
                    if (match.matches == 0)
                        return;
                } else {
                    return;
                }
            }

            // One match, or the match with the most required ancillary
            // interfaces when several apply.
            rateIndex = match.bestRateIndex;
            rxnIndex = oneRxn.relRxnIndex;
            ++totMatches;
            return;
        }
    }
}
