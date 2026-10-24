#include "reactions/shared_reaction_functions.hpp"

void find_which_reaction(int ifaceIndex1, int ifaceIndex2, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Interface::State& currState, const Molecule& reactMol1, const Molecule& reactMol2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<MolTemplate>& molTemplateList)
{
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

        // The implicit lipid keeps the reactant position the loop above found
        // for it.  It used to be forced to position 1 here, which is right only
        // when the reaction lists the lipid second: for `IL(m) + B(b!*,m)` the
        // protein is at position 1 as well, so both reactants were matched
        // against the same requirement list, and the protein's bound-b
        // requirement was looked for on the lipid's interfaces and never found.

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
                const BackRxn& backRxn = backRxns[oneRxn.conjBackRxnIndex];
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
                            best_matching_rate(backRxn, reactIndex1, reactIndex2, reactMol1, reactMol2)
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
            RateMatch match { best_matching_rate(oneRxn, reactIndex1, reactIndex2, reactMol1, reactMol2) };

            if (match.matches == 0) {
                // if there are no matching rates and the reaction symmetric (i.e. interface 1 binding to interface 1,
                // swap the reactants and check their ancillary interfaces
                if (oneRxn.rxnType == ReactionType::bimolecular
                    && (oneRxn.intReactantList[0] == oneRxn.intReactantList[1])) {
                    match = best_matching_rate(oneRxn, reactIndex2, reactIndex1, reactMol1, reactMol2);
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

namespace {

/*! \brief best_matching_rate()'s answer when it cannot depend on the molecules, else -1.
 *
 * hasIntangibles() reads a molecule only for the requirements listed under the
 * pair's two reactant positions.  When those lists are empty for every rate
 * state, every state matches, and best_matching_rate() keeps the first one
 * with the most ancillary interfaces -- which it counts over both lists, so
 * that count is still taken here.
 */
int fixed_rate_index(const std::vector<RxnBase::RateState>& rateList, int reactIndex1, int reactIndex2)
{
    if (rateList.empty())
        return -1;
    int best { -1 };
    std::size_t mostAncillary { 0 };
    for (std::size_t rateItr { 0 }; rateItr < rateList.size(); ++rateItr) {
        const RxnBase::RateState& oneRate = rateList[rateItr];
        const bool noneListed { oneRate.otherIfaceLists[0].empty() && oneRate.otherIfaceLists[1].empty() };
        if (!noneListed
            && (!oneRate.otherIfaceLists[reactIndex1].empty() || !oneRate.otherIfaceLists[reactIndex2].empty()))
            return -1;
        const std::size_t ancillaryCount { oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size() };
        if (best == -1 || ancillaryCount > mostAncillary) {
            best = static_cast<int>(rateItr);
            mostAncillary = ancillaryCount;
        }
    }
    return best;
}

/*! \brief What find_which_reaction()'s reactant and product scans reach for one
 * (state, partner interface) pair, where the first molecule is not the
 * implicit lipid.  With lipidIsSecond the second one is, and is matched as
 * find_which_reaction() matches it, by molecule type (molTemplateList is read
 * only then).
 *
 * The same loops as above, with the rate tests taken out: a reaction whose
 * reactants or products the pair matches becomes a candidate with the
 * positions found, in myForwardRxns order.  At run time a matched back state
 * change either returns or moves on to the next reaction, so the scan
 * continues past it; an ordinary forward match returns whatever its rates say,
 * so nothing after it is ever reached and the scan stops there.
 */
std::vector<Interface::State::RxnCandidate> scan_rxn_candidates(const Interface::State& currState, int absIface1,
    int absIface2, bool lipidIsSecond, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<MolTemplate>* molTemplateList)
{
    std::vector<Interface::State::RxnCandidate> candidates;
    for (auto rxnItr : currState.myForwardRxns) {
        const ForwardRxn& oneRxn = forwardRxns[rxnItr];
        int reactIndex1 { -1 };
        int reactIndex2 { -1 };
        for (std::size_t reactItr { 0 }; reactItr < oneRxn.reactantListNew.size(); ++reactItr) {
            const RxnIface& oneReactant = oneRxn.reactantListNew[reactItr];
            if (absIface1 == oneReactant.absIfaceIndex) {
                if (reactIndex1 == -1) {
                    reactIndex1 = static_cast<int>(reactItr);
                    continue;
                }
            }
            if (lipidIsSecond && (*molTemplateList)[oneReactant.molTypeIndex].isImplicitLipid) {
                if (reactIndex2 == -1) {
                    reactIndex2 = static_cast<int>(reactItr);
                    continue;
                }
            } else if (absIface2 == oneReactant.absIfaceIndex) {
                if (reactIndex2 == -1) {
                    reactIndex2 = static_cast<int>(reactItr);
                    continue;
                }
            }
        }

        if (oneRxn.conjBackRxnIndex != -1 && oneRxn.rxnType == ReactionType::biMolStateChange
            && (reactIndex1 == -1 || reactIndex2 == -1)) {
            bool matchedProducts { false };
            for (std::size_t prodItr { 0 }; prodItr < oneRxn.productListNew.size(); ++prodItr) {
                if (absIface1 == oneRxn.productListNew[prodItr].absIfaceIndex)
                    reactIndex1 = static_cast<int>(prodItr);
                if (absIface2 == oneRxn.productListNew[prodItr].absIfaceIndex)
                    reactIndex2 = static_cast<int>(prodItr);
                if (reactIndex1 != -1 && reactIndex2 != -1) {
                    candidates.push_back({ rxnItr, reactIndex1, reactIndex2, true, false,
                        fixed_rate_index(backRxns[oneRxn.conjBackRxnIndex].rateList, reactIndex1, reactIndex2) });
                    matchedProducts = true;
                    break;
                }
            }
            if (matchedProducts)
                continue;
        }

        if (reactIndex1 != -1 && reactIndex2 != -1) {
            const bool trySwap { oneRxn.rxnType == ReactionType::bimolecular
                && (oneRxn.intReactantList[0] == oneRxn.intReactantList[1]) };
            candidates.push_back({ rxnItr, reactIndex1, reactIndex2, false, trySwap,
                fixed_rate_index(oneRxn.rateList, reactIndex1, reactIndex2) });
            break;
        }
    }
    return candidates;
}

/*! \brief The rate tests find_which_reaction() makes for one (state, partner)
 * pair, in its order and with its outcomes, from the pair's candidates.
 */
inline void test_rxn_candidates(const std::vector<Interface::State::RxnCandidate>& candidates, int& rxnIndex,
    int& rateIndex, bool& isStateChangeBackRxn, const Molecule& reactMol1, const Molecule& reactMol2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns)
{
    for (const Interface::State::RxnCandidate& candidate : candidates) {
        const ForwardRxn& oneRxn = forwardRxns[candidate.rxnItr];
        if (candidate.fixedRateIndex != -1) {
            // Every rate state matches whatever the molecules are doing.
            rateIndex = candidate.fixedRateIndex;
            rxnIndex = oneRxn.relRxnIndex;
            isStateChangeBackRxn = candidate.backStateChange;
            if (!candidate.backStateChange)
                ++totMatches;
            return;
        }
        if (candidate.backStateChange) {
            const RateMatch match { best_matching_rate(backRxns[oneRxn.conjBackRxnIndex],
                candidate.reactIndex1, candidate.reactIndex2, reactMol1, reactMol2) };
            if (match.matches == 0)
                continue;
            rateIndex = match.bestRateIndex;
            rxnIndex = oneRxn.relRxnIndex;
            isStateChangeBackRxn = true;
            return;
        }

        RateMatch match { best_matching_rate(
            oneRxn, candidate.reactIndex1, candidate.reactIndex2, reactMol1, reactMol2) };
        if (match.matches == 0) {
            if (!candidate.trySwap)
                return;
            match = best_matching_rate(oneRxn, candidate.reactIndex2, candidate.reactIndex1, reactMol1, reactMol2);
            if (match.matches == 0)
                return;
        }
        rateIndex = match.bestRateIndex;
        rxnIndex = oneRxn.relRxnIndex;
        ++totMatches;
        return;
    }
}

} // namespace

void find_which_reaction_for_partner(int ifaceIndex1, std::size_t partnerItr, int& rxnIndex,
    int& rateIndex, bool& isStateChangeBackRxn, const Interface::State& currState, const Molecule& reactMol1,
    const Molecule& reactMol2, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns)
{
    if (!currState.rxnCandidatesBuilt) {
        // A state belongs to one template interface, so absIface1 is the same
        // for every molecule that reaches it.
        const int absIface1 { reactMol1.interfaceList[ifaceIndex1].index };
        currState.rxnCandidatesByPartner.clear();
        for (auto partner : currState.rxnPartners)
            currState.rxnCandidatesByPartner.push_back(scan_rxn_candidates(
                currState, absIface1, static_cast<int>(partner), false, forwardRxns, backRxns, nullptr));
        currState.rxnCandidatesBuilt = true;
    }

    test_rxn_candidates(currState.rxnCandidatesByPartner[partnerItr], rxnIndex, rateIndex, isStateChangeBackRxn,
        reactMol1, reactMol2, forwardRxns, backRxns);
}

void find_which_reaction_for_lipid_partner(int ifaceIndex1, std::size_t partnerItr, int& rxnIndex, int& rateIndex,
    bool& isStateChangeBackRxn, const Interface::State& currState, const Molecule& reactMol1, const Molecule& lipid,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<MolTemplate>& molTemplateList)
{
    if (!currState.lipidRxnCandidatesBuilt) {
        // The lipid's interface index is the partner itself: that is what the
        // state offset made it.
        const int absIface1 { reactMol1.interfaceList[ifaceIndex1].index };
        currState.lipidRxnCandidatesByPartner.clear();
        for (auto partner : currState.rxnPartners)
            currState.lipidRxnCandidatesByPartner.push_back(scan_rxn_candidates(
                currState, absIface1, static_cast<int>(partner), true, forwardRxns, backRxns, &molTemplateList));
        currState.lipidRxnCandidatesBuilt = true;
    }

    test_rxn_candidates(currState.lipidRxnCandidatesByPartner[partnerItr], rxnIndex, rateIndex,
        isStateChangeBackRxn, reactMol1, lipid, forwardRxns, backRxns);
}
