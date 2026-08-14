/*! \file reaction_bookkeeping.cpp
 * \brief The two chores every reaction performs once it has fired.
 *
 * Whatever a reaction did geometrically, it finishes the same way: it takes
 * itself out of the running for the rest of the timestep, and it moves whatever
 * observable counts it into place.  Both were written out at every site that
 * completes a reaction - eleven copies of the first between ten files, four of
 * the second - and the copies differ only in which molecule or which reaction
 * they name.
 */
#include "reactions/shared_reaction_functions.hpp"

void record_crossing_pair(int pro1, int pro2, int relIface1, int relIface2,
    const std::array<int, 3>& crossRxn, bool alsoInitProbvec, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList)
{
    // Both molecules record each other, the interface each of them presents, and
    // the reaction they might do; both complexes count one more crossing.  Every
    // copy of this pushed the same crossRxn triple onto both molecules.
    moleculeList[pro1].crossbase.push_back(pro2);
    moleculeList[pro2].crossbase.push_back(pro1);
    moleculeList[pro1].mycrossint.push_back(relIface1);
    moleculeList[pro2].mycrossint.push_back(relIface2);
    moleculeList[pro1].crossrxn.push_back(crossRxn);
    moleculeList[pro2].crossrxn.push_back(crossRxn);
    ++complexList[moleculeList[pro1].myComIndex].ncross;
    ++complexList[moleculeList[pro2].myComIndex].ncross;

    // The volume-exclusion sites also seed probvec, so that the entry they just
    // added stays index-aligned with crossbase.  get_distance() does not: its
    // caller pushes the probability it is about to compute.
    if (alsoInitProbvec) {
        moleculeList[pro1].probvec.push_back(0);
        moleculeList[pro2].probvec.push_back(0);
    }
}

void zero_partner_probvec(
    const Molecule& mol, std::vector<Molecule>& moleculeList, bool skipImplicitLipidPartners)
{
    // Every molecule that listed `mol` as a possible partner this timestep holds
    // a reaction probability for it.  Zeroing those stops `mol` reacting twice,
    // while leaving the partners free to keep avoiding overlap with it.
    for (unsigned crossItr { 0 }; crossItr < mol.crossbase.size(); ++crossItr) {
        int skipMol { mol.crossbase[crossItr] };
        if (skipImplicitLipidPartners && moleculeList[skipMol].isImplicitLipid)
            continue;
        for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
            if (moleculeList[skipMol].crossbase[crossItr2] == mol.index)
                moleculeList[skipMol].probvec[crossItr2] = 0;
        }
    }
}

void count_unimolecular_observable(bool isStateChangeBackRxn, int rxnIndex,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    std::map<std::string, int>& observablesList)
{
    // Unlike update_state_change_observable() this indexes *both* lists with the
    // same rxnIndex, and counts up whichever direction fired.  That is what the
    // four unimolecular sites did; whether it is intended is not recorded.
    bool isObserved { false };
    std::string observeLabel {};
    if (!isStateChangeBackRxn) {
        isObserved = forwardRxns[rxnIndex].isObserved;
        observeLabel = forwardRxns[rxnIndex].observeLabel;
    } else {
        isObserved = backRxns[rxnIndex].isObserved;
        observeLabel = backRxns[rxnIndex].observeLabel;
    }
    if (isObserved) {
        auto observeItr = observablesList.find(observeLabel);
        if (observeItr == observablesList.end()) {
            // std::cerr << "WARNING: Observable " << observeLabel << " not defined.\n";
        } else {
            ++observeItr->second;
        }
    }
}

void update_state_change_observable(bool isStateChangeBackRxn, int rxnIndex, int backRxnIndex,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    std::map<std::string, int>& observablesList)
{
    // TODO: Temporarily, if backRxn, iterate down, if forwardRxn, iterate up
    if (!isStateChangeBackRxn && forwardRxns[rxnIndex].isObserved) {
        auto observeItr = observablesList.find(forwardRxns[rxnIndex].observeLabel);
        if (observeItr == observablesList.end()) {
            // std::cerr << "WARNING: Observable " << forwardRxns[rxnIndex].observeLabel << " not defined.\n";
        } else {
            ++observeItr->second;
        }
    } else if (isStateChangeBackRxn && backRxns[backRxnIndex].isObserved) {
        auto observeItr = observablesList.find(backRxns[backRxnIndex].observeLabel);
        if (observeItr == observablesList.end()) {
            // std::cerr << "WARNING: Observable " << backRxns[rxnIndex].observeLabel << " not defined.\n";
        } else {
            --observeItr->second;
        }
    }
}
