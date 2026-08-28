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
    const std::array<int, 3>& crossRxn, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList)
{
    // Both molecules record each other, the interface each of them presents, and
    // the reaction they might do; both complexes count one more crossing.  Every
    // copy of this pushed the same crossRxn triple onto both molecules.
    //
    // The probability starts at zero and the caller overwrites it through
    // crossings.back().  This used to be an `alsoInitProbvec` flag, because the
    // four lists were separate and only some callers seeded probvec here while
    // get_distance()'s caller pushed it afterwards.  One entry carries all four
    // values, so there is nothing left to keep in step and the flag is gone.
    moleculeList[pro1].crossings.emplace_back(pro2, relIface1, crossRxn);
    moleculeList[pro2].crossings.emplace_back(pro1, relIface2, crossRxn);
    ++complexList[moleculeList[pro1].myComIndex].ncross;
    ++complexList[moleculeList[pro2].myComIndex].ncross;
}

void zero_partner_probvec(
    const Molecule& mol, std::vector<Molecule>& moleculeList, bool skipImplicitLipidPartners)
{
    // Every molecule that listed `mol` as a possible partner this timestep holds
    // a reaction probability for it.  Zeroing those stops `mol` reacting twice,
    // while leaving the partners free to keep avoiding overlap with it.
    for (unsigned crossItr { 0 }; crossItr < mol.crossings.size(); ++crossItr) {
        int skipMol { mol.crossings[crossItr].partner };
        if (skipImplicitLipidPartners && moleculeList[skipMol].isImplicitLipid)
            continue;
        for (auto& partnerCross : moleculeList[skipMol].crossings) {
            if (partnerCross.partner == mol.index)
                partnerCross.prob = 0;
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

size_t find_bond_slot(const Molecule& mol, int relIface)
{
    for (size_t slot { 0 }; slot < mol.bndlist.size(); ++slot) {
        if (mol.bndlist[slot] == relIface)
            return slot;
    }
    return mol.bndlist.size();
}

bool erase_bond(Molecule& mol, int relIface)
{
    const size_t slot { find_bond_slot(mol, relIface) };
    if (slot == mol.bndlist.size())
        return false;
    mol.bndlist.erase(mol.bndlist.begin() + slot);
    // bndpartner can be shorter than bndlist only if some other site has
    // already broken the invariant; bound the access rather than trust it.
    if (slot < mol.bndpartner.size())
        mol.bndpartner.erase(mol.bndpartner.begin() + slot);
    return true;
}
