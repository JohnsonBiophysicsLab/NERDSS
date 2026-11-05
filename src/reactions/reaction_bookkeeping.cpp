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
#include "io/io.hpp"
#include "reactions/shared_reaction_functions.hpp"

#include <fstream>

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
    //
    // Each entry also records where the other one sits, so that the decision
    // loop can reach the partner's entry without searching for it.
    int itr1 { static_cast<int>(moleculeList[pro1].crossings.size()) };
    int itr2 { static_cast<int>(moleculeList[pro2].crossings.size()) };
    moleculeList[pro1].crossings.emplace_back(pro2, relIface1, crossRxn);
    moleculeList[pro2].crossings.emplace_back(pro1, relIface2, crossRxn);
    moleculeList[pro1].crossings.back().partnerItr = itr2;
    moleculeList[pro2].crossings.back().partnerItr = itr1;
    ++complexList[moleculeList[pro1].myComIndex].ncross;
    ++complexList[moleculeList[pro2].myComIndex].ncross;
}

void remove_from_monomer_list(const Molecule& oneMol, std::vector<MolTemplate>& molTemplateList)
{
    MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
    bool isMonomer { oneMol.bndpartner.empty() };
    bool canDestroy { oneTemp.canDestroy };
    if (isMonomer && canDestroy) {
        std::vector<int>& oneList { oneTemp.monomerList };
        std::vector<int>::iterator result { std::find(std::begin(oneList), std::end(oneList), oneMol.index) };
        if (result != std::end(oneList)) {
            oneList.erase(result);
        }
    }
}

void add_to_monomer_list(const Molecule& oneMol, std::vector<MolTemplate>& molTemplateList)
{
    MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
    bool isMonomer { oneMol.bndpartner.empty() };
    bool canDestroy { oneTemp.canDestroy };
    if (isMonomer && canDestroy) {
        oneTemp.monomerList.emplace_back(oneMol.index);
    }
}

void perform_coupled_state_change(const ForwardRxn& coupledRxn, int pro1Index, int pro2Index,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, copyCounters& counterArrays,
    std::map<std::string, int>& observablesList)
{
    int stateChangeProIndex { -1 };
    int relIndex { -1 };
    // make sure the molecule iface have the same state with the uniMolStateChange reaction's reactant
    for (auto& tmpIface : moleculeList[pro1Index].interfaceList) {
        if (tmpIface.index == coupledRxn.reactantListNew[0].absIfaceIndex) {
            stateChangeProIndex = pro1Index; // Is it A or B?
            relIndex = tmpIface.relIndex;
        }
    }
    for (auto& tmpIface : moleculeList[pro2Index].interfaceList) {
        if (tmpIface.index == coupledRxn.reactantListNew[0].absIfaceIndex) {
            stateChangeProIndex = pro2Index; // Is it A or B?
            relIndex = tmpIface.relIndex;
        }
    }

    if (stateChangeProIndex == -1) {
        std::cerr << "The products of the disscociation do not match the corresponding uniMolStateChange reactant." << std::endl;
        exit(1);
    }

    const auto& stateList = molTemplateList[moleculeList[stateChangeProIndex].molTypeIndex].interfaceList[relIndex].stateList;
    const auto& newState = coupledRxn.productListNew[0];
    int relStateIndex { -1 };
    for (auto& state : stateList) {
        if (state.index == newState.absIfaceIndex) {
            relStateIndex = static_cast<int>(&state - &stateList[0]);
            break;
        }
    }

    // check observables.  parse_reaction() enters every observeLabel it reads
    // into observablesList, so the warning is for a label that arrived some
    // other way; one of the three copies kept it, two had commented it out.
    if (coupledRxn.isObserved) {
        auto observeItr = observablesList.find(coupledRxn.observeLabel);
        if (observeItr == observablesList.end()) {
            std::cerr << "WARNING: Observable " << coupledRxn.observeLabel << " not defined.\n";
        } else {
            ++observeItr->second;
        }
    }

    --counterArrays.copyNumSpecies[coupledRxn.reactantListNew[0].absIfaceIndex];
    ++counterArrays.copyNumSpecies[coupledRxn.productListNew[0].absIfaceIndex];

    // Change the state happens here
    moleculeList[stateChangeProIndex].interfaceList[relIndex].change_state(
        relStateIndex, newState.absIfaceIndex, newState.requiresState);
    moleculeList[stateChangeProIndex].refresh_iface_state_bits();

    moleculeList[stateChangeProIndex].trajStatus = TrajStatus::propagated;
    complexList[moleculeList[stateChangeProIndex].myComIndex].trajStatus = TrajStatus::propagated;
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

void zero_reciprocal_crossings(const Molecule& mol, std::vector<Molecule>& moleculeList)
{
    // The same entries zero_partner_probvec() finds by scanning each partner's
    // list: every entry naming `mol` was recorded in a pair with one of mol's
    // own, duplicates from several interface pairs included.  A partner that
    // has reacted has an empty list, and the bound skips it as the scan found
    // nothing in it.
    for (const auto& cross : mol.crossings) {
        if (cross.partnerItr < 0)
            continue;
        auto& partnerCrossings { moleculeList[cross.partner].crossings };
        if (cross.partnerItr < static_cast<int>(partnerCrossings.size()))
            partnerCrossings[cross.partnerItr].prob = 0;
    }
}

void remove_partner_crossings(
    const Molecule& mol, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList)
{
    // A partner that `mol` lists more than once, through several interface
    // pairs, loses all of its entries for `mol` on the first visit and finds
    // none on the others.  The implicit lipid never lists anyone, so it loses
    // nothing.
    //
    // The compaction keeps the survivors in order, as std::remove_if did, and
    // tells the entry on the other side of each moved one where it went, which
    // is what keeps zero_reciprocal_crossings() exact after a transmission.
    for (unsigned crossItr { 0 }; crossItr < mol.crossings.size(); ++crossItr) {
        Molecule& partner { moleculeList[mol.crossings[crossItr].partner] };
        auto& list { partner.crossings };
        unsigned kept { 0 };
        for (unsigned k { 0 }; k < list.size(); ++k) {
            if (list[k].partner == mol.index)
                continue;
            if (kept != k) {
                list[kept] = list[k];
                int j { list[kept].partnerItr };
                auto& otherCrossings { moleculeList[list[kept].partner].crossings };
                if (j >= 0 && j < static_cast<int>(otherCrossings.size()))
                    otherCrossings[j].partnerItr = static_cast<int>(kept);
            }
            ++kept;
        }
        complexList[partner.myComIndex].ncross -= static_cast<int>(list.size() - kept);
        list.resize(kept);
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

void finalize_bond(long long int iter, int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, const ForwardRxn& currRxn, const Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList,
    copyCounters& counterArrays, std::ofstream& assocDissocFile)
{
    remove_from_monomer_list(reactMol1, molTemplateList);
    remove_from_monomer_list(reactMol2, molTemplateList);

    // update Molecule interface statuses
    reactMol1.interfaceList[ifaceIndex1].interaction.partnerIndex = reactMol2.index;
    reactMol2.interfaceList[ifaceIndex2].interaction.partnerIndex = reactMol1.index;
    reactMol1.interfaceList[ifaceIndex1].interaction.partnerIfaceIndex = ifaceIndex2;
    reactMol2.interfaceList[ifaceIndex2].interaction.partnerIfaceIndex = ifaceIndex1;
    if (currRxn.isReversible) {
      reactMol1.interfaceList[ifaceIndex1].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
      reactMol2.interfaceList[ifaceIndex2].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
    }
    if (assocDissocFile.is_open()) {
        assocDissocFile << "ITR:" << iter << "," << "BOND," 
        << molTemplateList[reactMol1.molTypeIndex].molName << "," << assoc_dissoc_label(reactMol1) << "," << ifaceIndex1 << "," 
        << molTemplateList[reactMol2.molTypeIndex].molName << "," << assoc_dissoc_label(reactMol2) << "," << ifaceIndex2 << std::endl;
    }
    reactMol1.interfaceList[ifaceIndex1].isBound = true;
    reactMol1.interfaceList[ifaceIndex1].index = currRxn.productListNew[0].absIfaceIndex;
    reactMol2.interfaceList[ifaceIndex2].isBound = true;
    reactMol2.interfaceList[ifaceIndex2].index = currRxn.productListNew[0].absIfaceIndex;
    reactMol1.refresh_iface_state_bits();
    reactMol2.refresh_iface_state_bits();

    // add to the list of bound interfaces and remove from the list of free
    // interfaces
    reactMol1.bndlist.push_back(ifaceIndex1);
    reactMol2.bndlist.push_back(ifaceIndex2);
    reactMol1.bndpartner.push_back(reactMol2.index);
    reactMol2.bndpartner.push_back(reactMol1.index);

    {
      size_t tmpItr{0};
      for (unsigned i{0}; i < reactMol1.freelist.size(); ++i) {
        if (reactMol1.freelist[i] == ifaceIndex1)
          tmpItr = i;
      }
      reactMol1.freelist[tmpItr] = reactMol1.freelist.back();
      reactMol1.freelist.pop_back();
    }
    {
      size_t tmpItr{0};
      for (unsigned i{0}; i < reactMol2.freelist.size(); ++i) {
        if (reactMol2.freelist[i] == ifaceIndex2)
          tmpItr = i;
      }
      reactMol2.freelist[tmpItr] = reactMol2.freelist.back();
      reactMol2.freelist.pop_back();
    }

    // Set probability of this protein to zero in all reactions so it doesn't try
    // to react again but the partners still will avoid overlapping.
    zero_partner_probvec(reactMol1, moleculeList);

    // Set probability of this protein to zero in all reactions so it doesn't try
    // to react again but the partners still will avoid overlapping.
    zero_partner_probvec(reactMol2, moleculeList);

    // Update the crossed molecule lists so that the current molecules won't avoid
    // anything, but others will.
    reactCom1.ncross = -1;
    reactMol1.crossings.clear();
    reactMol2.crossings.clear();

    /*update the number of bound species*/
    update_Nboundpairs(reactMol1.molTypeIndex, reactMol2.molTypeIndex, 1, params, counterArrays);
    /*Update species copy numbers*/
    counterArrays.copyNumSpecies[currRxn.reactantListNew[0].absIfaceIndex]--; // decrement ifaceIndex1
    counterArrays.copyNumSpecies[currRxn.reactantListNew[1].absIfaceIndex]--; // decrement ifaceIndex2
    counterArrays.copyNumSpecies[currRxn.productListNew[0].absIfaceIndex]++; // increment product state

    //-------------------------UPDATE BINDPAIRLIST----------------------------
    if (counterArrays.canDissociate[currRxn.productListNew[0].absIfaceIndex]) {
      // int molIndex { reactMol1.index };
      // if (reactMol1.index > reactMol2.index)
      //     molIndex = reactMol2.index;
      int molIndex{-1};
      if (ifaceIndex1 == currRxn.reactantListNew[0].relIfaceIndex) {
        molIndex = reactMol1.index;
      } else {
        molIndex = reactMol2.index;
      }

      counterArrays.bindPairList[currRxn.productListNew[0].absIfaceIndex]
          .emplace_back(molIndex);
    }
    //-----------------------END UPDATE BINDPAIRLIST--------------------------

    // std::cout << " After ASSOCIATE, CHANGE COPY NUMBERS, interfaces: " <<
    // ifaceIndex1 << ' ' << ifaceIndex2
    //           << " add to product: " << currRxn.productListNew[0].absIfaceIndex
    //           << " sub from reactants: " <<
    //           currRxn.reactantListNew[0].absIfaceIndex << " "
    //           << currRxn.reactantListNew[1].absIfaceIndex << std::endl;

    // species tracking here
    if (currRxn.isObserved) {
      auto obsItr = observablesList.find(currRxn.observeLabel);
      if (obsItr != observablesList.end())
        ++obsItr->second;
    }
}
