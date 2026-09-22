#include <algorithm>
#include <iostream>

#include "parser/parser_functions.hpp"

/*
 * `counterArrays.bindPairList` is indexed by species, like `copyNumSpecies`:
 * every interface state first, in template order, then one product per
 * bimolecular reaction. Unlike `copyNumSpecies` it is not recounted on a
 * restart; `read_restart()` restores it exactly as written, so it arrives in the
 * species order the restart file was written with.
 *
 * An add file changes that order. Its interface states are numbered after the
 * existing ones, which moves every existing product up by `numStateAdd`, and its
 * own products come after all of those. The add block renumbers molecules and
 * reactions to match, so the lists have to move with them. Left in place, there
 * are fewer lists than species and each existing product's list sits at another
 * species' index: the first association into, or dissociation of, a product
 * past the end reads and writes outside `bindPairList`. That is silent heap
 * corruption, which surfaces thousands of steps later as a std::length_error or
 * a segfault.
 */
void reindex_bindPairList_for_add(copyCounters& counterArrays, const std::vector<ForwardRxn>& forwardRxns,
    int lastStateIndexBeforeAdd, int numStateAdd, int numDoubleBeforeAdd)
{
    std::vector<std::vector<int>>& bindPairList { counterArrays.bindPairList };
    const int numStateBeforeAdd { lastStateIndexBeforeAdd + 1 };

    if (static_cast<int>(bindPairList.size()) != numStateBeforeAdd + numDoubleBeforeAdd) {
        std::cerr << "ERROR: the restart file holds bound-pair lists for " << bindPairList.size()
                  << " species, but its molecule types and reactions define " << numStateBeforeAdd
                  << " interface states and " << numDoubleBeforeAdd << " products. Exiting...\n";
        exit(1);
    }

    const int numDoubleAfterAdd { static_cast<int>(std::count_if(forwardRxns.begin(), forwardRxns.end(),
        [](const ForwardRxn& oneRxn) { return oneRxn.rxnType == ReactionType::bimolecular; })) };

    // the added interface states start with no bound pairs, and move the
    // existing products up by numStateAdd
    bindPairList.insert(bindPairList.begin() + numStateBeforeAdd, numStateAdd, std::vector<int> {});
    // the added products start empty too, after all of them
    bindPairList.resize(bindPairList.size() + numDoubleAfterAdd - numDoubleBeforeAdd);
}
