#include "io/io.hpp"
#include "reactions/shared_reaction_functions.hpp"

#ifdef NERDSS_CHECK_RATE_STATE_CACHE
#include <cstdlib>
#endif

namespace {

#ifdef NERDSS_CHECK_RATE_STATE_CACHE
// How many selections the both-ways build compared, reported at exit.
struct SelectionTally {
    unsigned long long checked { 0 };
    ~SelectionTally() { std::cerr << "NERDSS_CHECK_RATE_STATE_CACHE: " << checked << " selections agree\n"; }
} selectionTally;
#endif

// The reference form: scans the reaction's reactants for the bond's two
// interfaces, then selects the rate state.
int find_reaction_rate_state_scan(int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1,
    const Molecule& reactMol2, const BackRxn& backRxn, const std::vector<MolTemplate>& molTemplateList)
{
    int reactIndex1 { -1 };
    int reactIndex2 { -1 };
    for (std::size_t reactItr { 0 }; reactItr < backRxn.reactantListNew.size(); ++reactItr) {
        const RxnIface& oneReactant = backRxn.reactantListNew[reactItr];
        if (reactMol1.molTypeIndex == oneReactant.molTypeIndex && relIfaceIndex1 == oneReactant.relIfaceIndex) {
            if (reactIndex1 == -1) {
                reactIndex1 = static_cast<int>(reactItr);
                continue;
            }
        }
        if (reactMol2.molTypeIndex == oneReactant.molTypeIndex && relIfaceIndex2 == oneReactant.relIfaceIndex) {
            if (reactIndex2 == -1) {
                reactIndex2 = static_cast<int>(reactItr);
                continue;
            }
        }
    }

    if (reactIndex1 == -1 || reactIndex2 == -1) {
        std::cerr << llinebreak;
        std::cerr << "ERROR: Association product has incorrect dissociation reaction.\n";
        reactMol1.display(molTemplateList[reactMol1.molTypeIndex]);
        std::cerr << linebreak;
        reactMol2.display(molTemplateList[reactMol2.molTypeIndex]);
        std::cerr << llinebreak;
        backRxn.display();
        std::cout << std::endl;
        exit(1);
    }

    // Selection rule: one match is used directly, and if several rate states
    // match, the one with the most required ancillary interfaces wins.  Both
    // are handled by the single streaming scan, so no temporary match list is
    // built and rescanned.
    const RateMatch match { best_matching_rate(backRxn, reactIndex1, reactIndex2, reactMol1, reactMol2) };
    if (match.matches == 0)
        return -1;

    return match.bestRateIndex;
}

// The entry for this bond kind, made on first use; nullptr when the reactant
// scan does not find both interfaces, which the scan reports.
const BackRxn::RatePositions* rate_positions(
    int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1, const Molecule& reactMol2, const BackRxn& backRxn)
{
    std::vector<BackRxn::RatePositions>& entries = backRxn.ratePositions.entries;
    for (const BackRxn::RatePositions& entry : entries) {
        if (entry.molType1 == reactMol1.molTypeIndex && entry.relIface1 == relIfaceIndex1
            && entry.molType2 == reactMol2.molTypeIndex && entry.relIface2 == relIfaceIndex2
            && entry.rateCount == backRxn.rateList.size())
            return &entry;
    }

    // The scan above, keeping the positions.
    BackRxn::RatePositions entry;
    entry.molType1 = reactMol1.molTypeIndex;
    entry.relIface1 = relIfaceIndex1;
    entry.molType2 = reactMol2.molTypeIndex;
    entry.relIface2 = relIfaceIndex2;
    for (std::size_t reactItr { 0 }; reactItr < backRxn.reactantListNew.size(); ++reactItr) {
        const RxnIface& oneReactant = backRxn.reactantListNew[reactItr];
        if (entry.molType1 == oneReactant.molTypeIndex && relIfaceIndex1 == oneReactant.relIfaceIndex) {
            if (entry.reactIndex1 == -1) {
                entry.reactIndex1 = static_cast<int>(reactItr);
                continue;
            }
        }
        if (entry.molType2 == oneReactant.molTypeIndex && relIfaceIndex2 == oneReactant.relIfaceIndex) {
            if (entry.reactIndex2 == -1) {
                entry.reactIndex2 = static_cast<int>(reactItr);
                continue;
            }
        }
    }
    if (entry.reactIndex1 == -1 || entry.reactIndex2 == -1)
        return nullptr;
    entry.fixedRateIndex = fixed_rate_index(backRxn.rateList, entry.reactIndex1, entry.reactIndex2);
    entry.rateCount = backRxn.rateList.size();
    // A stale entry (rateList changed length) is replaced rather than kept.
    for (BackRxn::RatePositions& old : entries) {
        if (old.molType1 == entry.molType1 && old.relIface1 == entry.relIface1 && old.molType2 == entry.molType2
            && old.relIface2 == entry.relIface2) {
            old = entry;
            return &old;
        }
    }
    entries.push_back(entry);
    return &entries.back();
}

} // namespace

int find_reaction_rate_state(int simItr, int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1,
    const Molecule& reactMol2, const BackRxn& backRxn, const std::vector<MolTemplate>& molTemplateList)
{
    // The reactant scan depends only on the bond kind, so its outcome is kept
    // in backRxn.ratePositions; only a rate list that reads the molecules'
    // other interfaces is still tested per bond.
    const BackRxn::RatePositions* positions { rate_positions(
        relIfaceIndex1, relIfaceIndex2, reactMol1, reactMol2, backRxn) };
    if (positions == nullptr)
        return find_reaction_rate_state_scan(
            relIfaceIndex1, relIfaceIndex2, reactMol1, reactMol2, backRxn, molTemplateList);

    int rateIndex { positions->fixedRateIndex };
    if (rateIndex == -1) {
        const RateMatch match { best_matching_rate(
            backRxn, positions->reactIndex1, positions->reactIndex2, reactMol1, reactMol2) };
        rateIndex = match.matches == 0 ? -1 : match.bestRateIndex;
    }

#ifdef NERDSS_CHECK_RATE_STATE_CACHE
    // Both-ways build: the scan must agree on every call.
    const int scanned { find_reaction_rate_state_scan(
        relIfaceIndex1, relIfaceIndex2, reactMol1, reactMol2, backRxn, molTemplateList) };
    if (scanned != rateIndex) {
        std::cerr << "NERDSS_CHECK_RATE_STATE_CACHE: mismatch after " << selectionTally.checked << " selections: scan "
                  << scanned << ", cached " << rateIndex << '\n';
        std::abort();
    }
    ++selectionTally.checked;
#endif
    return rateIndex;
}
