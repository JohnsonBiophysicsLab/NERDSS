#include "io/io.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

int find_reaction_rate_state(int simItr, int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1,
    const Molecule& reactMol2, const BackRxn& backRxn, const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
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
    const RateMatch match { best_matching_rate(backRxn.rateList, reactIndex1, reactIndex2, reactMol1, reactMol2) };
    if (match.matches == 0)
        return -1;

    return match.bestRateIndex;
}
