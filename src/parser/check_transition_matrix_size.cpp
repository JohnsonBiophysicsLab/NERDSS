#include <iostream>

#include "parser/parser_functions.hpp"

/*
 * `transitionMatrix` and `lifeTime` are sized once, to `transitionMatrixSize`,
 * and are then indexed by the number of copies of the molecule type held in a
 * single complex. Nothing bounds that cluster size, so a complex larger than
 * `transitionMatrixSize` writes outside both containers. Report the risk while
 * the user can still fix the .mol file, since the failure itself is silent heap
 * corruption followed by a segfault much later in the run.
 */
void check_transition_matrix_size(
    const std::vector<MolTemplate>& molTemplateList,
    const std::vector<CreateDestructRxn>& createDestructRxns)
{
    for (const auto& molTemp : molTemplateList) {
        if (molTemp.countTransition == false) continue;

        // the largest cluster possible is every copy of this type in one complex
        int copies { 0 };
        if (molTemp.molTypeIndex >= 0
            && molTemp.molTypeIndex < static_cast<int>(MolTemplate::numEachMolType.size()))
            copies = MolTemplate::numEachMolType[molTemp.molTypeIndex];

        // creation reactions raise the copy number at runtime, so `copies` is
        // not an upper bound for these types
        bool isCreated { false };
        for (const auto& oneRxn : createDestructRxns) {
            int asProduct { 0 };
            int asReactant { 0 };
            for (const auto& oneMol : oneRxn.productMolList)
                if (oneMol.molTypeIndex == molTemp.molTypeIndex) ++asProduct;
            for (const auto& oneMol : oneRxn.reactantMolList)
                if (oneMol.molTypeIndex == molTemp.molTypeIndex) ++asReactant;
            if (asProduct > asReactant) {
                isCreated = true;
                break;
            }
        }

        std::cout << "Molecule " << molTemp.molName << ": counting transitions with "
                  << "transitionMatrixSize = " << molTemp.transitionMatrixSize
                  << ", system holds " << copies << " copies." << std::endl;

        if (copies > molTemp.transitionMatrixSize) {
            std::cout << "WARNING: molecule " << molTemp.molName
                      << " has transitionMatrixSize = " << molTemp.transitionMatrixSize
                      << " but the system holds " << copies << " copies, so a complex "
                      << "large enough to exceed the matrix is possible. If any single "
                      << "complex ever holds more than " << molTemp.transitionMatrixSize
                      << " copies of " << molTemp.molName << ", NERDSS writes outside "
                      << "transitionMatrix and lifeTime: memory corruption, and a crash "
                      << "later in the run. Set transitionMatrixSize >= " << copies
                      << " in " << molTemp.molName << ".mol to be safe for any cluster "
                      << "size, or set countTransition = false." << std::endl;
        } else if (isCreated == true) {
            std::cout << "WARNING: molecule " << molTemp.molName << " is produced by a "
                      << "creation reaction, so its copy number is not bounded by the "
                      << copies << " copies present now. transitionMatrixSize = "
                      << molTemp.transitionMatrixSize << " must stay above the largest "
                      << "cluster reached during the run, or transitionMatrix and "
                      << "lifeTime will be written outside their bounds." << std::endl;
        }
    }
}
