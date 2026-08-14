/*! \file association_rejection_checks.cpp
 * \brief The two checks that can still cancel an association after the move.
 *
 * By the time these run the two complexes have already been rotated into
 * contact in their tmp coordinates, and the span check has passed.  Two things
 * can still reject the move: the new complex overlapping something else in the
 * system, and either complex having travelled further than one timestep of
 * diffusion should allow.  Each rejection is tallied, and the displacement one
 * is tallied three ways depending on the dimensionality of the reaction, which
 * is the part that made this worth writing once instead of four times.
 *
 * All four association routines - box, sphere, and the two implicit-lipid
 * variants - carried the same sixteen lines.
 */
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"

void run_association_rejection_checks(bool& cancelAssoc, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns, bool isOnMembrane,
    bool transitionToSurface, copyCounters& counterArrays)
{
    if (cancelAssoc == false) {
        check_for_structure_overlap_system(cancelAssoc, reactCom1, reactCom2, moleculeList, params, molTemplateList,
            complexList, forwardRxns, backRxns);
        if (cancelAssoc == true)
            counterArrays.nCancelOverlapSystem++;
    }
    if (cancelAssoc == false) {
        measure_complex_displacement(
            cancelAssoc, reactCom1, reactCom2, moleculeList, params, molTemplateList, complexList);
        if (cancelAssoc == true) {
            if (isOnMembrane)
                counterArrays.nCancelDisplace2D++;
            else if (transitionToSurface)
                counterArrays.nCancelDisplace3Dto2D++;
            else
                counterArrays.nCancelDisplace3D++;
        }
    }
}
