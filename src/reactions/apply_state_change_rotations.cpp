/*! \file apply_state_change_rotations.cpp
 * \brief The five association angles, applied to a state-change pair.
 *
 * A bimolecular state change orients the two molecules exactly the way an
 * association does: theta1, theta2, omega, phi1, phi2, in that order, each
 * skipped when the reaction leaves the angle unspecified.  All four
 * `perform_*_state_change_{box,sphere}` files carried their own copy of the
 * sequence, and the copies had not drifted - only the guard *around* them had,
 * which is why that guard stays at the call sites: the bimolecular version
 * skips the rotations only when both molecules are points, the implicit-lipid
 * version when the facilitator alone is.
 *
 * The `associate_*` routines run the same five rotations with the two molecules
 * in the other roles and with theta guarded by `isnan`; they are not folded in
 * here, because making one function serve both would take more parameters than
 * the duplication costs.
 */
#include "reactions/association/association.hpp"

void apply_state_change_rotations(Vec3D& reactIface1, Vec3D& reactIface2, int stateChangeIface, int facilitatorIface,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    const ForwardRxn::Angles& assocAngles, const ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList)
{
    /* THETA */
    theta_rotation(reactIface1, reactIface2, facilitatorMol, stateChangeMol, assocAngles.theta1, facilitatorCom,
        stateChangeCom, moleculeList);

    theta_rotation(reactIface2, reactIface1, stateChangeMol, facilitatorMol, assocAngles.theta2, stateChangeCom,
        facilitatorCom, moleculeList);

    /* OMEGA */
    // if protein has theta M_PI, uses protein norm instead of com_iface vector.
    // Guarded on currRxn rather than on the local copy, as it always has been.
    if (!std::isnan(currRxn.assocAngles.omega)) {
        omega_rotation(reactIface1, reactIface2, stateChangeIface, facilitatorMol, stateChangeMol, facilitatorCom,
            stateChangeCom, assocAngles.omega, currRxn, moleculeList, molTemplateList);
    } // else P1 or P2 is a rod-type protein, no dihedral for associated complex.

    /* PHI */
    if (!std::isnan(assocAngles.phi1)) {
        phi_rotation(reactIface1, reactIface2, stateChangeIface, facilitatorMol, stateChangeMol, facilitatorCom,
            stateChangeCom, currRxn.norm1, assocAngles.phi1, currRxn, moleculeList, molTemplateList);
    } // else P1 has no valid phi angle.

    if (!std::isnan(assocAngles.phi2)) {
        phi_rotation(reactIface2, reactIface1, facilitatorIface, stateChangeMol, facilitatorMol, stateChangeCom,
            facilitatorCom, currRxn.norm2, assocAngles.phi2, currRxn, moleculeList, molTemplateList);
    } // else P2 has no valid phi angle.
}
