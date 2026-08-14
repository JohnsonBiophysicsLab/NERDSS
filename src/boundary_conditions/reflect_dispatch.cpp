/*! \file reflect_dispatch.cpp
 * \brief Picks the box, sphere or compartment reflector for each boundary operation.
 *
 * These six functions were six files, each holding ten to twenty live lines
 * behind one to three hundred lines of its own former body, commented out.
 * There is nothing left in any of them but the geometry choice, and reading
 * them side by side is the only way to see that the choice is not made the same
 * way in all six: the two that take a compartment flag consult it, the other
 * four look only at `membraneObject.isSphere`.
 *
 * The commented-out bodies are gone.  Every one of them was an earlier version
 * of the `_box` routine that now lives in its own file, and `git log` has them.
 *
 * Created 2020-02 by Yiben Fu; see the individual `_box` and `_sphere` files
 * for the original headers.
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"

#include <iostream>

/* DISSOCIATION REFLECTION */

void reflect_complex_rad_rot(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput, bool isInsideCompartment)
{
    if (isInsideCompartment == false) {
        if (membraneObject.isSphere == true)
            reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, membraneObject.sphereR, RS3Dinput);
        else
            reflect_complex_rad_rot_box(membraneObject, targCom, moleculeList, RS3Dinput);
        if (moleculeList[targCom.memberList[0]].enforceCompartmentBC == true) {
            reflect_complex_compartment(membraneObject, targCom, moleculeList, RS3Dinput);
        }
    } else {
        reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, membraneObject.compartmentR, RS3Dinput);
    }
}

/* TRAJ PROPAGATION FUNCTIONS */

void reflect_traj_complex_rad_rot(
    const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput, bool isInsideCompartment)
{
    if (isInsideCompartment == false) {
        if (membraneObject.isSphere == true)
            reflect_traj_complex_rad_rot_sphere(params, moleculeList, targCom, membraneObject, membraneObject.sphereR, RS3Dinput);
        else
            reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, RS3Dinput);
        if (moleculeList[targCom.memberList[0]].enforceCompartmentBC == true) {
            reflect_traj_complex_compartment(params, moleculeList, targCom, membraneObject, RS3Dinput);
        }
    } else {
        reflect_traj_complex_rad_rot_sphere(params, moleculeList, targCom, membraneObject, membraneObject.compartmentR, RS3Dinput);
    }
}

void reflect_traj_check_span(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    if (membraneObject.isSphere)
        reflect_traj_check_span_sphere(params, targCom, moleculeList, membraneObject, membraneObject.sphereR, RS3Dinput);
    else
        reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, RS3Dinput);
}

void reflect_traj_complex_rad_rot_nocheck(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    if (membraneObject.isSphere == true)
        reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, RS3Dinput);
    else
        reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, RS3Dinput);
}

/*evaluates reflection out of box based on tmpCoords of all proteins. targCOM tmpCOM coords also must be consistent-updated.
 *Does not update molecule traj vectors, just updates the passed in vector traj, to collect for both complexes.
 */
void reflect_traj_tmp_crds(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput, bool isInsideCompartment)
{
    /*This routine updated February 2020 to test if a large complex that spans the box or sphere could extend out in both directions
    if so, it attempts to correct for this by resampling the complex's translational and rotational updates.
    */
    if (isInsideCompartment == false) {
        if (membraneObject.isSphere == true)
            reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, membraneObject.sphereR, RS3Dinput);
        else
            reflect_traj_tmp_crds_box(params, moleculeList, targCom, traj, membraneObject, RS3Dinput);
        if (moleculeList[targCom.memberList[0]].enforceCompartmentBC == true) {
            reflect_traj_tmp_crds_compartment(params, moleculeList, targCom, traj, membraneObject, RS3Dinput);
        }
    } else {
        reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, membraneObject.compartmentR, RS3Dinput);
    }
}

/* ASSOCIATION SPAN CHECK */

void check_if_spans(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Membrane& membraneObject)
{
    // Associating proteins have been moved to contact. Before assigning them to the same
    // complex, test to see if the complex is too big to fit in the box.
    if (membraneObject.isSphere == true)
        check_if_spans_sphere(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject, membraneObject.sphereR);
    else
        check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);
}
