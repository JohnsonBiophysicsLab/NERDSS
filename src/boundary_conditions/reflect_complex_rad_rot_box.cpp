/*! \file reflect_complex_rad_rot_box.cpp
 * ### Created on 11/8/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "tracing.hpp"

void reflect_complex_rad_rot_box(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput)
{
    // TRACE();
    // only works for the complex after association or diffusion, box system.
    double RS3D;
    // if (targCom.D.z < 1E-8) { // on surface
    if (targCom.OnSurface) {
        RS3D = 0.0;
    } else { // in solution
        RS3D = RS3Dinput;
    }

    // if moleculeList[memMol].enforceCompartmentBC is true, we need to prevent the molecule from crossing the compartment.
    // if its inside the compartment, it only needs to check the compartment boundaries, not the box boundaries, so we should skip the box boundary checks.
    // if its outside the compartment, we need to check both boundaries: box and compartment.

    // declare the six boundary sides of the system box;
    const Vec3D posSide { membraneObject.waterBox.x / 2.0, membraneObject.waterBox.y / 2.0,
        membraneObject.waterBox.z / 2.0 };
    const Vec3D negSide { -membraneObject.waterBox.x / 2.0, -membraneObject.waterBox.y / 2.0,
        -membraneObject.waterBox.z / 2.0 + RS3D };

    static const char* axisName[3] = { "X", "Y", "Z" };
    static const char* axisLabel[3] = { "x", "y", "z" };

    for (int axis { 0 }; axis < 3; ++axis) {
        const double posWallSide { axis_value(posSide, axis) };
        const double negWallSide { axis_value(negSide, axis) };
        const double base { axis_value(targCom.comCoord, axis) };

        if (!((base + targCom.radius) > posWallSide || (base - targCom.radius) < negWallSide))
            continue;

        // find the farthest point in +axis (posWall) and -axis (negWall)
        const WallExtent extent { scan_axis_extent_placed(
            targCom, moleculeList, axis, WallExtent { negWallSide, posWallSide }) };

        // check whether this complex is out of the box
        const double posd { extent.overshoot(posWallSide) };
        const double negd { extent.undershoot(negWallSide) };
        const bool outsidePos { posd != 0.0 };
        const bool outsideNeg { negd != 0.0 };

        if (outsideNeg && outsidePos) {
            // extends out both the front and back.
            std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in BOTH directions of " << axisName[axis]
                      << " . ALREADY UPDATED POSITIONS. EXITING..." << '\n';
            exit(1);
        }
        if (outsideNeg && !outsidePos) {
            if (extent.posWall - 2.0 * negd > posWallSide) {
                std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in NEGATIVE side of " << axisName[axis]
                          << ": try to put back in the box, " << axisLabel[axis] << ": " << -negd
                          << "BUT will EXTEND again in POSITIVE side of " << axisName[axis] << '\n';
                exit(1);
            }
            // just update positions. Put back inside the box
            translate_coords_along_axis(targCom, moleculeList, axis, 2.0 * negd);
        }
        if (!outsideNeg && outsidePos) {
            if (extent.negWall - 2.0 * posd < negWallSide) {
                std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in POSITIVE side of " << axisName[axis]
                          << ": try to put back in the box, " << axisLabel[axis] << ": " << -posd
                          << "BUT will EXTEND again in NEGATIVE side of " << axisName[axis] << '\n';
                exit(1);
            }
            // just update positions. Put back inside the box
            translate_coords_along_axis(targCom, moleculeList, axis, 2.0 * posd);
        }
    }
}
