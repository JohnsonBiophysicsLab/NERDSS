/*! \file reflect_complex_compartment.cpp
 * ### Created on 02/25/2020 by Yiben Fu
 * ### Purpose: pushes a complex that has drifted into the compartment back out
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

void reflect_complex_compartment(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput)
{
    // TRACE();
    // only works for the complex after association or diffusion, spherical system

    // declare the boundary
    double sphereR;
    if (targCom.D.z < 1E-8) {
        sphereR = membraneObject.compartmentR;
    } else {
        sphereR = membraneObject.compartmentR + RS3Dinput;
    }

    if (!((targCom.comCoord.length() + targCom.radius) < sphereR))
        return;

    // find the point that reaches deepest into the compartment.  Scored as a
    // negated radius so the same "largest score wins" accumulator applies;
    // negation is exact.  The `currR < sphereR` half of the original test is
    // implied by `currR < rtmp`, because rtmp starts at sphereR and only shrinks.
    ExtremePoint deepest { Vec3D { 0, 0, sphereR }, radial_signed_boundary(sphereR, RadialSide::Outside) };
    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        deepest.consider(point, radial_signed_radius(point, RadialSide::Outside));
    });

    // Applied unconditionally: when nothing was found inside, the seed makes the
    // shift zero.  Unlike reflect_complex_rad_rot_sphere this does not loop until
    // the complex is clear - the re-check was commented out long ago.
    const Vec3D dtrans { radial_reflection_shift(deepest.point, sphereR) };
    targCom.comCoord += dtrans;
    for (auto memMol : targCom.memberList) {
        moleculeList[memMol].comCoord += dtrans;
        for (auto& iface : moleculeList[memMol].interfaceList)
            iface.coord += dtrans;
    }
}
