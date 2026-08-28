/*! \file reflect_complex_rad_rot_sphere.cpp
 * ### Created on 02/25/2020 by Yiben Fu
 * ### Purpose: only works for complex already bound on the spherical surface
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

namespace {

/*! \brief Point of the complex farthest from the sphere centre, where it stands.
 *
 * Seeded at the membrane, so a complex entirely inside comes back with
 * `score == sphereR` and nothing to correct.  The original also tested
 * `currR > sphereR` alongside `currR > rtmp`; that is implied, because `rtmp`
 * starts at `sphereR` and only grows.
 */
ExtremePoint farthest_point(const Complex& targCom, const std::vector<Molecule>& moleculeList, double sphereR)
{
    ExtremePoint farthest { Vec3D { 0, 0, sphereR }, radial_signed_boundary(sphereR, RadialSide::Inside) };
    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        farthest.consider(point, radial_signed_radius(point, RadialSide::Inside));
    });
    return farthest;
}

} // namespace

void reflect_complex_rad_rot_sphere(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double radius, double RS3Dinput)
{
    // TRACE();
    // only works for the complex after association or diffusion, spherical system
    // For the sphere system, many times of reflections may need to move the complex back inside the sphere!!
    //
    // NOTE: reflect_complex_compartment() is this routine's excluding twin and
    // has no such loop - it reflects once and accepts the result.  The two are
    // not merged for that reason; see docs/debranching_plan.md.

    // declare the boundary
    double sphereR;
    if (targCom.OnSurface) {
        sphereR = membraneObject.sphereR;
    } else {
        sphereR = radius - RS3Dinput;
    }

    if ((targCom.comCoord.length() + targCom.radius) > sphereR) {
        ExtremePoint farthest { farthest_point(targCom, moleculeList, sphereR) };

        int times = 0; // to count the loop-times of 'while'
        while (farthest.score > radial_signed_boundary(sphereR, RadialSide::Inside)) {
            times++;
            const Vec3D dtrans { radial_reflection_shift(farthest.point, sphereR) };
            targCom.comCoord += dtrans;
            for (auto memMol : targCom.memberList) {
                moleculeList[memMol].comCoord += dtrans;
                for (auto& iface : moleculeList[memMol].interfaceList)
                    iface.coord += dtrans;
            }
            // reflecting may make the complex outside the sphere in other direction,
            // thus we need to recheck whether outside
            farthest = farthest_point(targCom, moleculeList, sphereR);

            if (times > 100) {
                // so many times reflection still cannot make the complex back inside the sphere, thus we may need report 'WRONG!!'
                std::cout << "ALREADY UPDATED POSITIONS.BUT, IN REFLECT COMPLEX_RAD_ROT_SPHERE, 100 TIMES REFLECTIONS CAN'T MOVE THE COMPLEX BACK INSIDE SPHERE." << '\n';
                std::cout << "COMPLEX " << targCom.index << ", COMPLEX COM " << targCom.comCoord.x << ", " << targCom.comCoord.y << ", " << targCom.comCoord.z << '\n';
                std::cout << "COMPLEX RADIUS " << targCom.radius << ", COMPLEX SIZE " << targCom.memberList.size() << '\n';
                std::cout << "EXITING..." << '\n';
                exit(1);
            }
        } // end of while-loop
    } // update reflection

    // ALSO, the position update may cause lipids off sphere, so adjust lipids back onto surface,
    int molnumber = -1;
    double dr = 0.0;
    Vec3D dtrans;
    for (auto& memMol : targCom.memberList) {
        if (moleculeList[memMol].isLipid == true) {
            Vec3D targ = moleculeList[memMol].comCoord;
            double drtmp = std::abs(targ.length() - membraneObject.sphereR);
            if (drtmp > 1E-4 && drtmp > dr) {
                dr = drtmp;
                molnumber = memMol;
                dtrans = (membraneObject.sphereR / targ.length()) * targ - targ;
            }
        }
    }
    if (molnumber != -1) {
        for (auto& mol : targCom.memberList) {
            moleculeList[mol].comCoord += dtrans;
            for (auto& iface : moleculeList[mol].interfaceList) {
                iface.coord += dtrans;
            }
        }
    }
}
