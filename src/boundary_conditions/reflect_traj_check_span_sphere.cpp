/*! \file reflect_traj_check_span_sphere.cpp
 * ### Created on 02/27/2020 by Yiben Fu
 * ### Purpose: to check whether the complex.trajTrans and complex.trajRot will make the complex outside the sphere boundary
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
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

namespace {

//! \brief Farthest point of the complex from the sphere centre, after its sampled move.
ExtremePoint farthest_point(const Complex& targCom, const std::vector<Molecule>& moleculeList,
    const std::array<double, 9>& M, double sphereR)
{
    const Vec3D base { targCom.comCoord + targCom.trajTrans };
    // Seeded at the membrane, so "nothing found" reads as "nothing is outside".
    ExtremePoint farthest { Vec3D { 0, 0, sphereR },
        radial_signed_boundary(sphereR, RadialSide::Inside) };

    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D rot { matrix_rotate(point - targCom.comCoord, M) };
        const Vec3D curr { base + rot };
        farthest.consider(curr, radial_signed_radius(curr, RadialSide::Inside));
    });
    return farthest;
}

} // namespace

void reflect_traj_check_span_sphere(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double radius, double RS3Dinput)
{
    // TRACE();
    bool needsRecheck { true };
    int maxItr { 50 };
    int checkItr { 0 };

    double RS3D { reflecting_surface_offset(targCom, RS3Dinput) };
    double sphereR = radius - RS3D;

    // if (targCom.D.z < 1E-14 || targCom.OnSurface) { // for the complex on the sphere surface
    if (targCom.OnSurface) { // for the complex on the sphere surface
        // in this case, the movement only involves theta and phi, and R doesn't change,
        // so it won't make the complex outside the sphere.
        return;
    }

    // for the complex inside the sphere
    while (checkItr < maxItr && needsRecheck) {
        needsRecheck = false;

        std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);
        ExtremePoint farthest { farthest_point(targCom, moleculeList, M, sphereR) };

        // check whether this complex is out of the box, if so, change trajTrans by considering the reflection
        if (farthest.score > sphereR + 1E-15) {
            targCom.trajTrans += radial_reflection_shift(farthest.point, sphereR);
            // check whether the reflection made the complex inside the sphere
            farthest = farthest_point(targCom, moleculeList, M, sphereR);
        }

        // recheck whether this complex is still out sphere, if so, regenerate trajTrans
        if (farthest.score > sphereR + 1E-15) {
            resample_complex_trajectory(targCom, params);

            reflect_traj_complex_rad_rot_nocheck(params, targCom, moleculeList, membraneObject, RS3Dinput);
            ++checkItr;
            needsRecheck = true; // will need to recheck after resampling traj and trajR
        }
    } // end of while-loop

    // A failure to converge used to walk the complex again to print the positions
    // it ended up with; every one of those prints is commented out, so the walk
    // did nothing.  The caller cannot see the outcome either way.
}
