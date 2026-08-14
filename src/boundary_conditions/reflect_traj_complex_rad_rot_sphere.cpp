/*! \file reflect_traj_complex_rad_rot_sphere.cpp
 * ### Created on 02/25/2020 by Yiben Fu
 * ### Purpose:  works for complex inside a sphere
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
#include "tracing.hpp"

void reflect_traj_complex_rad_rot_sphere(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double radius, double RS3Dinput)
{
    // TRACE();

    // for implicit-lipid model, the boundary surface must consider the reflecting-surface RS3D
    // for explicit-lipid model, membraneObject.RS3D = 0. And, for those proteins that are bound on surface, they are not allowed to reflect along Z-axis
    // for the complex on the surface, its propagation is different from those inside the sphere
    double RS3D;
    if (targCom.OnSurface) {
        RS3D = 0;
    } else {
        RS3D = RS3Dinput;
    }

    double sphereR = radius - RS3D;

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    bool recheck = false;
    /*calculate distance to the center of the sphere. */
    /*first just test Complex COM+ radius, if it fits inside sphereR.
      if yes, then test if all interfaces inside sphereR.
      then, if they are outside, what is max displacement beyond the sphereR.
      move it along the radial direction inside by the displacement.
      Also, check if it fits inside the sphere.
     */

    // if (targCom.D.z < 1E-14 || targCom.OnSurface) { // for the complex on the sphere surface
    if (!targCom.OnSurface) { // for the complex inside the sphere
        /*assume the origin of the sphere is at zero. */
        const Vec3D base { targCom.comCoord + targCom.trajTrans };

        if (base.length() + targCom.radius > sphereR) {
            /*Now evaluate all molecules and interfaces distance from boundaries.*/
            ExtremePoint farthest { Vec3D {}, sphereR };
            for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
                const Vec3D rot { matrix_rotate(point - targCom.comCoord, M) };
                const Vec3D curr { base + rot };
                farthest.consider(curr, curr.length());
            });

            if (farthest.score > sphereR) {
                recheck = true;
                double lamda = -2.0 * (farthest.score - sphereR) / farthest.score;
                targCom.trajTrans = Vec3D(targCom.trajTrans + lamda * farthest.point);
            }
        }
    }

    if (recheck) {
        // Test that new coordinates have not pushed you out of the sphere for a very
        // large complex, if so, resample rotation matrix.
        reflect_traj_check_span_sphere(params, targCom, moleculeList, membraneObject, radius, RS3Dinput);
    }
}
