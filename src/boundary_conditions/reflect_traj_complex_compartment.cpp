/*! \file reflect_traj_complex_compartment.cpp
 * ### Created on 02/25/2020 by Yiben Fu
 * ### Purpose:  works for complex outside a compartment
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

void reflect_traj_complex_compartment(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput)
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

    double sphereR = membraneObject.compartmentR + RS3D;

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    /*calculate distance to the center of the sphere. */
    /*first just test Complex COM+ radius, if it fits outside the compartment.
      if yes, then test if all interfaces outside sphereR.
      then, if they are inside, what is max displacement beyond the sphereR.
      move it along the radial direction outside by the displacement.
     */
    /*assume the origin of the sphere is at zero. */
    const Vec3D base { targCom.comCoord + targCom.trajTrans };

    if (base.length() + targCom.radius >= sphereR)
        return;

    /*Now evaluate all molecules and interfaces distance from boundaries.*/
    // The compartment excludes rather than contains, so the point of interest is
    // the one nearest the centre, not the farthest.  Scored as a negated radius
    // so that the same "largest score wins" accumulator applies; negation is exact.
    ExtremePoint nearest { Vec3D {}, -sphereR };
    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D rot { matrix_rotate(point - targCom.comCoord, M) };
        const Vec3D curr { base + rot };
        nearest.consider(curr, -curr.length());
    });

    if (nearest.score > -sphereR) {
        const double targR { -nearest.score };
        double lamda = -2.0 * (targR - sphereR) / targR;
        targCom.trajTrans = Vec3D(targCom.trajTrans + lamda * nearest.point);
    }

    // assume the complex is not huge enough, so the reflection from the compartment
    // won't make it outside of the box: unlike the sphere and box twins, this one
    // never re-checks the span.
}
