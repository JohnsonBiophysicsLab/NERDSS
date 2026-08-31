#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "tracing.hpp"

void reflect_traj_complex_rad_rot_nocheck_sphere(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    double RS3D { reflecting_surface_offset(targCom, RS3Dinput) };

    double sphereR = membraneObject.sphereR - RS3D;

    std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);

    // if (targCom.D.z < 1E-14 || targCom.OnSurface) { // for the complex on the sphere surface
    if (targCom.OnSurface) {
        // in this case, the movement only involves theta and phi, and R doesn't change,
        // so it won't make the complex outside the sphere.
        return;
    }

    // for the complex inside the sphere
    const Vec3D base { targCom.comCoord + targCom.trajTrans };
    if (base.length() + targCom.radius <= sphereR)
        return;

    // for the outside sphere situation: find the furthest point from the sphere centre
    ExtremePoint farthest { Vec3D { 0, 0, sphereR },
        radial_signed_boundary(sphereR, RadialSide::Inside) };
    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D rot { matrix_rotate(point - targCom.comCoord, M) };
        const Vec3D curr { base + rot };
        farthest.consider(curr, radial_signed_radius(curr, RadialSide::Inside));
    });

    // Applied even when nothing was found outside, where the seed makes it a no-op.
    // Bit-identical to the open-coded form this replaced: the score is
    // `curr.length()` of the same vector radial_reflection_shift() re-measures.
    targCom.trajTrans += radial_reflection_shift(farthest.point, sphereR);
}
