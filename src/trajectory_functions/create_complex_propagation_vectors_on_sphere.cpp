#include "math/rand_gsl.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"


// only works for complex on sphere surface
Vec3D create_complex_propagation_vectors_on_sphere(const Parameters& params, Complex& targCom)
{
    // TRACE();
    Vec3D trajTrans;

    const double R_fixed = params.sphereR;
    const Vec3D COM_original = targCom.comCoord;
    const double COM_norm = COM_original.length();

    // project COM onto the surface
    Vec3D COM = Vec3D {
        R_fixed * COM_original.x / COM_norm,
        R_fixed * COM_original.y / COM_norm,
        R_fixed * COM_original.z / COM_norm,
    };

    double dx = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    double dy = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    double dl = sqrt(dx * dx + dy * dy); // propagation length
    double dangle = acos(dx / dl);
    if (dy < 0.0) {
        dangle = 2.0 * M_PI - dangle;
    } // propagation direction
    double rotangle = dangle - M_PI / 2.0;
    // Step along the polar angle by the arc length dl, then rotate that step
    // into the drawn direction below.  The named members are what keeps the
    // spherical triple from being read as a cartesian one.
    const SphericalCoord COMsphere = find_spherical_coords(COM);
    double dtheta = dl / COMsphere.r;
    Vec3D COMnewtmp = find_cartesian_coords(
        SphericalCoord { COMsphere.theta - dtheta, COMsphere.phi, COMsphere.r });
    // define the inner-coords-set
    Vec3D i = Vec3D { COM.x, COM.y, COM.z };
    i.normalize();
    Vec3D temp = Vec3D { 0.0, 0.0, COM.length() };
    if (std::abs(std::abs(COM.z) - COM.length()) < 1E-8) { // inner_coord_set()'s onAxisTolerance
        // COM is on the z axis, so the seed above is parallel to i and the cross
        // product below is the zero vector.  Vec3D::normalize() leaves a zero
        // vector alone rather than making NaNs, so j and k would come back as
        // (0, 0, 0), rotate_on_sphere() would build a zero tangential component,
        // and the complex would lose the whole tangential half of its step:
        // stuck on the axis for good, sinking dl^2/2R per step, with no warning.
        // The same guard, seed and tolerance as inner_coord_set()'s no-motion
        // branch, which this block is otherwise a copy of.  Calling that function
        // instead would not be a copy: it builds k from a j this block has
        // normalized once more, so the two frames differ in the last ulp at 44%
        // of surface positions.
        temp = Vec3D { -1.0, 0.0, 0.0 };
    }
    temp.normalize();
    Vec3D j = temp.unit_cross(i);
    j.normalize();
    Vec3D k = i.unit_cross(j);
    k.normalize();
    //std::vector<Vec3D> crdset;
    std::array<double, 9> crdset;
    crdset[0] = i.x;
    crdset[1] = i.y;
    crdset[2] = i.z;
    crdset[3] = j.x;
    crdset[4] = j.y;
    crdset[5] = j.z;
    crdset[6] = k.x;
    crdset[7] = k.y;
    crdset[8] = k.z;
    Vec3D COMnew = rotate_on_sphere(COMnewtmp, COM, crdset, rotangle);

    trajTrans.x = COMnew.x - COM.x;
    trajTrans.y = COMnew.y - COM.y;
    trajTrans.z = COMnew.z - COM.z;
    return trajTrans;
}
