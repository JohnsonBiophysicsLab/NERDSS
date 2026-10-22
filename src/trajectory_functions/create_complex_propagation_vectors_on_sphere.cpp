#include "math/rand_gsl.hpp"
#include "math/sincos.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "trajectory_functions/trajectory_functions.hpp"


// only works for complex on sphere surface
Vec3D create_complex_propagation_vectors_on_sphere(const Parameters& params, Complex& targCom)
{
    const double R_fixed = params.sphereR;
    const Vec3D COM_original = targCom.comCoord;
    const double COM_norm = COM_original.length();

    // project COM onto the surface
    Vec3D COM = Vec3D {
        R_fixed * COM_original.x / COM_norm,
        R_fixed * COM_original.y / COM_norm,
        R_fixed * COM_original.z / COM_norm,
    };

    // A step of length dl in the tangent plane, dx along j and dy along k.
    double dx = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    double dy = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    double dl = sqrt(dx * dx + dy * dy); // propagation length
    if (dl == 0.0) {
        // D = 0: no direction to step in.
        return Vec3D { 0.0, 0.0, 0.0 };
    }

    // define the inner-coords-set: i outward, j = (z x i) east, k = i x j north
    const double r = COM.length();
    Vec3D i = Vec3D { COM.x, COM.y, COM.z };
    i.normalize();
    Vec3D temp = Vec3D { 0.0, 0.0, r };
    if (std::abs(std::abs(COM.z) - r) < 1E-8) { // inner_coord_set()'s onAxisTolerance
        // COM is on the z axis, so the seed above is parallel to i and the cross
        // product below is the zero vector.  Vec3D::normalize() leaves a zero
        // vector alone rather than making NaNs, so j and k would come back as
        // (0, 0, 0) and the complex would lose the whole tangential half of its
        // step: stuck on the axis for good, sinking dl^2/2R per step, with no
        // warning.  The same guard, seed and tolerance as inner_coord_set()'s
        // no-motion branch.
        temp = Vec3D { -1.0, 0.0, 0.0 };
    }
    temp.normalize();
    Vec3D j = temp.unit_cross(i);
    j.normalize();
    Vec3D k = i.unit_cross(j);
    k.normalize();

    // Walk the arc dl along the great circle through COM in the direction
    // (dx j + dy k) / dl: a rotation by a = dl / r, which moves COM by
    //     r (cos a - 1) i + r sin a (dx j + dy k) / dl.
    //
    // This used to be built from angles: the direction as acos(dx / dl), COM
    // to spherical coordinates (two acos and a sin), a step of dl / r towards
    // the north pole back to cartesian (two sincos), and that step rotated into
    // the drawn direction about i (an acos and a sincos).  Off the z axis, k
    // points due north, so that is the same move and the two agree to
    // rounding; within 1e-8 of the axis the angles were taken from the meridian
    // and the frame from the (-1, 0, 0) seed, and the two directions differed by
    // a rotation, which the isotropic step does not see.
    const double arc = dl / r;
    double sinArc, cosArc;
    sin_cos(arc, sinArc, cosArc);
    const double radial = r * (cosArc - 1.0);
    const double tangential = r * sinArc / dl;
    return Vec3D {
        radial * i.x + tangential * (dx * j.x + dy * k.x),
        radial * i.y + tangential * (dx * j.y + dy * k.y),
        radial * i.z + tangential * (dx * j.z + dy * k.z),
    };
}
