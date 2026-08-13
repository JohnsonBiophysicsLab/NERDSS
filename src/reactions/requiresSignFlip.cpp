#include "reactions/association/association.hpp"

bool requiresSignFlip(Vec3D axis, Vec3D v1, Vec3D v2)
{
    const Vec3D zAxis { 0, 0, 1 }; // length exactly one
    const Vec3D xAxis { 1, 0, 0 }; // length exactly one

    /* Measured once, on entry.  `axis` is rotated in place further down, and
     * the rotation does not maintain a length, so every angle in this function
     * - including the one taken *after* the rotation - is measured against the
     * length the axis arrived with.  That was previously implicit in a
     * `magnitude` cache the rotation left untouched.
     */
    const double axisLength { axis.length() };

    Vec3D u { zAxis.unit_cross(axis) };
    double theta { zAxis.angle_between(axis, 1.0, axisLength) };
    double useXAxis { false };
    if (std::abs(u.x) < 1E-8 && std::abs(u.y) < 1E-8 && std::abs(u.z) < 1E-8) {
        u  = xAxis.unit_cross(axis);
        theta = xAxis.angle_between(axis, 1.0, axisLength);
        useXAxis = true;
    }

    // rotate
    Quat rot(cos(theta / 2), sin(theta / 2) * u.x, sin(theta / 2) * u.y, sin(theta / 2) * u.z);
    rot.normalize();

    // Each of these rotations used to rebuild the inverse of the same
    // quaternion; it is built once per quaternion instead.
    QuatRotation prep { rot };
    prep.rotate(v1);
    prep.rotate(v2);
    prep.rotate(axis);

    // if we rotated wrong way, reverse and rotate the other way
    // TODO: check this 0.01 business
    if ((zAxis.angle_between(axis, 1.0, axisLength) > 0.01 && !useXAxis)
        || (useXAxis && xAxis.angle_between(axis, 1.0, axisLength) < 0.01)) {
        rot.invert();
        prep = QuatRotation { rot };
        prep.rotate(v1);
        prep.rotate(v2);
        rot = Quat(cos(-theta / 2), sin(-theta / 2) * u.x, sin(-theta / 2) * u.y, sin(-theta / 2) * u.z);
        prep = QuatRotation { rot };
        prep.rotate(v1);
        prep.rotate(v2);
    }

    Vec3D projectedVec1;
    Vec3D projectedVec2;
    if (!useXAxis) {
        projectedVec1 = Vec3D{v1.x, v1.y, 0};
        projectedVec2 = Vec3D{v2.x, v2.y, 0};
        return projectedVec1.unit_cross(projectedVec2).z > 0;
    } else {
        projectedVec1 = Vec3D(0, v1.y, v1.z);
        projectedVec2 = Vec3D(0, v2.y, v2.z);
        return projectedVec1.unit_cross(projectedVec2).x > 0;
    }

    // if the angle sign isn't correct, return true
//    return !angleSignIsCorrect(projectedVec1, projectedVec2);
}
