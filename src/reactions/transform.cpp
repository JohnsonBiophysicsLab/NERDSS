#include "reactions/association/association.hpp"

void transform(Vec3D& reactIface, Molecule& reactMol1, Molecule& reactMol2, const Vec3D& axis, double axisLength)
{
//     if (isOnMembrane)
//         const Vec3D alignAxis { 0, 1, 0 };
//     else
    const Vec3D alignAxis { 0, 0, 1 }; // z-axis, whose length is exactly one

    // return without transformation if axis is already on the z-axis.
    //
    // A caller that passes zero for axisLength gets an angle of 0 back and so
    // skips the whole transform; that is the "unmeasured axis" case, and
    // check_bases.cpp relies on it.  See Vec3D::angle_between.
    double ang { alignAxis.angle_between(axis, 1.0, axisLength) };
    if (ang == M_PI || ang == 0)
        return;

    Vec3D rotAxis { alignAxis.unit_cross(axis) };
    rotAxis.normalize();
    double theta { axis.angle_between(alignAxis, axisLength, 1.0) };

    // I cant remember why the angles are negative
    // TODO: might need to have a signs check
    Quat rotQuat(
        cos(-theta / 2), sin(-theta / 2) * rotAxis.x, sin(-theta / 2) * rotAxis.y, sin(-theta / 2) * rotAxis.z);

    // One inverse for both molecules and all of their interfaces, instead of one
    // per rotated vector.
    const QuatRotation rot { rotQuat };

    { // base1
        // rotate COM
        reactMol1.tmpComCoord = rot.rotate_about(reactMol1.tmpComCoord, reactIface);

        // rotate all the other interfaces
        for (auto& coord : reactMol1.tmpICoords)
            coord = rot.rotate_about(coord, reactIface);
    }
    { // base2
        reactMol2.tmpComCoord = rot.rotate_about(reactMol2.tmpComCoord, reactIface);

        // rotate the ifaces
        for (auto& coord : reactMol2.tmpICoords)
            coord = rot.rotate_about(coord, reactIface);
    }
}
