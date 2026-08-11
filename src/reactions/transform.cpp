#include "reactions/association/association.hpp"

void transform(Coord& reactIface, Molecule& reactMol1, Molecule& reactMol2, const Vector& axis)
{
    Vector alignAxis {};
//     if (isOnMembrane)
//         alignAxis = Vector { 0, 1, 0 };
//     else
    alignAxis = Vector { 0, 0, 1 };//z-axis

    alignAxis.magnitude = 1.0;

    // return without transformation if axis is already on the z-axis
    double ang { alignAxis.dot_theta(axis) };
    if (ang == M_PI || ang == 0)
        return;

    Vector rotAxis { alignAxis.cross(axis) };
    rotAxis.normalize();
    double theta { axis.dot_theta(alignAxis) };

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
