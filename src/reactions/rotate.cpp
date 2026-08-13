#include "reactions/association/association.hpp"

void rotate(const Vec3D& rotOrigin, const Quat& rotQuat, Complex& targCom,
    std::vector<Molecule>& moleculeList)
{
    // The inverse is built once for the whole complex.  Quat::rotate() rebuilt
    // it on every call, so the cost used to scale with the number of interfaces
    // rotated here; the per-vector arithmetic is unchanged.
    const QuatRotation rot { rotQuat };

    // First rotate all points in complex 1 around the reacting interface of
    // protein p1, and translate them by vector
    for (auto& mol : targCom.memberList) {
        Molecule& oneMol { moleculeList[mol] };
        oneMol.tmpComCoord = rot.rotate_about(oneMol.tmpComCoord, rotOrigin);

        // now rotate each member molecule of the complex
        for (auto& iface : oneMol.tmpICoords)
            iface = rot.rotate_about(iface, rotOrigin);
    }
}
