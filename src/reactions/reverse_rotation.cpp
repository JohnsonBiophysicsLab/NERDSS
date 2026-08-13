#include "reactions/association/association.hpp"

void reverse_rotation(Vec3D& reactIface1, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1, Complex& reactCom2, Quat rotQuatPos, Quat rotQuatNeg, std::vector<Molecule>& moleculeList)
{
    rotQuatPos.invert();
    rotQuatNeg.invert();
    rotate(reactIface1, rotQuatPos, reactCom1, moleculeList);
    rotate(reactIface1, rotQuatNeg, reactCom2, moleculeList);
}
