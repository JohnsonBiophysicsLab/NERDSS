#include "reactions/association/association.hpp"

bool conservedRigid(const Complex& targCom, const std::vector<Molecule>& moleculeList)
{
    for (auto& memMol : targCom.memberList) {
        Vec3D v1_0 { moleculeList[memMol].interfaceList[0].coord - moleculeList[memMol].comCoord };
        Vec3D v1_1 { moleculeList[memMol].tmpICoords[0] - moleculeList[memMol].tmpComCoord };
        for (unsigned i { 1 }; i < moleculeList[memMol].tmpICoords.size(); ++i) {
            Vec3D v2_0 { moleculeList[memMol].interfaceList[i].coord - moleculeList[memMol].comCoord };
            Vec3D v2_1 { moleculeList[memMol].tmpICoords[i] - moleculeList[memMol].tmpComCoord };
            if (roundv(v1_0.angle_between(v2_0)) != roundv(v1_1.angle_between(v2_1))) {
                // std::cerr << "IFACE-COM vector " << i << " in protein " << memMol << " of complex " << targCom.index
                //           << " did not conserve its angle to icoord[0]-COM. Exiting..." << std::endl;
                return false;
            }
        }
    }
    return true;
}
