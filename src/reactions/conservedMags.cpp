#include "reactions/association/association.hpp"

bool conservedMags(const Complex& targCom, const std::vector<Molecule>& moleculeList)
{
    for (auto& memMol : targCom.memberList) {
        for (unsigned i { 0 }; i < moleculeList[memMol].tmpICoords.size(); ++i) {
            Vec3D tmpVec1 { moleculeList[memMol].interfaceList[i].coord - moleculeList[memMol].comCoord };
            Vec3D tmpVec2 { moleculeList[memMol].tmpICoords[i] - moleculeList[memMol].tmpComCoord };
            if (roundv(tmpVec2.length()) != roundv(tmpVec1.length())) {
                // std::cerr << "IFACE-COM vector " << i << " in protein " << memMol << " of complex " << targCom.index
                //           << " did not conserve its magnitude. Exiting..." << std::endl;
                // std::cerr << "Mag tempVec2: " << roundv(tmpVec2.length())
                //           << " Mag tempVec1: " << roundv(tmpVec1.length()) << std::endl;
                // std::cout << "Before association: " << tmpVec1 << '\n';
                // std::cout << "After association: " << tmpVec2 << '\n';
                return false;
            }
        }
    }
    return true;
}
