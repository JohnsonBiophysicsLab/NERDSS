#include "reactions/association/association.hpp"
#include "tracing.hpp"
/*What is the angular displacement between the site-COM of each associating molecule
  calculate the dot-product between the two vectors, the tmpICoords and the Coords.
*/
/*What is the angular displacement between the site-complexCOM of each associating molecule*/

void calc_angular_displacement(int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList)
{

    Vec3D v1 { reactMol1.tmpICoords[ifaceIndex1] - reactMol1.tmpComCoord }; //temporary vector due to association moves.
    Vec3D v2 { reactMol1.interfaceList[ifaceIndex1].coord - reactMol1.comCoord }; //original vector
    if (v1.length() < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }

    double currTheta = v2.angle_between(v1);
    // std::cout << "Angle between Molecule1 interface to COM, temporary vs original: " << currTheta << std::endl;

    /*molecule 2*/
    Vec3D v3 { reactMol2.tmpICoords[ifaceIndex2] - reactMol2.tmpComCoord }; //temporary vector due to association moves.
    Vec3D v4 { reactMol2.interfaceList[ifaceIndex2].coord - reactMol2.comCoord }; //original vector
    if (v3.length() < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }

    currTheta = v4.angle_between(v3);
    // std::cout << "Angle between Molecule2 interface to COM, temporary vs original: " << currTheta << std::endl;

    /*calculate angle between interface and the whole complex COM.*/
    Vec3D v5 { reactMol1.tmpICoords[ifaceIndex1] - reactCom1.tmpComCoord }; //temporary vector due to association moves.
    Vec3D v6 { reactMol1.interfaceList[ifaceIndex1].coord - reactCom1.comCoord }; //original vector
    if (v5.length() < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }

    currTheta = v6.angle_between(v5);
    // std::cout << "Angle between Molecule1 interface to ComplexCOM, temporary vs original: " << currTheta << std::endl;

    Vec3D v7 { reactMol2.tmpICoords[ifaceIndex2] - reactCom2.tmpComCoord }; //temporary vector due to association moves.
    Vec3D v8 { reactMol2.interfaceList[ifaceIndex2].coord - reactCom2.comCoord }; //original vector
    if (v7.length() < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }

    currTheta = v8.angle_between(v7);
    // std::cout << "Angle between Molecule2 interface to ComplexCOM, temporary vs original: " << currTheta << std::endl;
}
