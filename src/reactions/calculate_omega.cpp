#include "reactions/association/association.hpp"
#include "tracing.hpp"

double calculate_omega(Vec3D reactIface1, int reactIface2, const Vec3D& sigma, double sigmaLength,
    const ForwardRxn& currRxn, Molecule reactMol1, Molecule reactMol2, const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    /*Re-aligns the molecules so that Sigma faces purely along the z-axis. 
     */
    transform(reactIface1, reactMol1, reactMol2, sigma, sigmaLength);

    Vec3D v1 {};
    Vec3D v2 {};

    if (areSameAngle(currRxn.assocAngles.theta1, M_PI) || areSameAngle(currRxn.assocAngles.theta2, M_PI)) {
        v1 = determine_normal(currRxn.norm1, molTemplateList[reactMol1.molTypeIndex], reactMol1);
        v2 = determine_normal(currRxn.norm2, molTemplateList[reactMol2.molTypeIndex], reactMol2);
    } else {
        v1 = Vec3D(reactIface1 - reactMol1.tmpComCoord);
        v2 = Vec3D(reactMol2.tmpICoords[reactIface2] - reactMol2.tmpComCoord);
    }

    Vec3D projVec1 { v1.x, v1.y, 0 };
    Vec3D projVec2 { v2.x, v2.y, 0 };
    Vec3D test = projVec1.unit_cross(projVec2);
    double omega = projVec1.angle_between(projVec2);
    double tol = 1E-11;
    /*DO NOT FLIP SIGN IF IT IS PI OR ZERO, DUE TO PRECISION, z COULD BE >0! */
    if (test.z > 0 and std::abs(omega) > tol and (M_PI - std::abs(omega)) > tol) //positive z, flip theta sign
        omega = -omega;
    return omega;

    //(projVec1.unit_cross(projVec2).z < 0) ? projVec1.angle_between(projVec2) : -projVec1.angle_between(projVec2);

    //    projVec1 = Vec3D { v1.x, v1.y, 0 };
    //    projVec2 = Vec3D { v2.x, v2.y, 0 };
    //
    //    projVec1.calc_magnitude();
    //    projVec2.calc_magnitude();
    //
    //    return (angleSignIsCorrect(projVec1, projVec2)) ? projVec1.angle_between(projVec2) : -projVec1.angle_between(projVec2);
}
