/*! \file association_predicates.cpp
 * \brief The yes/no geometry questions the association code asks.
 *
 * Six one-function files, none longer than a screen, all called from the same
 * association routines and all answering the same kind of question: is this
 * angle the one we think it is, and did the move we just made preserve the
 * shape of the complex?  Read together, `areParallel` and `areSameAngle` show
 * that "the same angle" is a tolerance and "parallel" is exact, and
 * `conservedMags` and `conservedRigid` show what "rigid" is checked against.
 */
#include "reactions/association/association.hpp"

/* ------------------------------------------------------------------ angles */

bool areParallel(const double& angle)
{
    return angle == M_PI || angle == 0;
}

bool areSameAngle(double ang1, double ang2, const ComparisonTolerance& tolerance)
{
    //    return std::abs(ang1) == std::abs(ang2);
    return approximately_equal(ang1, ang2, tolerance);
}

bool angleSignIsCorrect(const Vec3D& vec1, const Vec3D& vec2)
{
    /*if their normal points in -z, keep theta, otherwise flip sign
    ADDED: Do not flip the sign if it is PI or Zero.
    double *test=new double[3];
    crossproduct(sigma, normal, test);
    if(test[2]>0  and abs(theta) > 1E-12 and (M_PI - abs(theta)) > 1E-12)//positive z, flip theta.
	theta=-theta;
    */
    if (std::abs(vec1.z) > 1E-12 && std::abs(vec2.z) > 1E-12) {
        Vec3D newVec1 { vec1.x, 0, vec1.z };
        Vec3D newVec2 { vec2.x, 0, vec2.z };
        return newVec1.unit_cross(newVec2).y < 0;
    } else {
        Vec3D newVec1 { vec1.x, vec1.y, 0 };
        Vec3D newVec2 { vec2.x, vec2.y, 0 };
        return newVec1.unit_cross(newVec2).z < 0;
    }
}

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
        u = xAxis.unit_cross(axis);
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
        projectedVec1 = Vec3D { v1.x, v1.y, 0 };
        projectedVec2 = Vec3D { v2.x, v2.y, 0 };
        return projectedVec1.unit_cross(projectedVec2).z > 0;
    } else {
        projectedVec1 = Vec3D(0, v1.y, v1.z);
        projectedVec2 = Vec3D(0, v2.y, v2.z);
        return projectedVec1.unit_cross(projectedVec2).x > 0;
    }

    // if the angle sign isn't correct, return true
    //    return !angleSignIsCorrect(projectedVec1, projectedVec2);
}

/* ------------------------------------------- did the move keep the shape? */

bool conservedMags(const Complex& targCom, const std::vector<Molecule>& moleculeList)
{
    for (auto& memMol : targCom.memberList) {
        for (unsigned i { 0 }; i < moleculeList[memMol].tmpICoords.size(); ++i) {
            Vec3D tmpVec1 { moleculeList[memMol].interfaceList[i].coord - moleculeList[memMol].comCoord };
            Vec3D tmpVec2 { moleculeList[memMol].tmpICoords[i] - moleculeList[memMol].tmpComCoord };
            if (roundv(tmpVec2.length()) != roundv(tmpVec1.length())) {
                // std::cerr << "IFACE-COM vector " << i << " in protein " << memMol << " of complex " << targCom.index
                //           << " did not conserve its magnitude. Exiting..." << std::endl;
                return false;
            }
        }
    }
    return true;
}

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
