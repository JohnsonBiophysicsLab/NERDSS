#include "reactions/association/association.hpp"

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
