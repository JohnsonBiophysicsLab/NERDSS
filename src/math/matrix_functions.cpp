/*! \file matrix_functions.cpp

 * ### Created on 11/2/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "math/matrix.hpp"
#include <cmath>

std::array<double, 9> create_euler_rotation_matrix(double x, double y, double z)
{
    double sx{ sin(x) };
    double cx{ cos(x) };
    double sy{ sin(y) };
    double cy{ cos(y) };
    double sz{ sin(z) };
    double cz{ cos(z) };

    std::array<double, 9> M{};
    M[0] = cz * cy;
    M[1] = cz * sx * sy - sz * cx;
    M[2] = cz * sy * cx + sz * sx;
    M[3] = sz * cy;
    M[4] = sz * sx * sy + cz * cx;
    M[5] = sz * sy * cx - cz * sx;
    M[6] = -sy;
    M[7] = cy * sx;
    M[8] = cy * cx;

    return M;
}

// The Vec3D overload used to be a second, character-for-character copy of the
// body above.
std::array<double, 9> create_euler_rotation_matrix(const Vec3D& angles)
{
    return create_euler_rotation_matrix(angles.x, angles.y, angles.z);
}
