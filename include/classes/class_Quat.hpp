/*! \file class_Quat.hpp
 * \brief Header file for Quat class.
 */


#pragma once
#include "classes/class_Vector.hpp"
#include <iostream>

struct Quat {
    double w{ 0 };
    double x{ 0 };
    double y{ 0 };
    double z{ 0 };

    Quat operator*(const Quat& q);
    friend std::ostream& operator<<(std::ostream& os, const Quat& q);
    double norm();
    double mag();
    Quat scale(double scal);
    Quat inverse();

    /*!
     * \brief Returns a Quat with magnitude unity.
     */
    Quat unit();

    /*!
     * \brief Takes the conjugate of the Quat, i.e. \$f Q^* \$f.
     */
    Quat conjugate();

    /*!
     * \brief Performs a vector rotation with a quaternion. See the \ref association page.
     */
    void rotate(Vector& vec);

    Quat() = default;
    Quat(double _w, double _x, double _y, double _z)
        : w(_w)
        , x(_x)
        , y(_y)
        , z(_z)
    {
    }
};

/*!
 * \brief Returns a unit quaternion drawn uniformly over all orientations.
 *
 * Normalizing four independent U(-1,1) components, which is what the random
 * orientation code used to do, does not sample rotations uniformly: it samples
 * uniformly inside a 4-cube and then projects onto the unit 3-sphere, so
 * directions towards the corners of the cube receive more probability mass than
 * directions towards its face centers (issue #10).
 *
 * This uses Shoemake's subgroup algorithm, which is exactly uniform on the
 * 3-sphere and therefore on the rotation group.  It draws exactly three uniform
 * variates with no rejection, and the result is a unit quaternion by
 * construction: r1^2 + r2^2 = (1 - u1) + u1 = 1.
 */
Quat rand_unit_quat();
