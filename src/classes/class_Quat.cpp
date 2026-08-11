/*! \file class_Quat.cpp
 * \brief Functions related to association
 *
 * The quaternion arithmetic itself lives in class_Quat.hpp so that it inlines
 * into the propagation and association loops.  What is left here is the printing
 * and the sampler, neither of which is on a hot path.
 */

#include "classes/class_Quat.hpp"
#include "math/rand_gsl.hpp"
#include <cmath>

std::ostream& operator<<(std::ostream& os, const Quat& q)
{
    os << '[' << q.w << ", " << q.x << "i, " << q.y << "j, " << q.z << "k]";
    return os;
}

Quat rand_unit_quat()
{
    const double u1 { rand_gsl() };
    const double u2 { rand_gsl() };
    const double u3 { rand_gsl() };

    const double r1 { std::sqrt(1.0 - u1) };
    const double r2 { std::sqrt(u1) };
    const double theta1 { 2.0 * M_PI * u2 };
    const double theta2 { 2.0 * M_PI * u3 };

    return { r2 * std::cos(theta2), r1 * std::sin(theta1), r1 * std::cos(theta1), r2 * std::sin(theta2) };
}
