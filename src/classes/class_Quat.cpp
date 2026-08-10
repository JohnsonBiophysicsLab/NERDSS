/*! \file class_Quat.cpp
 * \brief Functions related to association
 */

#include "classes/class_Quat.hpp"
#include "math/rand_gsl.hpp"
#include <cmath>

Quat Quat::operator*(const Quat& q)
{
    return { w * q.w - x * q.x - y * q.y - z * q.z, w * q.x + x * q.w + y * q.z - z * q.y,
        w * q.y + y * q.w + z * q.x - x * q.z, w * q.z + z * q.w + x * q.y - y * q.x };
}

std::ostream& operator<<(std::ostream& os, const Quat& q)
{
    os << '[' << q.w << ", " << q.x << "i, " << q.y << "j, " << q.z << "k]";
    return os;
}

double Quat::norm() { return (w * w + x * x + y * y + z * z); }

double Quat::mag() { return sqrt((*this).norm()); }

Quat Quat::scale(double scal) { return { w * scal, x * scal, y * scal, z * scal }; }

Quat Quat::conjugate() { return { w, -x, -y, -z }; }

Quat Quat::inverse() { return conjugate().scale(1 / norm()); }

Quat Quat::unit() { return (*this).scale(1 / (*this).mag()); }

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

void Quat::rotate(Vector& vec)
{
    // need to revisit this
    //    vec.calc_magnitude();
    //    if (vec.magnitude == 0) {
    //        std::cout << "Cannot rotate vector, magnitude is 0.";
    //        exit(1);
    //    }

    Quat qv { 0, vec.x, vec.y, vec.z };
    Quat qm = (*this) * qv * this->inverse();

    vec.x = qm.x;
    vec.y = qm.y;
    vec.z = qm.z;
}
