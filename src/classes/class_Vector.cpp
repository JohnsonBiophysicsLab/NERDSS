/*! \file class_Vector.cpp
* \ingroup Associate
 * Created on 5/20/18 by Matthew Varga
 * Purpose:
 * Notes:
 */

#include "classes/class_Vector.hpp"

#include <cmath>

/* CONSTRUCTORS */
// The trivial constructors are defined in class_Vector.hpp so they inline into
// the association and propagation loops.  Only the two that are cold, or that
// need error handling, are left out of line here.
Vector::Vector(std::array<double, 3>& arr)
    : Coord(arr)
{
}

Vector::Vector(const std::vector<double>& arr)
{
    if (arr.size() != 3) {
        std::cerr << "ERROR: Attempting to create Vector from an array with size > 3.";
        exit(1);
    }

    x = arr[0];
    y = arr[1];
    z = arr[2];
}

/* OPERATORS */
std::ostream& operator<<(std::ostream& os, Vector& vec)
{
    return os << '[' << vec.x << "i + " << vec.y << "j + " << vec.z << "k]";
}

/* MEMBER FUNCTIONS */
double Vector::dot_theta(const Vector& vec) const
{
    double dp { this->dot(vec) };
    double cTheta = dp / (this->magnitude * vec.magnitude);

    if (std::abs(cTheta - 1) < 1E-12)
        cTheta = 1;
    if (std::abs(-1 - cTheta) < 1E-12)
        cTheta = -1;

    if (this->magnitude < 1E-8 || vec.magnitude < 1E-8) {
        std::cout << "WARNING: Attempted to find angle between vectors with at least one of magnitude 0.\n";
        return 0.0;
    } else
        return acos(cTheta);
}

Vector Vector::vector_projection(Vector normal)
{
    double coefficient { this->dot(normal) / normal.dot(normal) };
    Vector sTerm { normal * coefficient };
    sTerm.calc_magnitude();
    return { *this - sTerm };
}
