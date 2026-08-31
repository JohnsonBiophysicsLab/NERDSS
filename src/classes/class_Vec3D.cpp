/*! \file class_Vec3D.cpp
 *
 * The cold half of \ref Vec3D.  Everything that runs inside the propagation and
 * association loops is inline in class_Vec3D.hpp; what is left here is the
 * validating constructor, the two output formats, the angle function - which
 * calls `acos` and can print a warning, so there is nothing for inlining to
 * save - and the co-linearity test.
 *
 * Replaces class_Coord.cpp and class_Vector.cpp.
 */

#include "classes/class_Vec3D.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <vector>

int NumericalSettings::Vec3D::coordinateEqualityPrecision = 10000;

/* ---------------------------------------------------------------- construction */

Vec3D::Vec3D(const std::vector<double>& vals)
{
    // Validate before reading, rather than initializing from vals[0..2] and
    // checking the size afterwards: a shorter vector used to be indexed out of
    // range before the check could reject it.
    if (vals.size() != 3) {
        std::cout << "Coordinate must have exactly 3 points, got " << vals.size() << ". Exiting." << '\n';
        exit(1);
    }
    x = vals[0];
    y = vals[1];
    z = vals[2];
}

/* ---------------------------------------------------------------------- angles */

double Vec3D::angle_between(const Vec3D& vec, double selfNorm, double vecNorm) const
{
    double dp { this->dot(vec) };
    double cTheta = dp / (selfNorm * vecNorm);

    if (std::abs(cTheta - 1) < 1E-12)
        cTheta = 1;
    if (std::abs(-1 - cTheta) < 1E-12)
        cTheta = -1;

    if (selfNorm < 1E-8 || vecNorm < 1E-8) {
        std::cout << "WARNING: Attempted to find angle between vectors with at least one of magnitude 0.\n";
        return 0.0;
    } else
        return acos(cTheta);
}

/* ---------------------------------------------------------------------- output */

std::ostream& operator<<(std::ostream& os, const Vec3D& v)
{
    os << std::setprecision(6) << std::setw(12) << v.x << " " << std::setw(12) << v.y << " " << std::setw(12) << v.z;
    return os;
}

std::ostream& write_ijk(std::ostream& os, const Vec3D& v)
{
    return os << '[' << v.x << "i + " << v.y << "j + " << v.z << "k]";
}

/* -------------------------------------------------------------------- geometry */

bool is_co_linear(const Vec3D& c1, const Vec3D& c2, const Vec3D& c3)
{
    // The three points are co-linear when the triangle they span has zero area.
    // Heron's formula needs three square roots to get the side lengths and is
    // badly conditioned for exactly the sliver triangles this test cares about,
    // so use the cross product instead: |u x v| is twice the triangle area, and
    // comparing squared quantities keeps the same 1E-8 area threshold without
    // taking any square root.
    const Vec3D u { c2 - c1 };
    const Vec3D v { c3 - c1 };

    const Vec3D crossProduct { (u.y * v.z) - (u.z * v.y), (u.z * v.x) - (u.x * v.z), (u.x * v.y) - (u.y * v.x) };

    // area = |u x v| / 2, so area < 1E-8 is |u x v|^2 < 4E-16.
    return crossProduct.length_squared() < 4E-16;
}
