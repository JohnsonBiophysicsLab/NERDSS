/*! \file class_Coord.cpp

 * Created on 6/1/18 by Matthew Varga
 * Purpose:
 * Notes:
 */

#include "classes/class_Coord.hpp"
#include <classes/class_Parameters.hpp>

#include <array>
#include <cmath>
#include <iomanip>
#include <vector>

// CONSTRUCTORS //
Coord::Coord(const std::vector<double>& vals)
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

bool is_co_linear(const Coord& c1, const Coord& c2, const Coord& c3)
{
    // The three points are co-linear when the triangle they span has zero area.
    // Heron's formula needs three square roots to get the side lengths and is
    // badly conditioned for exactly the sliver triangles this test cares about,
    // so use the cross product instead: |u x v| is twice the triangle area, and
    // comparing squared quantities keeps the same 1E-8 area threshold without
    // taking any square root.
    const Coord u { c2 - c1 };
    const Coord v { c3 - c1 };

    const Coord crossProduct { (u.y * v.z) - (u.z * v.y), (u.z * v.x) - (u.x * v.z), (u.x * v.y) - (u.y * v.x) };

    // area = |u x v| / 2, so area < 1E-8 is |u x v|^2 < 4E-16.
    return crossProduct.magnitude_squared() < 4E-16;
}

// OPERATORS //
std::ostream& operator<<(std::ostream& os, const Coord& c)
{
    os << std::setprecision(6) << std::setw(12) << c.x << " " << std::setw(12) << c.y << " " << std::setw(12) << c.z;
    return os;
}
