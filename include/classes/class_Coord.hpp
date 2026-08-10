/*! \file class_coord.hpp

 * Created on 6/1/18 by Matthew Varga
 * Purpose:
 * Notes:
 */
#pragma once

#include "classes/class_Membrane.hpp"
#include "classes/class_Parameters.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

/*!
 * \brief Rounds to four decimal places.
 *
 * Coordinate comparisons are made on rounded values, so this is what defines
 * "the same point" for Coord::operator== and for the geometry checks in
 * check_bases.cpp and conservedMags.cpp.  Defined here rather than in the
 * translation unit so those comparisons inline.
 */
inline double roundv(double var) noexcept
{
    // if-else is because neg and pos values will round differently
    double val = (int)(var > 0 ? var * 10000 + 0.5 : var * 10000 - 0.5);
    return val / 10000;
}

/*! \struct Coord
 * \brief Class to hold xyz coordinates
 *
 * The arithmetic below is defined inline because it is used inside the
 * propagation, reflection and association loops, which run once per molecule
 * per timestep.  With the definitions in class_Coord.cpp every one of those
 * uses was an out-of-line call that the compiler could not inline across
 * translation units.
 */
struct Coord {
public:
    double x { 0 };
    double y { 0 };
    double z { 0 };

    // operator overloads
    Coord& operator+=(const Coord& coord) noexcept
    {
        x += coord.x;
        y += coord.y;
        z += coord.z;
        return *this;
    }

    friend bool operator==(const Coord& c1, const Coord& c2) noexcept
    {
        return roundv(c1.x) == roundv(c2.x) && roundv(c1.y) == roundv(c2.y) && roundv(c1.z) == roundv(c2.z);
    }

    friend bool operator!=(const Coord& c1, const Coord& c2) noexcept { return !(c1 == c2); }

    friend std::ostream& operator<<(std::ostream& os, const Coord& c);

    friend Coord operator+(const std::array<double, 3>& arr, const Coord& c) noexcept
    {
        return { arr[0] + c.x, arr[1] + c.y, arr[2] + c.z };
    }

    friend Coord operator+(const Coord& c1, const Coord& c2) noexcept
    {
        return { c1.x + c2.x, c1.y + c2.y, c1.z + c2.z };
    }

    friend Coord operator-(const Coord& c1, const double val) noexcept
    {
        return { c1.x - val, c1.y - val, c1.z - val };
    }

    //! Adds an offset in place.  This used to return a new Coord and leave the
    //! operand untouched, which is not what `+=` means.
    friend Coord& operator+=(Coord& c, const std::array<double, 3>& arr) noexcept
    {
        c.x += arr[0];
        c.y += arr[1];
        c.z += arr[2];
        return c;
    }

    friend Coord& operator/=(Coord& c, double scal) noexcept
    {
        c.x = c.x / scal;
        c.y = c.y / scal;
        c.z = c.z / scal;
        return c;
    }

    Coord& operator=(const std::array<double, 3>& arr) noexcept
    {
        this->x = arr[0];
        this->y = arr[1];
        this->z = arr[2];
        return *this;
    }

    Coord& operator-=(const Coord& c) noexcept
    {
        this->x -= c.x;
        this->y -= c.y;
        this->z -= c.z;
        return *this;
    }

    Coord operator-(const Coord& coord2) const noexcept
    {
        return { this->x - coord2.x, this->y - coord2.y, this->z - coord2.z };
    };

    void zero_crds() noexcept
    {
        /// Just sets all values to zero
        x = 0;
        y = 0;
        z = 0;
    }

    /*!
     * \brief Checks if the coordinate is outside the waterbox. Only used in reflecting boundary conditions.
     */
    bool isOutOfBox(const Membrane& membraneObject) const noexcept
    {
        if ((x > (membraneObject.waterBox.x / 2.0)) || (x < -(membraneObject.waterBox.x / 2.0)))
            return true;
        if ((y > (membraneObject.waterBox.y / 2.0)) || (y < -(membraneObject.waterBox.y / 2.0)))
            return true;
        if ((z > (membraneObject.waterBox.z / 2.0)) || (z < -(membraneObject.waterBox.z / 2.0)))
            return true;

        return false;
    }

    //! Squared distance from the origin.  Use this instead of get_magnitude()
    //! whenever the result is only compared against another length, so the
    //! square root is never taken.
    double magnitude_squared() const noexcept { return x * x + y * y + z * z; }

    double get_magnitude() const noexcept { return sqrt(x * x + y * y + z * z); }

    Coord() = default;
    Coord(double x, double y, double z) noexcept
        : x(x)
        , y(y)
        , z(z)
    {
    }

    // TODO: include this in association 2.0
    explicit Coord(const std::array<double, 3>& arr) noexcept
        : x(arr[0])
        , y(arr[1])
        , z(arr[2])
    {
    }

    explicit Coord(const std::vector<double>& vals);

    /*
    Function serialize serializes the Coord
    into array of bytes.
    */
    void serialize(unsigned char* arrayRank, int& nArrayRank) const
    {
        // memcpy rather than a typed store through arrayRank + nArrayRank: that
        // pointer carries no alignment guarantee, so writing a double through it
        // is undefined behavior even where it happens to work.  The bytes
        // written are the same, so serialized buffers stay compatible.
        std::memcpy(arrayRank + nArrayRank, &x, sizeof(x));
        nArrayRank += sizeof(x);
        std::memcpy(arrayRank + nArrayRank, &y, sizeof(y));
        nArrayRank += sizeof(y);
        std::memcpy(arrayRank + nArrayRank, &z, sizeof(z));
        nArrayRank += sizeof(z);
    }
    /*
    Function deserialize deserializes the Coord
    from array of bytes.
    */
    void deserialize(unsigned char* arrayRank, int& nArrayRank)
    {
        std::memcpy(&x, arrayRank + nArrayRank, sizeof(x));
        nArrayRank += sizeof(x);
        std::memcpy(&y, arrayRank + nArrayRank, sizeof(y));
        nArrayRank += sizeof(y);
        std::memcpy(&z, arrayRank + nArrayRank, sizeof(z));
        nArrayRank += sizeof(z);
    }
};

inline Coord round(const Coord& c) noexcept { return { roundv(c.x), roundv(c.y), roundv(c.z) }; }

template <typename Scal>
Coord operator*(Scal scal, const Coord& coord)
{
    return { scal * coord.x, scal * coord.y, scal * coord.z };
}

bool is_co_linear(const Coord& c1, const Coord& c2, const Coord& c3);
