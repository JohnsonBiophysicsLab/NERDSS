/*! \file class_Vec3D.hpp
 *
 * \author Matthew Varga
 * \author Yue Ying
 * 
 * \brief The one three-dimensional vector type.
 *
 * Replaces `Coord` (created 6/1/18 by Matthew Varga) and `Vector` (created
 * 5/20/18 by Matthew Varga), which were the same three doubles twice over:
 * `Vector` derived from `Coord` and added a cached `magnitude`.  Having both
 * meant every expression had a type that depended on which of the two headers
 * the operand came from, and the two halves of the API had drifted apart:
 *
 *   * `Coord` compared with rounding, `Vector` inherited that but printed in a
 *     different format;
 *   * `Coord::get_magnitude()`, `Coord::magnitude_squared()` and
 *     `Vector::calc_magnitude()` / `Vector::magnitude` were four spellings of
 *     two operations;
 *   * `operator+(Vector, Coord)` existed twice, once returning `Vector` and
 *     once returning `Coord`, and which one a call site got was decided by
 *     whether its left operand happened to be const;
 *   * a conversion `Coord` -> `Vector` was needed to reach `dot`, `cross` and
 *     `normalize`, so the association code is full of `Vector { a - b }`
 *     round-trips that exist only to change the type name.
 *
 * \section vec3d_magnitude The removed magnitude cache
 *
 * `Vector::magnitude` was a *stale-able* cache: it was written by
 * `calc_magnitude()` and then not maintained by any operation that changed
 * x, y or z.  Nearly every read was preceded by a `calc_magnitude()` call one
 * or two lines above, so for those the cache was pure bookkeeping and `length()`
 * gives the same bits.  A handful of call sites, however, depended on the
 * cache holding something *other* than the current length - a magnitude of
 * zero on a vector that was never measured, or a pre-rotation length on a
 * vector that has since been rotated - and those decided real branches.  Those
 * sites now pass the length they mean as an explicit argument to
 * \ref Vec3D::angle_between, so the dependency is visible in the call instead
 * of hidden in an object's history.  See requiresSignFlip.cpp,
 * create_arbitrary_vector.cpp, transform.cpp and calculate_phi.cpp.
 *
 * Dropping the member takes the type from 32 back to 24 bytes and removes one
 * of the two `sqrt` calls that `normalize()` used to make (the second only
 * refreshed the cache).
 *
 * \section vec3d_inline Why everything is inline
 *
 * This arithmetic runs inside the propagation, reflection and association
 * loops, once per interface per molecule per timestep.  With the definitions in
 * a translation unit of their own none of it could be inlined into those loops:
 * neither build file enables link-time optimization.  Only the cold paths - the
 * validating constructor, stream output, and the angle function with its
 * warning - are left out of line in class_Vec3D.cpp.
 *
 * Every expression below is character for character what the corresponding
 * `Coord` or `Vector` member computed, in the same order, so each one rounds
 * and contracts exactly as it did before.
 */
#pragma once

#include "attributes.hpp"
#include "classes/class_Membrane.hpp"
#include "classes/class_Parameters.hpp"
#include "numerics/numerical_settings.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

/*!
 * \brief Rounds according to the configured Vec3D coordinate precision.
 *
 * Coordinate comparisons are made on rounded values, so this defines what
 * constitutes "the same point" for Vec3D::operator== and related geometry
 * checks.
 *
 * \param var The value to round.
 *
 * \return The value rounded to the precision specified by
 *         NumericalSettings::Vec3D::coordinateEqualityPrecision.
 *
 * \note A precision of 10000 corresponds to four decimal places, 100000
 *       corresponds to five decimal places, etc.
 *
 * \note The arithmetic is intentionally kept equivalent to the original
 *       implementation. In particular, the multiplication, addition or
 *       subtraction, integer conversion, and division are kept as separate
 *       operations to preserve floating-point behavior.
 */
inline double roundv(double var) noexcept
{
    const int precision =
        NumericalSettings::Vec3D::coordinateEqualityPrecision;

    // if-else is because neg and pos values will round differently
    double val = (int)(var > 0
        ? var * precision + 0.5
        : var * precision - 0.5);

    return val / precision;
}

/*! \struct Vec3D
 * \brief Three doubles: a point, a displacement or a direction.
 *
 * Trivially copyable and 24 bytes, so it is a homogeneous floating-point
 * aggregate and passes and returns in three FP registers.
 */
struct Vec3D {
    double x { 0 };
    double y { 0 };
    double z { 0 };

    /* ----------------------------------------------------------- construction */

    Vec3D() = default;

    Vec3D(double x, double y, double z) noexcept
        : x(x)
        , y(y)
        , z(z)
    {
    }

    explicit Vec3D(const std::array<double, 3>& arr) noexcept
        : x(arr[0])
        , y(arr[1])
        , z(arr[2])
    {
    }

    //! \brief Exits unless `vals` holds exactly three numbers.  Parse paths only.
    explicit Vec3D(const std::vector<double>& vals);

    /* ------------------------------------------------------ compound assignment */

    Vec3D& operator+=(const Vec3D& v) noexcept
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vec3D& operator-=(const Vec3D& v) noexcept
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    Vec3D& operator*=(double scal) noexcept
    {
        x *= scal;
        y *= scal;
        z *= scal;
        return *this;
    }

    Vec3D& operator/=(double scal) noexcept
    {
        x /= scal;
        y /= scal;
        z /= scal;
        return *this;
    }

    /* ------------------------------------------------- std::array interop */

    Vec3D& operator=(const std::array<double, 3>& arr) noexcept
    {
        this->x = arr[0];
        this->y = arr[1];
        this->z = arr[2];
        return *this;
    }

    Vec3D& operator+=(const std::array<double, 3>& arr) noexcept
    {
        x += arr[0];
        y += arr[1];
        z += arr[2];
        return *this;
    }

    /* ---------------------------------------------------------------- products */

    NERDSS_NODISCARD double dot(const Vec3D& v) const noexcept
    {
        return ((this->x * v.x) + (this->y * v.y) + (this->z * v.z));
    }

    /*!
     * \brief The cross product \f$ this \times v \f$, unnormalized.
     *
     * `Vector::cross()` normalized its result, which is not what a cross
     * product is and cost a `sqrt` that several callers threw away by
     * normalizing again.  That behaviour is \ref unit_cross; every former
     * `.cross()` call site was moved to whichever of the two it actually meant.
     */
    NERDSS_NODISCARD Vec3D cross(const Vec3D& v) const noexcept
    {
        return { this->y * v.z - this->z * v.y, this->z * v.x - this->x * v.z, this->x * v.y - this->y * v.x };
    }

    //! \brief The cross product, normalized to unit length.  Was `Vector::cross()`.
    NERDSS_NODISCARD Vec3D unit_cross(const Vec3D& v) const noexcept
    {
        Vec3D result { cross(v) };
        result.normalize();
        return result;
    }

    /* --------------------------------------------------------------- magnitude */

    /*!
     * \brief \f$ |v|^2 \f$.  Use this instead of length() whenever the result is
     * only compared against another length, so the square root is never taken.
     */
    NERDSS_NODISCARD double length_squared() const noexcept { return x * x + y * y + z * z; }

    //! \brief \f$ |v| \f$.  Replaces `Coord::get_magnitude()` and the
    //! `calc_magnitude()` / `magnitude` pair.
    NERDSS_NODISCARD double length() const noexcept { return sqrt(x * x + y * y + z * z); }

    /*!
     * \brief Scales to unit length in place.
     *
     * A zero-length vector is left alone rather than turned into NaNs, which is
     * what dividing by its length would do; this is the guard `Vector::normalize()`
     * expressed as `if (magnitude == 0) magnitude = 1`.
     */
    void normalize() noexcept
    {
        double magnitude { length() };
        if (magnitude == 0) {
            // dividing by zero results in NANs
            magnitude = 1;
        }
        *this /= magnitude;
    }

    //! \brief Out-of-place counterpart of normalize().
    NERDSS_NODISCARD Vec3D normalized() const noexcept
    {
        Vec3D result { *this };
        result.normalize();
        return result;
    }

    /* ------------------------------------------------------ angles, projections */

    /*!
     * \brief The angle between two vectors whose lengths are already known.
     *
     * This is the primitive because the caller, not the vector, is what knows
     * which length is meant: `Vector::dot_theta()` read a cached `magnitude`
     * that several call sites had deliberately left holding something other
     * than the current length.  Passing zero for either length reproduces the
     * "undefined angle" case, which warns and returns 0.
     */
    NERDSS_NODISCARD double angle_between(const Vec3D& vec, double selfNorm, double vecNorm) const;

    //! \brief The angle between two vectors, measuring both lengths now.
    NERDSS_NODISCARD double angle_between(const Vec3D& vec) const
    {
        return angle_between(vec, length(), vec.length());
    }

    /*!
     * \brief The component of this vector perpendicular to `normal`, i.e. the
     * vector rejection \f$ v - \frac{v \cdot n}{n \cdot n} n \f$.
     *
     * Named `vector_projection()` before, which is the name of the other half
     * of the decomposition - the part along `normal`, which this subtracts off.
     *
     * Defined below, after the free operators its body uses.
     */
    NERDSS_NODISCARD Vec3D rejection_from(const Vec3D& normal) const noexcept;

    /* ------------------------------------------------------------------- misc */

    //! \brief Sets all three components to zero.
    void zero() noexcept
    {
        x = 0;
        y = 0;
        z = 0;
    }

    /*!
     * \brief Checks if the coordinate is outside the waterbox. Only used in reflecting boundary conditions.
     */
    NERDSS_NODISCARD bool isOutOfBox(const Membrane& membraneObject) const noexcept
    {
        if ((x > (membraneObject.waterBox.x / 2.0)) || (x < -(membraneObject.waterBox.x / 2.0)))
            return true;
        if ((y > (membraneObject.waterBox.y / 2.0)) || (y < -(membraneObject.waterBox.y / 2.0)))
            return true;
        if ((z > (membraneObject.waterBox.z / 2.0)) || (z < -(membraneObject.waterBox.z / 2.0)))
            return true;

        return false;
    }

    /* -------------------------------------------------------------- serialization */

    /*!
     * \brief Writes the three components to `arrayRank` and advances `nArrayRank`.
     *
     * memcpy rather than a typed store through `arrayRank + nArrayRank`: that
     * pointer carries no alignment guarantee, so writing a double through it
     * is undefined behavior even where it happens to work.  The bytes written
     * are the same.
     *
     * `Vector::serialize()` used to append a fourth double for the magnitude
     * cache; with the cache gone every vector on the wire is three doubles, and
     * both ends of every message are built from this one definition.
     */
    void serialize(unsigned char* arrayRank, int& nArrayRank) const
    {
        std::memcpy(arrayRank + nArrayRank, &x, sizeof(x));
        nArrayRank += sizeof(x);
        std::memcpy(arrayRank + nArrayRank, &y, sizeof(y));
        nArrayRank += sizeof(y);
        std::memcpy(arrayRank + nArrayRank, &z, sizeof(z));
        nArrayRank += sizeof(z);
    }

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

/* ------------------------------------------------------------- free operators */

inline Vec3D operator+(const Vec3D& v1, const Vec3D& v2) noexcept
{
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}

inline Vec3D operator-(const Vec3D& v1, const Vec3D& v2) noexcept
{
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}

inline Vec3D operator-(const Vec3D& v) noexcept { return { -v.x, -v.y, -v.z }; }

//! \brief Subtracts `val` from each component.
inline Vec3D operator-(const Vec3D& v, double val) noexcept { return { v.x - val, v.y - val, v.z - val }; }

inline Vec3D operator*(const Vec3D& v, double scal) noexcept { return { v.x * scal, v.y * scal, v.z * scal }; }

inline Vec3D operator*(double scal, const Vec3D& v) noexcept { return { scal * v.x, scal * v.y, scal * v.z }; }

inline Vec3D operator/(const Vec3D& v, double scal) noexcept { return { v.x / scal, v.y / scal, v.z / scal }; }

inline Vec3D operator+(const std::array<double, 3>& arr, const Vec3D& v) noexcept
{
    return { arr[0] + v.x, arr[1] + v.y, arr[2] + v.z };
}

//! \brief Equality on values rounded to the configured precision; see roundv().
inline bool operator==(const Vec3D& v1, const Vec3D& v2) noexcept
{
    return roundv(v1.x) == roundv(v2.x) && roundv(v1.y) == roundv(v2.y) && roundv(v1.z) == roundv(v2.z);
}

inline bool operator!=(const Vec3D& v1, const Vec3D& v2) noexcept { return !(v1 == v2); }

inline Vec3D round(const Vec3D& v) noexcept { return { roundv(v.x), roundv(v.y), roundv(v.z) }; }

/**
 * @brief Computes the vector rejection of this vector from a given normal.
 * The rejection is the component of this vector that is perpendicular to
 * @p normal. It is computed by subtracting the projection onto @p normal:
 *
 * rejection = *this - normal * ((*this · normal) / (normal · normal))
 *
 * Equivalently, the result is the component of this vector lying in the
 * hyperplane perpendicular to @p normal.
 * 
 * @param normal The vector defining the direction from which this vector
 * is rejected. It must be non-zero.
 * 
 * @return The component of this vector perpendicular to @p normal.
 * 
 * @note The intermediate projection term is intentionally computed as a
 * separate statement. Keeping the multiplication and subtraction
 * separate prevents the compiler from combining them into a fused
 * multiply-add (FMA), which can produce slightly different rounding
 * results on platforms that support FMA.
 * 
 * @note This function is @c noexcept and does not modify either vector.
 * 
 * @warning Passing a zero vector as @p normal results in division by zero
 * and produces an undefined/invalid result.
*/
inline Vec3D Vec3D::rejection_from(const Vec3D& normal) const noexcept
{
    double coefficient { this->dot(normal) / normal.dot(normal) };
    // Kept as two statements, as it was: written as one expression the multiply
    // and the subtract may become a single fused multiply-add on targets that
    // have one, which rounds differently.
    Vec3D sTerm { normal * coefficient };
    return *this - sTerm;
}

/* -------------------------------------------------------------------- output */

//! \brief The fixed-width form every coordinate dump uses.
std::ostream& operator<<(std::ostream& os, const Vec3D& v);

/*!
 * \brief Streams the `[xi + yj + zk]` form the reaction parser echoes normals in.
 *
 * A named function rather than a second `operator<<`: with one vector type
 * there is one `operator<<`, and the two formats used to be selected by which
 * of `Coord` and `Vector` a value was declared as.
 */
std::ostream& write_ijk(std::ostream& os, const Vec3D& v);

/* ------------------------------------------------------------------ geometry */

NERDSS_NODISCARD bool is_co_linear(const Vec3D& c1, const Vec3D& c2, const Vec3D& c3);
