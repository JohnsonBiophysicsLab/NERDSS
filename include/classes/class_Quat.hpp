/*! \file class_Quat.hpp
 * \brief Header file for Quat class.
 */


#pragma once
#include "attributes.hpp"
#include "classes/class_Vector.hpp"

#include <cmath>
#include <iostream>

/*! \struct Quat
 * \brief A quaternion, used for every rotation in the association and
 * propagation code.
 *
 * Every operation that only reads the quaternion is `const`, `noexcept` and
 * marked \ref NERDSS_NODISCARD, so a call whose result is dropped is a compile
 * warning rather than a silent no-op.  That distinction is not academic here:
 * `unit()` returned the normalized quaternion and left the object untouched, and
 * two call sites wrote `rotQuat.unit();` expecting it to normalize in place
 * (issue #9).  The copy-returning and in-place forms are now named apart -
 * `normalized()` / `normalize()`, `inverse()` / `invert()`, `scaled()` /
 * `operator*=` - so the mistake no longer compiles to a no-op.
 *
 * The definitions are inline because the rotation path runs once per interface
 * per molecule per timestep, and with them in class_Quat.cpp none of it could be
 * inlined into the loops that use it: there is no link-time optimization in
 * either build file.  This follows what Coord does for the same reason.
 *
 * To rotate more than one vector with the same quaternion, use \ref
 * QuatRotation, which builds the inverse once instead of once per vector.
 */
struct Quat {
    double w{ 0 };
    double x{ 0 };
    double y{ 0 };
    double z{ 0 };

    /* ---------------------------------------------------------------- products */

    /*!
     * \brief Hamilton product.  `a * b` is the rotation that applies `b` first,
     * then `a`.
     */
    NERDSS_NODISCARD Quat operator*(const Quat& q) const noexcept
    {
        return { w * q.w - x * q.x - y * q.y - z * q.z, w * q.x + x * q.w + y * q.z - z * q.y,
            w * q.y + y * q.w + z * q.x - x * q.z, w * q.z + z * q.w + x * q.y - y * q.x };
    }

    Quat& operator*=(const Quat& q) noexcept
    {
        *this = *this * q;
        return *this;
    }

    //! \brief Scales all four components, which leaves the rotation unchanged.
    NERDSS_NODISCARD Quat scaled(double scal) const noexcept { return { w * scal, x * scal, y * scal, z * scal }; }

    //! \brief In-place counterpart of scaled().
    Quat& operator*=(double scal) noexcept
    {
        w *= scal;
        x *= scal;
        y *= scal;
        z *= scal;
        return *this;
    }

    /* ------------------------------------------------------------- magnitudes */

    /*!
     * \brief The *squared* magnitude \f$ |Q|^2 \f$, not the magnitude.
     *
     * The name is kept from the original class.  mag() is the magnitude.
     */
    NERDSS_NODISCARD double norm() const noexcept { return (w * w + x * x + y * y + z * z); }

    NERDSS_NODISCARD double mag() const noexcept { return sqrt((*this).norm()); }

    /* --------------------------------------------- conjugate, inverse, unit */

    /*!
     * \brief Takes the conjugate of the Quat, i.e. \$f Q^* \$f.
     */
    NERDSS_NODISCARD Quat conjugate() const noexcept { return { w, -x, -y, -z }; }

    /*!
     * \brief Returns \f$ Q^{-1} = Q^* / |Q|^2 \f$.
     *
     * A zero quaternion has no inverse.  It used to yield infinities here and
     * NaN coordinates one rotation later, so the identity is returned instead,
     * which makes the rotation that follows a no-op.  Issue #7 suggested
     * throwing on this case; returning the identity is preferred because the
     * only callers are inside the propagation and association loops, which have
     * no handler and run once per interface per timestep.  Nothing else in the
     * code throws from that depth, and Vector::normalize() already degrades the
     * same way for a zero-length vector.
     *
     * The guard costs one comparison and cannot change any result that was
     * previously finite: for `norm() != 0` the arithmetic below is what the
     * original `conjugate().scale(1 / norm())` performed, term for term.
     */
    NERDSS_NODISCARD Quat inverse() const noexcept
    {
        const double sqNorm{ norm() };
        if (sqNorm == 0.0)
            return { 1.0, 0.0, 0.0, 0.0 };
        return conjugate().scaled(1 / sqNorm);
    }

    /*!
     * \brief Returns a Quat with magnitude unity.
     *
     * Replaces `unit()`, whose name did not say that the object is left alone.
     * As with inverse(), a zero quaternion yields the identity rather than NaNs.
     */
    NERDSS_NODISCARD Quat normalized() const noexcept
    {
        const double magnitude{ mag() };
        if (magnitude == 0.0)
            return { 1.0, 0.0, 0.0, 0.0 };
        return (*this).scaled(1 / magnitude);
    }

    //! \brief In-place counterpart of normalized().
    void normalize() noexcept { *this = normalized(); }

    //! \brief In-place counterpart of inverse().
    void invert() noexcept { *this = inverse(); }

    /* ---------------------------------------------------------------- rotation */

    /*!
     * \brief Performs a vector rotation with a quaternion. See the \ref association page.
     *
     * This builds the inverse on every call.  When one quaternion rotates
     * several vectors, build a \ref QuatRotation once and use that instead.
     */
    void rotate(Vector& vec) const noexcept;

    friend std::ostream& operator<<(std::ostream& os, const Quat& q);

    Quat() = default;
    Quat(double _w, double _x, double _y, double _z)
        : w(_w)
        , x(_x)
        , y(_y)
        , z(_z)
    {
    }
};

/*! \struct QuatRotation
 * \brief A quaternion together with its precomputed inverse, for rotating many
 * vectors with one rotation.
 *
 * `Quat::rotate()` rebuilds `inverse()` - four multiplies, three adds and a
 * divide - on every call, and it is called once per member molecule and once per
 * interface of each of them inside loops that reuse a single quaternion
 * (Complex::propagate(), rotate(), transform(), orient_crds_to_template()).  The
 * redundant work therefore scales with the number of interfaces in the complex.
 * Hoisting it by hand is what is needed rather than trusting the optimizer:
 * those loops write through references that the compiler cannot prove disjoint
 * from the quaternion, so it has to assume the quaternion may have changed and
 * recompute the inverse every iteration.
 *
 * Results are bit-for-bit identical to repeated `Quat::rotate()` calls: the
 * inverse is built by the same expression from the same values, and the
 * per-vector arithmetic below is `Quat::operator*` term for term.
 */
struct QuatRotation {
    Quat rot{};
    Quat inv{ 1.0, 0.0, 0.0, 0.0 };

    QuatRotation() = default;

    explicit QuatRotation(const Quat& q) noexcept
        : rot(q)
        , inv(q.inverse())
    {
    }

    /*!
     * \brief Rotates `vec` in place, \f$ v \mapsto Q v Q^{-1} \f$.
     */
    void rotate(Vector& vec) const noexcept
    {
        const Quat qv{ 0, vec.x, vec.y, vec.z };
        const Quat qm{ rot * qv };

        // Only the vector part of `qm * inv` is ever read, so the scalar part -
        // zero up to rounding, and discarded by every caller - is not formed.
        // The three expressions are the x, y and z components of
        // Quat::operator* with `qm` on the left and `inv` on the right, in the
        // same order, so each one rounds exactly as it did before.
        vec.x = qm.w * inv.x + qm.x * inv.w + qm.y * inv.z - qm.z * inv.y;
        vec.y = qm.w * inv.y + qm.y * inv.w + qm.z * inv.x - qm.x * inv.z;
        vec.z = qm.w * inv.z + qm.z * inv.w + qm.x * inv.y - qm.y * inv.x;
    }

    /*!
     * \brief Rotates `point` about `origin` and returns the rotated point.
     *
     * Replaces the four-line "subtract the origin, rotate, add the origin back"
     * block that every batch-rotation loop in the association and propagation
     * code writes out by hand.
     */
    NERDSS_NODISCARD Coord rotate_about(const Coord& point, const Coord& origin) const noexcept
    {
        Vector vec{ point - origin };
        rotate(vec);
        return Coord{ vec.x, vec.y, vec.z } + origin;
    }
};

inline void Quat::rotate(Vector& vec) const noexcept { QuatRotation{ *this }.rotate(vec); }

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
NERDSS_NODISCARD Quat rand_unit_quat();
