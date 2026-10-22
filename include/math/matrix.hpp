/*! \file matrix.hpp

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
#pragma once

#include "attributes.hpp"
#include "classes/class_Vec3D.hpp"

#include <cstring>

/*!
 * \brief Rotate a vector using a rotation matrix (LEGACY).
 *
 * Defined inline because it is called once per interface inside the
 * reflect_traj_* and sweep_separation_* loops, and with the definition in
 * matrix_functions.cpp every one of those was an out-of-line call with the nine
 * matrix entries reloaded from memory; there is no link-time optimization in
 * either build file.  The arithmetic is unchanged.
 *
 * Both arguments used to be taken by non-const reference even though neither is
 * written, which meant a caller with a `const` matrix could not use it.
 */
NERDSS_NODISCARD inline Vec3D matrix_rotate(const Vec3D& vec, const std::array<double, 9>& M)
{
    return { M[0] * vec.x + M[1] * vec.y + M[2] * vec.z,
             M[3] * vec.x + M[4] * vec.y + M[5] * vec.z,
             M[6] * vec.x + M[7] * vec.y + M[8] * vec.z };
}

/*!
 * \brief Create an Euler (Tait-Bryan angles) rotation matrix from x, y, z values.
 */
NERDSS_NODISCARD std::array<double, 9> create_euler_rotation_matrix(double x, double y, double z);

/*!
 * \brief Create an Euler (Tait–Bryan angles) rotation matrix from a Vec3D containing angles of rotation.
 */
NERDSS_NODISCARD std::array<double, 9> create_euler_rotation_matrix(const Vec3D& angles);

/*!
 * \brief An Euler rotation matrix that is built the first time it is needed.
 *
 * The reflect_traj_* boundary routines need the matrix only for a complex close
 * enough to a wall to be scanned interface by interface; on the in vivo
 * clathrin model that is 3.4% of the calls to reflect_traj_complex_rad_rot_box,
 * and building it eagerly cost six trig calls for each of the other 96.6%.
 *
 * The angles are copied when the object is made, which is where the eager code
 * built the matrix, so a later change to the complex's trajRot cannot leak into
 * it.  The matrix itself comes from the same create_euler_rotation_matrix().
 */
class LazyEulerRotationMatrix {
public:
    explicit LazyEulerRotationMatrix(const Vec3D& angles)
        : angles_(angles)
    {
    }

    const std::array<double, 9>& get()
    {
        if (!built_) {
            matrix_ = create_euler_rotation_matrix(angles_);
            built_ = true;
        }
        return matrix_;
    }

private:
    Vec3D angles_;
    bool built_ { false };
    std::array<double, 9> matrix_ {};
};

/*!
 * \brief The Euler rotation matrix for the angles it was last asked for.
 *
 * The overlap sweeps rotate one interface per pair, and built a fresh matrix for
 * each, although consecutive pairs mostly belong to the same complex and a
 * complex's trajRot changes only when it is resampled.  On the in vivo clathrin
 * model 86% of the matrices the box sweep built were ones it had already built
 * in the same call.
 *
 * The cache is keyed on the bits of the three angles, not on which complex they
 * came from, so a reused matrix is exactly the one create_euler_rotation_matrix()
 * would build again, whatever happened to the complex in between.
 */
class EulerRotationMatrixCache {
public:
    const std::array<double, 9>& get(const Vec3D& angles)
    {
        if (!valid_ || !same_bits(angles.x, angles_.x) || !same_bits(angles.y, angles_.y)
            || !same_bits(angles.z, angles_.z)) {
            matrix_ = create_euler_rotation_matrix(angles);
            angles_ = angles;
            valid_ = true;
        }
        return matrix_;
    }

private:
    static bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

    Vec3D angles_;
    bool valid_ { false };
    std::array<double, 9> matrix_ {};
};
