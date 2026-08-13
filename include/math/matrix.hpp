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
