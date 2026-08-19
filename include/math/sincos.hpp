/*! \file sincos.hpp
 * \brief The sine and cosine of one angle, from one argument reduction.
 *
 * Every libm computes sin and cos of the same angle from the same reduced
 * argument, so asking for both together costs barely more than asking for
 * either alone.  Two separate calls do the reduction twice.
 *
 * Both compilers this project is built with fuse an adjacent `sin(x)`/`cos(x)`
 * pair into the platform's sincos entry point on their own -- Apple clang at
 * -O3 emits `__sincos_stret`, and that fusion is what the disassembly in
 * commit 64417c9 was read against.  This helper exists because that fusion is
 * not something to rely on: it needs the two calls to be recognizably adjacent
 * with the same argument, it does not survive `-fno-builtin`, and it never
 * happened at all on older MSVC.  Spelling the pairing out costs nothing where
 * the compiler would have done it anyway and keeps it where the compiler would
 * not.
 *
 * `__builtin_sincos` rather than POSIX `sincos`: the latter is a GNU extension
 * that Apple's libm spells `__sincos` and that a strict-conformance libc need
 * not declare at all, whereas the builtin is understood by both GCC and clang
 * on every target and lowers to whichever of the two the platform has (or to
 * two plain calls where it has neither).  On this project's reference host the
 * results are bit-for-bit those of separate `sin` and `cos` calls, which is
 * what lets it be introduced without moving any simulation output.
 */
#pragma once

#include <cmath>

/*!
 * \brief Writes \f$ \sin(angle) \f$ and \f$ \cos(angle) \f$ in one evaluation.
 *
 * \param angle     radians
 * \param sinAngle  set to the sine
 * \param cosAngle  set to the cosine
 */
inline void sin_cos(double angle, double& sinAngle, double& cosAngle) noexcept
{
#if defined(__GNUC__) || defined(__clang__)
    __builtin_sincos(angle, &sinAngle, &cosAngle);
#else
    sinAngle = std::sin(angle);
    cosAngle = std::cos(angle);
#endif
}
