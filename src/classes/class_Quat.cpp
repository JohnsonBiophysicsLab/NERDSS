/*! \file class_Quat.cpp
 * \brief Functions related to association
 *
 * The quaternion arithmetic itself lives in class_Quat.hpp so that it inlines
 * into the propagation and association loops.  What is left here is the
 * sampler, which is not on a hot path.
 */

#include "classes/class_Quat.hpp"
#include "math/rand_gsl.hpp"
#include <cmath>

Quat rand_unit_quat()
{
#ifdef NERDSS_LEGACY_STREAM
    // master: Quat rotQuat{rand_gsl()*2-1 x4}; rotQuat = rotQuat.unit();
    // unit() = scale(1 / mag()), mag() = sqrt(w*w + x*x + y*y + z*z)
    const double lw { rand_gsl() * 2 - 1 };
    const double lx { rand_gsl() * 2 - 1 };
    const double ly { rand_gsl() * 2 - 1 };
    const double lz { rand_gsl() * 2 - 1 };
    const double lnorm { lw * lw + lx * lx + ly * ly + lz * lz };
    const double lmag { std::sqrt(lnorm) };
    const double lscal { 1 / lmag };
    return { lw * lscal, lx * lscal, ly * lscal, lz * lscal };
#endif
    const double u1 { rand_gsl() };
    const double u2 { rand_gsl() };
    const double u3 { rand_gsl() };

    const double r1 { std::sqrt(1.0 - u1) };
    const double r2 { std::sqrt(u1) };
    const double theta1 { 2.0 * M_PI * u2 };
    const double theta2 { 2.0 * M_PI * u3 };

    return { r2 * std::cos(theta2), r1 * std::sin(theta1), r1 * std::cos(theta1), r2 * std::sin(theta2) };
}
