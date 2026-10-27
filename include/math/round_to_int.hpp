/*! \file round_to_int.hpp
 * \brief `int(std::round(x))` without the libm call.
 *
 * `std::round` rounds half away from zero.  Without SSE4.1 (the baseline
 * x86-64 target the cluster builds use) it has no instruction to lower to and
 * stays a call into libm, which made `int(round(dtmp))` in
 * discretize_2D_Dtot() 8% of a 2D-heavy model (secretion kd1) on Rockfish.
 * arm64 rounds in one instruction, so on the Mac it was cheap already.
 *
 * Inside the range below the result is exactly `int(std::round(x))`: the
 * truncation toward zero is exact, and so is `x - t`, because t and x have the
 * same sign and |t| <= |x| < |t| + 1, so either t is 0 or the two are within a
 * factor of two of each other (Sterbenz).  `x - t >= 0.5` is then exactly the
 * condition for rounding away from zero, and unlike `floor(x + 0.5)` nothing is
 * rounded on the way, so 0.49999999999999994 stays 0.  Outside the range
 * (and for NaN) the call is kept, which leaves the conversion as it was.
 *
 * Where `std::round` is one instruction (arm64's `frinta`, or x86 with SSE4.1)
 * it is kept: there the truncate-and-compare form costs a few instructions more (+2.5% of
 * secretion kd1's instructions on the M5).
 */
#pragma once

#include <cmath>

//! The truncate-and-compare form; exact as described above.  Tested against
//! std::round on every platform by run_code_tests/round_to_int_check.cpp.
inline int round_half_away_to_int(double x) noexcept
{
    if (!(x > -2147483647.0 && x < 2147483647.0))
        return int(std::round(x));
    const int t { static_cast<int>(x) };
    const double frac { x - t };
    if (frac >= 0.5)
        return t + 1;
    if (frac <= -0.5)
        return t - 1;
    return t;
}

//! int(std::round(x)), by whichever of the two forms is cheaper here.
inline int round_to_int(double x) noexcept
{
#if defined(__aarch64__) || defined(__arm64__) || defined(__SSE4_1__)
    return int(std::round(x));
#else
    return round_half_away_to_int(x);
#endif
}
