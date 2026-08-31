/*! \file comparison.hpp
 * \brief Scale-aware floating-point comparison helpers.
 */
#pragma once

#include <algorithm>
#include <cmath>

struct ComparisonTolerance {
    double absolute;
    double relative;

    ComparisonTolerance(double absoluteValue = 0.0, double relativeValue = 0.0)
        : absolute(absoluteValue)
        , relative(relativeValue)
    {
    }
};

/*! \brief Tests equality using the larger of an absolute and relative bound.
 *
 * Callers must choose a tolerance for the quantity being compared.  There is
 * intentionally no universal default: angles, distances, rates, and diffusion
 * coefficients do not share a meaningful epsilon.
 */
inline bool approximately_equal(double lhs, double rhs, const ComparisonTolerance& tolerance)
{
    const double scale { std::max(std::abs(lhs), std::abs(rhs)) };
    const double bound { std::max(tolerance.absolute, tolerance.relative * scale) };
    return bound == 0.0 ? lhs == rhs : std::abs(lhs - rhs) < bound;
}

inline bool approximately_zero(double value, double absoluteTolerance)
{
    return value == 0.0 || std::abs(value) < absoluteTolerance;
}
