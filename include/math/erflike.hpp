/*
 * ERFLIKE interface for NERDSS.
 *
 * Based on Federico Maria Guercilena's ERFLIKE implementation:
 * https://doi.org/10.5281/zenodo.11261631
 */

#ifndef NERDSS_MATH_ERFLIKE_HPP
#define NERDSS_MATH_ERFLIKE_HPP

#include <complex>

namespace erflike {

// w(z) = exp(-z^2) erfc(-i z), the Faddeeva function.
std::complex<double> w(std::complex<double> z);
double w_im(double x);

// Scaled complementary error function.
std::complex<double> erfcx(std::complex<double> z);
double erfcx(double x);

std::complex<double> erfc(std::complex<double> z);
double erfc(double x);

std::complex<double> erf(std::complex<double> z);
double erf(double x);

std::complex<double> erfi(std::complex<double> z);
double erfi(double x);

std::complex<double> Dawson(std::complex<double> z);
double Dawson(double x);

} // namespace erflike

#endif // NERDSS_MATH_ERFLIKE_HPP
