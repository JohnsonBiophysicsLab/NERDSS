#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>
#include <math.h>

#include <algorithm>
#include <complex>

// unbinding probability
// h is the time-step; sigma is the bind_radius, Na is the number of proteins in solution,
// Nlipid is the number of lipids on the membrane surface, A is the area of membrane surface
double dissociate2D(paramsIL& parameters2D)
{
    double h = parameters2D.dt;
    double D = parameters2D.Dtot;
    double sigma = parameters2D.sigma;
    double ka = parameters2D.ka;
    double kb = parameters2D.kb / 1.0e6;
    int Na = parameters2D.Na;
    int Nlipid = parameters2D.Nlipid;
    double A = parameters2D.area;
    if (kb < 1E-15) {
        return 0.0;
    }

    double KD = kb / ka;
    double maxNaNlipid = 0;
    if (Na > Nlipid) {
        maxNaNlipid = Na;
    } else {
        maxNaNlipid = Nlipid;
    }

    double b = 2.0 * sqrt(A / M_PI / maxNaNlipid + sigma * sigma);
    double kon = 1.0 / (1.0 / ka + 1.0 / (8.0 * M_PI * D) * (4.0 * log(b / sigma) / pow(1.0 - pow(sigma / b, 2.0), 2.0) - 2.0 / (1.0 - pow(sigma / b, 2.0)) - 1.0));
    double koff = kon * KD;
    double out = 1.0 - exp(-koff * h);
    return out;
}

// The implicit-lipid 2D binding probability.
//
// A protein on the membrane meets lipids spread at density rho.  Its chance of
// binding one in a step is rho times the 2D association probability p_rxn(r0)
// of an isolated pair (what create_survMatrix tabulates), integrated over the
// annulus of lipid positions R2D < r0 < rmax:
//   p = int_R2D^rmax 2 pi r0 p_rxn(r0) dr0 = 4 ka I(R2D),
//   I(R) = int_0^inf (exp(-D dt u^2) - 1) / u^2 / (alpha^2 + beta^2)
//          * (alpha [rmax J1(rmax u) - R J1(R u)] - beta [rmax Y1(rmax u) - R Y1(R u)]) du,
//   alpha = h u Y1(sigma u) + ka Y0(sigma u),  beta = h u J1(sigma u) + ka J0(sigma u),
// with h = 2 pi sigma D, Rmax = sigma + 3 sqrt(4 D dt) and rmax = 5 Rmax.  The
// bracketed differences are r0 J1(r0 u)/u and r0 Y1(r0 u)/u between the
// limits, the antiderivatives of r0 J0(r0 u) and r0 Y0(r0 u).
//
// R2D, the blocking distance, removes the lipids nearest the protein so that
// binding and unbinding balance: it is the root of 4 kb I(R2D) = pOff, the
// one-step unbinding probability dissociate2D(), searched on [sigma, Rmax].
// At an interior root the probability is therefore 4 ka I(R2D) = ka pOff / kb
// exactly, whatever I is; I only decides which of three cases holds:
//   4 kb I(sigma) <= pOff:  R2D = sigma, p = 4 ka I(sigma),
//   4 kb I(Rmax)  >  pOff:  R2D = Rmax,  p = 4 ka I(Rmax),
//   otherwise:              p = ka pOff / kb.
// These are the limits the bisection that stood here converged to (it kept
// sigma when the condition failed everywhere, and Rmax when it held
// everywhere).  It halved [sigma, Rmax] down to 1e-5, about 20 steps, each a
// qagiu integral of I with a fresh 1e6-interval workspace, and returned
// 4 ka I at its last midpoint: up to ~1e-5 from the root, so ~1e-5 relative
// off ka pOff / kb, and qagiu's error estimate never settled on the
// oscillatory integrand (as integrate_survival() describes for the 2D tables).
//
// I(sigma) and I(Rmax) share alpha, beta and the rmax terms, so both come
// from one pass over one quadrature rule.  Up to X it is composite 20-point
// Gauss-Legendre; beyond X the Bessel functions are in their large-argument
// regime and the integral is taken in closed form (function2D_tail).
//   o [1/rmax 2^-40, 1/rmax] in 40 panels halving toward 0.  There the
//     integrand is smooth in ln u rather than u (beta ~ ka, alpha ~ ln u),
//     and it vanishes as u / ln^2 u, so what lies below is negligible.
//   o [1/rmax, X] in panels two periods long of the fastest factor,
//     J1(rmax u) Y1(sigma u) at frequency rmax + sigma: 4 pi radians, which
//     the 20-point rule integrates to far below 1e-20.
//   o X = 40 max(1/sigma, ka/h), as for the 2D tables, and no less than
//     sqrt(46/(D dt)), so that exp(-D dt u^2) < 1e-20 beyond it.  The tail is
//     a series in 1/(sigma u), 1/(r u) and ka/(h u); taken to third order its
//     remainder at X is of relative size 40^-4.  On the keys
//     run_code_tests/implicit_lipid_2d_accuracy_check.cpp covers, I is within
//     2e-10 of I(sigma) of its reference at both ends; with 20 max(...) in
//     place of 40, within 3e-9.
IL2DEndIntegrals integrate_function2D_ends(const paramsIL& params)
{
    const double sigma { params.sigma };
    const double D { params.Dtot };
    const double ka { params.ka };
    const double dt { params.dt };
    const double Rmax { sigma + 3.0 * sqrt(4.0 * D * dt) };
    const double rmax { 5.0 * Rmax };
    const double h { 2.0 * M_PI * sigma * D };
    // A pair that does not move never meets a lipid beyond contact.  (Dtot
    // rounds to 0 below 5e-6 in discretize_2D_Dtot.)
    if (!(D * dt > 0.0))
        return { 0.0, 0.0 };

    const double uGraded { 1.0 / rmax };
    const double X { std::max(40.0 * std::max(1.0 / sigma, ka / h), sqrt(46.0 / (D * dt))) };
    const int nUniform { std::max(1, static_cast<int>(std::ceil((X - uGraded) * (rmax + sigma) / (4.0 * M_PI)))) };
    const double width { (X - uGraded) / nUniform };
    const int nGraded { 40 };
    const size_t order { 20 };
    gsl_integration_glfixed_table* gl { gsl_integration_glfixed_table_alloc(order) };

    double atSigma {};
    double atRmax {};
    auto add_panel = [&](double lo, double hi) {
        for (size_t i = 0; i < order; ++i) {
            double u {};
            double w {};
            gsl_integration_glfixed_point(lo, hi, i, &u, &w, gl);
            const double J0s { j0(sigma * u) };
            const double J1s { j1(sigma * u) };
            const double Y0s { y0(sigma * u) };
            const double Y1s { y1(sigma * u) };
            const double alpha { h * u * Y1s + ka * Y0s };
            const double beta { h * u * J1s + ka * J0s };
            const double aOut { rmax * j1(rmax * u) };
            const double bOut { rmax * y1(rmax * u) };
            const double c { w * expm1(-D * dt * u * u) / (u * u) / (alpha * alpha + beta * beta) };
            atSigma += c * (alpha * (aOut - sigma * J1s) - beta * (bOut - sigma * Y1s));
            atRmax += c * (alpha * (aOut - Rmax * j1(Rmax * u)) - beta * (bOut - Rmax * y1(Rmax * u)));
        }
    };
    double hi { uGraded };
    for (int j = 0; j < nGraded; ++j) {
        add_panel(hi / 2.0, hi);
        hi /= 2.0;
    }
    for (int m = 0; m < nUniform; ++m)
        add_panel(uGraded + m * width, uGraded + (m + 1) * width);
    gsl_integration_glfixed_table_free(gl);

    const double Xend { uGraded + nUniform * width };
    const double tailOut { function2D_tail(params, rmax, Xend) };
    return { atSigma + function2D_tail(params, sigma, Xend) - tailOut,
        atRmax + function2D_tail(params, Rmax, Xend) - tailOut };
}

namespace {

// E_n(z) = int_1^inf exp(-z t) t^-n dt, by its continued fraction (modified
// Lentz), which converges for |z| >= 1 off the negative real axis.
std::complex<double> expint_En(int n, std::complex<double> z)
{
    const double tiny { 1e-300 };
    std::complex<double> b { z + static_cast<double>(n) };
    std::complex<double> c { 1.0 / tiny };
    std::complex<double> d { 1.0 / b };
    std::complex<double> f { d };
    for (int i = 1; i < 100000; ++i) {
        const double an { -static_cast<double>(i) * (n - 1 + i) };
        b += 2.0;
        d = 1.0 / (an * d + b);
        c = b + an / c;
        const std::complex<double> delta { c * d };
        f *= delta;
        if (std::abs(delta - 1.0) < 1e-16)
            break;
    }
    return f * std::exp(-z);
}

} // namespace

double function2D_tail(const paramsIL& params, double s, double X)
{
    // Beyond X, exp(-D dt u^2) = 0 and the integrand of I is the difference
    // of two terms of the same form, one per limit s of the r0 integral:
    //   -Q(u; R)/u^2 + Q(u; rmax)/u^2,  with
    //   Q(u; s) = -s Im[ H1(s u) / (h u H1(sigma u) + ka H0(sigma u)) ]
    // and H_n = J_n + i Y_n.  This returns int_X^inf Q(u; s)/u^2 du, so the
    // tail of I(R) is function2D_tail(R) - function2D_tail(rmax).  Hankel's
    // expansions H_n(z) ~ sqrt(2/(pi z)) e^{i(z - n pi/2 - pi/4)} sum_k i^k a_k(n) z^-k,
    //   a_1(n) = (4n^2 - 1)/8,  a_2(n) = (4n^2 - 1)(4n^2 - 9)/128,
    //   a_3(n) = (4n^2 - 1)(4n^2 - 9)(4n^2 - 25)/3072,
    // carried to third order, together with ka/(h u), give
    //   Q(u; s)/u^2 = -sqrt(s sigma)/h [sin(w u)/u^3 + E cos(w u)/u^4
    //                 + F sin(w u)/u^5 + G cos(w u)/u^6],  w = s - sigma,
    // where, with b1 = 3/(8 sigma) + ka/h, c2 = 15/(128 sigma^2) + ka/(8 sigma h),
    // b3 = -105/(1024 sigma^3) - 9 ka/(128 sigma^2 h), p1 = 3/(8 s),
    // p2 = 15/(128 s^2), p3 = -105/(1024 s^3):
    //   E = p1 - b1,  F = p2 - c2 - b1^2 + p1 b1,
    //   G = p3 - b3 + 2 b1 c2 + b1^3 - p1 c2 - p1 b1^2 - p2 b1.
    // The integrals of sin(w u)/u^3 ... cos(w u)/u^6 from X are Im and Re of
    // X^(1-n) E_n(-i w X).  For w X > 1 they come from the continued fraction:
    // the integration-by-parts recursion from the sine integral, used below
    // that, multiplies the rounding error of pi/2 - Si(w X) by w at each step
    // and loses everything by the u^-5 term once w X is in the thousands.
    const double sigma { params.sigma };
    const double ka { params.ka };
    const double h { 2.0 * M_PI * sigma * params.Dtot };
    const double w { s - sigma };
    const double b1 { 3.0 / (8.0 * sigma) + ka / h };
    const double c2 { 15.0 / (128.0 * sigma * sigma) + ka / (8.0 * sigma * h) };
    const double b3 { -105.0 / (1024.0 * sigma * sigma * sigma) - 9.0 * ka / (128.0 * sigma * sigma * h) };
    const double p1 { 3.0 / (8.0 * s) };
    const double p2 { 15.0 / (128.0 * s * s) };
    const double p3 { -105.0 / (1024.0 * s * s * s) };
    const double E { p1 - b1 };
    const double F { p2 - c2 - b1 * b1 + p1 * b1 };
    const double G { p3 - b3 + 2.0 * b1 * c2 + b1 * b1 * b1 - p1 * c2 - p1 * b1 * b1 - p2 * b1 };

    double S3 {};
    double C4 {};
    double S5 {};
    double C6 {};
    const double wX { w * X };
    if (wX > 1.0) {
        const std::complex<double> z { 0.0, -wX };
        S3 = std::imag(expint_En(3, z)) / (X * X);
        C4 = std::real(expint_En(4, z)) / (X * X * X);
        S5 = std::imag(expint_En(5, z)) / (X * X * X * X);
        C6 = std::real(expint_En(6, z)) / (X * X * X * X * X);
    } else {
        const double C2 { cos(wX) / X - w * (M_PI / 2.0 - gsl_sf_Si(wX)) };
        S3 = sin(wX) / (2.0 * X * X) + w * C2 / 2.0;
        C4 = cos(wX) / (3.0 * X * X * X) - w * S3 / 3.0;
        S5 = sin(wX) / (4.0 * X * X * X * X) + w * C4 / 4.0;
        C6 = cos(wX) / (5.0 * X * X * X * X * X) - w * S5 / 5.0;
    }
    return -sqrt(s * sigma) / h * (S3 + E * C4 + F * S5 + G * C6);
}

double pimplicitlipid_2D(const paramsIL& parameters2D)
{
    // The binding probability per unit lipid density; see
    // integrate_function2D_ends() for the three cases.
    const double ka { parameters2D.ka };
    if (ka < 1E-15) {
        return 0.0;
    }
    paramsIL params { parameters2D };
    const double pOff { dissociate2D(params) };
    const double kb { parameters2D.kb / 1.0e6 };
    const IL2DEndIntegrals I { integrate_function2D_ends(parameters2D) };
    if (!(4.0 * kb * I.atSigma > pOff))
        return 4.0 * ka * I.atSigma;
    if (4.0 * kb * I.atRmax > pOff)
        return 4.0 * ka * I.atRmax;
    return ka * pOff / kb;
}

///////////////////////////////////////////////////////////////
// 3D
double dissociate3D(double h, double D, double sigma, double ka, double kbsecond)
{
    //double h = parameters3D.dt;
    //double D = parameters3D.Dtot;
    //double sigma = parameters3D.sigma;
    //double ka = parameters3D.ka;
    //double kb = parameters3D.kb;
    double kb = kbsecond / 1.0e6; // change the unit S into us.
    if (kb < 1E-15) {
        return 0.0;
    }
    double KD = 2.0 * kb / ka;
    double kon = 0.5 / (1.0 / ka + 1.0 / (4.0 * M_PI * D * sigma));
    double koff = kon * KD;
    double out = 1.0 - exp(-koff * h);
    return out;
}

// binding probability, but must time the lipid density
double pimplicitlipid_3D(double z, paramsIL& parameters3D)
{
    double h = parameters3D.dt;
    double D = parameters3D.Dtot;
    double sigma = parameters3D.sigma;
    double ka = parameters3D.ka;
    if (ka < 1E-15) {
        return 0.0;
    }

    double out;
    if (z > sigma) {
        double alpha = sqrt(D) / sigma * (1.0 + ka / (4.0 * M_PI * sigma * D));
        double conf = 2.0 * M_PI * sigma * sigma * ka * (4.0 * M_PI * sigma * D) / (ka + 4.0 * M_PI * sigma * D) / (ka + 4.0 * M_PI * sigma * D);
        double a = (z - sigma) / sqrt(4.0 * D * h);
        double b = alpha * sqrt(h);
        if (std::isinf(exp(2.0 * a * b + b * b))) {
            out = conf * (exp(-a * a) / sqrt(M_PI) / (a + b) - (2.0 * a * b + 1.0) * erfc(a) + 2.0 * alpha * sqrt(h / M_PI) * exp(-a * a));
        } else {
            out = conf * (exp(2.0 * a * b + b * b) * erfc(a + b) - (2.0 * a * b + 1.0) * erfc(a) + 2.0 * alpha * sqrt(h / M_PI) * exp(-a * a));
        }
    } else {
        z = sigma;
        double alpha = sqrt(D) / sigma * (1.0 + ka / (4.0 * M_PI * sigma * D));
        double conf = 2.0 * M_PI * sigma * ka * sqrt(D) / alpha / (ka + 4.0 * M_PI * sigma * D);
        double a = 0;
        double b = alpha * sqrt(h);
        if (std::isinf(exp(b * b))) {
            out = conf * (1.0 / sqrt(M_PI) / b - 1.0 + 2.0 * alpha * sqrt(h / M_PI));
        } else {
            out = conf * (exp(b * b) * erfc(b) - 1.0 + 2.0 * alpha * sqrt(h / M_PI));
        }
    }
    return out;
}

// for the droplet (compartment)
// binding probability, no need to time the lipid density
double prob_entering_compartment(double dr, paramsIL& parameters)
{
    double h = parameters.dt;
    double D = parameters.Dtot;
    double sigma = parameters.sigma;
    double ka = parameters.ka;
    double kd = 4.0 * M_PI * sigma * D;
    double R = parameters.compartmentR;
    double r = dr + R;
    double rho = parameters.compartSiteRho;
    if (ka < 1E-15) {
        return 0.0;
    }
    if ( dr < sigma ){
        dr = sigma;
        r = dr + R;
    }
    double alpha = sqrt(D*h) * (ka + kd) / ( sigma*kd );
    double conf = R/r * 2.0 * M_PI * rho * sigma * sigma * ka * kd / (ka + kd) / (ka + kd);
    double x1 = (r+R-sigma)/sqrt(4.0*D*h);
    double x2 = (r-R-sigma)/sqrt(4.0*D*h);
    double func1{0};
    double func2{0};
    if (std::isinf(exp(alpha*alpha+2.0*alpha*x2)) == true) { // so for x1, it is also inf
        func1 = - ( 2.0*alpha/sqrt(M_PI) + 1.0/(alpha+x1)/sqrt(M_PI) )*exp(-x1*x1) + (2.0*alpha*x1+1.0)*erfc(x1);
        func2 = - ( 2.0*alpha/sqrt(M_PI) + 1.0/(alpha+x2)/sqrt(M_PI) )*exp(-x2*x2) + (2.0*alpha*x2+1.0)*erfc(x2);
    }else if (std::isinf(exp(alpha*alpha+2.0*alpha*x1)) == true && std::isinf(exp(alpha*alpha+2.0*alpha*x2)) == false) {
        func1 = - ( 2.0*alpha/sqrt(M_PI) + 1.0/(alpha+x1)/sqrt(M_PI) )*exp(-x1*x1) + (2.0*alpha*x1+1.0)*erfc(x1);
        func2 = - exp(alpha*alpha+2.0*alpha*x2)*erfc(alpha+x2) + (2.0*alpha*x2+1.0)*erfc(x2) - 2.0*alpha/sqrt(M_PI)*exp(-x2*x2);
    }else{
        func1 = - exp(alpha*alpha+2.0*alpha*x1)*erfc(alpha+x1) + (2.0*alpha*x1+1.0)*erfc(x1) - 2.0*alpha/sqrt(M_PI)*exp(-x1*x1);
        func2 = - exp(alpha*alpha+2.0*alpha*x2)*erfc(alpha+x2) + (2.0*alpha*x2+1.0)*erfc(x2) - 2.0*alpha/sqrt(M_PI)*exp(-x2*x2);
    }
    double out = conf * ( func1 - func2 );
    return out;
}

double prob_exiting_compartment(double dr, paramsIL& parameters)
{
    double h = parameters.dt;
    double D = parameters.Dtot;
    double sigma = parameters.sigma;
    double ka = parameters.ka;
    double kd = 4.0 * M_PI * sigma * D;
    double R = parameters.compartmentR;
    double r = - dr + R;
    double rho = parameters.compartSiteRho;

    if (ka < 1E-15) {
        return 0.0;
    }
    if ( dr < sigma ){
        dr = sigma;
        r = - dr + R;
    }
    double alpha = sqrt(D*h) * (ka + kd) / ( sigma*kd );
    double conf = R/r * 2.0 * M_PI * rho * sigma * sigma * ka * kd / (ka + kd) / (ka + kd);
    double x1 = (R-sigma+r)/sqrt(4.0*D*h);
    double x2 = (R-sigma-r)/sqrt(4.0*D*h);
    double func1 = 0.0;
    double func2 = 0.0;
    if (std::isinf(exp(alpha*alpha+2.0*alpha*x2)) == true) { // so for x1, it is also inf
        func1 = - ( 2.0*alpha/sqrt(M_PI) + 1.0/(alpha+x1)/sqrt(M_PI) )*exp(-x1*x1) + (2.0*alpha*x1+1.0)*erfc(x1);
        func2 = - ( 2.0*alpha/sqrt(M_PI) + 1.0/(alpha+x2)/sqrt(M_PI) )*exp(-x2*x2) + (2.0*alpha*x2+1.0)*erfc(x2);
    }else if (std::isinf(exp(alpha*alpha+2.0*alpha*x1)) == true && std::isinf(exp(alpha*alpha+2.0*alpha*x2)) == false) {
        func1 = - ( 2.0*alpha/sqrt(M_PI) + 1.0/(alpha+x1)/sqrt(M_PI) )*exp(-x1*x1) + (2.0*alpha*x1+1.0)*erfc(x1);
        func2 = - exp(alpha*alpha+2.0*alpha*x2)*erfc(alpha+x2) + (2.0*alpha*x2+1.0)*erfc(x2) - 2.0*alpha/sqrt(M_PI)*exp(-x2*x2);
    }else{
        func1 = - exp(alpha*alpha+2.0*alpha*x1)*erfc(alpha+x1) + (2.0*alpha*x1+1.0)*erfc(x1) - 2.0*alpha/sqrt(M_PI)*exp(-x1*x1);
        func2 = - exp(alpha*alpha+2.0*alpha*x2)*erfc(alpha+x2) + (2.0*alpha*x2+1.0)*erfc(x2) - 2.0*alpha/sqrt(M_PI)*exp(-x2*x2);
    }
    double out = conf * ( func1 - func2 );

    return out;
}
