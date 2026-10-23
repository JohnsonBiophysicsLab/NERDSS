// Accuracy and cost check for the implicit-lipid 2D binding probability,
// pimplicitlipid_2D().
//
// The probability is 4 ka I(R2D), where I(R) integrates the binding
// probability of a lipid at r0 over R < r0 < rmax, and the blocking distance
// R2D solves 4 kb I(R2D) = pOff on [sigma, Rmax] (functions_implicitlipid.cpp
// has the details).  At an interior root that makes it ka pOff / kb exactly;
// otherwise it is 4 ka I at sigma or at Rmax.  So what has to be right is
// I(sigma) and I(Rmax), and this checks the simulator's
// integrate_function2D_ends() and pimplicitlipid_2D() against references
// computed here, independently of the code under test:
//
//   refA: the integrand as written in the old function2D(), integrated
//         panel by panel (GSL qag, GK61, one panel per half period of the
//         fastest Bessel factor) on [0, X2], plus the asymptotic tail beyond
//         X2, for X2 ten times the production X (and, from the command line,
//         also thirty times, which shows the tail is consistent).
//   refB: 2D Brownian motion is recurrent, so the binding probability of
//         every lipid tends to 1 at long times, and
//           I(R) = pi (rmax^2 - R^2) / (4 ka) + int_0^inf exp(-D dt u^2) g(u) du,
//         with g the integrand without its (exp - 1) factor.  That integral is
//         Gaussian-damped and needs no tail.  Near u = 0, g goes as
//         1/(u ln^2 u): below u_s it is integrated in closed form (an
//         arctangent), and from u_s to 1/rmax in ln u, with Y1's 1/u pole
//         removed analytically so the rmax and R terms do not cancel.
//     The two share nothing but the integrand.  refB subtracts its damped
//     integral from pi (rmax^2 - R^2) / (4 ka), which is huge when ka is
//     small (6e7 times I on the clathrin in-vitro keys), so it is compared
//     only to the precision that cancellation leaves.
//
// Without arguments it runs the keys below and fails if p, I(sigma) or I(Rmax)
// is off by more than 1e-8 of its value at sigma (4 ka I(sigma) for p, which
// bounds it), or the two references disagree beyond what refB's cancellation
// explains.
// With arguments it reports on one key, including the reference R2D:
//   implicit_lipid_2d_accuracy_check sigma Dtot ka kb dt area Na Nlipid
//   implicit_lipid_2d_accuracy_check --bench sigma Dtot ka kb dt area Na Nlipid
// kb is in 1/s, as in the input files; ka is the 2D rate the simulator passes
// (the 3D rate divided by length3Dto2D).  --bench prints the cost of one
// pimplicitlipid_2D() call, with no reference.
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_expint.h>
#ifdef __APPLE__
#include <libproc.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;

namespace {

double instructions()
{
#ifdef __APPLE__
    rusage_info_v4 ri {};
    proc_pid_rusage(getpid(), RUSAGE_INFO_V4, reinterpret_cast<rusage_info_t*>(&ri));
    return static_cast<double>(ri.ri_instructions);
#else
    return 0.0;
#endif
}

struct Key {
    double a, D, k, t; // sigma, Dtot, ka, dt
    double Rmax() const { return a + 3.0 * sqrt(4.0 * D * t); }
    double rmax() const { return 5.0 * Rmax(); }
    double h() const { return 2.0 * M_PI * a * D; }
    double X() const { return std::max(40.0 * std::max(1.0 / a, k / h()), sqrt(46.0 / (D * t))); }
};

struct KR {
    Key p;
    double R;
};

// s Y1(s u) + 2/(pi u): Y1's pole removed, by its power series for s u < 2,
// so that the rmax and R terms can be subtracted at small u without
// cancellation (DLMF 10.8.1).
double sY1_regular(double s, double u)
{
    double z { s * u };
    if (z > 2.0)
        return s * y1(z) + 2.0 / (M_PI * u);
    double half { z / 2 };
    double term { half };
    double psi1 { -M_EULER };
    double psi2 { 1.0 - M_EULER };
    double sum { (psi1 + psi2) * term };
    for (int k = 1; k < 60; ++k) {
        term *= -half * half / (k * (k + 1.0));
        psi1 += 1.0 / k;
        psi2 += 1.0 / (k + 1.0);
        sum += (psi1 + psi2) * term;
        if (std::abs(term) < 1e-18 * std::abs(sum))
            break;
    }
    return s * ((2.0 / M_PI) * log(half) * j1(z) - sum / M_PI);
}

// The old function2D(), copied so the reference does not move with the code
// under test, without its (exp(-D dt u^2) - 1) factor.
double g(double u, const Key& p, double R, bool regular)
{
    double H { p.h() }, rm { p.rmax() };
    double alpha { H * u * y1(p.a * u) + p.k * y0(p.a * u) };
    double peta { H * u * j1(p.a * u) + p.k * j0(p.a * u) };
    double a { u * rm * j1(rm * u) - u * R * j1(R * u) };
    double b { regular ? u * (sY1_regular(rm, u) - sY1_regular(R, u)) : u * rm * y1(rm * u) - u * R * y1(R * u) };
    return (alpha * a - peta * b) / (alpha * alpha + peta * peta) / (u * u * u);
}

double f_full(double u, void* v)
{
    const KR& q { *static_cast<KR*>(v) };
    return expm1(-q.p.D * q.p.t * u * u) * g(u, q.p, q.R, false);
}

double f_damped(double u, void* v)
{
    const KR& q { *static_cast<KR*>(v) };
    return exp(-q.p.D * q.p.t * u * u) * g(u, q.p, q.R, true);
}

double f_damped_log(double s, void* v) // in s = ln u
{
    double u { exp(s) };
    return f_damped(u, v) * u;
}

gsl_integration_workspace* W;
bool panelFailed { false };

double panels(double (*f)(double, void*), KR& q, double lo, double hi, double width)
{
    gsl_function F { f, &q };
    double sum {};
    int n { std::max(1, static_cast<int>(std::ceil((hi - lo) / width))) };
    double dx { (hi - lo) / n };
    for (int i = 0; i < n; ++i) {
        double res {}, e {};
        int st { gsl_integration_qag(&F, lo + i * dx, lo + (i + 1) * dx, 1e-18, 1e-13, 100000, 6, W, &res, &e) };
        // qag reports a roundoff stop when it cannot reach 1e-13; accept the
        // panel if its error estimate is still under 1e-11 of its value.
        if (st != GSL_SUCCESS && std::abs(e) > std::max(1e-14, 1e-11 * std::abs(res))) {
            panelFailed = true;
            std::printf("    reference panel [%.6g, %.6g]: status %d, error %.3g on %.3g\n", lo + i * dx,
                lo + (i + 1) * dx, st, e, res);
        }
        sum += res;
    }
    return sum;
}

// E_n(z) = int_1^inf exp(-z t) t^-n dt by its continued fraction, |z| >= 1.
std::complex<double> En(int n, std::complex<double> z)
{
    std::complex<double> b { z + double(n) }, c { 1e300 }, d { 1.0 / b }, f { d };
    for (int i = 1; i < 100000; ++i) {
        double an { -double(i) * (n - 1 + i) };
        b += 2.0;
        d = 1.0 / (an * d + b);
        c = b + an / c;
        std::complex<double> del { c * d };
        f *= del;
        if (std::abs(del - 1.0) < 1e-16)
            break;
    }
    return f * std::exp(-z);
}

// int_X^inf of the R or rmax term (s) of the integrand, from the third-order
// Hankel expansion, as in function2D_tail(); only used at X2, where it is
// many orders below the tolerance.
double tail_s(const Key& p, double s, double X)
{
    double a { p.a }, k { p.k }, h { p.h() }, w { s - a };
    double b1 { 3 / (8 * a) + k / h }, c2 { 15 / (128 * a * a) + k / (8 * a * h) };
    double b3 { -105 / (1024 * a * a * a) - 9 * k / (128 * a * a * h) };
    double p1 { 3 / (8 * s) }, p2 { 15 / (128 * s * s) }, p3 { -105 / (1024 * s * s * s) };
    double E { p1 - b1 }, F { p2 - c2 - b1 * b1 + p1 * b1 };
    double G { p3 - b3 + 2 * b1 * c2 + b1 * b1 * b1 - p1 * c2 - p1 * b1 * b1 - p2 * b1 };
    double S3, C4, S5, C6;
    if (w * X > 1) {
        std::complex<double> z { 0, -w * X };
        S3 = std::imag(En(3, z)) / pow(X, 2);
        C4 = std::real(En(4, z)) / pow(X, 3);
        S5 = std::imag(En(5, z)) / pow(X, 4);
        C6 = std::real(En(6, z)) / pow(X, 5);
    } else {
        double C2 { cos(w * X) / X - w * (M_PI / 2 - gsl_sf_Si(w * X)) };
        S3 = sin(w * X) / (2 * X * X) + w * C2 / 2;
        C4 = cos(w * X) / (3 * pow(X, 3)) - w * S3 / 3;
        S5 = sin(w * X) / (4 * pow(X, 4)) + w * C4 / 4;
        C6 = cos(w * X) / (5 * pow(X, 5)) - w * S5 / 5;
    }
    return -sqrt(s * a) / h * (S3 + E * C4 + F * S5 + G * C6);
}

double refA(const Key& p, double R, double X2)
{
    KR q { p, R };
    double rm { p.rmax() };
    double u1 { 1.0 / rm };
    return panels(f_full, q, 0.0, u1, u1 / 8) + panels(f_full, q, u1, X2, M_PI / (rm + p.a)) + tail_s(p, R, X2)
        - tail_s(p, rm, X2);
}

// refB and the size of the term it cancels against
double refB(const Key& p, double R, double& cancelled)
{
    KR q { p, R };
    double rm { p.rmax() }, h { p.h() }, k { p.k };
    double us { 1e-6 / rm };
    double u1 { 1.0 / (rm + sqrt(p.D * p.t)) };
    double ug { sqrt(90.0 / (p.D * p.t)) };
    // Below us: H1(z) = -2i/(pi z) + (z/2)(1 + (2i/pi)(ln(z/2) + gamma - 1/2)),
    // H0(z) = 1 + (2i/pi)(ln(z/2) + gamma), so with beta0 the imaginary part
    // of h u H1(a u) + k H0(a u), g = -(M/2) / (u (k^2 + beta0^2)) and its
    // integral is an arctangent in beta0.
    double bb { 2.0 * k / M_PI };
    auto m = [&](double s) { return bb * (log(s / p.a) - 0.5) + 2.0 * h / (M_PI * p.a); };
    double M { rm * rm * m(rm) - R * R * m(R) };
    double beta0s { bb * (log(us * p.a / 2.0) + M_EULER) - 2.0 * h / (M_PI * p.a) };
    double piece0 { -(M / 2.0) / (bb * k) * (atan(beta0s / k) + M_PI / 2.0) };
    double ilog { panels(f_damped_log, q, log(us), log(u1), 0.25) };
    double ilin { panels(f_damped, q, u1, ug, M_PI / (rm + p.a)) };
    cancelled = M_PI * (rm * rm - R * R) / (4.0 * k);
    return cancelled + piece0 + ilog + ilin;
}

// pimplicitlipid_2D()'s three cases, from the two end integrals
double p_from(const Key& p, double kb, double pOff, double Isig, double IRmax, int& branch)
{
    if (!(4 * kb * Isig > pOff)) {
        branch = 0;
        return 4 * p.k * Isig;
    }
    if (4 * kb * IRmax > pOff) {
        branch = 2;
        return 4 * p.k * IRmax;
    }
    branch = 1;
    return p.k * pOff / kb;
}

const char* branchName[] { "R2D = sigma", "interior R2D", "R2D = Rmax" };

double rootf(double x, void* v) { return (*static_cast<std::function<double(double)>*>(v))(x); }

struct Result {
    double pErr {}, endErr {}, refSpread {};
    bool fail {};
};

Result check_key(paramsIL q, bool bench, bool verbose)
{
    Result res;
    Key p { q.sigma, q.Dtot, q.ka, q.dt };
    double kb { q.kb / 1.0e6 };
    paramsIL qq { q };
    double pOff { dissociate2D(qq) };
    double density { q.Nlipid / q.area };
    std::printf("key sigma=%g Dtot=%g ka=%.8g kb=%.8g dt=%g area=%.8g Na=%d Nlipid=%d: Rmax=%.6g, rmax=%.6g, ka/h=%.4g, "
                "pOff=%.6g\n",
        p.a, p.D, p.k, q.kb, p.t, q.area, q.Na, q.Nlipid, p.Rmax(), p.rmax(), p.k / p.h(), pOff);

    using clk = std::chrono::steady_clock;
    double i0 { instructions() };
    auto t0 { clk::now() };
    double pSim { pimplicitlipid_2D(q) };
    double i1 { instructions() };
    auto t1 { clk::now() };
    std::printf("  cost: %.4g instructions, %.4f s\n", i1 - i0, std::chrono::duration<double>(t1 - t0).count());
    if (bench)
        return res;
    gsl_set_error_handler_off();
    panelFailed = false;

    IL2DEndIntegrals ends { integrate_function2D_ends(q) };
    double X2 { 10.0 * p.X() };
    double Isig { refA(p, p.a, X2) };
    double IRmax { refA(p, p.Rmax(), X2) };
    double cs {}, cr {};
    double IsigB { refB(p, p.a, cs) };
    double IRmaxB { refB(p, p.Rmax(), cr) };
    int branch {};
    double pRef { p_from(p, kb, pOff, Isig, IRmax, branch) };
    double scale { 4 * p.k };

    // p on the scale of 4 ka I(sigma), the probability with no blocking,
    // which bounds it: at R2D = Rmax p can be 1e-6 of that.
    res.pErr = std::abs(pSim - pRef) / (scale * Isig);
    // Both ends on the scale of the larger, 4 ka I(sigma): I(Rmax) is often
    // orders below it and matters only next to that.
    res.endErr = std::max(std::abs(ends.atSigma - Isig), std::abs(ends.atRmax - IRmax)) / Isig;
    // refB carries the rounding of the term it cancels against: measured at
    // up to 1.2e-11 of it, on the keys with the smallest ka.
    double allowedB { 5e-11 * std::max(cs, cr) * scale + 1e-10 * pRef };
    res.refSpread = std::max(std::abs(Isig - IsigB), std::abs(IRmax - IRmaxB)) * scale;
    res.fail = res.pErr > 1e-8 || res.endErr > 1e-8 || res.refSpread > allowedB || panelFailed;
    if (panelFailed)
        std::printf("  a reference panel did not converge\n");

    std::printf("  reference: %s, p = %.15g (rho p = %.6g at density %.6g)\n", branchName[branch], pRef, density * pRef,
        density);
    std::printf("  pimplicitlipid_2D: p = %.15g, relative error %.3g (%.3g of 4 ka I(sigma))\n", pSim, (pSim - pRef) / pRef,
        res.pErr);
    std::printf("  4 ka I(sigma): sim %.15g ref %.15g | 4 ka I(Rmax): sim %.6g ref %.6g | max error / I(sigma) %.3g\n",
        scale * ends.atSigma, scale * Isig, scale * ends.atRmax, scale * IRmax, res.endErr);
    std::printf("  references: max 4 ka |refA - refB| %.3g, allowed %.3g (refB cancels %.3g)\n", res.refSpread, allowedB,
        scale * std::max(cs, cr));
    if (verbose) {
        double far { std::max(std::abs(refA(p, p.a, 3 * X2) - Isig), std::abs(refA(p, p.Rmax(), 3 * X2) - IRmax)) };
        std::printf("  refA with X2 three times further: max 4 ka change %.3g\n", scale * far);
        if (branch == 1) {
            std::function<double(double)> f { [&](double R) { return 4 * kb * refA(p, R, X2) - pOff; } };
            gsl_function F { rootf, &f };
            gsl_root_fsolver* s { gsl_root_fsolver_alloc(gsl_root_fsolver_brent) };
            gsl_root_fsolver_set(s, &F, p.a, p.Rmax());
            for (int it = 0; it < 100; ++it) {
                gsl_root_fsolver_iterate(s);
                if (gsl_root_test_interval(gsl_root_fsolver_x_lower(s), gsl_root_fsolver_x_upper(s), 1e-13, 0)
                    == GSL_SUCCESS)
                    break;
            }
            std::printf("  reference R2D = %.13g\n", gsl_root_fsolver_root(s));
            gsl_root_fsolver_free(s);
        }
    }
    return res;
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    gsl_set_error_handler_off();
    W = gsl_integration_workspace_alloc(100000);

    auto parse = [](char** v) {
        paramsIL q {};
        q.sigma = atof(v[0]);
        q.Dtot = atof(v[1]);
        q.ka = atof(v[2]);
        q.kb = atof(v[3]);
        q.dt = atof(v[4]);
        q.area = atof(v[5]);
        q.Na = atoi(v[6]);
        q.Nlipid = atoi(v[7]);
        return q;
    };
    if (argc > 1) {
        bool bench { std::strcmp(argv[1], "--bench") == 0 };
        int o { bench ? 2 : 1 };
        if (argc < o + 8) {
            std::fprintf(stderr,
                "usage: implicit_lipid_2d_accuracy_check [--bench] sigma Dtot ka kb dt area Na Nlipid\n");
            return 2;
        }
        Result res { check_key(parse(argv + o), bench, !bench) };
        return res.fail ? 1 : 0;
    }

    // One key from each implicit-lipid model that builds these values
    // (clathrin in vitro, pucadyil, the clathrin pioneer family, gag on a
    // sphere, membrane localization, secretion with kb = 0), then keys that
    // push the asymptotics: a large sigma, a large ka/h, D dt far from 1.
    // The last has kb below the 1e-9/s at which dissociate2D() returns 0,
    // which is the only way found to reach R2D = Rmax: otherwise that needs
    // the binding probability beyond Rmax, a few 1e-5 of the total, to
    // exceed kon dt, with kon >= 2 pi D / ln(b / sigma).
    const char* keys[][8] = {
        { "1", "1.09", "6.633333333333334e-05", "0.03", "3", "1000000", "3", "8999" },
        { "1", "0.51", "0.00022133333333333334", "0.2", "3", "1000000", "3", "6999" },
        { "1", "1.12", "13.611790946289451", "11.710315261218806", "1", "490000", "13", "9988" },
        { "1", "0.58", "0.018496238291903243", "1.0034870905463473", "1", "490000", "72", "9947" },
        { "1", "0.53", "2.5", "1.002", "3", "122500", "18", "9983" },
        { "1", "0.99", "0.277", "1.002", "3", "122500", "72", "9980" },
        { "1", "0.95", "86.795", "0.01", "0.1", "61575.216010359945", "200", "458" },
        { "1", "0.99", "1.6693", "1.00526", "0.1", "220900", "1330", "3622" },
        { "1", "0.11", "0.08527769693136615", "0", "0.1", "110224", "200", "4984" },
        { "5", "1", "3", "1", "0.5", "1000000", "10", "10000" },
        { "1", "0.2", "300", "100", "0.1", "100000", "10", "1000" },
        { "1", "10", "1", "10", "10", "1000000", "10", "10000" },
        { "1", "0.01", "1", "10", "0.01", "100000", "10", "10000" },
        { "1", "1", "1", "1e-10", "0.1", "100000", "10", "1000" },
    };
    int failures { 0 };
    for (const auto& key : keys) {
        char* v[8];
        for (int i = 0; i < 8; ++i)
            v[i] = const_cast<char*>(key[i]);
        Result res { check_key(parse(v), false, false) };
        if (res.fail) {
            std::printf("  FAIL\n");
            ++failures;
        }
    }
    if (failures) {
        std::printf("implicit_lipid_2d_accuracy_check: %d key(s) failed\n", failures);
        return 1;
    }
    std::printf("implicit_lipid_2d_accuracy_check: all keys within 1e-8 of the reference\n");
    return 0;
}
