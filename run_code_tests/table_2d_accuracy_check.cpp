// Accuracy and cost check for the 2D reaction tables (survMatrix, pirMatrix).
//
// Builds the tables for a (sigma, Dtot, ka, deltaT) key with the simulator's
// own create_normMatrix / create_survMatrix / create_pirMatrix and compares
// the entries against references computed here, independently of the code
// under test:
//
//   survival, finite ka.  The table holds the association probability
//       S(r0) = int_0^inf a k T P (1 - exp(-D t x^2)) dx.
//     refA: the integrand as written, integrated panel by panel (GSL qag,
//           GK61, one panel per half period of the fastest Bessel factor) on
//           [0, X2], plus the asymptotic tail beyond X2, for X2 an order of
//           magnitude past the production X (and, from the command line, two
//           orders).
//     refB: 2D Brownian motion is recurrent, so S -> 1 at long times and
//           S(r0) = 1 - int_0^inf a k T P exp(-D t x^2) dx,
//           an integrand that is Gaussian-damped and needs no tail.  Near
//           x = 0 it goes as 1/(x ln^2 x): below x_s it is integrated in closed
//           form (an arctangent), and from x_s to 1/r0 in ln x.
//     The two share nothing but the integrand; their agreement is what makes
//     them a reference.
//   pir.  int_0^inf x exp(-D t x^2) P(r) T(r0) / (2 pi) dx, Gaussian-damped:
//     panel-by-panel qag on [0, sqrt(90/(D t))].  For ka = inf (absorbing) only
//     pir is checked: its survival table still comes from the old integrator.
//
// Without arguments it runs the keys below, sampling the tables, and fails if
// an entry is off by more than 1e-6 or the two survival references disagree.
// With arguments it reports on one key:
//   table_2d_accuracy_check <sigma> <Dtot> <ka> <deltaT> [surv-stride] [pir-stride]
//   table_2d_accuracy_check --bench <sigma> <Dtot> <ka> <deltaT>
// --bench builds the three tables and prints the cost of each, with no
// reference.  RMax is computed the way determine_2D_bimolecular_reaction_probability
// does; ka may be inf.
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_expint.h>
#ifdef __APPLE__
#include <libproc.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;

namespace {

struct Key {
    double a, D, k, t, r0, r;
};

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

// The production integrands, copied so the reference does not move with the
// code under test.
void alpha_beta(const Key& p, double x, double& alp, double& bet)
{
    if (p.k < 1.0 / 0.0) {
        double h { 2.0 * M_PI * p.a * p.D };
        alp = h * x * j1(x * p.a) + p.k * j0(x * p.a);
        bet = h * x * y1(x * p.a) + p.k * y0(x * p.a);
    } else {
        alp = j0(x * p.a);
        bet = y0(x * p.a);
    }
}

double surv_noexp(double x, const Key& p) // a k T P
{
    double alp, bet;
    alpha_beta(p, x, alp, bet);
    double P { j0(x * p.a) * y1(x * p.a) - j1(x * p.a) * y0(x * p.a) };
    double T { (j0(x * p.r0) * bet - y0(x * p.r0) * alp) / (alp * alp + bet * bet) };
    return T * P * p.a * p.k;
}

double f_surv(double x, void* v)
{
    const Key& p { *static_cast<Key*>(v) };
    return surv_noexp(x, p) * (1.0 - exp(-p.D * p.t * x * x));
}

double f_surv_damped(double x, void* v)
{
    const Key& p { *static_cast<Key*>(v) };
    return surv_noexp(x, p) * exp(-p.D * p.t * x * x);
}

double f_surv_damped_log(double s, void* v) // in s = ln x
{
    double x { exp(s) };
    return f_surv_damped(x, v) * x;
}

double f_pir(double x, void* v)
{
    const Key& p { *static_cast<Key*>(v) };
    double alp, bet;
    alpha_beta(p, x, alp, bet);
    double tet { sqrt(alp * alp + bet * bet) };
    double P { (j0(x * p.r) * bet - y0(x * p.r) * alp) / tet };
    double T { (j0(x * p.r0) * bet - y0(x * p.r0) * alp) / tet };
    return x * exp(-p.D * p.t * x * x) * P * T / (2.0 * M_PI);
}

gsl_integration_workspace* W;
bool panelFailed { false };

double panels(double (*f)(double, void*), Key& p, double lo, double hi, double width)
{
    gsl_function F { f, &p };
    double sum {};
    int n { std::max(1, static_cast<int>(std::ceil((hi - lo) / width))) };
    double dx { (hi - lo) / n };
    for (int i = 0; i < n; ++i) {
        double r {}, e {};
        int st { gsl_integration_qag(&F, lo + i * dx, lo + (i + 1) * dx, 1e-17, 1e-13, 100000, 6, W, &r, &e) };
        if (st != GSL_SUCCESS && std::abs(e) > 1e-14)
            panelFailed = true;
        sum += r;
    }
    return sum;
}

// a k T P beyond X: C [cos(wx)/x^2 + e1 sin(wx)/x^3 - d2 cos(wx)/x^4], as in
// survival_tail(), but the reference only uses it where it is below 1e-12.
double tail(const Key& p, double X)
{
    double h { 2.0 * M_PI * p.a * p.D };
    double C { 2.0 * p.k * sqrt(p.a / p.r0) / (M_PI * h) };
    double w { p.r0 - p.a };
    double b { 3.0 / (8.0 * p.a) + p.k / h };
    double e1 { 1.0 / (8.0 * p.r0) + b };
    double d2 { 9.0 / (128.0 * p.r0 * p.r0) + 15.0 / (128.0 * p.a * p.a) + p.k / (8.0 * p.a * h) + b * b
        + b / (8.0 * p.r0) };
    double I2 { cos(w * X) / X - w * (M_PI / 2.0 - gsl_sf_Si(w * X)) };
    double I3 { sin(w * X) / (2.0 * X * X) + 0.5 * w * I2 };
    double I4 { cos(w * X) / (3.0 * X * X * X) - w * I3 / 3.0 };
    return C * (I2 + e1 * I3 - d2 * I4);
}

double refA(Key p, double X2)
{
    return panels(f_surv, p, 0.0, X2, M_PI / (p.r0 + p.a)) + tail(p, X2);
}

double refB(Key p)
{
    double h { 2.0 * M_PI * p.a * p.D };
    double xs { 1e-7 / p.r0 };
    double x1 { 1.0 / (p.r0 + sqrt(p.D * p.t)) };
    double xg { sqrt(90.0 / (p.D * p.t)) };
    // Below xs: J0(x r0) = 1, alpha = k, beta = b (ln(x a/2) + gamma) - 2h/(pi a)
    // and the numerator of T is the constant N0, so the integrand is
    // -(2k/pi) N0 / (x (k^2 + beta^2)), whose integral is an arctangent in beta.
    double b { 2.0 * p.k / M_PI };
    double N0 { b * log(p.a / p.r0) - 2.0 * h / (M_PI * p.a) };
    double betas { b * (log(xs * p.a / 2.0) + M_EULER) - 2.0 * h / (M_PI * p.a) };
    double piece0 { -(N0 / p.k) * (atan(betas / p.k) + M_PI / 2.0) };
    double ilog { panels(f_surv_damped_log, p, log(xs), log(x1), 0.25) };
    double ilin { panels(f_surv_damped, p, x1, xg, M_PI / (p.r0 + p.a)) };
    return 1.0 - (piece0 + ilog + ilin);
}

double refPir(Key p)
{
    double xg { sqrt(90.0 / (p.D * p.t)) };
    double x1 { 1.0 / (p.r0 + p.r + sqrt(p.D * p.t)) };
    return panels(f_pir, p, 0.0, x1, x1) + panels(f_pir, p, x1, xg, M_PI / (p.r0 + p.r + 2 * p.a));
}

struct Result {
    double survErr {}, survRefSpread {}, pirErr {}, pirMax {}, pirAsym {};
    int survOver {}, pirOver {}, survChecked {}, pirChecked {};
};

// verbose: also refA with X2 another order of magnitude out, and per-entry lines
Result check_key(double a, double D, double k, double t, int survStride, int pirStride, bool bench, bool verbose)
{
    Result res;
    Parameters params;
    params.timeStep = t;
    const double RMax { 3.5 * sqrt(4.0 * D * t) + a };
    const size_t n { size_lookup(a, D, params, RMax) };
    const double step { sqrt(D * t) / 50 };
    gsl_matrix* surv { gsl_matrix_alloc(2, n) };
    gsl_matrix* norm { gsl_matrix_alloc(2, n) };
    gsl_matrix* pir { gsl_matrix_alloc(n, n) };
    std::printf("key sigma=%g Dtot=%g ka=%g deltaT=%g: RMax=%.6g, %zu radii, k/(2 pi a D)=%.4g\n", a, D, k, t, RMax,
        n, k / (2 * M_PI * a * D));

    using clk = std::chrono::steady_clock;
    auto secs = [](clk::time_point x, clk::time_point y) { return std::chrono::duration<double>(y - x).count(); };
    double i0 { instructions() };
    auto t0 { clk::now() };
    create_normMatrix(norm, a, D, k, RMax, step, params);
    double i1 { instructions() };
    auto t1 { clk::now() };
    if (k < 1.0 / 0.0) {
        create_survMatrix(surv, a, D, k, RMax, step, params);
    } else {
        // ka = inf keeps the old qagiu integrator, unchanged and slow; only its
        // pir table is checked here, so lay out the abscissae the same way.
        double RIndex { a };
        for (size_t j = 0; j < n; ++j, RIndex += step)
            gsl_matrix_set(surv, 0, j, RIndex);
    }
    double i2 { instructions() };
    auto t2 { clk::now() };
    create_pirMatrix(pir, a, D, k, RMax, step, params);
    double i3 { instructions() };
    auto t3 { clk::now() };
    std::printf("  cost: norm %.4g instr %.3f s | surv %.4g instr %.3f s | pir %.4g instr %.3f s\n", i1 - i0,
        secs(t0, t1), i2 - i1, secs(t1, t2), i3 - i2, secs(t2, t3));
    if (bench)
        return res;
    gsl_set_error_handler_off(); // create_survMatrix restores the aborting default

    const double h { 2.0 * M_PI * a * D };
    const double X2 { 400.0 * std::max(1.0 / a, k / h) }; // ten times the production X
    if (k < 1.0 / 0.0) {
        double spreadFar {}, worstR0 {};
        for (size_t j = 0; j < n; j += survStride) {
            Key p { a, D, k, t, gsl_matrix_get(surv, 0, j), 0 };
            double B { refB(p) };
            double A { refA(p, X2) };
            double v { gsl_matrix_get(surv, 1, j) };
            double e { std::abs(v - B) };
            if (e > res.survErr) {
                res.survErr = e;
                worstR0 = p.r0;
            }
            res.survRefSpread = std::max(res.survRefSpread, std::abs(A - B));
            res.survOver += e > 1e-6;
            ++res.survChecked;
            if (verbose) {
                double A10 { refA(p, 10 * X2) };
                spreadFar = std::max(spreadFar, std::abs(A10 - B));
                if (getenv("TABLE_CHECK_VERBOSE"))
                    std::printf("    r0=%.6f table=%.12g refB=%.12g refA=%.12g diff=%.3g\n", p.r0, v, B, A, v - B);
            }
        }
        std::printf("  surv: %d entries, max |table - ref| %.3g at r0=%.6g, %d over 1e-6\n", res.survChecked,
            res.survErr, worstR0, res.survOver);
        std::printf("  surv references: max |refA - refB| %.3g", res.survRefSpread);
        if (verbose)
            std::printf(", with X2 ten times further %.3g", spreadFar);
        std::printf("\n");
    }

    for (size_t i = 0; i < n; i += pirStride)
        for (size_t j = 0; j < n; j += pirStride) {
            // pirMatrix(r0 index, r index), on the table's own abscissae
            Key p { a, D, k, t, gsl_matrix_get(surv, 0, j), gsl_matrix_get(surv, 0, i) };
            double ref { refPir(p) };
            double v { gsl_matrix_get(pir, j, i) };
            res.pirErr = std::max(res.pirErr, std::abs(v - ref));
            res.pirMax = std::max(res.pirMax, std::abs(ref));
            res.pirAsym = std::max(res.pirAsym, std::abs(v - gsl_matrix_get(pir, i, j)));
            res.pirOver += std::abs(v - ref) > 1e-6;
            ++res.pirChecked;
        }
    std::printf("  pir: %d entries, max |table - ref| %.3g (max |pir| %.4g), %d over 1e-6, max asymmetry %.3g\n",
        res.pirChecked, res.pirErr, res.pirMax, res.pirOver, res.pirAsym);
    gsl_matrix_free(surv);
    gsl_matrix_free(norm);
    gsl_matrix_free(pir);
    return res;
}

} // namespace

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    gsl_set_error_handler_off();
    W = gsl_integration_workspace_alloc(100000);

    if (argc > 1) {
        bool bench { std::strcmp(argv[1], "--bench") == 0 };
        int o { bench ? 2 : 1 };
        if (argc < o + 4) {
            std::fprintf(stderr,
                "usage: table_2d_accuracy_check [--bench] sigma Dtot ka deltaT [surv-stride] [pir-stride]\n");
            return 2;
        }
        check_key(atof(argv[o]), atof(argv[o + 1]), atof(argv[o + 2]), atof(argv[o + 3]),
            argc > o + 4 ? atoi(argv[o + 4]) : 1, argc > o + 5 ? atoi(argv[o + 5]) : 8, bench, true);
        if (panelFailed)
            std::printf("warning: a reference panel did not converge\n");
        return 0;
    }

    // The rev_2D key, the two gagsphere rates, an ap2clafchoeps key, and keys
    // that push the asymptotics: a large sigma, a large k/(2 pi a D), Dt far
    // from 1, and absorbing boundaries (pir only).
    const double keys[][4] = {
        { 2, 4, 100, 0.1 },
        { 1, 0.96, 10, 0.1 },
        { 1, 0.48, 0.25, 0.1 },
        { 1, 1.09, 0.44874013741642377, 3 },
        { 20, 1, 4.6, 3 },
        { 0.5, 0.3, 10, 0.1 },
        { 1, 0.01, 1, 0.01 },
        { 1, 10, 1, 10 },
        { 1, 1, 1.0 / 0.0, 0.1 },
    };
    int failures { 0 };
    for (const auto& key : keys) {
        Result res { check_key(key[0], key[1], key[2], key[3], 7, 23, false, false) };
        if (res.survOver || res.pirOver || res.survRefSpread > 1e-9 || res.pirAsym != 0.0) {
            std::printf("  FAIL\n");
            ++failures;
        }
    }
    if (panelFailed) {
        std::printf("FAIL: a reference panel did not converge\n");
        ++failures;
    }
    if (failures) {
        std::printf("table_2d_accuracy_check: %d key(s) failed\n", failures);
        return 1;
    }
    std::printf("table_2d_accuracy_check: all keys within 1e-6\n");
    return 0;
}
