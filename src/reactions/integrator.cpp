#include "io/io.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>

double integrator(gsl_function F, IntegrandParams params, gsl_integration_workspace* w, double r0, double bindrad,
    double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*),
    const NumericalSettings::Integration& settings)
{
    double result {};
    double error {};
    const double xlow { 0 };
    const double epsabs { settings.tableAbsoluteError };
    const double epsrel { settings.tableRelativeError };
    const double publicationcrit { settings.fallbackError };

    /*For successful qagiu, takes ~0.01:0.3s, depending on R0.
      using qags even just twice takes ~4s-10s. qags seems necessary for large ka values. .*/
    std::chrono::steady_clock::time_point start;
    int status { gsl_integration_qagiu(&F, xlow, epsabs, epsrel, 1e6, w, &result, &error) };

    int it { 1 };
    /*Try again with a weaker criterion*/
    if (status != GSL_SUCCESS)
        status = gsl_integration_qagiu(&F, xlow, publicationcrit, publicationcrit, 1e6, w, &result, &error);

    /*Now try with a different upper integration bound. Issues can arise due to oscillations of Bessel functions,
      convergence can vary with choice of upper integration.
     */
    if (status != GSL_SUCCESS) {
        start = std::chrono::steady_clock::now();
        /*This value of umax is not optimized*/
        double umax { 10000.0 }; // sqrt(-log(DBL_MIN) / (Dtot * deltat));
        while (std::abs((*f)(umax, &params)) > settings.tailCutoff)
            umax = umax * 1.2;

        // std::cout << linebreak;
        // std::cout << std::setprecision(12);
        // std::cout << funcID << " qagiu failed with status: " << status << '\n';
        // std::cout << "No solution found with rel/abs error smaller than: " << publicationcrit << '\n';
        // std::cout << "   umax: " << umax << '\n';
        // std::cout << "   @time: " << deltat << '\n';
        // std::cout << "   ka: " << kr << '\n';
        // std::cout << "   D: " << Dtot << '\n';
        // std::cout << "   r0: " << r0 << '\n';
        // std::cout << "Truncation will be performed on the semi-infinite domain..." << '\n';
        // std::cout << std::flush;

        it = 0;
        while (status != GSL_SUCCESS) {
            /*Can we use qag, instead of qags?*/
            int key { 2 };
            status = gsl_integration_qag(&F, xlow, umax, epsabs, publicationcrit, 1e6, key, w, &result, &error);
            umax = umax * 0.9;
            it += 1;
        }

        // std::cout << "New integration upper bound " << umax / 0.9 << " found after " << it << " iterations." << '\n';
        // std::cout << "   err: " << error << '\n';
        // std::cout << linebreak;
        // std::cout << "Time of 2D table generation: "
        //           << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    }
    return result;
}

double integrate_survival(IntegrandParams params, gsl_integration_workspace* w, size_t limit,
    const NumericalSettings::Integration& settings)
{
    // The association probability for a finite rate.  Its integrand decays
    // only as cos((r0 - a) x)/x^2, which defeats qagiu: its error estimate
    // never settles, and the qag retries of integrator() then truncate the
    // range wherever the estimate happens to pass, which left entries off by
    // up to 1e-4.  Instead integrate [0, X] with GK61 and add the tail beyond X
    // in closed form (survival_tail).  X sits well inside the asymptotic
    // regime of that tail: many times past 1/a, where the Bessel expansions
    // start to hold, and past k/(2 pi a D), where the h x term takes over
    // from k in alpha and beta.  At 40 times the larger of the two, the tail
    // is off by 2e-8 or less on every key run_code_tests/table_2d_accuracy_check
    // covers; at 10 times, by up to 1e-5.
    const double h { 2.0 * M_PI * params.a * params.D };
    double X { 40.0 * std::max(1.0 / params.a, params.k / h) };
    X = std::max(X, sqrt(46.0 / (params.D * params.t))); // exp(-D t x^2) < 1e-20 beyond X

    // qag over [0, X] in one piece is not safe: with hundreds of periods in
    // the range, its first bisections leave intervals of tens of periods,
    // where the Gauss and Kronrod estimates can agree by accident (r0 = 25.8,
    // D t = 100: estimated 7e-8, actual 4e-6).  Start it instead on panels of
    // eight periods of J0(x r0), the fastest factor, which the 30-point Gauss
    // rule alone integrates to 1e-14, and share the tolerance among them.
    const double panelWidth { 16.0 * M_PI / (params.r0 + params.a) };
    const int nPanels { std::max(1, static_cast<int>(std::ceil(X / panelWidth))) };
    const double dx { X / nPanels };

    gsl_function F;
    F.function = &survival_function;
    F.params = reinterpret_cast<void*>(&params);
    double result {};
    for (int i = 0; i < nPanels; ++i) {
        double panel {};
        double error {};
        int status { gsl_integration_qag(&F, i * dx, (i + 1) * dx, settings.tableAbsoluteError / nPanels,
            settings.tableRelativeError, limit, GSL_INTEG_GAUSS61, w, &panel, &error) };
        if (status != GSL_SUCCESS)
            gsl_integration_qag(&F, i * dx, (i + 1) * dx, settings.fallbackError / nPanels, settings.fallbackError,
                limit, GSL_INTEG_GAUSS61, w, &panel, &error);
        result += panel;
    }
    return result + survival_tail(params, X);
}
