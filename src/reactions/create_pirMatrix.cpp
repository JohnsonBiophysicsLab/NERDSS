#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_integration.h>

#include <cstdlib>
#include <vector>

void create_pirMatrix(gsl_matrix*& pirMatrix, double bindRadius, double Dtot, double kr, double comRMax,
    double RStepSize, const Parameters& params)
{
    // pir(r, r0) = int_0^inf x exp(-D t x^2) g(x; r) g(x; r0) / (2 pi) dx, with
    //   g(x; r) = (J0(x r) beta - Y0(x r) alpha) / sqrt(alpha^2 + beta^2),
    // where alpha and beta depend on x alone.  The integrand factors, so one
    // quadrature rule serves every entry: evaluate g once per (node, radius)
    // and the table is sum_n w_n g_n(r) g_n(r0), symmetric by construction.
    // Integrating each entry on its own, with qagiu, took 124 000 integrals of
    // about 300 evaluations each.
    //
    // The rule is composite 10-point Gauss-Legendre in x, scaled by
    // s = sqrt(D t).  Up to 8/s (exp(-64) beyond) the panels are 0.25/s wide.
    // Since r, r0 <= RMax = a + 7 s, g oscillates no faster than
    // cos((r - a) x), so a panel spans at most 3.5 radians of the product.
    // Below 0.25/s the panels halve toward zero, 40 times: g goes as 1/ln(x)
    // there (alpha -> k, beta -> (2k/pi) ln x), smooth in ln x but not in x,
    // and the integrand below the last panel is O(x^2), under 1e-25 / s^2.
    std::vector<double> radii;
    double RIndex { bindRadius };
    while (RIndex <= comRMax + RStepSize) {
        radii.push_back(RIndex);
        RIndex += RStepSize;
    }
    const size_t nRadii { radii.size() };

    const double s { sqrt(Dtot * params.timeStep) };
    const double width { 0.25 / s };
    const int nUniform { 31 }; // [width, 32 width] = [0.25/s, 8/s]
    const int nGraded { 40 }; // [width 2^-40, width]
    const size_t order { 10 };
    gsl_integration_glfixed_table* gl { gsl_integration_glfixed_table_alloc(order) };

    std::vector<double> nodes;
    std::vector<double> weights;
    auto add_panel = [&](double lo, double hi) {
        for (size_t i = 0; i < order; ++i) {
            double xi {};
            double wi {};
            gsl_integration_glfixed_point(lo, hi, i, &xi, &wi, gl);
            nodes.push_back(xi);
            weights.push_back(wi);
        }
    };
    double hi { width };
    for (int j = 0; j < nGraded; ++j) {
        add_panel(hi / 2.0, hi);
        hi /= 2.0;
    }
    for (int m = 1; m <= nUniform; ++m)
        add_panel(m * width, (m + 1) * width);
    gsl_integration_glfixed_table_free(gl);
    const size_t nNodes { nodes.size() };

    // g[i * nNodes + n] = g(x_n; r_i); gw additionally carries the weight and
    // the r-independent factor x exp(-D t x^2) / (2 pi).
    std::vector<double> g(nRadii * nNodes);
    std::vector<double> gw(nRadii * nNodes);
    const double h { 2.0 * M_PI * bindRadius * Dtot };
    for (size_t n = 0; n < nNodes; ++n) {
        const double x { nodes[n] };
        double alp {};
        double bet {};
        if (kr < 1.0 / 0.0) {
            alp = h * x * j1(x * bindRadius) + kr * j0(x * bindRadius);
            bet = h * x * y1(x * bindRadius) + kr * y0(x * bindRadius);
        } else { // absorbing boundary conditions
            alp = j0(x * bindRadius);
            bet = y0(x * bindRadius);
        }
        const double invTet { 1.0 / sqrt(alp * alp + bet * bet) };
        const double kernel { weights[n] * x * exp(-Dtot * params.timeStep * x * x) / (2.0 * M_PI) };
        for (size_t i = 0; i < nRadii; ++i) {
            const double gi { (j0(x * radii[i]) * bet - y0(x * radii[i]) * alp) * invTet };
            g[i * nNodes + n] = gi;
            gw[i * nNodes + n] = gi * kernel;
        }
    }

    // pirMatrix(r0 index, r index), as the per-entry loop filled it.
    for (size_t i = 0; i < nRadii; ++i) {
        const double* gwi { &gw[i * nNodes] };
        for (size_t j = i; j < nRadii; ++j) {
            const double* gj { &g[j * nNodes] };
            double sum {};
            for (size_t n = 0; n < nNodes; ++n)
                sum += gwi[n] * gj[n];
            gsl_matrix_set(pirMatrix, i, j, sum);
            gsl_matrix_set(pirMatrix, j, i, sum);
        }
    }
}
