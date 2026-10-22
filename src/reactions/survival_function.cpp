#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_sf_expint.h>

double survival_function(double x, void* p)
{
    double f{};
    IntegrandParams& params = *reinterpret_cast<IntegrandParams*>(p);

    if (params.k < 1.0 / 0.0) { // This gives the association probability
        double h = (2.0 * M_PI * params.a * params.D);

        double alp = h * x * j1(x * params.a) + params.k * j0(x * params.a);
        double bet = h * x * y1(x * params.a) + params.k * y0(x * params.a);
        double tet = alp * alp + bet * bet;
        double P = j0(x * params.a) * y1(x * params.a) - j1(x * params.a) * y0(x * params.a);
        double T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;

        f = T * P * (1 - exp(-params.D * params.t * x * x));
        f = f * params.a * params.k;
    } else { // absorbing boundary conditions... this gives the survival probability
        double alp = j0(x * params.a);
        double bet = y0(x * params.a);
        double tet = alp * alp + bet * bet;
        double P = j0(x * params.a) * y0(x * params.r0) - j0(x * params.r0) * y0(x * params.a);
        double T = 1.0 / tet;

        f = (2.0 / M_PI) * T * P * exp(-params.D * params.t * x * x) / x;
    }

    return f;
}

double survival_tail(const IntegrandParams& params, double X)
{
    // Beyond X the exponential has died (D t X^2 >= 46) and the Bessel
    // functions are in their large-argument regime.  There P is the Wronskian
    // -2/(pi a x), and with H_n = J_n + i Y_n
    //   T = -Im[ H0(x r0) / (h x H1(x a) + k H0(x a)) ].
    // Taking the Hankel expansions H_n(z) ~ sqrt(2/(pi z)) e^{i(z - n pi/2 - pi/4)}
    // (1 + i (4n^2 - 1)/(8z) - (4n^2 - 1)(4n^2 - 9)/(128 z^2)) through second
    // order, and k/(h x) along with them,
    //   T = -sqrt(a/r0)/(h x) [(1 - d2/x^2) cos(w x) + (e1/x) sin(w x) + O(x^-3)],
    //   w = r0 - a,  b = 3/(8a) + k/h,  e1 = 1/(8 r0) + b,
    //   d2 = 9/(128 r0^2) + 15/(128 a^2) + k/(8 a h) + b^2 + b/(8 r0),
    // so the integrand is C [cos(w x)/x^2 + e1 sin(w x)/x^3 - d2 cos(w x)/x^4]
    // with C = 2 k sqrt(a/r0)/(pi h).  Its integrals from X to infinity follow
    // by parts from the sine integral:
    //   I2 = cos(w X)/X - w (pi/2 - Si(w X)),
    //   I3 = sin(w X)/(2 X^2) + w I2/2,
    //   I4 = cos(w X)/(3 X^3) - w I3/3.
    // The remainder is O(C / X^4).
    const double h { 2.0 * M_PI * params.a * params.D };
    const double C { 2.0 * params.k * sqrt(params.a / params.r0) / (M_PI * h) };
    const double w { params.r0 - params.a };
    const double b { 3.0 / (8.0 * params.a) + params.k / h };
    const double e1 { 1.0 / (8.0 * params.r0) + b };
    const double d2 { 9.0 / (128.0 * params.r0 * params.r0) + 15.0 / (128.0 * params.a * params.a)
        + params.k / (8.0 * params.a * h) + b * b + b / (8.0 * params.r0) };
    const double I2 { cos(w * X) / X - w * (M_PI / 2.0 - gsl_sf_Si(w * X)) };
    const double I3 { sin(w * X) / (2.0 * X * X) + 0.5 * w * I2 };
    const double I4 { cos(w * X) / (3.0 * X * X * X) - w * I3 / 3.0 };
    return C * (I2 + e1 * I3 - d2 * I4);
}
