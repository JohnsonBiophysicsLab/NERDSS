#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

double passocF(double r0, double tCurr, double Dtot, double bindRadius, double alpha, double cof)
{
    const double fDt { 4.0 * Dtot * tCurr };
    const double sqrtfDt { sqrt(fDt) };

    const double f1 { cof * bindRadius / r0 };

    const double sqrttCurr { sqrt(tCurr) };
    const double a2 { alpha * alpha };
    const double sep { (r0 - bindRadius) / sqrtfDt }; // a

    const double e1 { 2.0 * sep * sqrttCurr * alpha + a2 * tCurr };
    const double ef1 { sep + alpha * sqrttCurr }; // a+b
    const double ep1 { exp(e1) };

    double term1 { erfc(sep) };
    double term2 {};
    // erfc(ef1) can underflow before the complete product does. Use the
    // scaled form if that happens, as well as when exp(e1) overflows.
    if (!std::isinf(ep1)) {
        const double erfc1 { erfc(ef1) };
        if (erfc1 != 0.0)
            term2 = ep1 * erfc1;
        else
            term2 = exp(-sep * sep) * Faddeeva::erfcx(ef1);
    } else {
        term2 = exp(-sep * sep) * Faddeeva::erfcx(ef1);
    }

    return (term1 - term2) * f1;
}
