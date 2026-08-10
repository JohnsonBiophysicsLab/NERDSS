/*
 * ERFLIKE: evaluation of complex-valued error-like functions via an
 * exponentially convergent trapezoidal rule.
 *
 * Original implementation by Federico Maria Guercilena (2024):
 * https://doi.org/10.5281/zenodo.11261631
 * The Zenodo source archive is distributed under CC BY 4.0:
 * https://creativecommons.org/licenses/by/4.0/
 *
 * NERDSS adaptation: the upstream writable static scratch variables are
 * automatic locals so evaluations from multiple threads do not race.
 */

#ifndef ERFLIKE_H
#define ERFLIKE_H


/*
h = 5.022020510055717734e-1
N = 12
*/


/*
    If this is C, check for complex arithmetic support of C99
    and include the headers. Use C complex numbers ("double complex z").
*/
#ifndef __cplusplus

#ifdef __STDC_NO_COMPLEX__
#error "This implementation does not support C complex arithmetic"
#endif

#include <math.h>
#include <complex.h>
#include <stddef.h>

#define COMPLEX double complex
#define J I
#define NAMING(NAME) erflike_##NAME
#define NAMING_RE(NAME) erflike_##NAME##_re
#define USE_STD

/*
    If this is C++, load the complex header and put everything in a namespace.
    Use C++ complex numbers ("std::complex<double> z").
*/
#else

#include <cmath>
#include <complex>

#define COMPLEX std::complex<double>
#define J std::complex<double>(0., 1.)
#define NAMING(NAME) NAME
#define NAMING_RE(NAME) NAME
#define USE_STD using namespace std;
#define cexp std::exp
#define cabs std::abs
#define creal std::real
#define cimag std::imag

namespace erflike
{

#endif

/* Real scaled complementary error function */
extern double NAMING_RE(erfcx)(const double x)
{
    USE_STD

    const double x2 = x * x;

    if (x < 0 && x2 > 4.1024025e+1)
    {
        return 2. * exp(x2);
    }

    const double ax = fabs(x);
    double r;

#ifndef PURE_TRAPEZOIDAL
    if (ax > 5.e7)
    {
        return 5.641895835477562869e-1 / ax;
    }
    else if (ax > 30.)
    {
        const double ix  = 1. / x;
        const double ix2 = ix * ix;

        return ix *
               (5.641895835477562869e-1 +
                ix2 *
                    (-2.820947917738781435e-1 +
                     ix2 * (4.231421876608172152e-1 +
                            ix2 * (-1.057855469152043038 +
                                   ix2 * (3.702494142032150633 +
                                          ix2 * (-1.666122363914467785e+1))))));
    }
    else if (ax < 0.08)
    {
        return exp(x2) *
               (1. -
                (x *
                 (1.128379167095512574 +
                  x2 *
                      (-3.761263890318375246e-1 +
                       x2 * (1.128379167095512574e-1 +
                             x2 * (-2.686617064513125176e-2 +
                                   x2 * (5.223977625442187842e-3 +
                                         x2 * (-8.548327023450852833e-4))))))));
    }
    else
    {
#endif

        r = ax * (3.001757391943707186e-1 / (x2 + 6.305172500855072826e-2) +
                  1.812639579936974154e-1 / (x2 + 5.674655250769565543e-1) +
                  6.609722077788082798e-2 / (x2 + 1.576293125213768206) +
                  1.455428278875993587e-2 / (x2 + 3.089534525418985685) +
                  1.935237580905366322e-3 / (x2 + 5.107189725692608989) +
                  1.553866269261963503e-4 / (x2 + 7.629258726034638119) +
                  7.534056269654093622e-6 / (x2 + 1.065574152644507308e+1) +
                  2.205870180031759912e-7 / (x2 + 1.418663812692391386e+1) +
                  3.900020681070635733e-9 / (x2 + 1.822194852747116047e+1) +
                  4.163798844433521215e-11 / (x2 + 2.27616727280868129e+1) +
                  2.684407482331670116e-13 / (x2 + 2.780581072877087116e+1) +
                  1.045064297123954134e-15 / (x2 + 3.335436252952333525e+1));

        if (ax <= 6.255634853141884637)
        {
            r += ((ax < 6.255634853141884637) + 1.) * exp(x2) /
                 (1.0 + exp(1.251126970628376927e+1 * ax));
        }

#ifndef PURE_TRAPEZOIDAL
    }
#endif

    if (x < 0)
    {
        r = 2. * exp(x2) - r;
    }

    return r;
}

/* Imaginary part of Faddeeva function of real argument */
extern double NAMING(w_im)(const double x)
{
    USE_STD

    const double x2 = x * x;
    double r;

#ifndef PURE_TRAPEZOIDAL
    if (x2 > 0.25e16)
    {
        return 5.641895835477562869e-1 / x;
    }
    if (x2 > 900.)
    {
        const double ix  = 1. / x;
        const double ix2 = ix * ix;

        r = ix *
            (5.641895835477562869e-1 +
             ix2 * (2.820947917738781435e-1 +
                    ix2 * (4.231421876608172152e-1 +
                           ix2 * (1.057855469152043038 +
                                  ix2 * (3.702494142032150633 +
                                         ix2 * (1.666122363914467785e+1))))));
    }
    else if (x2 < 0.0064)
    {
        return (x *
                (1.128379167095512574 +
                 x2 * (-7.522527780636750493e-1 +
                       x2 * (3.009011112254700197e-1 +
                             x2 * (-8.597174606442000563e-2 +
                                   x2 * (1.910483245876000125e-2 +
                                         x2 * (-3.4736059015927275e-3)))))));
    }
    else
    {
#endif

        double fuffa;
        const double x_frac_h = modf(fabs(x) * 1.991230418110947371, &fuffa);

        if (0.25 <= x_frac_h && x_frac_h <= 0.75)
        {
            r = x * (2.484428714219530481e-1 / (x2 - 2.52206900034202913e-1) +
                     1.165816895727056702e-1 / (x2 - 1.008827600136811652) +
                     3.303467360515457584e-2 / (x2 - 2.269862100307826217) +
                     5.652565002709598381e-3 / (x2 - 4.035310400547246609) +
                     5.840593495428009957e-4 / (x2 - 6.305172500855072826) +
                     3.644217301209265032e-5 / (x2 - 9.079448401231304869) +
                     1.373053377768266636e-6 / (x2 - 1.235813810167594274e+1) +
                     3.123967061563795513e-8 / (x2 - 1.614124160218898643e+1) +
                     4.292016090250947202e-10 / (x2 - 2.042875890277043596e+1) +
                     3.560837104703782992e-12 / (x2 - 2.52206900034202913e+1) +
                     1.783933832854722189e-14 / (x2 - 3.051703490413855248e+1) +
                     5.396861184587406183e-17 /
                         (x2 - 3.631779360492521948e+1)) +
                1.598558776968497887e-1 / x;

            if (x2 <= 3.80689e+1)
            {
                r += exp(-x2) * sin(1.251126970628376927e+1 * x) /
                     (cos(1.251126970628376927e+1 * x) - 1.);
            }
        }
        else
        {
            r = x * (3.001757391943707186e-1 / (x2 - 6.305172500855072826e-2) +
                     1.812639579936974154e-1 / (x2 - 5.674655250769565543e-1) +
                     6.609722077788082798e-2 / (x2 - 1.576293125213768206) +
                     1.455428278875993587e-2 / (x2 - 3.089534525418985685) +
                     1.935237580905366322e-3 / (x2 - 5.107189725692608989) +
                     1.553866269261963503e-4 / (x2 - 7.629258726034638119) +
                     7.534056269654093622e-6 / (x2 - 1.065574152644507308e+1) +
                     2.205870180031759912e-7 / (x2 - 1.418663812692391386e+1) +
                     3.900020681070635733e-9 / (x2 - 1.822194852747116047e+1) +
                     4.163798844433521215e-11 / (x2 - 2.27616727280868129e+1) +
                     2.684407482331670116e-13 / (x2 - 2.780581072877087116e+1) +
                     1.045064297123954134e-15 / (x2 - 3.335436252952333525e+1));

            if (x2 <= 3.80689e+1)
            {
                r += exp(-x2) * sin(1.251126970628376927e+1 * x) /
                     (cos(1.251126970628376927e+1 * x) + 1.);
            }
        }

#ifndef PURE_TRAPEZOIDAL
    }
#endif

    return r;
}

/* Complex Faddeeva function */
extern COMPLEX NAMING(w)(const COMPLEX z_in)
{
    USE_STD

    double RE        = creal(z_in);
    double IM        = cimag(z_in);
    const double RE2 = RE * RE;
    const double IM2 = IM * IM;

    if (RE == 0)
    {
        return NAMING_RE(erfcx)(IM) + 0. * J;
    }
    if (IM == 0)
    {
        return exp(-RE2) + J * NAMING(w_im)(RE);
    }

    COMPLEX z        = z_in;
    const COMPLEX z2 = z_in * z_in;
    COMPLEX r;
    int flip = 0;

    if (IM < 0)
    {
        z  = -z_in;
        RE = -RE;
        IM = -IM;

        if (IM2 - RE2 > 4.1024025e+1)
        {
            return 2. * cexp(-z2);
        }
        else if (IM2 - RE2 < -4.1024025e+1)
        {
            flip = 2;
        }
        else
        {
            flip = 1;
        }
    }

#ifndef PURE_TRAPEZOIDAL
    if (RE2 + IM2 > 0.25e16)
    {
        r = J * 5.641895835477562869e-1 / z;
    }
    if (RE2 + IM2 > 900.)
    {
        const COMPLEX iz  = 1. / z;
        const COMPLEX iz2 = iz * iz;

        r = J * iz *
            (5.641895835477562869e-1 +
             iz2 * (2.820947917738781435e-1 +
                    iz2 * (4.231421876608172152e-1 +
                           iz2 * (1.057855469152043038 +
                                  iz2 * (3.702494142032150633 +
                                         iz2 * (1.666122363914467785e+1))))));
    }
    else
    {
#endif

        double fuffa;
        const double RE_frac_h = modf(fabs(RE) * 1.991230418110947371, &fuffa);

        if (0.25 <= RE_frac_h && RE_frac_h <= 0.75)
        {
            r = J * z *
                    (2.484428714219530481e-1 / (z2 - 2.52206900034202913e-1) +
                     1.165816895727056702e-1 / (z2 - 1.008827600136811652) +
                     3.303467360515457584e-2 / (z2 - 2.269862100307826217) +
                     5.652565002709598381e-3 / (z2 - 4.035310400547246609) +
                     5.840593495428009957e-4 / (z2 - 6.305172500855072826) +
                     3.644217301209265032e-5 / (z2 - 9.079448401231304869) +
                     1.373053377768266636e-6 / (z2 - 1.235813810167594274e+1) +
                     3.123967061563795513e-8 / (z2 - 1.614124160218898643e+1) +
                     4.292016090250947202e-10 / (z2 - 2.042875890277043596e+1) +
                     3.560837104703782992e-12 / (z2 - 2.52206900034202913e+1) +
                     1.783933832854722189e-14 / (z2 - 3.051703490413855248e+1) +
                     5.396861184587406183e-17 /
                         (z2 - 3.631779360492521948e+1)) +
                1.598558776968497887e-1 * J / z;

            if ((IM <= 6.255634853141884637) && (RE2 <= 3.80689e+1))
            {
                r += ((IM < 6.255634853141884637) + 1.) * cexp(-z2) /
                     (1.0 - cexp(-1.251126970628376927e+1 * J * z));
            }
        }
        else
        {
            r = J * z *
                (3.001757391943707186e-1 / (z2 - 6.305172500855072826e-2) +
                 1.812639579936974154e-1 / (z2 - 5.674655250769565543e-1) +
                 6.609722077788082798e-2 / (z2 - 1.576293125213768206) +
                 1.455428278875993587e-2 / (z2 - 3.089534525418985685) +
                 1.935237580905366322e-3 / (z2 - 5.107189725692608989) +
                 1.553866269261963503e-4 / (z2 - 7.629258726034638119) +
                 7.534056269654093622e-6 / (z2 - 1.065574152644507308e+1) +
                 2.205870180031759912e-7 / (z2 - 1.418663812692391386e+1) +
                 3.900020681070635733e-9 / (z2 - 1.822194852747116047e+1) +
                 4.163798844433521215e-11 / (z2 - 2.27616727280868129e+1) +
                 2.684407482331670116e-13 / (z2 - 2.780581072877087116e+1) +
                 1.045064297123954134e-15 / (z2 - 3.335436252952333525e+1));

            if ((IM <= 6.255634853141884637) && (RE2 <= 3.80689e+1))
            {
                r += ((IM < 6.255634853141884637) + 1.) * cexp(-z2) /
                     (1.0 + cexp(-1.251126970628376927e+1 * J * z));
            }
        }
#ifndef PURE_TRAPEZOIDAL
    }
#endif

    if (flip == 1)
    {
        r = 2. * cexp(-z2) - r;
    }
    else if (flip == 2)
    {
        r = -r;
    }

    return r;
}

/* Complex scaled complementary error function */
extern COMPLEX NAMING(erfcx)(const COMPLEX z)
{
    USE_STD

    return NAMING(w)(J * z);
}

/* Complex complementary error function */
extern COMPLEX NAMING(erfc)(const COMPLEX z)
{
    USE_STD

    if (creal(z) >= 0.)
    {
        return cexp(-z * z) * NAMING(w)(J * z);
    }
    else
    {
        return 2. - cexp(-z * z) * NAMING(w)(-J * z);
    }
}

/* Real complementary error function */
extern double NAMING_RE(erfc)(const double x)
{
    USE_STD

#ifndef PURE_TRAPEZOIDAL
    // #ifndef __cplusplus
    //     return erfc(x); // C99 supplies erfc in math.h
    // #elif (__cplusplus >= 201103L)
    //     return std::erfc(x); // C++11 supplies std::erfc in cmath
    // #else
    if (x < -7.)
    {
        return 2.0;
    }
    else if (fabs(x) < 0.08)
    {
        const double x2 = x * x;
        return 1. -
               (x *
                (1.128379167095512574 +
                 x2 * (-3.761263890318375246e-1 +
                       x2 * (1.128379167095512574e-1 +
                             x2 * (-2.686617064513125176e-2 +
                                   x2 * (5.223977625442187842e-3 +
                                         x2 * (-8.548327023450852833e-4)))))));
    }
// #endif
#endif

    if (x >= 0.)
    {
        return exp(-x * x) * NAMING_RE(erfcx)(x);
    }
    else
    {
        return 2. - exp(-x * x) * NAMING_RE(erfcx)(-x);
    }
}

/* Complex error function */
extern COMPLEX NAMING(erf)(const COMPLEX z)
{
    USE_STD

    if (creal(z) >= 0.)
    {
        return 1. - cexp(-z * z) * NAMING(w)(J * z);
    }
    else
    {
        return cexp(-z * z) * NAMING(w)(-J * z) - 1.;
    }
}

/* Real error function */
extern double NAMING_RE(erf)(const double x)
{
    USE_STD

#ifndef PURE_TRAPEZOIDAL
    // #ifndef __cplusplus
    //    return erf(x); // C99 supplies erf in math.h
    // #elif (__cplusplus >= 201103L)
    //    return std::erf(x); // C++11 supplies std::erf in cmath
    // #else
    if (fabs(x) > 7.)
    {
        return copysign(1.0, x);
    }
    else if (fabs(x) < 0.08)
    {
        const double x2 = x * x;
        return (x *
                (1.128379167095512574 +
                 x2 * (-3.761263890318375246e-1 +
                       x2 * (1.128379167095512574e-1 +
                             x2 * (-2.686617064513125176e-2 +
                                   x2 * (5.223977625442187842e-3 +
                                         x2 * (-8.548327023450852833e-4)))))));
    }
// #endif
#endif

    if (x >= 0.)
    {
        return 1. - exp(-x * x) * NAMING_RE(erfcx)(x);
    }
    else
    {
        return exp(-x * x) * NAMING_RE(erfcx)(-x) - 1.;
    }
}

/* Complex imaginary error function */
extern COMPLEX NAMING(erfi)(const COMPLEX z)
{
    USE_STD

    if (cimag(z) <= 0.)
    {
        return -J + J * cexp(z * z) * NAMING(w)(-z);
    }
    else
    {
        return -J * cexp(z * z) * NAMING(w)(z) + J;
    }
}

/* Real imaginary error function */
extern double NAMING_RE(erfi)(const double x)
{
    USE_STD

    return exp(x * x) * NAMING(w_im)(x);
}

/* Complex Dawson integral */
extern COMPLEX NAMING(Dawson)(const COMPLEX z)
{
    USE_STD

    if (creal(z) >= 0.)
    {
        return J * 8.862269254527580136e-1 * (cexp(-z * z) - NAMING(w)(z));
    }
    else
    {
        return J * 8.862269254527580136e-1 * (NAMING(w)(-z) - cexp(-z * z));
    }
}

/* Real Dawson integral */
extern double NAMING_RE(Dawson)(const double x)
{
    USE_STD

    return 8.862269254527580136e-1 * NAMING(w_im)(x);
}


#ifdef __cplusplus
}

#undef cexp
#undef cabs
#undef creal
#undef cimag

#endif

#undef COMPLEX
#undef J
#undef NAMING
#undef NAMING_RE
#undef USE_STD

#endif
