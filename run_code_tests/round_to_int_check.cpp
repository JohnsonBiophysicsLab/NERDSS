// Proves round_half_away_to_int(x) == int(std::round(x)), and so
// round_to_int(x) on the targets that use it, for every x the two 2D
// diffusion-snapping sites can pass and the cases rounding gets wrong most
// often: every half-integer and its two neighbouring doubles up to 2^20, both
// signs, zeros and denormals, the largest values below the range limit, and
// 2e7 random doubles spread over all exponents up to 2^31.
#include "math/round_to_int.hpp"
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); every check links all objects
unsigned long totMatches = 0;

namespace {
long long compared { 0 }, mismatches { 0 };
void test(double x)
{
    ++compared;
    const int expected { int(std::round(x)) };
    const int got { round_half_away_to_int(x) };
    if (expected != got && ++mismatches < 4)
        std::printf("  x = %.17g: round %d, round_half_away_to_int %d\n", x, expected, got);
    if (round_to_int(x) != expected && ++mismatches < 4)
        std::printf("  x = %.17g: round %d, round_to_int %d\n", x, expected, round_to_int(x));
}
} // namespace

int main()
{
    for (double sign : { 1.0, -1.0 }) {
        test(sign * 0.0);
        test(sign * std::numeric_limits<double>::denorm_min());
        test(sign * std::numeric_limits<double>::min());
        test(sign * 0.49999999999999994);
        for (int k = 0; k < (1 << 20); ++k) {
            const double h { sign * (k + 0.5) };
            test(h);
            test(std::nextafter(h, 0.0));
            test(std::nextafter(h, 2.0 * h));
            test(sign * k);
        }
        for (double big : { 2147483646.0, 2147483646.5, 2147483646.9999998, 1073741823.5 }) {
            test(sign * big);
            test(std::nextafter(sign * big, 0.0));
        }
    }
    std::mt19937_64 rng(20261026);
    std::uniform_real_distribution<double> exponent(-40.0, 31.0), unit(-1.0, 1.0);
    for (int i = 0; i < 20000000; ++i) {
        const double x { unit(rng) * std::exp2(exponent(rng)) };
        if (std::fabs(x) < 2147483646.0)
            test(x);
    }
    std::printf("round_to_int_check: %lld values compared, %lld mismatches\n", compared, mismatches);
    if (mismatches != 0) {
        std::printf("round_to_int_check: FAILED\n");
        return 1;
    }
    std::printf("round_to_int_check: all passed\n");
    return 0;
}
