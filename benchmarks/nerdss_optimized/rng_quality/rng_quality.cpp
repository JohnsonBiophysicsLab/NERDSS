/*! \file rng_quality.cpp
 *
 * Standalone evidence for the two changes that cannot be validated by bitwise
 * comparison, because they deliberately alter the random stream:
 *
 *   issue #10 -- random orientations were sampled by normalizing four
 *                independent U(-1,1) components, which is not uniform over
 *                rotations.  Replaced by Shoemake's subgroup algorithm.
 *   issue #12 -- GaussV() used a hand-written Marsaglia polar sampler.
 *                Replaced by gsl_ran_gaussian_ziggurat().
 *
 * Both the old and the new sampler are reproduced here so they can be measured
 * side by side with the same generator NERDSS uses (gsl_rng_default, seeded the
 * same way).
 *
 * Build and run:
 *   ./benchmarks/nerdss_optimized/rng_quality/run.sh
 */

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

gsl_rng* rng = nullptr;

double uniform() { return gsl_rng_uniform(rng); }

struct Quat {
    double w, x, y, z;
};

//! The sampler NERDSS used before issue #10: four U(-1,1) components, then
//! normalize.  Uniform inside a 4-cube, not uniform on the 3-sphere.
Quat naive_quat()
{
    Quat q { uniform() * 2 - 1, uniform() * 2 - 1, uniform() * 2 - 1, uniform() * 2 - 1 };
    const double mag { std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z) };
    return { q.w / mag, q.x / mag, q.y / mag, q.z / mag };
}

//! The replacement: Shoemake's subgroup algorithm, exactly uniform on the
//! 3-sphere, three uniform variates, no rejection.
Quat shoemake_quat()
{
    const double u1 { uniform() };
    const double u2 { uniform() };
    const double u3 { uniform() };

    const double r1 { std::sqrt(1.0 - u1) };
    const double r2 { std::sqrt(u1) };
    const double theta1 { 2.0 * M_PI * u2 };
    const double theta2 { 2.0 * M_PI * u3 };

    return { r2 * std::cos(theta2), r1 * std::sin(theta1), r1 * std::cos(theta1), r2 * std::sin(theta2) };
}

//! z component of the rotation of (0,0,1) by q, which is what a random molecular
//! orientation actually produces.  For a uniformly random rotation this is
//! uniform on [-1,1]; NERDSS relies on that when it orients molecules.
double rotated_z(const Quat& q)
{
    // Third column of the rotation matrix built from q.
    return q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z;
}

struct Deviation {
    double chiSquare;
    int bins;
    double maxRelDeviation;
};

//! Chi-square of a histogram against a uniform expectation.
Deviation uniform_deviation(const std::vector<long>& counts, long total)
{
    const int bins { static_cast<int>(counts.size()) };
    const double expected { static_cast<double>(total) / bins };
    double chiSquare { 0.0 };
    double maxRel { 0.0 };
    for (long c : counts) {
        const double diff { c - expected };
        chiSquare += (diff * diff) / expected;
        const double rel { std::fabs(diff) / expected };
        if (rel > maxRel)
            maxRel = rel;
    }
    return { chiSquare, bins, maxRel };
}

void report_orientation(long samples, int bins)
{
    std::printf("\n== issue #10: uniformity of random orientations ==\n");
    std::printf("Statistic: z component of (0,0,1) rotated by the sampled quaternion.\n");
    std::printf("For uniformly distributed rotations this is exactly uniform on [-1,1].\n");
    std::printf("Samples per sampler: %ld, histogram bins: %d\n\n", samples, bins);

    struct Row {
        const char* name;
        Quat (*sample)();
    };
    const Row rows[2] = { { "old: normalized U(-1,1)^4", naive_quat }, { "new: Shoemake subgroup", shoemake_quat } };

    std::printf("%-28s %14s %8s %14s %18s\n", "sampler", "chi^2", "dof", "chi^2/dof", "max bin deviation");
    for (const Row& row : rows) {
        gsl_rng_set(rng, 20260810);
        std::vector<long> counts(bins, 0);
        for (long i = 0; i < samples; ++i) {
            double z { rotated_z(row.sample()) };
            if (z < -1.0)
                z = -1.0;
            if (z > 1.0)
                z = 1.0;
            int bin { static_cast<int>((z + 1.0) * 0.5 * bins) };
            if (bin >= bins)
                bin = bins - 1;
            ++counts[bin];
        }
        const Deviation dev { uniform_deviation(counts, samples) };
        std::printf("%-28s %14.1f %8d %14.2f %17.2f%%\n", row.name, dev.chiSquare, dev.bins - 1,
            dev.chiSquare / (dev.bins - 1), dev.maxRelDeviation * 100.0);
    }
    std::printf("\nchi^2/dof near 1 means the sampler matches the uniform expectation.\n");
}

//! The sampler GaussV() used before issue #12.
double polar_gauss()
{
    double R { 2.0 };
    double V1 {};
    while (R >= 1.0) {
        V1 = 2.0 * uniform() - 1.0;
        const double V2 { 2.0 * uniform() - 1.0 };
        R = (V1 * V1) + (V2 * V2);
    }
    return (V1 * std::sqrt(-2.0 * std::log(R) / R));
}

double ziggurat_gauss() { return gsl_ran_gaussian_ziggurat(rng, 1.0); }

void report_gauss(long samples, int bins)
{
    std::printf("\n== issue #12: GaussV() distribution and speed ==\n");
    std::printf("Samples per sampler: %ld\n\n", samples);

    struct Row {
        const char* name;
        double (*sample)();
    };
    const Row rows[2] = { { "old: Marsaglia polar", polar_gauss }, { "new: GSL ziggurat", ziggurat_gauss } };

    // Bin over [-5,5]; the tails outside that hold ~6e-7 of the mass and are
    // counted separately so nothing is silently dropped.
    const double lo { -5.0 };
    const double hi { 5.0 };

    std::printf("%-24s %10s %10s %10s %10s %12s %12s %10s\n", "sampler", "mean", "variance", "skewness", "kurtosis",
        "chi^2/dof", "outside +-5", "seconds");
    for (const Row& row : rows) {
        gsl_rng_set(rng, 20260810);
        std::vector<long> counts(bins, 0);
        long outside { 0 };
        double m1 { 0 }, m2 { 0 }, m3 { 0 }, m4 { 0 };

        const auto start = std::chrono::steady_clock::now();
        for (long i = 0; i < samples; ++i) {
            const double v { row.sample() };
            m1 += v;
            m2 += v * v;
            m3 += v * v * v;
            m4 += v * v * v * v;
            if (v < lo || v >= hi) {
                ++outside;
            } else {
                int bin { static_cast<int>((v - lo) / (hi - lo) * bins) };
                if (bin >= bins)
                    bin = bins - 1;
                ++counts[bin];
            }
        }
        const auto finish = std::chrono::steady_clock::now();
        const double seconds { std::chrono::duration<double>(finish - start).count() };

        const double n { static_cast<double>(samples) };
        const double mean { m1 / n };
        const double var { m2 / n - mean * mean };
        const double sd { std::sqrt(var) };
        const double skew { (m3 / n - 3 * mean * var - mean * mean * mean) / (sd * sd * sd) };
        const double kurt { (m4 / n - 4 * mean * (m3 / n) + 6 * mean * mean * (m2 / n) - 3 * mean * mean * mean * mean)
            / (var * var) };

        // Expected counts from the standard normal CDF over each bin.
        double chiSquare { 0.0 };
        int usedBins { 0 };
        for (int b = 0; b < bins; ++b) {
            const double left { lo + (hi - lo) * b / bins };
            const double right { lo + (hi - lo) * (b + 1) / bins };
            const double p { 0.5 * (std::erf(right / std::sqrt(2.0)) - std::erf(left / std::sqrt(2.0))) };
            const double expected { p * n };
            if (expected < 10.0)
                continue; // chi-square is unreliable on nearly empty bins
            const double diff { counts[b] - expected };
            chiSquare += (diff * diff) / expected;
            ++usedBins;
        }

        std::printf("%-24s %10.5f %10.5f %10.5f %10.5f %12.2f %12ld %10.3f\n", row.name, mean, var, skew, kurt,
            chiSquare / (usedBins - 1), outside, seconds);
    }
    std::printf("\nA standard normal has mean 0, variance 1, skewness 0, kurtosis 3.\n");
}

} // namespace

int main(int argc, char** argv)
{
    const long orientationSamples { argc > 1 ? std::atol(argv[1]) : 4000000L };
    const long gaussSamples { argc > 2 ? std::atol(argv[2]) : 20000000L };

    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);

    std::printf("generator: %s\n", gsl_rng_name(rng));

    report_orientation(orientationSamples, 50);
    report_gauss(gaussSamples, 100);

    gsl_rng_free(rng);
    return 0;
}
