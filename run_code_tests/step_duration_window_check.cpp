// Proves StepDurationWindow gives the "Avg timestep duration" sum the shifted
// std::vector gave: after every push, sum() equals, bit for bit,
// std::accumulate() over a vector that dropped its first entry and appended
// the new one, for window sizes 1 to 1000 and random durations.
#include "io/step_duration_window.hpp"
#include <cstdio>
#include <cstring>
#include <numeric>
#include <random>
#include <vector>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); every check links all objects
unsigned long totMatches = 0;

int main()
{
    using duration = std::chrono::duration<double>;
    std::mt19937_64 rng(20261026);
    std::uniform_real_distribution<double> step(1e-7, 3e-3);
    long long compared { 0 }, mismatches { 0 };
    for (std::size_t size : { std::size_t(1), std::size_t(2), std::size_t(7), std::size_t(1000) }) {
        const duration initial { 0.123456789 };
        std::vector<duration> shifted(size, initial);
        StepDurationWindow window(size, initial);
        for (int i = 0; i < 5000; ++i) {
            const duration d { step(rng) };
            shifted.erase(shifted.begin());
            shifted.emplace_back(d);
            window.push(d);
            const double a { std::accumulate(shifted.begin(), shifted.end(), duration { 0 }).count() };
            const double b { window.sum().count() };
            ++compared;
            if (std::memcmp(&a, &b, sizeof a) != 0 && ++mismatches < 4)
                std::printf("  size %zu step %d: vector %.17g, window %.17g\n", size, i, a, b);
        }
    }
    std::printf("step_duration_window_check: %lld sums compared, %lld mismatches\n", compared, mismatches);
    if (mismatches != 0) {
        std::printf("step_duration_window_check: FAILED\n");
        return 1;
    }
    std::printf("step_duration_window_check: all passed\n");
    return 0;
}
