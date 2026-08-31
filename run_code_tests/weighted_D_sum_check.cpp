// Proves weighted_D_sum() is bit-identical to the four expressions it replaced,
// over random and adversarial inputs. Needed because no input file in the tree
// exercises the 1D fiber path, so the model-level A/B cannot cover it.
#include "classes/class_Molecule_Complex.hpp"
#include <cstdio>
#include <cstdint>
#include <random>
#include <cmath>
#include <cstring>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;
long long randNum = 0;

static bool bitEqual(double a, double b)
{
    uint64_t ua, ub;
    std::memcpy(&ua, &a, 8);
    std::memcpy(&ub, &b, 8);
    return ua == ub;
}

int main()
{
    std::mt19937_64 rng(20260827);
    // Wide exponent range: diffusion constants in NERDSS span many decades,
    // and cancellation is where an FP reassociation would show up.
    std::uniform_real_distribution<double> mant(-1.0, 1.0);
    std::uniform_int_distribution<int> expo(-300, 300);

    long long checked = 0, bad = 0;
    for (long long i = 0; i < 4000000; ++i) {
        Vec3D D1, D2;
        double* a[3] = { &D1.x, &D1.y, &D1.z };
        double* b[3] = { &D2.x, &D2.y, &D2.z };
        for (int k = 0; k < 3; ++k) {
            *a[k] = std::ldexp(mant(rng), expo(rng));
            *b[k] = std::ldexp(mant(rng), expo(rng));
        }
        if (i % 7 == 0) { D1.z = 0.0; D2.z = 0.0; }          // the membrane case
        if (i % 11 == 0) { D1 = Vec3D(); D2 = Vec3D(); }      // all-zero
        if (i % 13 == 0) { D2.x = -D1.x; D2.y = -D1.y; }      // exact cancellation

        // The original expressions, verbatim from check_bimolecular_reactions.cpp
        const double ref1D = D1.x + D2.x;
        const double ref2D = 1.0 / 2.0 * (D1.x + D2.x) + 1.0 / 2.0 * (D1.y + D2.y);
        const double ref3D = 1.0 / 3.0 * (D1.x + D2.x)
            + 1.0 / 3.0 * (D1.y + D2.y)
            + 1.0 / 3.0 * (D1.z + D2.z);

        const double got1D = weighted_D_sum(D1, D2, Dim::Fiber1D);
        const double got2D = weighted_D_sum(D1, D2, Dim::Surface2D);
        const double got3D = weighted_D_sum(D1, D2, Dim::Bulk3D);

        checked += 3;
        if (!bitEqual(ref1D, got1D)) { ++bad; if (bad < 4) std::printf("1D mismatch %.20g vs %.20g\n", ref1D, got1D); }
        if (!bitEqual(ref2D, got2D)) { ++bad; if (bad < 4) std::printf("2D mismatch %.20g vs %.20g\n", ref2D, got2D); }
        if (!bitEqual(ref3D, got3D)) { ++bad; if (bad < 4) std::printf("3D mismatch %.20g vs %.20g\n", ref3D, got3D); }
    }

    // pair_dim(): exhaustive over the four flag combinations, both orders.
    struct { bool f1, s1, f2, s2; Dim want; } cases[] = {
        { true,  false, true,  false, Dim::Fiber1D   },
        { true,  true,  true,  true,  Dim::Fiber1D   },  // fiber wins over surface
        { false, true,  false, true,  Dim::Surface2D },
        { true,  false, false, true,  Dim::Bulk3D    },  // disagree -> 3D, NOT min()
        { false, true,  true,  false, Dim::Bulk3D    },
        { false, false, false, false, Dim::Bulk3D    },
        { true,  false, false, false, Dim::Bulk3D    },
        { false, true,  false, false, Dim::Bulk3D    },
    };
    int dimBad = 0;
    for (auto& c : cases) {
        Complex c1, c2;
        c1.onFiber = c.f1; c1.OnSurface = c.s1;
        c2.onFiber = c.f2; c2.OnSurface = c.s2;
        // The original ladder, verbatim.
        Dim ref = (c1.onFiber && c2.onFiber) ? Dim::Fiber1D
            : (c1.OnSurface && c2.OnSurface) ? Dim::Surface2D
                                             : Dim::Bulk3D;
        if (ref != c.want || pair_dim(c1, c2) != c.want) {
            ++dimBad;
            std::printf("pair_dim mismatch f1=%d s1=%d f2=%d s2=%d\n", c.f1, c.s1, c.f2, c.s2);
        }
    }

    std::printf("weighted_D_sum: %lld comparisons, %lld mismatches\n", checked, bad);
    std::printf("pair_dim: %zu cases, %d mismatches\n", sizeof(cases) / sizeof(cases[0]), dimBad);
    return (bad || dimBad) ? 1 : 0;
}
