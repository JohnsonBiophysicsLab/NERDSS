// Proves Interface::for_each_excluded_iface() reports the same hits, in the
// same order, as the loops it replaced in check_bimolecular_reactions(): every
// interface of the other molecule in index order, and for each one the
// exclusion entries naming it in list order.  The crossings are recorded in that
// order.  No input file in the tree has an interface whose exclusion entries are
// both out of sorted order and hit together in one pair, so the model-level A/B
// cannot cover the ordering.
#include "classes/class_MolTemplate.hpp"
#include <cstdio>
#include <random>
#include <tuple>
#include <vector>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;

using Hit = std::tuple<int, double, int>;

int main()
{
    std::mt19937_64 rng(20261020);
    std::uniform_int_distribution<int> listLen(0, 12);
    std::uniform_int_distribution<int> molType(0, 3);
    std::uniform_int_distribution<int> iface(-1, 8); // -1 and >= numIfaceOther are out of range
    std::uniform_int_distribution<int> numIface(0, 8);
    std::uniform_int_distribution<int> smallInt(0, 3); // few values, so duplicates are common

    long long cases = 0, bad = 0, hits = 0;
    for (int trial = 0; trial < 200000; ++trial) {
        Interface boundIface;
        const int n { listLen(rng) };
        for (int i = 0; i < n; ++i) {
            boundIface.excludeVolumeBoundList.push_back(molType(rng));
            boundIface.excludeVolumeBoundIfaceList.push_back(iface(rng));
            boundIface.excludeRadiusList.push_back(1.0 + smallInt(rng));
            boundIface.excludeVolumeBoundReactList.push_back(i); // identifies the entry
        }
        boundIface.build_excluded_ifaces();

        for (int otherType = 0; otherType <= 3; ++otherType) {
            const int numIfaceOther { numIface(rng) };

            // The original loops, verbatim apart from the names.
            std::vector<Hit> want;
            for (int relIfaceOther { 0 }; relIfaceOther < numIfaceOther; ++relIfaceOther) {
                for (unsigned indexItr { 0 }; indexItr < boundIface.excludeVolumeBoundList.size(); ++indexItr) {
                    if (boundIface.excludeVolumeBoundList[indexItr] == otherType
                        && boundIface.excludeVolumeBoundIfaceList[indexItr] == relIfaceOther) {
                        want.emplace_back(relIfaceOther, boundIface.excludeRadiusList[indexItr],
                            boundIface.excludeVolumeBoundReactList[indexItr]);
                    }
                }
            }

            std::vector<Hit> got;
            boundIface.for_each_excluded_iface(otherType, numIfaceOther,
                [&](int relIfaceOther, double radius, int rxnIndex) { got.emplace_back(relIfaceOther, radius, rxnIndex); });

            ++cases;
            hits += static_cast<long long>(want.size());
            if (got != want) {
                ++bad;
                if (bad < 4)
                    std::printf("mismatch: trial %d, other type %d: %zu hits expected, %zu reported\n", trial,
                        otherType, want.size(), got.size());
            }
        }
    }

    std::printf("excluded_ifaces_check: %lld cases, %lld hits, %lld mismatches\n", cases, hits, bad);
    if (bad != 0)
        return 1;
    std::printf("excluded_ifaces_check: all passed\n");
    return 0;
}
