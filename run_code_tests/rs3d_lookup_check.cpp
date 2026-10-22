/*! \file rs3d_lookup_check.cpp
 * \brief Checks Membrane::RS3D_of_molType() against the scan it replaced.
 *
 * Every overlap sweep, create_complex_propagation_vectors(), check_overlap()
 * and SimulVolume's membrane check used to find a molecule's RS3D entry by
 * scanning RS3Dvect[400..499] for its type and reading RS3Dvect[300 + i] at
 * the first match, leaving RS3Dinput alone when nothing matched.  The lookup
 * keeps the first match per type instead.  This compares the two, member by
 * member, on random tables that hold repeated types, near misses on both sides
 * of the 1E-2 tolerance and the -1 the table starts from; on negative types,
 * which match those -1 entries; and after the table changes, through
 * forget_RS3D_matches() and through deserialize().
 */
#include "classes/class_Membrane.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <random>
#include <vector>

gsl_rng* r; // defined in main() by the simulator; unused here
unsigned long totMatches = 0;

namespace {

int failures = 0;
long compared = 0;

// The block each reader held before the lookup, verbatim apart from its input.
double scan(const Membrane& m, int molTypeIndex, double RS3Dinput)
{
    for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
        if (std::abs(m.RS3Dvect[RS3Dindex + 400] - molTypeIndex) < 1E-2) {
            RS3Dinput = m.RS3Dvect[RS3Dindex + 300];
            break;
        }
    }
    return RS3Dinput;
}

void fill(Membrane& m, std::mt19937_64& gen, int nTypes)
{
    std::uniform_real_distribution<double> value(0.0, 5.0);
    std::uniform_int_distribution<int> type(0, nTypes - 1);
    std::uniform_int_distribution<int> used(0, 100);
    std::uniform_int_distribution<int> kind(0, 9);
    m.RS3Dvect.assign(500, -1.0);
    const int nUsed { used(gen) };
    for (int i = 0; i < nUsed; ++i) {
        m.RS3Dvect[i + 300] = value(gen);
        double t = type(gen);
        switch (kind(gen)) {
        case 0: t += 0.0099; break; // still a match
        case 1: t -= 0.0099; break;
        case 2: t += 0.0101; break; // a near miss
        case 3: t -= 0.0101; break;
        default: break;
        }
        m.RS3Dvect[i + 400] = t;
    }
}

// Walk random "complexes" through both, the way the sweeps do: the last
// member with a match sets the value, starting from 0.0.
void compare(const Membrane& m, std::mt19937_64& gen, int nTypes, const char* label)
{
    std::uniform_int_distribution<int> type(-2, nTypes + 2);
    std::uniform_int_distribution<int> size(1, 12);
    for (int complex = 0; complex < 2000; ++complex) {
        double expected { 0.0 };
        double got { 0.0 };
        const int n { size(gen) };
        for (int member = 0; member < n; ++member) {
            const int t { type(gen) };
            expected = scan(m, t, expected);
            got = m.RS3D_of_molType(t, got);
            ++compared;
            if (std::memcmp(&expected, &got, sizeof got) != 0) {
                std::printf("FAIL: %s: type %d gives %.17g, scan gives %.17g\n", label, t, got, expected);
                ++failures;
                return;
            }
        }
    }
}

} // namespace

int main()
{
    std::mt19937_64 gen(20261022);
    for (int round = 0; round < 200; ++round) {
        const int nTypes { 1 + round % 40 };
        Membrane m;
        fill(m, gen, nTypes);
        compare(m, gen, nTypes, "fresh table");

        // A table that changes after lookups must be forgotten first.
        fill(m, gen, nTypes);
        m.forget_RS3D_matches();
        compare(m, gen, nTypes, "after forget_RS3D_matches()");

        // deserialize() replaces the table, so it forgets for itself.
        Membrane other;
        fill(other, gen, nTypes);
        std::vector<unsigned char> buf(1 << 16, 0);
        int written { 0 };
        other.serialize(buf.data(), written);
        int read { 0 };
        m.deserialize(buf.data(), read);
        compare(m, gen, nTypes, "after deserialize()");
    }

    if (failures == 0)
        std::printf("rs3d lookup: %ld lookups match the scan\n", compared);
    return failures ? 1 : 0;
}
