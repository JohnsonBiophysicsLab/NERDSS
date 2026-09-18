// Proves that splitting calc_interface_distance() out of get_distance() changed
// nothing: the refactored get_distance() and the new calc_interface_distance()
// are compared bit for bit against a verbatim copy of the body they replaced.
//
// The split also rewrote the fiber promoter override of `sep` from
// `(a && !b) || (!a && b)` inside the 1D branch into `a != b` guarded by the
// same onFiber test, outside it.  No input file in the tree sets isPromoter, so
// nothing sets onFiber, so the model-level A/B cannot reach that arm at all;
// this covers it along with the sphere, flat-membrane and bulk arms.
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <random>
#include <vector>

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

// The original body, verbatim from src/reactions/get_distance.cpp at 8025fd7,
// minus the record_crossing_pair() call, which the split did not touch and
// which is checked separately below.
//
// noinline, like the real get_distance(): GCC contracts a*b + c into an FMA by
// default (-ffp-contract=fast), and inlined into main()'s loop it contracted
// this body differently from the standalone function, which moved the last bit
// of 1 in 100 flat-membrane distances.  That was the reference drifting, not the
// refactor: the pre-refactor get_distance.cpp compiled by the same GCC agrees
// with the refactored one bit for bit.
__attribute__((noinline)) static bool reference_get_distance(int pro1, int pro2, int iface1, int iface2, double& sep, double& R1, double Rmax,
    std::vector<Complex>& complexList, const ForwardRxn& currRxn, std::vector<Molecule>& moleculeList, bool isSphere)
{
    bool is2D = false;
    bool is1D = false;
    if (complexList[moleculeList[pro1].myComIndex].onFiber && complexList[moleculeList[pro2].myComIndex].onFiber) {
        is1D = true;
    } else if (complexList[moleculeList[pro1].myComIndex].OnSurface &&
               complexList[moleculeList[pro2].myComIndex].OnSurface) {
      is2D = true;
    }

    if (isSphere == true && is2D == true) {
      Vec3D iface11 = moleculeList[pro1].interfaceList[iface1].coord;
      Vec3D iface22 = moleculeList[pro2].interfaceList[iface2].coord;
      double r1 = iface11.length();
      double r2 = iface22.length();
      double r = (r1 + r2) / 2.0;
      double theta = acos((iface11.x * iface22.x + iface11.y * iface22.y + iface11.z * iface22.z) / r1 / r2);
      R1 = r * theta;
      sep = R1 - currRxn.bindRadius;
    } else if (is1D) {
      double coordx1{moleculeList[pro1].interfaceList[iface1].coord.x};
      double coordx2{moleculeList[pro2].interfaceList[iface2].coord.x};
      R1 = abs(coordx1 - coordx2);
      sep = R1 - currRxn.bindRadius;
      if (moleculeList[pro1].isPromoter && !moleculeList[pro2].isPromoter) {
        sep = R1;
      } else if (!moleculeList[pro1].isPromoter && moleculeList[pro2].isPromoter) {
        sep = R1;
      }
    } else {
      double dx = moleculeList[pro1].interfaceList[iface1].coord.x -
                  moleculeList[pro2].interfaceList[iface2].coord.x;
      double dy = moleculeList[pro1].interfaceList[iface1].coord.y -
                  moleculeList[pro2].interfaceList[iface2].coord.y;
      double dz{};
      if (is2D == true) {
          dz = 0;
      } else {
          dz = moleculeList[pro1].interfaceList[iface1].coord.z - moleculeList[pro2].interfaceList[iface2].coord.z;
      }
      R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
      sep = R1 - currRxn.bindRadius;
    }
    return R1 < Rmax;
}

int main()
{
    std::mt19937_64 rng(20260917);
    std::uniform_real_distribution<double> coord(-40.0, 40.0);
    std::uniform_real_distribution<double> radius(0.5, 5.0);

    long long checked = 0, bad = 0;
    long long arm[4] = { 0, 0, 0, 0 }; // sphere 2D, fiber 1D, flat 2D, bulk 3D
    long long promoterOverride = 0;

    // Exhaustive over every flag that selects an arm, both orders, with and
    // without a sphere; random coordinates within each combination.
    for (int fiber1 = 0; fiber1 < 2; ++fiber1)
    for (int fiber2 = 0; fiber2 < 2; ++fiber2)
    for (int surf1 = 0; surf1 < 2; ++surf1)
    for (int surf2 = 0; surf2 < 2; ++surf2)
    for (int prom1 = 0; prom1 < 2; ++prom1)
    for (int prom2 = 0; prom2 < 2; ++prom2)
    for (int sphere = 0; sphere < 2; ++sphere)
    for (int trial = 0; trial < 2000; ++trial) {
        std::vector<Complex> complexList(2);
        std::vector<Molecule> moleculeList(2);
        for (int m = 0; m < 2; ++m) {
            moleculeList[m].index = m;
            moleculeList[m].myComIndex = m;
            moleculeList[m].interfaceList.resize(1);
            moleculeList[m].interfaceList[0].coord = Vec3D { coord(rng), coord(rng), coord(rng) };
            complexList[m].index = m;
        }
        if (trial % 50 == 0) // coincident interfaces: a zero separation
            moleculeList[1].interfaceList[0].coord = moleculeList[0].interfaceList[0].coord;
        moleculeList[0].isPromoter = prom1;
        moleculeList[1].isPromoter = prom2;
        complexList[0].onFiber = fiber1;
        complexList[1].onFiber = fiber2;
        complexList[0].OnSurface = surf1;
        complexList[1].OnSurface = surf2;

        ForwardRxn rxn {};
        rxn.bindRadius = radius(rng);
        const double Rmax = 3.0 * rxn.bindRadius;

        const bool is1D = fiber1 && fiber2;
        const bool is2D = !is1D && surf1 && surf2;
        if (sphere && is2D) ++arm[0];
        else if (is1D) { ++arm[1]; if (prom1 != prom2) ++promoterOverride; }
        else if (is2D) ++arm[2];
        else ++arm[3];

        double refSep = 0.0, refR1 = 0.0;
        const bool refRet = reference_get_distance(0, 1, 0, 0, refSep, refR1, Rmax, complexList, rxn, moleculeList, sphere);

        const double gotDistance = calc_interface_distance(0, 1, 0, 0, complexList, moleculeList, sphere);

        double gotSep = 0.0, gotR1 = 0.0;
        const bool gotRet = get_distance(0, 1, 0, 0, /*rxnIndex*/ 3, /*rateIndex*/ 1, /*isStateChangeBackRxn*/ false,
            gotSep, gotR1, Rmax, complexList, rxn, moleculeList, sphere);

        ++checked;
        const bool same = refRet == gotRet && bitEqual(refSep, gotSep) && bitEqual(refR1, gotR1)
            && bitEqual(refR1, gotDistance);
        // The crossing get_distance() records is outside the split, but a
        // change there would be silent, so it is held to the same standard.
        const size_t wantCrossings = gotRet ? 1 : 0;
        const bool crossingsOk = moleculeList[0].crossings.size() == wantCrossings
            && moleculeList[1].crossings.size() == wantCrossings
            && complexList[0].ncross == static_cast<int>(wantCrossings)
            && complexList[1].ncross == static_cast<int>(wantCrossings);
        if (!same || !crossingsOk) {
            if (++bad < 6)
                std::printf("mismatch fiber=%d%d surf=%d%d prom=%d%d sphere=%d: ret %d/%d sep %.17g/%.17g "
                            "R1 %.17g/%.17g/%.17g crossings %s\n",
                    fiber1, fiber2, surf1, surf2, prom1, prom2, sphere, refRet, gotRet, refSep, gotSep, refR1, gotR1,
                    gotDistance, crossingsOk ? "ok" : "WRONG");
        }
    }

    std::printf("interface distance: %lld cases (sphere-2D %lld, fiber-1D %lld of which promoter override %lld, "
                "flat-2D %lld, bulk-3D %lld), %lld mismatches\n",
        checked, arm[0], arm[1], promoterOverride, arm[2], arm[3], bad);
    return bad ? 1 : 0;
}
