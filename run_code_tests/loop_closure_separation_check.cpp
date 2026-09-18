// Calls associate() on two molecules that are already in one complex and checks
// that it bonds them only when their interfaces are within
// bindRadSameCom * bindRadius -- the same criterion evaluate_binding_within_complex()
// applies before it will record a same-complex crossing at all.
//
// A loop closure bonds in place, moving nothing, so a stale crossing that
// reaches the loop-closure branch after an earlier association in the same
// timestep merged the pair's complexes used to write a bond between interfaces
// many sigma apart.  run_code_tests/BondSiteSeparation reproduces that
// end-to-end, but only on seeds whose trajectories happen to line up; this
// drives associate() directly, so it does not depend on any trajectory.
#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "reactions/association/association.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main(); unused here
unsigned long totMatches = 0;
long long randNum = 0;

namespace {

// Species indices for the one reaction A(a) + A(b) -> A(a!1).A(b!1).
constexpr int kReactantA = 0;
constexpr int kReactantB = 1;
constexpr int kProduct = 2;

struct Fixture {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> templates;
    std::vector<Molecule> molecules;
    std::vector<Complex> complexes;
    ForwardRxn rxn {};
    copyCounters counters {};
    std::map<std::string, int> observables;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::ofstream assocDissocFile; // never opened: associate() skips the log

    // Two molecules of one type, each with interfaces `a` (0) and `b` (1), both
    // members of complex 0.  Molecule 0's `a` sits at the origin and molecule
    // 1's `b` at (separation, 0, 0); nothing else about the geometry matters to
    // a loop closure, which moves nothing.
    Fixture(double separation, double bindRadius, double bindRadSameCom)
    {
        params.numMolTypes = 1;

        MolTemplate oneTemplate {};
        oneTemplate.molName = "A";
        templates.push_back(oneTemplate);

        for (int m = 0; m < 2; ++m) {
            Molecule mol {};
            mol.index = m;
            mol.myComIndex = 0;
            mol.molTypeIndex = 0;
            mol.interfaceList.resize(2);
            mol.freelist = { 0, 1 };
            molecules.push_back(mol);
        }
        molecules[0].interfaceList[0].coord = Vec3D { 0.0, 0.0, 0.0 };
        molecules[1].interfaceList[1].coord = Vec3D { separation, 0.0, 0.0 };

        Complex com {};
        com.index = 0;
        com.memberList = { 0, 1 };
        complexes.push_back(com);

        rxn.bindRadius = bindRadius;
        rxn.bindRadSameCom = bindRadSameCom;
        rxn.reactantListNew.resize(2);
        rxn.reactantListNew[0].absIfaceIndex = kReactantA;
        rxn.reactantListNew[0].relIfaceIndex = 0;
        rxn.reactantListNew[1].absIfaceIndex = kReactantB;
        rxn.reactantListNew[1].relIfaceIndex = 1;
        rxn.productListNew.resize(1);
        rxn.productListNew[0].absIfaceIndex = kProduct;

        counters.nBoundPairs.assign(1, 0);
        counters.copyNumSpecies = { 10, 10, 0 };
        counters.canDissociate.assign(3, false);
        counters.bindPairList.resize(3);
    }

    void associate_pair()
    {
        associate(/*iter*/ 1, /*ifaceIndex1*/ 0, /*ifaceIndex2*/ 1, molecules[0], molecules[1], complexes[0],
            complexes[0], params, rxn, molecules, templates, observables, counters, complexes, membrane, forwardRxns,
            backRxns, assocDissocFile);
    }

    bool bonded() const
    {
        return molecules[0].interfaceList[0].isBound && molecules[1].interfaceList[1].isBound
            && molecules[0].interfaceList[0].interaction.partnerIndex == 1
            && molecules[1].interfaceList[1].interaction.partnerIndex == 0 && molecules[0].bndpartner.size() == 1
            && molecules[1].bndpartner.size() == 1 && counters.nLoops == 1 && counters.nBoundPairs[0] == 1
            && counters.copyNumSpecies[kReactantA] == 9 && counters.copyNumSpecies[kReactantB] == 9
            && counters.copyNumSpecies[kProduct] == 1;
    }

    // Refused means untouched: no flag, no partner, no counter, no species move.
    bool untouched() const
    {
        return !molecules[0].interfaceList[0].isBound && !molecules[1].interfaceList[1].isBound
            && molecules[0].bndpartner.empty() && molecules[1].bndpartner.empty() && molecules[0].bndlist.empty()
            && molecules[1].bndlist.empty() && molecules[0].freelist.size() == 2 && molecules[1].freelist.size() == 2
            && counters.nLoops == 0 && counters.nBoundPairs[0] == 0
            && counters.copyNumSpecies[kReactantA] == 10 && counters.copyNumSpecies[kReactantB] == 10
            && counters.copyNumSpecies[kProduct] == 0;
    }
};

int failures = 0;

void expect(bool ok, const char* what)
{
    std::printf("  %-66s %s\n", what, ok ? "ok" : "FAILED");
    if (!ok)
        ++failures;
}

} // namespace

int main()
{
    // bindRadius 1 and bindRadSameCom 1.25 are exact in binary, so the
    // boundary case below compares 1.25 against exactly 1.25.
    const double sigma = 1.0;
    const double sameCom = 1.25;

    {
        Fixture f(0.9 * sigma, sigma, sameCom);
        expect(!f.membrane.isSphere(), "fixture uses the box geometry");
        f.associate_pair();
        expect(f.bonded(), "loop closure inside the limit (0.9 sigma) bonds");
    }
    {
        Fixture f(1.2 * sigma, sigma, sameCom);
        f.associate_pair();
        expect(f.bonded(), "loop closure just inside the limit (1.2 sigma) bonds");
    }
    {
        Fixture f(1.25 * sigma, sigma, sameCom);
        f.associate_pair();
        expect(f.untouched(), "loop closure exactly at the limit is refused, as at admission");
    }
    {
        Fixture f(5.2 * sigma, sigma, sameCom);
        f.associate_pair();
        expect(f.untouched(), "loop closure 5.2 sigma apart (the 6BNO case) is refused");
    }
    {
        Fixture f(27.8 * sigma, sigma, sameCom);
        f.associate_pair();
        expect(f.untouched(), "loop closure 27.8 sigma apart is refused");
    }
    {
        Fixture f(std::numeric_limits<double>::quiet_NaN(), sigma, sameCom);
        f.associate_pair();
        expect(f.untouched(), "loop closure at a NaN separation is refused");
    }

    std::printf("loop closure separation: %s\n", failures ? "FAILED" : "all cases OK");
    return failures ? 1 : 0;
}
