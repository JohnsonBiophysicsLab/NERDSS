// Proves find_reaction_rate_state(), which keeps each bond kind's reactant
// positions -- and, when the rate list cannot depend on the molecules, its
// answer -- in BackRxn::ratePositions, selects what the reactant scan it
// replaced selects, for both orders of every bond a back reaction breaks, over
// random states of every other interface of both molecules.  Both
// check_dissociation() and check_dissociation_implicitlipid() ask it once per
// bound pair per step.
#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <map>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main()
unsigned long totMatches = 0;

namespace {

struct Model {
    std::string dir; // relative to the repository root
    std::string parm;
};

// A molecule of template t with interface states drawn at random.
Molecule random_molecule(const MolTemplate& tmpl, int molType, std::mt19937_64& rng)
{
    Molecule mol;
    mol.molTypeIndex = molType;
    mol.isImplicitLipid = tmpl.isImplicitLipid;
    for (std::size_t i = 0; i < tmpl.interfaceList.size(); ++i) {
        const Interface& tIface = tmpl.interfaceList[i];
        Molecule::Iface iface;
        iface.relIndex = static_cast<int>(i);
        iface.molTypeIndex = molType;
        const int s { static_cast<int>(rng() % tIface.stateList.size()) };
        iface.stateIndex = s;
        iface.stateIden = tIface.stateList[s].iden;
        iface.index = tIface.stateList[s].index;
        iface.isBound = (rng() % 2) == 0;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

// find_reaction_rate_state() as it was before the bond-kind entries.
int reference_rate_state(int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1, const Molecule& reactMol2,
    const BackRxn& backRxn)
{
    int reactIndex1 { -1 };
    int reactIndex2 { -1 };
    for (std::size_t reactItr { 0 }; reactItr < backRxn.reactantListNew.size(); ++reactItr) {
        const RxnIface& oneReactant = backRxn.reactantListNew[reactItr];
        if (reactMol1.molTypeIndex == oneReactant.molTypeIndex && relIfaceIndex1 == oneReactant.relIfaceIndex) {
            if (reactIndex1 == -1) {
                reactIndex1 = static_cast<int>(reactItr);
                continue;
            }
        }
        if (reactMol2.molTypeIndex == oneReactant.molTypeIndex && relIfaceIndex2 == oneReactant.relIfaceIndex) {
            if (reactIndex2 == -1) {
                reactIndex2 = static_cast<int>(reactItr);
                continue;
            }
        }
    }
    if (reactIndex1 == -1 || reactIndex2 == -1)
        return -2; // the simulator exits here; never asked below
    const RateMatch match { best_matching_rate(backRxn, reactIndex1, reactIndex2, reactMol1, reactMol2) };
    return match.matches == 0 ? -1 : match.bestRateIndex;
}

struct Totals {
    long long bondKinds { 0 }, fixed { 0 }, tested { 0 }, selections { 0 }, none { 0 }, mismatches { 0 };
};

void check_model(const std::string& root, const Model& model, std::mt19937_64& rng, Totals& totals)
{
    const std::string dir { root + "/" + model.dir };
    if (chdir(dir.c_str()) != 0) {
        std::printf("  cannot enter %s\n", dir.c_str());
        std::exit(1);
    }
    Parameters params;
    std::map<std::string, int> observables;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    std::vector<TransmissionRxn> transmissionRxns;
    std::vector<MolTemplate> templates;
    Membrane membrane;
    std::string parm { model.parm };
    // The parser reports on stdout; keep the check's own output readable.
    std::fflush(stdout);
    const int savedStdout { dup(1) };
    const int devNull { open("/dev/null", O_WRONLY) };
    dup2(devNull, 1);
    parse_input(parm, params, observables, forwardRxns, backRxns, createDestructRxns, transmissionRxns, templates,
        membrane);
    std::fflush(stdout);
    dup2(savedStdout, 1);
    close(devNull);
    close(savedStdout);

    Totals here;
    for (const BackRxn& backRxn : backRxns) {
        if (backRxn.reactantListNew.size() != 2)
            continue;
        for (int order = 0; order < 2; ++order) {
            const RxnIface& iface1 { backRxn.reactantListNew[order] };
            const RxnIface& iface2 { backRxn.reactantListNew[1 - order] };
            ++here.bondKinds;
            for (int trial = 0; trial < 2000; ++trial) {
                Molecule m1 { random_molecule(templates[iface1.molTypeIndex], iface1.molTypeIndex, rng) };
                Molecule m2 { random_molecule(templates[iface2.molTypeIndex], iface2.molTypeIndex, rng) };
                m1.interfaceList[iface1.relIfaceIndex].isBound = true;
                m2.interfaceList[iface2.relIfaceIndex].isBound = true;
                const int expected { reference_rate_state(iface1.relIfaceIndex, iface2.relIfaceIndex, m1, m2, backRxn) };
                const int got { find_reaction_rate_state(
                    0, iface1.relIfaceIndex, iface2.relIfaceIndex, m1, m2, backRxn, templates) };
                ++here.selections;
                if (expected != got) {
                    if (++here.mismatches < 4)
                        std::printf("  mismatch in %s: old %d, new %d\n", model.dir.c_str(), expected, got);
                }
                if (got == -1)
                    ++here.none;
            }
            for (const BackRxn::RatePositions& entry : backRxn.ratePositions.entries) {
                if (entry.molType1 == iface1.molTypeIndex && entry.relIface1 == iface1.relIfaceIndex
                    && entry.molType2 == iface2.molTypeIndex && entry.relIface2 == iface2.relIfaceIndex) {
                    if (entry.fixedRateIndex != -1)
                        ++here.fixed;
                    else
                        ++here.tested;
                }
            }
        }
    }
    std::printf("  %-62s %3lld bond kinds (%lld fixed, %lld tested), %7lld selections, %lld none, %lld mismatches\n",
        model.dir.c_str(), here.bondKinds, here.fixed, here.tested, here.selections, here.none, here.mismatches);
    totals.bondKinds += here.bondKinds;
    totals.fixed += here.fixed;
    totals.tested += here.tested;
    totals.selections += here.selections;
    totals.none += here.none;
    totals.mismatches += here.mismatches;
}

} // namespace

int main()
{
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 20261024);
    std::mt19937_64 rng(20261024);

    char cwd[4096];
    if (getcwd(cwd, sizeof cwd) == nullptr)
        return 1;
    const std::string root { cwd };

    const std::vector<Model> models {
        { "sample_inputs/VALIDATE_SUITE/implicit_lipid", "parms.inp" },
        { "sample_inputs/VALIDATE_SUITE/sphere", "parms_sphere.inp" },
        { "sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL", "parms.inp" },
        { "sample_inputs/VALIDATE_SUITE/trimer", "parm_hetTri.inp" },
        { "sample_inputs/VALIDATE_SUITE/hetTrimer", "parm_autodiff_hetTri.inp" },
        { "sample_inputs/VALIDATE_SUITE/clathrin", "parms_clath_kon1uM.inp" },
        { "sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible", "uni_state_rev.inp" },
        { "benchmarks/nerdss_optimized/regression/bimol_state_change", "parms.inp" },
        { "sample_inputs/clathrin_invitro_invivo/invitro", "parms.inp" },
        { "sample_inputs/clathrin_invitro_invivo/invivo", "parms.inp" },
        { "sample_inputs/clathrin_pioneer/ffl", "parmsffl.inp" },
        { "sample_inputs/clathrin_pioneer/test", "parmscross.inp" },
        { "sample_inputs/clathrin_pioneer/ap2clafchoeps", "parms.inp" },
        { "sample_inputs/testAdd/box", "parms.inp" },
        { "sample_inputs/gagsphere", "parms.inp" },
        { "sample_inputs/gagLatticeRemodeling/1-initialAssemblyInSolution", "parm.inp" },
        { "sample_inputs/auto_phos", "autophos_D10.inp" },
        { "sample_inputs/enzyme", "parms_clat_enzyme.inp" },
    };

    Totals totals;
    for (const Model& model : models)
        check_model(root, model, rng, totals);
    gsl_rng_free(r);

    std::printf("rate_state_cache_check: %lld bond kinds (%lld fixed, %lld tested), %lld selections, %lld none, "
                "%lld mismatches\n",
        totals.bondKinds, totals.fixed, totals.tested, totals.selections, totals.none, totals.mismatches);
    // Every path must have been taken, or the check proves less than it says.
    if (totals.mismatches != 0 || totals.fixed == 0 || totals.tested == 0 || totals.none == 0
        || totals.none == totals.selections) {
        std::printf("rate_state_cache_check: FAILED\n");
        return 1;
    }
    std::printf("rate_state_cache_check: all passed\n");
    return 0;
}
