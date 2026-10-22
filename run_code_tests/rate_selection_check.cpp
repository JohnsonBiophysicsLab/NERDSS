// Proves best_matching_rate_compiled(), which selects from a reaction's
// compiled requirements, returns the RateMatch best_matching_rate_scan()
// returns -- the match count, the selected rate state and its ancillary count
// -- for every rate list of several parsed models, the cooperative clathrin
// ones included, over random interface states of both molecules.  So does
// best_matching_rate(), which hands single-state lists to the scan.
//
// The molecules are drawn to land near the rate states, not uniformly: each
// trial picks a rate state, gives both molecules the interfaces it asks for,
// then disturbs every interface with probability 1/4, so the lists' many-match,
// tie and no-match cases all occur.  Some molecules also break the
// relIndex == position invariant or lose interfaces, which sends the probe to
// hasAncillaryIface()'s linear scan; some calls pass reactant positions other
// than (0, 1) and (1, 0), which go to the scan whole.
//
// Synthetic lists then cover what no model has: more than 64 distinct
// requirements (left to the scan), a list not shaped [2][n], a list that grows
// after it was compiled, and a copied reaction, which must compile its own.
//
// It also prints, per model, what the compiled lists look like: rate states,
// interfaces read and distinct requirements, against the requirement tests the
// scan makes when it visits every state.
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

struct Totals {
    long long lists { 0 }, calls { 0 }, mismatches { 0 };
    long long noMatch { 0 }, oneMatch { 0 }, severalMatches { 0 }, tieKeptFirst { 0 };
    long long compiledCalls { 0 }, scanCalls { 0 }, fallbackProbes { 0 };
};

const char states[] = { '\0', 'A', 'B', 'P', 'U' };

// A molecule of template t with every interface's state and bound status drawn at random.
Molecule random_molecule(const MolTemplate& tmpl, int molType, std::mt19937_64& rng)
{
    Molecule mol;
    mol.molTypeIndex = molType;
    for (std::size_t i = 0; i < tmpl.interfaceList.size(); ++i) {
        const Interface& tIface = tmpl.interfaceList[i];
        Molecule::Iface iface;
        iface.relIndex = static_cast<int>(i);
        iface.molTypeIndex = molType;
        if (!tIface.stateList.empty()) {
            const int s { static_cast<int>(rng() % tIface.stateList.size()) };
            iface.stateIndex = s;
            iface.stateIden = tIface.stateList[s].iden;
            iface.index = tIface.stateList[s].index;
        }
        iface.isBound = (rng() % 2) == 0;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

// Gives mol what one requirement list asks of it, then disturbs it.
void steer(Molecule& mol, const std::vector<RxnIface>& reqs, std::mt19937_64& rng)
{
    for (const RxnIface& req : reqs) {
        for (auto& iface : mol.interfaceList) {
            if (iface.relIndex == req.relIfaceIndex) {
                iface.isBound = req.requiresInteraction;
                iface.stateIden = req.requiresState;
                if (rng() % 8 == 0)
                    iface.molTypeIndex = req.molTypeIndex;
            }
        }
    }
    for (auto& iface : mol.interfaceList) {
        if (rng() % 4 != 0)
            continue;
        if (rng() % 2 == 0)
            iface.isBound = !iface.isBound;
        else
            iface.stateIden = states[rng() % sizeof states];
    }
    // Break the relIndex == position invariant, or drop interfaces, now and then.
    if (mol.interfaceList.size() >= 2 && rng() % 16 == 0)
        std::swap(mol.interfaceList[0], mol.interfaceList[1 + rng() % (mol.interfaceList.size() - 1)]);
    if (!mol.interfaceList.empty() && rng() % 32 == 0)
        mol.interfaceList.resize(rng() % mol.interfaceList.size());
}

bool same(const RateMatch& a, const RateMatch& b)
{
    return a.matches == b.matches && a.bestRateIndex == b.bestRateIndex
        && a.mostAncillaryIfaces == b.mostAncillaryIfaces;
}

// Compares both selections on one pair, at all four reactant-position pairs.
void compare(const RxnBase& rxn, const Molecule& m1, const Molecule& m2, const char* where, Totals& totals)
{
    static const int positions[4][2] { { 0, 1 }, { 1, 0 }, { 0, 0 }, { 1, 1 } };
    for (const auto& pos : positions) {
        const RateMatch expected { best_matching_rate_scan(rxn.rateList, pos[0], pos[1], m1, m2) };
        const RateMatch got { best_matching_rate_compiled(rxn, pos[0], pos[1], m1, m2) };
        const RateMatch dispatched { best_matching_rate(rxn, pos[0], pos[1], m1, m2) };
        ++totals.calls;
        if (!same(expected, got) || !same(expected, dispatched)) {
            if (++totals.mismatches < 6)
                std::printf("  mismatch in %s at (%d, %d): scan (%d, %d, %zu), compiled (%d, %d, %zu), dispatched "
                            "(%d, %d, %zu)\n",
                    where, pos[0], pos[1], expected.matches, expected.bestRateIndex, expected.mostAncillaryIfaces,
                    got.matches, got.bestRateIndex, got.mostAncillaryIfaces, dispatched.matches,
                    dispatched.bestRateIndex, dispatched.mostAncillaryIfaces);
        }
        const bool compiled { rxn.compiledRates.usable && pos[0] != pos[1] };
        if (compiled)
            ++totals.compiledCalls;
        else
            ++totals.scanCalls;
        if (compiled) {
            for (const auto& probe : rxn.compiledRates.probes) {
                const Molecule& mol { (probe.side == pos[0]) ? m1 : m2 };
                if (probe.relIfaceIndex < 0 || static_cast<std::size_t>(probe.relIfaceIndex) >= mol.interfaceList.size()
                    || mol.interfaceList[probe.relIfaceIndex].relIndex != probe.relIfaceIndex)
                    ++totals.fallbackProbes;
            }
        }
        if (expected.matches == 0)
            ++totals.noMatch;
        else if (expected.matches == 1)
            ++totals.oneMatch;
        else
            ++totals.severalMatches;
        // A later state with as many ancillary interfaces as the selected one matched too.
        if (expected.matches > 1) {
            int equallySpecific { 0 };
            for (std::size_t k = 0; k < rxn.rateList.size(); ++k) {
                const auto& lists { rxn.rateList[k].otherIfaceLists };
                if (lists.size() == 2 && lists[0].size() + lists[1].size() == expected.mostAncillaryIfaces
                    && hasIntangibles(pos[0], pos[1], m1, m2, rxn.rateList[k]))
                    ++equallySpecific;
            }
            if (equallySpecific > 1)
                ++totals.tieKeptFirst;
        }
    }
}

int side_template(const RxnBase& rxn, int side, int fallback)
{
    for (const auto& rate : rxn.rateList)
        if (rate.otherIfaceLists.size() == 2 && !rate.otherIfaceLists[side].empty())
            return rate.otherIfaceLists[side][0].molTypeIndex;
    if (static_cast<std::size_t>(side) < rxn.reactantListNew.size())
        return rxn.reactantListNew[side].molTypeIndex;
    return fallback;
}

void check_list(const RxnBase& rxn, const std::vector<MolTemplate>& templates, const std::string& where,
    std::mt19937_64& rng, Totals& totals, std::map<std::string, long long>& shape)
{
    if (rxn.rateList.empty())
        return;
    // Unimolecular lists hold one requirement list per state; the bimolecular
    // selection never sees them, and the scan would read past them.
    for (const auto& rate : rxn.rateList)
        if (rate.otherIfaceLists.size() < 2)
            return;
    ++totals.lists;
    const int t0 { side_template(rxn, 0, 0) };
    const int t1 { side_template(rxn, 1, t0) };
    for (int trial = 0; trial < 3000; ++trial) {
        // Mostly the templates the requirements name; sometimes any template.
        const int a { (rng() % 16 == 0) ? static_cast<int>(rng() % templates.size()) : t0 };
        const int b { (rng() % 16 == 0) ? static_cast<int>(rng() % templates.size()) : t1 };
        Molecule m1 { random_molecule(templates[a], a, rng) };
        Molecule m2 { random_molecule(templates[b], b, rng) };
        if (rng() % 8 != 0) {
            const auto& target { rxn.rateList[rng() % rxn.rateList.size()] };
            if (target.otherIfaceLists.size() == 2) {
                steer(m1, target.otherIfaceLists[0], rng);
                steer(m2, target.otherIfaceLists[1], rng);
            }
        }
        compare(rxn, m1, m2, where.c_str(), totals);
        compare(rxn, m2, m1, where.c_str(), totals);
    }

    // What selection costs: the scan's requirement tests if it visited every
    // state to the end, against the compiled form's interface reads.
    long long scanTests { 0 };
    for (const auto& rate : rxn.rateList)
        for (const auto& list : rate.otherIfaceLists)
            scanTests += static_cast<long long>(list.size());
    const auto& c { rxn.compiledRates };
    char key[128];
    std::snprintf(key, sizeof key, "%2zu states, %2lld tests -> %2zu ifaces, %2zu distinct%s", rxn.rateList.size(),
        scanTests, c.probes.size(), c.tests.size(), c.usable ? "" : " (scan)");
    ++shape[key];
}

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
    std::map<std::string, long long> shape;
    for (const auto& rxn : forwardRxns)
        check_list(rxn, templates, model.dir + " forward " + rxn.rxnLabel, rng, here, shape);
    for (const auto& rxn : backRxns)
        check_list(rxn, templates, model.dir + " back " + rxn.rxnLabel, rng, here, shape);

    std::printf("  %-55s %3lld lists %8lld calls, %lld mismatches\n", model.dir.c_str(), here.lists, here.calls,
        here.mismatches);
    for (const auto& entry : shape)
        if (entry.first.compare(0, 2, " 1") != 0) // single-state lists have nothing to select
            std::printf("      %3lld x %s\n", entry.second, entry.first.c_str());
    totals.lists += here.lists;
    totals.calls += here.calls;
    totals.mismatches += here.mismatches;
    totals.noMatch += here.noMatch;
    totals.oneMatch += here.oneMatch;
    totals.severalMatches += here.severalMatches;
    totals.tieKeptFirst += here.tieKeptFirst;
    totals.compiledCalls += here.compiledCalls;
    totals.scanCalls += here.scanCalls;
    totals.fallbackProbes += here.fallbackProbes;
}

// A random requirement on one of n interfaces of molecule type `type`.
RxnIface random_requirement(int type, int n, std::mt19937_64& rng)
{
    RxnIface req;
    req.molTypeIndex = type;
    req.relIfaceIndex = static_cast<int>(rng() % (n + 1)) - (rng() % 32 == 0 ? n + 1 : 0); // a few out of range
    req.requiresState = states[rng() % sizeof states];
    req.requiresInteraction = rng() % 2 == 0;
    return req;
}

MolTemplate synthetic_template(int n)
{
    MolTemplate tmpl;
    for (int i = 0; i < n; ++i) {
        Interface iface;
        for (char iden : states) {
            Interface::State state;
            state.iden = iden;
            iface.stateList.push_back(state);
        }
        tmpl.interfaceList.push_back(iface);
    }
    return tmpl;
}

// Random lists over synthetic templates: small ones that compile, large ones
// that do not, and the structural cases the parser never produces.
void check_synthetic(std::mt19937_64& rng, Totals& totals)
{
    const std::vector<MolTemplate> templates { synthetic_template(4), synthetic_template(40) };
    Totals here;
    std::map<std::string, long long> shape;
    for (int list = 0; list < 300; ++list) {
        const int n { list % 3 == 2 ? 40 : 4 };
        const int tmpl { n == 40 ? 1 : 0 };
        ForwardRxn rxn;
        const int stateCount { 1 + static_cast<int>(rng() % 9) };
        for (int k = 0; k < stateCount; ++k) {
            RxnBase::RateState rate;
            rate.otherIfaceLists.resize(2);
            for (int side = 0; side < 2; ++side) {
                const int count { static_cast<int>(rng() % (n == 40 ? 24 : 4)) };
                for (int q = 0; q < count; ++q)
                    rate.otherIfaceLists[side].push_back(random_requirement(tmpl, n, rng));
            }
            rxn.rateList.push_back(rate);
        }
        check_list(rxn, templates, "synthetic", rng, here, shape);

        // A copy compiles its own, and a list that grows is recompiled.
        ForwardRxn copy { rxn };
        if (copy.compiledRates.built)
            ++here.mismatches;
        RxnBase::RateState extra;
        extra.otherIfaceLists.resize(2);
        extra.otherIfaceLists[0].push_back(random_requirement(tmpl, n, rng));
        copy.rateList.push_back(extra);
        check_list(copy, templates, "synthetic, grown", rng, here, shape);
        if (copy.compiledRates.stateCount != copy.rateList.size())
            ++here.mismatches;
    }
    // Not shaped [2][n]: selection is left to the scan.  (The scan reads
    // otherIfaceLists[1], so the short list has a third, unread, list instead.)
    ForwardRxn odd;
    odd.rateList.resize(2);
    odd.rateList[0].otherIfaceLists.resize(2);
    odd.rateList[1].otherIfaceLists.resize(3);
    odd.rateList[1].otherIfaceLists[0].push_back(random_requirement(0, 4, rng));
    check_list(odd, templates, "synthetic, [3][n]", rng, here, shape);
    if (odd.compiledRates.usable)
        ++here.mismatches;

    std::printf("  %-55s %3lld lists %8lld calls, %lld mismatches\n", "synthetic", here.lists, here.calls,
        here.mismatches);
    long long unusable { 0 };
    for (const auto& entry : shape)
        if (entry.first.find("(scan)") != std::string::npos)
            unusable += entry.second;
    std::printf("      %lld of them left to the scan\n", unusable);
    if (unusable == 0)
        ++here.mismatches;
    totals.lists += here.lists;
    totals.calls += here.calls;
    totals.mismatches += here.mismatches;
    totals.noMatch += here.noMatch;
    totals.oneMatch += here.oneMatch;
    totals.severalMatches += here.severalMatches;
    totals.tieKeptFirst += here.tieKeptFirst;
    totals.compiledCalls += here.compiledCalls;
    totals.scanCalls += here.scanCalls;
    totals.fallbackProbes += here.fallbackProbes;
}

} // namespace

int main()
{
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 20261022);
    std::mt19937_64 rng(20261022);

    char cwd[4096];
    if (getcwd(cwd, sizeof cwd) == nullptr)
        return 1;
    const std::string root { cwd };

    const std::vector<Model> models {
        { "sample_inputs/clathrin_pioneer/ap2clafchoeps", "parms.inp" },
        { "sample_inputs/clathrin_pioneer/ap2clafchoeps4site", "parms.inp" },
        { "sample_inputs/clathrin_pioneer/test", "parmscross.inp" },
        { "sample_inputs/clathrin_invitro_invivo/invivo", "parms.inp" },
        { "sample_inputs/clathrin_coat/cooperative_flat_clat.dir", "parms.inp" },
        { "sample_inputs/clathrin_coat/flat_clat-ap2-pip2.dir", "parms.inp" },
        { "sample_inputs/VALIDATE_SUITE/clathrin", "parms_clath_kon1uM.inp" },
        { "sample_inputs/VALIDATE_SUITE/trimer", "parm_hetTri.inp" },
        { "sample_inputs/VALIDATE_SUITE/hetTrimer", "parm_autodiff_hetTri.inp" },
        { "sample_inputs/VALIDATE_SUITE/implicit_lipid", "parms.inp" },
        { "sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible", "uni_state_rev.inp" },
        { "sample_inputs/gagsphere", "gagOriginalModelSolution.inp" },
        { "sample_inputs/auto_phos", "autophos_D10.inp" },
        { "benchmarks/nerdss_optimized/regression/bimol_state_change", "parms.inp" },
    };

    Totals totals;
    for (const Model& model : models)
        check_model(root, model, rng, totals);
    check_synthetic(rng, totals);
    gsl_rng_free(r);

    std::printf("rate_selection_check: %lld lists, %lld selections (%lld compiled, %lld scanned): %lld no match, %lld "
                "one, %lld several (%lld with a tie for most specific), %lld probes fell back to the interface scan, "
                "%lld mismatches\n",
        totals.lists, totals.calls, totals.compiledCalls, totals.scanCalls, totals.noMatch, totals.oneMatch,
        totals.severalMatches, totals.tieKeptFirst, totals.fallbackProbes, totals.mismatches);
    // Every case must have occurred, or the check proves less than it says.
    if (totals.mismatches != 0 || totals.noMatch == 0 || totals.oneMatch == 0 || totals.severalMatches == 0
        || totals.tieKeptFirst == 0 || totals.compiledCalls == 0 || totals.scanCalls == 0
        || totals.fallbackProbes == 0) {
        std::printf("rate_selection_check: FAILED\n");
        return 1;
    }
    std::printf("rate_selection_check: all passed\n");
    return 0;
}
