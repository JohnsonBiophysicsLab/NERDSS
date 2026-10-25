// Proves find_which_reaction_for_partner() returns what find_which_reaction()
// returns, and counts totMatches the same way, for every (state, partner) pair a
// model can offer the pair search, over random states of every other interface
// of both molecules.
//
// The model runs compare the two on every lookup they make (both suites and the
// production models, ~1.5e9 lookups, no difference), but they only reach two of
// the per-partner lookup's paths: a rate that cannot depend on the molecules,
// and an ordinary forward rate tested per pair.  The back state change, the
// symmetric-reaction swap and "no rate state matches" need states the models
// never put next to each other, so this check builds them: each model's input
// is parsed as the simulator parses it, and the molecules' interfaces are set at
// random.  Two of the models are written here to reach those paths.
//
// It also proves rate selection from Molecule::ifaceStateBits: for every rate
// list a lookup tests, in both reactant orders, the selection made from the two
// molecules' words alone must equal best_matching_rate_scan()'s from their
// interfaces, and with one molecule's word invalidated the compiled selection
// must fall back to the interfaces and still agree.
#include "classes/class_Membrane.hpp"
#include "classes/iface_state_bits.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main()
unsigned long totMatches = 0;

namespace {

struct Model {
    std::string dir; // relative to the repository root, or a scratch directory
    std::string parm;
};

// A homodimer with an ancillary requirement on one side only: for a pair where
// just the second molecule has b bound, only the swapped positions match, and
// where neither has it bound, no rate state matches.
const char* symmetricParm = R"(start parameters
    nItr = 10
    timeStep = 0.1
end parameters
start boundaries
    WaterBox = [100,100,100]
end boundaries
start molecules
    S : 10
end molecules
start reactions
    S(b) + S(b) <-> S(b!1).S(b!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [M_PI,M_PI,nan,nan,0]

    S(a,b!*) + S(a) <-> S(a!1,b!*).S(a!1)
    onRate3Dka = 2.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [M_PI,M_PI,nan,nan,0]
end reactions
)";
const char* symmetricMol = R"(Name = S
checkOverlap = true
D = [10.0, 10.0, 10.0]
Dr = [0.1, 0.1, 0.1]
COM 0.0 0.0 0.0
a 1.0 0.0 0.0
b -1.0 0.0 0.0
bonds = 2
com a
com b
)";

// A bimolecular state change that asks for c to be bound on A, in both
// directions, so the back state change's rates are tested per pair and fail
// (the lookup moves on) whenever c is free.  A(c) + A(c) lets c be bound.
const char* stateChangeParm = R"(start parameters
    nItr = 10
    timeStep = 0.1
end parameters
start boundaries
    WaterBox = [100,100,100]
end boundaries
start molecules
    K : 10
    A : 10 (a~X)
end molecules
start reactions
    K(k) + A(a~X,c!*) <-> K(k) + A(a~Y,c!*)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    A(c) + A(c) <-> A(c!1).A(c!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [M_PI,M_PI,nan,nan,0]
end reactions
)";
const char* stateChangeMolA = R"(Name = A
checkOverlap = true
D = [10.0, 10.0, 10.0]
Dr = [0.1, 0.1, 0.1]
COM 0.0 0.0 0.0
a 1.0 0.0 0.0
c -1.0 0.0 0.0
bonds = 2
com a
com c
state = a~X~Y
)";
const char* stateChangeMolK = R"(Name = K
checkOverlap = true
D = [10.0, 10.0, 10.0]
Dr = [0.1, 0.1, 0.1]
COM 0.0 0.0 0.0
k 1.0 0.0 0.0
bonds = 1
com k
)";

void write_file(const std::string& path, const char* text)
{
    std::ofstream out(path);
    out << text;
}

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

void set_state(Molecule& mol, const MolTemplate& tmpl, int relIface, int stateIndex)
{
    Molecule::Iface& iface = mol.interfaceList[relIface];
    iface.stateIndex = stateIndex;
    iface.stateIden = tmpl.interfaceList[relIface].stateList[stateIndex].iden;
    iface.index = tmpl.interfaceList[relIface].stateList[stateIndex].index;
    iface.isBound = false;
}

struct Totals {
    long long pairs { 0 }, lookups { 0 }, mismatches { 0 }, found { 0 }, none { 0 }, backStateChange { 0 };
    // Which path the per-partner lookup took, read back from its candidate list.
    long long fixedRate { 0 }, testedRate { 0 }, backSkipped { 0 }, swapped { 0 };
    // Selections from the state words: made, made from a list of several rate
    // states, made with several states matching, declined for an invalid word.
    long long words { 0 }, wordsMultiState { 0 }, wordsMultiMatch { 0 }, wordsDeclined { 0 }, wordMismatches { 0 };
};

bool same(const RateMatch& a, const RateMatch& b)
{
    return a.matches == b.matches && a.bestRateIndex == b.bestRateIndex
        && a.mostAncillaryIfaces == b.mostAncillaryIfaces;
}

// The selection from the words, and from the interfaces with m1's word
// invalidated, against the scan's.
void compare_words(const RxnBase& rxn, int reactIndex1, int reactIndex2, const Molecule& m1, const Molecule& m2,
    const std::string& where, Totals& totals)
{
    const RateMatch scan { best_matching_rate_scan(rxn.rateList, reactIndex1, reactIndex2, m1, m2) };
    RateMatch words;
    if (best_matching_rate_from_state_bits(rxn, reactIndex1, reactIndex2, m1, m2, words)) {
        ++totals.words;
        if (rxn.rateList.size() > 1)
            ++totals.wordsMultiState;
        if (words.matches > 1)
            ++totals.wordsMultiMatch;
        if (!same(words, scan) && ++totals.wordMismatches < 4)
            std::printf("  word mismatch in %s: rxn %d (%d, %d): words (%d, %d), scan (%d, %d)\n", where.c_str(),
                rxn.absRxnIndex, reactIndex1, reactIndex2, words.matches, words.bestRateIndex, scan.matches,
                scan.bestRateIndex);
    }

    Molecule stale { m1 };
    stale.ifaceStateBits = iface_state_bits::invalid;
    RateMatch declined;
    if (!best_matching_rate_from_state_bits(rxn, reactIndex1, reactIndex2, stale, m2, declined))
        ++totals.wordsDeclined;
    if (!same(best_matching_rate(rxn, reactIndex1, reactIndex2, stale, m2), scan) && ++totals.wordMismatches < 4)
        std::printf("  fallback mismatch in %s: rxn %d (%d, %d)\n", where.c_str(), rxn.absRxnIndex, reactIndex1,
            reactIndex2);
}

// Walks the candidates the way find_which_reaction_for_partner() does, to tell
// which of its paths produced the answer.  Only for the tally above.
void classify(const std::vector<Interface::State::RxnCandidate>& candidates, const Molecule& m1, const Molecule& m2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns, Totals& totals)
{
    for (const auto& c : candidates) {
        const ForwardRxn& rxn = forwardRxns[c.rxnItr];
        if (c.fixedRateIndex != -1) {
            ++totals.fixedRate;
            return;
        }
        ++totals.testedRate;
        if (c.backStateChange) {
            if (best_matching_rate(backRxns[rxn.conjBackRxnIndex], c.reactIndex1, c.reactIndex2, m1, m2).matches
                == 0) {
                ++totals.backSkipped;
                continue;
            }
            return;
        }
        if (best_matching_rate(rxn, c.reactIndex1, c.reactIndex2, m1, m2).matches == 0 && c.trySwap
            && best_matching_rate(rxn, c.reactIndex2, c.reactIndex1, m1, m2).matches != 0)
            ++totals.swapped;
        return;
    }
}

void check_model(const std::string& root, const Model& model, std::mt19937_64& rng, Totals& totals)
{
    const std::string dir { model.dir[0] == '/' ? model.dir : root + "/" + model.dir };
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
    iface_state_bits::build_layouts(templates);

    Totals here;
    for (std::size_t t1 = 0; t1 < templates.size(); ++t1) {
        if (templates[t1].isImplicitLipid)
            continue; // the pair search hands those to find_which_reaction()
        for (std::size_t i1 = 0; i1 < templates[t1].interfaceList.size(); ++i1) {
            for (std::size_t s1 = 0; s1 < templates[t1].interfaceList[i1].stateList.size(); ++s1) {
                const Interface::State& state { templates[t1].interfaceList[i1].stateList[s1] };
                for (std::size_t p = 0; p < state.rxnPartners.size(); ++p) {
                    // The partner is a state index: find its template, interface and state.
                    int t2 { -1 }, i2 { -1 }, s2 { -1 };
                    for (std::size_t t = 0; t < templates.size() && t2 == -1; ++t)
                        for (std::size_t i = 0; i < templates[t].interfaceList.size() && t2 == -1; ++i)
                            for (std::size_t s = 0; s < templates[t].interfaceList[i].stateList.size(); ++s)
                                if (templates[t].interfaceList[i].stateList[s].index
                                    == static_cast<int>(state.rxnPartners[p])) {
                                    t2 = static_cast<int>(t);
                                    i2 = static_cast<int>(i);
                                    s2 = static_cast<int>(s);
                                    break;
                                }
                    if (t2 == -1 || templates[t2].isImplicitLipid)
                        continue;
                    ++here.pairs;
                    for (int trial = 0; trial < 2000; ++trial) {
                        Molecule m1 { random_molecule(templates[t1], static_cast<int>(t1), rng) };
                        Molecule m2 { random_molecule(templates[t2], t2, rng) };
                        set_state(m1, templates[t1], static_cast<int>(i1), static_cast<int>(s1));
                        set_state(m2, templates[t2], i2, s2);
                        m1.refresh_iface_state_bits();
                        m2.refresh_iface_state_bits();

                        int rxnA { -1 }, rateA { -1 }, rxnB { -1 }, rateB { -1 };
                        bool scA { false }, scB { false };
                        const unsigned long before { totMatches };
                        find_which_reaction(static_cast<int>(i1), i2, rxnA, rateA, scA, state, m1, m2, forwardRxns,
                            backRxns);
                        const unsigned long countA { totMatches - before };
                        find_which_reaction_for_partner(
                            static_cast<int>(i1), p, rxnB, rateB, scB, state, m1, m2, forwardRxns, backRxns);
                        const unsigned long countB { totMatches - before - countA };

                        ++here.lookups;
                        if (rxnA != rxnB || rateA != rateB || scA != scB || countA != countB) {
                            if (++here.mismatches < 4)
                                std::printf("  mismatch in %s: old (%d, %d, %d, %lu), new (%d, %d, %d, %lu)\n",
                                    model.dir.c_str(), rxnA, rateA, scA, countA, rxnB, rateB, scB, countB);
                        }
                        if (rxnA != -1 && rateA != -1)
                            ++here.found;
                        else
                            ++here.none;
                        if (scA)
                            ++here.backStateChange;
                        classify(state.rxnCandidatesByPartner[p], m1, m2, forwardRxns, backRxns, here);
                        for (const auto& c : state.rxnCandidatesByPartner[p]) {
                            const ForwardRxn& rxn = forwardRxns[c.rxnItr];
                            const RxnBase& tested = c.backStateChange
                                ? static_cast<const RxnBase&>(backRxns[rxn.conjBackRxnIndex])
                                : static_cast<const RxnBase&>(rxn);
                            compare_words(tested, c.reactIndex1, c.reactIndex2, m1, m2, model.dir, here);
                            compare_words(tested, c.reactIndex2, c.reactIndex1, m1, m2, model.dir, here);
                        }
                    }
                }
            }
        }
    }
    std::printf("  %-62s %3lld partners %7lld lookups: %lld found (%lld back state change), %lld none, %lld mismatches\n"
                "  %-62s paths: fixed rate %lld, tested rate %lld, back state change skipped %lld, swapped %lld\n"
                "  %-62s state words: %lld selections (%lld from several states, %lld with several matches), "
                "%lld declined, %lld mismatches\n",
        model.dir.c_str(), here.pairs, here.lookups, here.found, here.backStateChange, here.none, here.mismatches, "",
        here.fixedRate, here.testedRate, here.backSkipped, here.swapped, "", here.words, here.wordsMultiState,
        here.wordsMultiMatch, here.wordsDeclined, here.wordMismatches);
    totals.words += here.words;
    totals.wordsMultiState += here.wordsMultiState;
    totals.wordsMultiMatch += here.wordsMultiMatch;
    totals.wordsDeclined += here.wordsDeclined;
    totals.wordMismatches += here.wordMismatches;
    totals.fixedRate += here.fixedRate;
    totals.testedRate += here.testedRate;
    totals.backSkipped += here.backSkipped;
    totals.swapped += here.swapped;
    totals.pairs += here.pairs;
    totals.lookups += here.lookups;
    totals.mismatches += here.mismatches;
    totals.found += here.found;
    totals.none += here.none;
    totals.backStateChange += here.backStateChange;
}

} // namespace

int main()
{
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 20261020);
    std::mt19937_64 rng(20261020);

    char cwd[4096];
    if (getcwd(cwd, sizeof cwd) == nullptr)
        return 1;
    const std::string root { cwd };

    char scratchTemplate[] = "/tmp/reaction_lookup_check.XXXXXX";
    const char* scratch { mkdtemp(scratchTemplate) };
    if (scratch == nullptr)
        return 1;
    const std::string symmetricDir { std::string(scratch) + "/symmetric" };
    const std::string stateChangeDir { std::string(scratch) + "/state_change" };
    mkdir(symmetricDir.c_str(), 0755);
    mkdir(stateChangeDir.c_str(), 0755);
    write_file(symmetricDir + "/parms.inp", symmetricParm);
    write_file(symmetricDir + "/S.mol", symmetricMol);
    write_file(stateChangeDir + "/parms.inp", stateChangeParm);
    write_file(stateChangeDir + "/A.mol", stateChangeMolA);
    write_file(stateChangeDir + "/K.mol", stateChangeMolK);

    const std::vector<Model> models {
        { symmetricDir, "parms.inp" },
        { stateChangeDir, "parms.inp" },
        { "benchmarks/nerdss_optimized/regression/bimol_state_change", "parms.inp" },
        { "sample_inputs/VALIDATE_SUITE/trimer", "parm_hetTri.inp" },
        { "sample_inputs/VALIDATE_SUITE/hetTrimer", "parm_autodiff_hetTri.inp" },
        { "sample_inputs/VALIDATE_SUITE/clathrin", "parms_clath_kon1uM.inp" },
        { "sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible", "uni_state_rev.inp" },
        { "sample_inputs/clathrin_invitro_invivo/invivo", "parms.inp" },
        { "sample_inputs/clathrin_pioneer/test", "parmscross.inp" },
        { "sample_inputs/clathrin_pioneer/ap2clafchoeps", "parms.inp" },
        { "sample_inputs/gagLatticeRemodeling/1-initialAssemblyInSolution", "parm.inp" },
        { "sample_inputs/auto_phos", "autophos_D10.inp" },
        { "sample_inputs/enzyme", "parms_clat_enzyme.inp" },
    };

    Totals totals;
    for (const Model& model : models)
        check_model(root, model, rng, totals);
    gsl_rng_free(r);

    std::printf("reaction_lookup_check: %lld partners, %lld lookups: %lld found (%lld back state change), %lld none, "
                "%lld mismatches\n"
                "reaction_lookup_check: paths: fixed rate %lld, tested rate %lld, back state change skipped %lld, "
                "swapped %lld\n",
        totals.pairs, totals.lookups, totals.found, totals.backStateChange, totals.none, totals.mismatches,
        totals.fixedRate, totals.testedRate, totals.backSkipped, totals.swapped);
    std::printf("reaction_lookup_check: state words: %lld selections (%lld from several states, %lld with several "
                "matches), %lld declined, %lld mismatches\n",
        totals.words, totals.wordsMultiState, totals.wordsMultiMatch, totals.wordsDeclined, totals.wordMismatches);
    // Every path must have been taken, or the check proves less than it says.
    if (totals.mismatches != 0 || totals.backStateChange == 0 || totals.none == 0 || totals.fixedRate == 0
        || totals.testedRate == 0 || totals.swapped == 0 || totals.backSkipped == 0 || totals.wordMismatches != 0
        || totals.words == 0 || totals.wordsMultiState == 0 || totals.wordsMultiMatch == 0
        || totals.wordsDeclined == 0) {
        std::printf("reaction_lookup_check: FAILED\n");
        return 1;
    }
    std::printf("reaction_lookup_check: all passed\n");
    return 0;
}
