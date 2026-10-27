// Proves find_which_reaction_for_partner(), called with the implicit lipid as
// the second molecule, returns what check_implicit_reactions() used to get from
// find_which_reaction() with the lipid's interface index moved by the lipid
// state's offset, and counts totMatches the same way, for every (state,
// implicit-lipid partner) pair a model offers, over random states of every
// other interface of the protein and of the lipid.  It also proves the reaction
// found is one with the lipid in the state the pair check is testing.
//
// No sample model gives the lipid more than one state, so the lipid offset is
// always 0 there, and no sample has a rate list that reads the protein's other
// interfaces.  The model written here does: a lipid with three states; B, whose
// rate depends on whether its a is bound (two rate states, so the compiled
// selection); C, whose only rate state needs a bound, so a free a matches
// nothing; and K, which changes the lipid from U to P and back.  The back
// state change starts from K and IL~P, so K's state must list IL~P as a
// partner and the lookup must return the reaction as a back state change.
// The parser used to file that partner only when the facilitator's own state
// changed too, which no reaction that parses as one state change does, so
// the pair check never tested K against IL~P and the lipid never went from P
// back to U.  The check fails unless some lookup returns a back state change
// (only K's can), with the lipid in the forward product's state.
//
// find_which_reaction() used to match the lipid by molecule type, not by
// state, so a protein was handed its first reaction with the lipid whatever
// lipid state the pair check was testing: A(m) + IL(m~P) and A(m,a!*) +
// IL(m~Q) below were never selected, and A bound every lipid state with the
// rates of A(m) + IL(m~U).  The state check below fails on that lookup: A's m
// must reach the lipid through three different reactions, one per state.
#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <random>
#include <set>
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

const char* lipidStatesParm = R"(start parameters
    nItr = 10
    timeStep = 0.1
end parameters
start boundaries
    WaterBox = [100,100,100]
    implicitLipid = true
end boundaries
start molecules
    IL : 100
    A : 10
    B : 10
    C : 10
    K : 10 (k~X)
end molecules
start reactions
    A(m) + IL(m~U) <-> A(m!1).IL(m!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    A(m,a) + IL(m~P) <-> A(m!1,a).IL(m!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    A(m,a!*) + IL(m~P) <-> A(m!1,a!*).IL(m!1)
    onRate3Dka = 2.0
    offRatekb = 0.5
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    B(m,a) + IL(m~U) <-> B(m!1,a).IL(m!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    B(m,a!*) + IL(m~U) <-> B(m!1,a!*).IL(m!1)
    onRate3Dka = 2.0
    offRatekb = 0.5
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    C(m,a!*) + IL(m~P) <-> C(m!1,a!*).IL(m!1)
    onRate3Dka = 2.0
    offRatekb = 0.5
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    B(a) + C(a) <-> B(a!1).C(a!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [M_PI,M_PI,nan,nan,0]

    A(a) + A(a) <-> A(a!1).A(a!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [M_PI,M_PI,nan,nan,0]

    A(m,a!*) + IL(m~Q) <-> A(m!1,a!*).IL(m!1)
    onRate3Dka = 3.0
    offRatekb = 0.5
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    K(k~X,c!*) + IL(m~U) <-> K(k~X,c!*) + IL(m~P)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]

    K(c) + K(c) <-> K(c!1).K(c!1)
    onRate3Dka = 1.0
    offRatekb = 1.0
    sigma = 1.0
    norm1 = [0,0,1]
    norm2 = [0,0,1]
    assocAngles = [nan,nan,nan,nan,nan]
end reactions
)";
const char* lipidStatesMolIL = R"(Name = IL
isImplicitLipid = true
checkOverlap = false
D = [0.5, 0.5, 0.0]
Dr = [0, 0, 0]
COM 0.0 0.0 0.0
m 0.0 0.0 0.0
bonds = 1
com m
state = m~U~P~Q
)";
const char* lipidStatesMolA = R"(Name = A
checkOverlap = true
D = [10.0, 10.0, 10.0]
Dr = [0.1, 0.1, 0.1]
COM 0.0 0.0 0.0
m 0.0 0.0 -1.0
a 1.0 0.0 0.0
bonds = 2
com m
com a
)";
const char* lipidStatesMolB = R"(Name = B
checkOverlap = true
D = [10.0, 10.0, 10.0]
Dr = [0.1, 0.1, 0.1]
COM 0.0 0.0 0.0
m 0.0 0.0 -1.0
a 1.0 0.0 0.0
bonds = 2
com m
com a
)";
const char* lipidStatesMolC = R"(Name = C
checkOverlap = true
D = [10.0, 10.0, 10.0]
Dr = [0.1, 0.1, 0.1]
COM 0.0 0.0 0.0
m 0.0 0.0 -1.0
a 1.0 0.0 0.0
bonds = 2
com m
com a
)";
const char* lipidStatesMolK = R"(Name = K
checkOverlap = true
D = [10.0, 10.0, 10.0]
Dr = [0.1, 0.1, 0.1]
COM 0.0 0.0 0.0
k 1.0 0.0 0.0
c -1.0 0.0 0.0
bonds = 2
com k
com c
state = k~X~Y
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
    long long partners { 0 }, lookups { 0 }, mismatches { 0 }, found { 0 }, none { 0 }, backStateChange { 0 },
        lipidOffset { 0 }, wrongLipidState { 0 }, foundWithOffset { 0 };
    // The most reactions one protein interface reaches the lipid through.
    std::size_t mostLipidRxns { 0 };
    // Which path the lipid lookup took, read back from its candidate list.
    long long fixedRate { 0 }, testedRate { 0 }, testedFound { 0 }, backSkipped { 0 };

    void add(const Totals& o)
    {
        partners += o.partners;
        lookups += o.lookups;
        mismatches += o.mismatches;
        found += o.found;
        none += o.none;
        backStateChange += o.backStateChange;
        lipidOffset += o.lipidOffset;
        wrongLipidState += o.wrongLipidState;
        foundWithOffset += o.foundWithOffset;
        mostLipidRxns = std::max(mostLipidRxns, o.mostLipidRxns);
        fixedRate += o.fixedRate;
        testedRate += o.testedRate;
        testedFound += o.testedFound;
        backSkipped += o.backSkipped;
    }
};

// Walks the candidates the way find_which_reaction_for_partner() does,
// to tell which of its paths produced the answer.  Only for the tally above;
// true when the answer came from a rate test.
bool classify(const std::vector<Interface::State::RxnCandidate>& candidates, const Molecule& m1, const Molecule& m2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns, Totals& totals)
{
    for (const auto& c : candidates) {
        const ForwardRxn& rxn = forwardRxns[c.rxnItr];
        if (c.fixedRateIndex != -1) {
            ++totals.fixedRate;
            return false;
        }
        ++totals.testedRate;
        if (c.backStateChange) {
            if (best_matching_rate(backRxns[rxn.conjBackRxnIndex], c.reactIndex1, c.reactIndex2, m1, m2).matches
                == 0) {
                ++totals.backSkipped;
                continue;
            }
        }
        return true;
    }
    return false;
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

    int lipidType { -1 };
    for (std::size_t t = 0; t < templates.size(); ++t)
        if (templates[t].isImplicitLipid)
            lipidType = static_cast<int>(t);
    if (lipidType == -1) {
        std::printf("  %s has no implicit lipid\n", model.dir.c_str());
        std::exit(1);
    }
    const MolTemplate& lipidTemplate { templates[lipidType] };
    // The simulator's lipid molecule carries its first state's index; the
    // pair check adds the offset of the state it tests.
    const int lipidBaseIndex { lipidTemplate.interfaceList[0].stateList[0].index };
    const int nStates { static_cast<int>(lipidTemplate.interfaceList[0].stateList.size()) };

    Totals here;
    for (std::size_t t1 = 0; t1 < templates.size(); ++t1) {
        if (templates[t1].isImplicitLipid)
            continue;
        for (std::size_t i1 = 0; i1 < templates[t1].interfaceList.size(); ++i1) {
            // Reactions found per lipid state, over all of this interface's states.
            std::map<unsigned, std::set<int>> lipidRxns;
            for (std::size_t s1 = 0; s1 < templates[t1].interfaceList[i1].stateList.size(); ++s1) {
                const Interface::State& state { templates[t1].interfaceList[i1].stateList[s1] };
                for (std::size_t p = 0; p < state.rxnPartners.size(); ++p) {
                    // The pair check tests a partner against the lipid when this offset is a lipid state.
                    const unsigned offset { state.rxnPartners[p] - static_cast<unsigned>(lipidBaseIndex) };
                    if (offset >= static_cast<unsigned>(nStates))
                        continue;
                    ++here.partners;
                    for (int trial = 0; trial < 2000; ++trial) {
                        Molecule m1 { random_molecule(templates[t1], static_cast<int>(t1), rng) };
                        Molecule lipid { random_molecule(lipidTemplate, lipidType, rng) };
                        set_state(m1, templates[t1], static_cast<int>(i1), static_cast<int>(s1));
                        lipid.interfaceList[0].index = lipidBaseIndex;

                        int rxnA { -1 }, rateA { -1 }, rxnB { -1 }, rateB { -1 };
                        bool scA { false }, scB { false };
                        const unsigned long before { totMatches };
                        lipid.interfaceList[0].index += static_cast<int>(offset);
                        find_which_reaction(static_cast<int>(i1), 0, rxnA, rateA, scA, state, m1, lipid, forwardRxns,
                            backRxns);
                        lipid.interfaceList[0].index -= static_cast<int>(offset);
                        const unsigned long countA { totMatches - before };
                        find_which_reaction_for_partner(static_cast<int>(i1), p, rxnB, rateB, scB, state, m1,
                            lipid, forwardRxns, backRxns);
                        const unsigned long countB { totMatches - before - countA };

                        ++here.lookups;
                        if (offset != 0)
                            ++here.lipidOffset;
                        if (rxnA != rxnB || rateA != rateB || scA != scB || countA != countB
                            || lipid.interfaceList[0].index != lipidBaseIndex) {
                            if (++here.mismatches < 4)
                                std::printf("  mismatch in %s: old (%d, %d, %d, %lu), new (%d, %d, %d, %lu)\n",
                                    model.dir.c_str(), rxnA, rateA, scA, countA, rxnB, rateB, scB, countB);
                        }
                        if (rxnA != -1 && rateA != -1) {
                            ++here.found;
                            if (offset != 0)
                                ++here.foundWithOffset;
                            lipidRxns[offset].insert(rxnA);
                            // The lipid state the reaction consumes: the
                            // forward reactant, or for a back state change the
                            // forward product, must be the state tested.
                            const std::vector<RxnIface>& sides { scA ? forwardRxns[rxnA].productListNew
                                                                     : forwardRxns[rxnA].reactantListNew };
                            const RxnIface& lipidSide { templates[sides[1].molTypeIndex].isImplicitLipid ? sides[1]
                                                                                                         : sides[0] };
                            if (lipidSide.absIfaceIndex != static_cast<int>(state.rxnPartners[p])) {
                                if (++here.wrongLipidState < 4)
                                    std::printf("  %s: %s state %zu meets lipid state %u through rxn %d, whose lipid "
                                                "is index %d, not %u\n",
                                        model.dir.c_str(), templates[t1].molName.c_str(), s1, offset, rxnA,
                                        lipidSide.absIfaceIndex, state.rxnPartners[p]);
                            }
                        } else
                            ++here.none;
                        if (scA)
                            ++here.backStateChange;
                        if (classify(state.rxnCandidatesByPartner[p], m1, lipid, forwardRxns, backRxns, here)
                            && rxnA != -1 && rateA != -1)
                            ++here.testedFound;
                    }
                }
            }
            std::set<int> allRxns;
            for (const auto& entry : lipidRxns)
                allRxns.insert(entry.second.begin(), entry.second.end());
            here.mostLipidRxns = std::max(here.mostLipidRxns, allRxns.size());
            if (nStates > 1 && !lipidRxns.empty()) {
                std::printf("  %s(%s) meets", templates[t1].molName.c_str(),
                    templates[t1].interfaceList[i1].name.c_str());
                for (const auto& entry : lipidRxns) {
                    std::printf(" IL~%c through rxn", lipidTemplate.interfaceList[0].stateList[entry.first].iden);
                    for (int rxn : entry.second)
                        std::printf(" %d", rxn);
                    std::printf(";");
                }
                std::printf("\n");
            }
        }
    }

    std::printf("  %-48s %3lld partners %7lld lookups: %lld found (%lld back state change, %lld with a lipid "
                "offset), %lld none, %lld with a lipid offset, %lld mismatches, %lld in the wrong lipid state\n"
                "  %-48s paths: fixed rate %lld, tested rate %lld (%lld found), back state change skipped %lld\n",
        model.dir.c_str(), here.partners, here.lookups, here.found, here.backStateChange, here.foundWithOffset,
        here.none, here.lipidOffset, here.mismatches, here.wrongLipidState, "", here.fixedRate, here.testedRate,
        here.testedFound, here.backSkipped);
    totals.add(here);
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

    char scratchTemplate[] = "/tmp/lipid_lookup_check.XXXXXX";
    const char* scratch { mkdtemp(scratchTemplate) };
    if (scratch == nullptr)
        return 1;
    const std::string lipidStatesDir { std::string(scratch) + "/lipid_states" };
    mkdir(lipidStatesDir.c_str(), 0755);
    write_file(lipidStatesDir + "/parms.inp", lipidStatesParm);
    write_file(lipidStatesDir + "/IL.mol", lipidStatesMolIL);
    write_file(lipidStatesDir + "/A.mol", lipidStatesMolA);
    write_file(lipidStatesDir + "/B.mol", lipidStatesMolB);
    write_file(lipidStatesDir + "/C.mol", lipidStatesMolC);
    write_file(lipidStatesDir + "/K.mol", lipidStatesMolK);

    const std::vector<Model> models {
        { lipidStatesDir, "parms.inp" },
        { "sample_inputs/VALIDATE_SUITE/implicit_lipid", "parms.inp" },
        { "sample_inputs/VALIDATE_SUITE/sphere", "parms_sphere.inp" },
        { "sample_inputs/VALIDATE_SUITE/mem_localization/LargeBox/IL", "parms.inp" },
        { "sample_inputs/clathrin_invitro_invivo/invitro", "parms.inp" },
        { "sample_inputs/clathrin_invitro_invivo/invivo", "parms.inp" },
        { "sample_inputs/clathrin_pioneer/ffl", "parmsffl.inp" },
        { "sample_inputs/clathrin_pioneer/ap2clafchoeps", "parms.inp" },
        { "sample_inputs/testAdd/box", "parms.inp" },
        { "sample_inputs/gagsphere", "parms.inp" },
        { "sample_inputs/secretion_catalyzed_assembly/membrane/kd1", "parms.inp" },
    };

    Totals totals;
    for (const Model& model : models)
        check_model(root, model, rng, totals);
    gsl_rng_free(r);

    std::printf("lipid_lookup_check: %lld partners, %lld lookups: %lld found (%lld back state change), %lld none, "
                "%lld with a lipid offset, %lld mismatches, %lld in the wrong lipid state, at most %zu reactions "
                "with the lipid per protein interface\n"
                "lipid_lookup_check: paths: fixed rate %lld, tested rate %lld (%lld found), back state change skipped %lld\n",
        totals.partners, totals.lookups, totals.found, totals.backStateChange, totals.none, totals.lipidOffset,
        totals.mismatches, totals.wrongLipidState, totals.mostLipidRxns, totals.fixedRate, totals.testedRate,
        totals.testedFound, totals.backSkipped);
    // Every path must have been taken, or the check proves less than it says.
    // A's m meets the three lipid states through three reactions; the
    // type-matched lookup found only the first.  K's k meets IL~P only through
    // the back state change.
    if (totals.mismatches != 0 || totals.wrongLipidState != 0 || totals.mostLipidRxns < 3
        || totals.backStateChange == 0
        || totals.foundWithOffset == 0 || totals.none == 0 || totals.found == 0
        || totals.lipidOffset == 0 || totals.fixedRate == 0 || totals.testedRate == 0 || totals.testedFound == 0
        || totals.testedFound == totals.testedRate) {
        std::printf("lipid_lookup_check: FAILED\n");
        return 1;
    }
    std::printf("lipid_lookup_check: all passed\n");
    return 0;
}
