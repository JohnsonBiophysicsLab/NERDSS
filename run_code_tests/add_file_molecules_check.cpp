// Calls generate_coordinates_for_restart(), which places the molecules an add
// file brings into a restarted simulation (`nerdss -r restart.dat -a add.inp`),
// and checks what it leaves behind.
//
// Each added molecule has to reach its template's monomerList: that list is the
// pool check_for_unimolecular_reactions_population() destroys from.  The
// function used to fill a copy of the template, so an added type with a
// destruction reaction was never destroyed.
//
// No added molecule may start with an interface closer to a reaction partner's
// than their reaction's bindRadius.  moleculeOverlapsForRestart() is the test;
// it used to stop at the first molecule out of range, report pairs *beyond*
// bindRadius as overlapping, read a second reactant from unimolecular
// reactions and a complex from emptied molecule slots.
#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "system_setup/system_setup.hpp"
#include <cstdio>
#include <string>
#include <vector>
#include <gsl/gsl_rng.h>

gsl_rng* r; // the globals the simulator defines in main()
unsigned long totMatches = 0;
long long randNum = 0;

namespace {

// One interface `name` at `offset` from the center, in a single state with
// absolute index `stateIndex`.
MolTemplate make_template(const std::string& name, int molTypeIndex, int copies, int stateIndex, double offset)
{
    MolTemplate oneTemp {};
    oneTemp.molName = name;
    oneTemp.molTypeIndex = molTypeIndex;
    oneTemp.copies = copies;
    oneTemp.mass = 1.0;
    oneTemp.radius = offset;
    oneTemp.D = Vec3D { 1.0, 1.0, 1.0 };
    Interface iface { name + "1", Vec3D { offset, 0.0, 0.0 } };
    iface.index = 0;
    iface.stateList.emplace_back(stateIndex);
    oneTemp.interfaceList.push_back(iface);
    return oneTemp;
}

// A molecule of template `molTypeIndex` alone in complex `index`, its one
// interface free, in state `stateIndex`, at `com + offset`.
Molecule make_molecule(int index, int molTypeIndex, int stateIndex, const Vec3D& com, const Vec3D& offset)
{
    Molecule mol {};
    mol.index = index;
    mol.myComIndex = index;
    mol.molTypeIndex = molTypeIndex;
    mol.comCoord = com;
    mol.interfaceList.emplace_back(com + offset, '\0', stateIndex, molTypeIndex, false);
    mol.interfaceList[0].relIndex = 0;
    mol.interfaceList[0].stateIndex = 0;
    mol.freelist = { 0 };
    return mol;
}

// Two free interfaces, of the given templates and states, reacting at
// bindRadius.
ForwardRxn two_reactant_rxn(ReactionType rxnType, int molType1, int state1, int molType2, int state2, double bindRadius)
{
    ForwardRxn rxn {};
    rxn.rxnType = rxnType;
    rxn.bindRadius = bindRadius;
    rxn.reactantListNew.emplace_back("x", molType1, state1, 0, '\0', false);
    rxn.reactantListNew.emplace_back("y", molType2, state2, 0, '\0', false);
    return rxn;
}

// A(a) + B(b) at bindRadius 2; each template's interface sits 1 nm from its
// center.  overlaps() asks moleculeOverlapsForRestart() about a B being placed
// with its center at `com` and its interface at `com + offset`.
struct Placement {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> templates;
    std::vector<Molecule> molecules;
    std::vector<ForwardRxn> forwardRxns;

    Placement()
    {
        templates.push_back(make_template("A", 0, 0, 0, 1.0));
        templates.push_back(make_template("B", 1, 0, 1, 1.0));
        forwardRxns.push_back(two_reactant_rxn(ReactionType::bimolecular, 0, 0, 1, 1, 2.0));
    }

    void add_A(const Vec3D& com, const Vec3D& offset)
    {
        molecules.push_back(make_molecule(molecules.size(), 0, 0, com, offset));
    }

    bool overlaps(const Vec3D& com, const Vec3D& offset)
    {
        molecules.push_back(make_molecule(molecules.size(), 1, 1, com, offset));
        bool result = moleculeOverlapsForRestart(params, molecules.back(), molecules, forwardRxns, templates, membrane);
        molecules.pop_back();
        return result;
    }
};

// A restarted system holding `existing` molecules of template A at random
// places in a box of side `boxSide`, about to receive the add file's template
// C (index 1) with `addedCopies` copies.
struct AddFixture {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> templates;
    std::vector<Molecule> molecules;
    std::vector<Complex> complexes;
    std::vector<ForwardRxn> forwardRxns;

    AddFixture(int existing, int addedCopies, double boxSide)
    {
        Molecule::numberOfMolecules = 0;
        Molecule::emptyMolList.clear();
        Complex::numberOfComplexes = 0;
        Complex::emptyComList.clear();

        membrane.waterBox = Membrane::WaterBox { std::vector<double> { boxSide, boxSide, boxSide } };

        templates.push_back(make_template("A", 0, existing, 0, 1.0));
        templates.push_back(make_template("C", 1, addedCopies, 1, 1.0));
        MolTemplate::numMolTypes = templates.size();
        MolTemplate::numEachMolType.assign(templates.size(), 0);

        for (int molIndex = 0; molIndex < existing; ++molIndex) {
            molecules.push_back(make_molecule(molIndex, 0, 0, Vec3D { 0.0, 0.0, 0.0 }, Vec3D { 1.0, 0.0, 0.0 }));
            molecules.back().create_random_coords(templates[0], membrane);
            complexes.emplace_back(molIndex, molecules.back(), templates[0]);
            ++Molecule::numberOfMolecules;
            ++Complex::numberOfComplexes;
            ++MolTemplate::numEachMolType[0];
        }
    }

    void add() { generate_coordinates_for_restart(params, molecules, complexes, templates, forwardRxns, membrane, 1, 0); }
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
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 20260922);

    const Vec3D origin { 0.0, 0.0, 0.0 };
    const Vec3D plusX { 1.0, 0.0, 0.0 };
    const Vec3D minusX { -1.0, 0.0, 0.0 };
    const Vec3D plusY { 0.0, 1.0, 0.0 };

    {
        AddFixture f(1, 3, 100.0);
        f.add();
        expect(f.molecules.size() == 4, "the add file's 3 molecules are created");
        expect(f.templates[1].monomerList == std::vector<int> { 1, 2, 3 },
            "the added template's monomerList holds each of them");
        expect(f.templates[0].monomerList.empty(), "the restart file's template is left alone");
        bool ownComplexes = true;
        for (int molIndex : { 1, 2, 3 }) {
            const Molecule& mol = f.molecules[molIndex];
            ownComplexes = ownComplexes && mol.molTypeIndex == 1 && !mol.isEmpty
                && f.complexes[mol.myComIndex].memberList == std::vector<int> { molIndex };
        }
        expect(ownComplexes, "each is a monomer of template C in its own complex");
    }

    // A's interface at (1, 0, 0) throughout, unless it points away.
    {
        Placement p;
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "reactive interfaces 1.41 nm apart (bindRadius 2) overlap");
        expect(!p.overlaps(Vec3D { 3.0, 0.0, 0.0 }, plusX), "reactive interfaces 3 nm apart do not");
    }
    {
        Placement p;
        p.add_A(origin, minusX);
        expect(!p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusX), "interfaces 4 nm apart, centers 2 nm apart, do not overlap");
    }
    {
        Placement p;
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 3.5, 0.0, 0.0 }, minusX),
            "interfaces 1.5 nm apart, centers 3.5 nm (past both radii), overlap");
    }
    {
        Placement p;
        p.add_A(Vec3D { 40.0, 40.0, 40.0 }, plusX);
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "a partner in range is found after one out of range");
    }
    {
        Placement p;
        p.molecules.emplace_back();
        p.molecules.back().index = 0;
        p.molecules.back().isEmpty = true; // destroyed: myComIndex -1, no interfaces
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "an emptied molecule slot is passed over");
    }
    {
        Placement p;
        p.add_A(origin, plusX);
        p.molecules.back().interfaceList[0].isBound = true;
        expect(!p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "a bound interface does not overlap");
    }
    {
        // A(s~U) -> A(s~P) ahead of the binding reaction: one reactant only
        Placement p;
        ForwardRxn stateChange {};
        stateChange.rxnType = ReactionType::uniMolStateChange;
        stateChange.bindRadius = 50.0;
        stateChange.reactantListNew.emplace_back("s", 0, 0, 0, '\0', false);
        p.forwardRxns.insert(p.forwardRxns.begin(), stateChange);
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "a unimolecular reaction is passed over");
    }
    {
        Placement p;
        p.forwardRxns[0].rxnType = ReactionType::biMolStateChange;
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "a bimolecular state change overlaps like a binding");
    }
    {
        // B(b) + M(m), M an implicit lipid whose molecule sits at the origin
        Placement p;
        p.templates.push_back(make_template("M", 2, 0, 2, 1.0));
        p.templates.back().isImplicitLipid = true;
        p.forwardRxns = { two_reactant_rxn(ReactionType::bimolecular, 1, 1, 2, 2, 2.0) };
        p.molecules.push_back(make_molecule(0, 2, 2, origin, plusX));
        p.molecules.back().isImplicitLipid = true;
        expect(!p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "an implicit lipid has no position to overlap");
        p.molecules.back().isImplicitLipid = false;
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "the same molecule as an explicit one does");
    }

    {
        // 400 A in a 20 nm box, then 100 C that bind A at bindRadius 2: about
        // four in five places drawn at random are too close to some A.
        AddFixture f(400, 100, 20.0);
        f.forwardRxns.push_back(two_reactant_rxn(ReactionType::bimolecular, 0, 0, 1, 1, 2.0));
        const int unitsBefore = f.params.numTotalUnits;
        f.add();
        int tooClose = 0;
        for (int added = 400; added < 500; ++added)
            for (int existing = 0; existing < 400; ++existing)
                if (Vec3D { f.molecules[added].interfaceList[0].coord - f.molecules[existing].interfaceList[0].coord }
                        .length()
                    < 2.0)
                    ++tooClose;
        expect(f.molecules.size() == 500 && tooClose == 0, "100 C added among 400 A: none within bindRadius of an A");
        expect(Molecule::numberOfMolecules == 500 && MolTemplate::numEachMolType[1] == 100
                && f.params.numTotalUnits == unitsBefore + 100 * 2,
            "each placed C is counted once, however often it was moved");
    }

    gsl_rng_free(r);
    if (failures) {
        std::printf("add_file_molecules_check: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("add_file_molecules_check: all passed\n");
    return 0;
}
