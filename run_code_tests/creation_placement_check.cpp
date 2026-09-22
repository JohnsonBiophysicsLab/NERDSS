// Calls moleculeOverlaps() and create_molecule_and_complex_from_rxn(), which
// place the molecules a creation reaction makes during a run -- zeroth order
// (`NULL -> A(a)`) and unimolecular (`A(a) -> A(a) + B(b)`).
//
// A created molecule may not land with an interface closer to a reaction
// partner's than their reaction's bindRadius.  moleculeOverlaps() is the test,
// and it was broken in three ways: it reported "no overlap" for the whole
// molecule at the first SubBox member whose complex lay out of range, so it
// rarely looked past that member; its range cut left out the bindRadius, so it
// skipped partners whose interfaces were in reach whenever the bounding spheres
// were less than a bindRadius apart; and it read reactantListNew[1] of every
// forward reaction, one past the end of a unimolecular state change's list.
//
// A working test resamples, which exposes the loop around it: each attempt used
// to re-run initialize_molecule_after_zeroth_reaction() /
// initialize_molecule_after_uni_reaction(), and those count the molecule in
// Molecule::numberOfMolecules, MolTemplate::numEachMolType,
// Parameters::numTotalUnits and Molecule::maxID.
#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
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

// One SubBox over the whole volume, so every molecule already placed is a
// member of the one a new molecule lands in.
void make_one_cell_volume(SimulVolume& simulVolume, double boxSide)
{
    simulVolume.numSubCells.x = 1;
    simulVolume.numSubCells.y = 1;
    simulVolume.numSubCells.z = 1;
    simulVolume.numSubCells.tot = 1;
    simulVolume.subCellSize = Vec3D { boxSide, boxSide, boxSide };
    simulVolume.subCellList.resize(1);
    simulVolume.subCellList[0].absIndex = 0;
    simulVolume.occupancyMask.assign(1, 0);
}

// A(a) + B(b) at bindRadius 2; each template's interface sits 1 nm from its
// center.  overlaps() asks moleculeOverlaps() about a B being placed with its
// center at `com` and its interface at `com + offset`.
struct Placement {
    Parameters params {};
    Membrane membrane {};
    SimulVolume simulVolume {};
    std::vector<MolTemplate> templates;
    std::vector<Molecule> molecules;
    std::vector<Complex> complexes;
    std::vector<ForwardRxn> forwardRxns;

    explicit Placement(double boxSide = 100.0)
    {
        membrane.waterBox = Membrane::WaterBox { std::vector<double> { boxSide, boxSide, boxSide } };
        make_one_cell_volume(simulVolume, boxSide);

        templates.push_back(make_template("A", 0, 0, 0, 1.0));
        templates.push_back(make_template("B", 1, 0, 1, 1.0));
        MolTemplate::numMolTypes = templates.size();
        MolTemplate::numEachMolType.assign(templates.size(), 0);
        forwardRxns.push_back(two_reactant_rxn(ReactionType::bimolecular, 0, 0, 1, 1, 2.0));
    }

    // A molecule already in place: in its own complex, and a member of the one
    // SubBox, which is where moleculeOverlaps() looks for it.
    int add_mol(int molTypeIndex, int stateIndex, const Vec3D& com, const Vec3D& offset)
    {
        int index = molecules.size();
        molecules.push_back(make_molecule(index, molTypeIndex, stateIndex, com, offset));
        molecules.back().myComIndex = complexes.size();
        complexes.emplace_back(molecules.back().myComIndex, molecules.back(), templates[molTypeIndex]);
        simulVolume.add_member(0, index, molTypeIndex);
        return index;
    }

    int add_A(const Vec3D& com, const Vec3D& offset) { return add_mol(0, 0, com, offset); }

    bool overlaps(const Vec3D& com, const Vec3D& offset)
    {
        molecules.push_back(make_molecule(molecules.size(), 1, 1, com, offset));
        bool result = moleculeOverlaps(params, simulVolume, molecules.back(), molecules, complexes, forwardRxns,
            templates, membrane);
        // moleculeOverlaps() makes the molecule a member of the SubBox when it
        // reports no overlap; this one is not staying.
        if (!result)
            simulVolume.subCellList[0].memberMolList.pop_back();
        molecules.pop_back();
        return result;
    }
};

// `NULL -> B(b~<productState>)`, creation from concentration.  Every draw leaves
// an interface in its template's first state, so a molecule moved off an overlap
// has to have the product's state put back on; pass '\0' for a template whose
// one state is the first.
CreateDestructRxn zeroth_order_rxn(int molTypeIndex, int relIfaceIndex, char productState)
{
    CreateDestructRxn rxn {};
    rxn.rxnType = ReactionType::zerothOrderCreation;
    rxn.productMolList.emplace_back();
    rxn.productMolList.back().molTypeIndex = molTypeIndex;
    rxn.productMolList.back().interfaceList.emplace_back(
        "b", molTypeIndex, 0, relIfaceIndex, productState, false);
    return rxn;
}

// `existing` A at random places in a box of side `boxSide`, and the reaction
// `NULL -> B(b)` that create() then runs.  A(a) + B(b) bind at `bindRadius`, so
// with the box crowded a place drawn at random is usually too close to some A.
struct CreationFixture {
    Parameters params {};
    Membrane membrane {};
    SimulVolume simulVolume {};
    std::vector<MolTemplate> templates;
    std::vector<Molecule> molecules;
    std::vector<Complex> complexes;
    std::vector<ForwardRxn> forwardRxns;
    CreateDestructRxn createRxn {};

    CreationFixture(int existing, double boxSide, double bindRadius)
    {
        Molecule::numberOfMolecules = 0;
        Molecule::maxID = 0;
        Molecule::emptyMolList.clear();
        Complex::numberOfComplexes = 0;
        Complex::emptyComList.clear();

        membrane.waterBox = Membrane::WaterBox { std::vector<double> { boxSide, boxSide, boxSide } };
        make_one_cell_volume(simulVolume, boxSide);

        templates.push_back(make_template("A", 0, existing, 0, 1.0));
        templates.push_back(make_template("B", 1, 0, 1, 1.0));
        MolTemplate::numMolTypes = templates.size();
        MolTemplate::numEachMolType.assign(templates.size(), 0);
        forwardRxns.push_back(two_reactant_rxn(ReactionType::bimolecular, 0, 0, 1, 1, bindRadius));
        createRxn = zeroth_order_rxn(1, 0, '\0');

        for (int molIndex = 0; molIndex < existing; ++molIndex) {
            molecules.push_back(make_molecule(molIndex, 0, 0, Vec3D { 0.0, 0.0, 0.0 }, Vec3D { 1.0, 0.0, 0.0 }));
            molecules.back().create_random_coords(templates[0], membrane);
            complexes.emplace_back(molIndex, molecules.back(), templates[0]);
            simulVolume.add_member(0, molIndex, 0);
            ++Molecule::numberOfMolecules;
            ++Complex::numberOfComplexes;
            ++MolTemplate::numEachMolType[0];
            ++Molecule::maxID;
        }
    }

    void create(int copies)
    {
        for (int copyItr = 0; copyItr < copies; ++copyItr) {
            int newMolIndex { 0 };
            int newComIndex { 0 };
            create_molecule_and_complex_from_rxn(0, newMolIndex, newComIndex, false, templates[1], params, createRxn,
                simulVolume, molecules, complexes, templates, forwardRxns, membrane);
        }
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
    r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, 20260922);

    const Vec3D origin { 0.0, 0.0, 0.0 };
    const Vec3D plusX { 1.0, 0.0, 0.0 };
    const Vec3D minusX { -1.0, 0.0, 0.0 };
    const Vec3D plusY { 0.0, 1.0, 0.0 };

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
        p.simulVolume.add_member(0, 0, 0);
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "an emptied molecule slot is passed over");
    }
    {
        // The new molecule can hold the slot an emptied one left behind, and
        // that slot can still be a member of the SubBox.  Compared with itself
        // it would overlap at every place it was ever drawn.
        Placement p;
        p.forwardRxns = { two_reactant_rxn(ReactionType::bimolecular, 1, 1, 1, 1, 2.0) };
        p.simulVolume.add_member(0, 0, 1);
        p.molecules.push_back(make_molecule(0, 1, 1, origin, plusX));
        p.complexes.emplace_back(0, p.molecules.back(), p.templates[1]);
        bool self = moleculeOverlaps(p.params, p.simulVolume, p.molecules[0], p.molecules, p.complexes, p.forwardRxns,
            p.templates, p.membrane);
        expect(!self, "the molecule being placed is not compared with itself");
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
        MolTemplate::numMolTypes = p.templates.size();
        MolTemplate::numEachMolType.assign(p.templates.size(), 0);
        p.forwardRxns = { two_reactant_rxn(ReactionType::bimolecular, 1, 1, 2, 2, 2.0) };
        p.add_mol(2, 2, origin, plusX);
        p.molecules.back().isImplicitLipid = true;
        expect(!p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "an implicit lipid has no position to overlap");
        p.molecules.back().isImplicitLipid = false;
        expect(p.overlaps(Vec3D { 2.0, 0.0, 0.0 }, plusY), "the same molecule as an explicit one does");
    }
    {
        Placement p;
        p.add_A(origin, plusX);
        expect(p.overlaps(Vec3D { 60.0, 0.0, 0.0 }, plusX), "a place outside the box is rejected");
    }

    {
        // 400 A in a 20 nm box, then 100 B that bind A at bindRadius 2: about
        // four in five places drawn at random are too close to some A.
        CreationFixture f(400, 20.0, 2.0);
        const int unitsBefore = f.params.numTotalUnits;
        const int idBefore = Molecule::maxID;
        f.create(100);

        int tooClose = 0;
        for (int created = 400; created < int(f.molecules.size()); ++created)
            for (int existing = 0; existing < 400; ++existing)
                if (Vec3D { f.molecules[created].interfaceList[0].coord
                        - f.molecules[existing].interfaceList[0].coord }
                        .length()
                    < 2.0)
                    ++tooClose;
        expect(f.molecules.size() == 500 && tooClose == 0, "100 B created among 400 A: none within bindRadius of an A");
        expect(Molecule::numberOfMolecules == 500 && MolTemplate::numEachMolType[1] == 100
                && f.params.numTotalUnits == unitsBefore + 100 * 2 && Molecule::maxID == idBefore + 100,
            "each created B is counted once, however often it was moved");

        int members = 0;
        for (int memMol : f.simulVolume.subCellList[0].memberMolList)
            if (memMol >= 400)
                ++members;
        expect(members == 100, "each created B is a member of its SubBox exactly once");
    }

    {
        // The product's interface state has to survive being moved: the draw
        // leaves every interface in its template's first state.
        MolTemplate twoState = make_template("B", 0, 0, 1, 1.0);
        twoState.interfaceList[0].stateList.emplace_back('P', 2);
        std::vector<MolTemplate> templates { twoState };
        MolTemplate::numMolTypes = 1;
        MolTemplate::numEachMolType.assign(1, 0);
        Molecule::numberOfMolecules = 0;

        Parameters params {};
        Membrane membrane {};
        membrane.waterBox = Membrane::WaterBox { std::vector<double> { 100.0, 100.0, 100.0 } };
        CreateDestructRxn rxn = zeroth_order_rxn(0, 0, 'P');

        Molecule mol = initialize_molecule_after_zeroth_reaction(0, params, templates[0], rxn, membrane);
        bool stateKept { mol.interfaceList[0].stateIden == 'P' && mol.interfaceList[0].stateIndex == 1 };
        bool moved { false };
        for (int drawItr = 0; drawItr < 20; ++drawItr) {
            Vec3D before { mol.comCoord };
            draw_coords_after_zeroth_reaction(mol, templates[0], rxn, membrane);
            moved = moved || Vec3D { mol.comCoord - before }.length() > 0.0;
            stateKept = stateKept && mol.interfaceList[0].stateIden == 'P' && mol.interfaceList[0].stateIndex == 1
                && Vec3D { mol.interfaceList[0].coord - mol.comCoord }.length() > 0.0;
        }
        expect(moved && stateKept, "a molecule drawn again moves and keeps the reaction's product state");
        expect(MolTemplate::numEachMolType[0] == 1 && Molecule::numberOfMolecules == 1,
            "drawing it again does not count it again");
    }

    gsl_rng_free(r);
    if (failures) {
        std::printf("creation_placement_check: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("creation_placement_check: all passed\n");
    return 0;
}
