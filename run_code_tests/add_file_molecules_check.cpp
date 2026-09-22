// Calls generate_coordinates_for_restart(), which places the molecules an add
// file brings into a restarted simulation (`nerdss -r restart.dat -a add.inp`),
// and checks what it leaves behind.
//
// Each added molecule has to reach its template's monomerList: that list is the
// pool check_for_unimolecular_reactions_population() destroys from.  The
// function used to fill a copy of the template, so an added type with a
// destruction reaction was never destroyed.
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

// A restarted system holding one molecule of template A, about to receive
// the add file's template C.
struct AddFixture {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> templates;
    std::vector<Molecule> molecules;
    std::vector<Complex> complexes;
    std::vector<ForwardRxn> forwardRxns;

    explicit AddFixture(int addedCopies)
    {
        Molecule::numberOfMolecules = 0;
        Molecule::emptyMolList.clear();
        Complex::numberOfComplexes = 0;
        Complex::emptyComList.clear();

        membrane.waterBox = Membrane::WaterBox { std::vector<double> { 100.0, 100.0, 100.0 } };

        templates.push_back(make_template("A", 0, 1, 0, 1.0));
        templates.push_back(make_template("C", 1, addedCopies, 1, 1.0));
        MolTemplate::numMolTypes = templates.size();
        MolTemplate::numEachMolType.assign(templates.size(), 0);

        Molecule mol {};
        mol.index = 0;
        mol.myComIndex = 0;
        mol.molTypeIndex = 0;
        mol.comCoord = Vec3D { 0.0, 0.0, 0.0 };
        mol.interfaceList.emplace_back(Vec3D { 1.0, 0.0, 0.0 }, '\0', 0, 0, false);
        mol.interfaceList[0].relIndex = 0;
        mol.freelist = { 0 };
        molecules.push_back(mol);
        complexes.emplace_back(0, molecules[0], templates[0]);
        ++Molecule::numberOfMolecules;
        ++Complex::numberOfComplexes;
        ++MolTemplate::numEachMolType[0];
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

    {
        AddFixture f(3);
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

    gsl_rng_free(r);
    if (failures) {
        std::printf("add_file_molecules_check: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("add_file_molecules_check: all passed\n");
    return 0;
}
