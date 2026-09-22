#include "io/io.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "system_setup/system_setup.hpp"
#include <algorithm>
#include <numeric>

//this is used to generate the added molecules and complexes
void generate_coordinates_for_restart(Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList,
    const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject, int numMolTemplateBeforeAdd, int numForwardRxnBdeforeAdd)
{
    for (int molTemplateListIndex = numMolTemplateBeforeAdd; molTemplateListIndex < molTemplateList.size(); molTemplateListIndex++) {
        // The template itself, not a copy: each new molecule goes onto its
        // monomerList, which is where destruction picks its victims.  With a
        // copy, the list was thrown away with it and no added molecule could be
        // destroyed.
        MolTemplate& oneTemp { molTemplateList[molTemplateListIndex] };
        for (unsigned itr { 0 }; itr < oneTemp.copies; ++itr) {
            create_molecule_and_complex_for_restart(oneTemp, params, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);
        }
    }
}

void create_molecule_and_complex_for_restart(MolTemplate& createdMolTemp, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject)
{
    int newMolIndex = 0;
    int newComIndex = 0;
    if (Molecule::emptyMolList.size() > 0) {
        // Check to make sure the object pointed to by the last element in the each list is actually empty
        // If it isn't delete it from the list and move on until you find one that is
        try {
            while (!moleculeList[Molecule::emptyMolList.back()].isEmpty)
                Molecule::emptyMolList.pop_back();

            // if there's an available empty Molecule spot, make the new Molecule in place
            newMolIndex = Molecule::emptyMolList.back();
            Molecule::emptyMolList.pop_back(); // remove the empty Molecule spot index from the list
        } catch (std::out_of_range& e) {
            newMolIndex = moleculeList.size();
            moleculeList.emplace_back(); // create new empty Molecule spot
        }
    } else {
        newMolIndex = moleculeList.size();
        moleculeList.emplace_back(); // create new empty Molecule spot
    }

    if (Complex::emptyComList.size() > 0) {
        // Check to make sure the object pointed to by the last element in the each list is actually empty
        // If it isn't delete it from the list and move on until you find one that is
        try {
            while (!complexList[Complex::emptyComList.back()].isEmpty)
                Complex::emptyComList.pop_back();
            // if there's an available empty Complex spot, make the new Complex in place
            newComIndex = Complex::emptyComList.back();
            Complex::emptyComList.pop_back(); // remove the empty Complex spot index from the list
        } catch (std::out_of_range) {
            newComIndex = complexList.size();
            complexList.emplace_back(); // create new empty Complex spot
        }

    } else {
        newComIndex = complexList.size();
        complexList.emplace_back(); // create new empty Complex spot
    }

    // Now create the new species, and move it for as long as it overlaps a
    // molecule already in place.  Only the coordinates are drawn again:
    // initialize_molecule_for_restart() also counts the molecule, in
    // numberOfMolecules, numEachMolType and numTotalUnits.
    moleculeList[newMolIndex] = initialize_molecule_for_restart(newMolIndex, params, createdMolTemp, membraneObject);
    while (moleculeOverlapsForRestart(
        params, moleculeList[newMolIndex], moleculeList, forwardRxns, molTemplateList, membraneObject))
        moleculeList[newMolIndex].create_random_coords(createdMolTemp, membraneObject);

    moleculeList[newMolIndex].myComIndex = newComIndex;
    complexList[newComIndex] = Complex { newComIndex, moleculeList.at(newMolIndex), createdMolTemp };
    ++Complex::numberOfComplexes;

    createdMolTemp.monomerList.emplace_back(newMolIndex); //add this new molecule to the monomerList
}

Molecule initialize_molecule_for_restart(
    int index, Parameters& params, MolTemplate& molTemplate, const Membrane& membraneObject)
{
    /*!
     * \brief Creates a molecule according to add.inp, assigns
     * coords.
     *
     * params[in] molTemplate MolTemplate of the molecule to be created
     */

    Molecule tmp {};
    // TODO: Need to finish this, state should be set to reaction's product's state
    tmp.molTypeIndex = molTemplate.molTypeIndex;
    tmp.mass = molTemplate.mass;
    tmp.isLipid = molTemplate.isLipid;
    tmp.isPromoter = molTemplate.isPromoter;

    // Set up interface state vectors
    tmp.freelist = std::vector<int>(molTemplate.interfaceList.size());
    std::iota(tmp.freelist.begin(), tmp.freelist.end(), 0);
    tmp.interfaceList = std::vector<Molecule::Iface>(molTemplate.interfaceList.size());

    // Create center of mass
    tmp.create_random_coords(molTemplate, membraneObject);

    // clean up
    tmp.isEmpty = false;

    // iterate number of molecules in the system and set index
    tmp.index = index;
    ++Molecule::numberOfMolecules;
    params.numTotalUnits = params.numTotalUnits + molTemplate.interfaceList.size() + 1;

    // Molecule() leaves id uninitialized; take the next one, as
    // initialize_molecule() does for a new simulation
    tmp.id = Molecule::maxID++;

    // keep track of molecule types
    ++MolTemplate::numEachMolType[molTemplate.molTypeIndex];

    return tmp;
}

bool moleculeOverlapsForRestart(const Parameters& params, const Molecule& createdMol,
    const std::vector<Molecule>& moleculeList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // The new molecule overlaps when one of its interfaces and one of another
    // molecule's could react with each other and sit closer than that
    // reaction's bindRadius, the criterion generate_coordinates() uses for a new
    // simulation.  Only a reaction between two interfaces can overlap: a
    // unimolecular state change has a single reactant.
    double maxBindRadius { 0.0 };
    for (const auto& oneRxn : forwardRxns) {
        if (oneRxn.reactantListNew.size() == 2)
            maxBindRadius = std::max(maxBindRadius, oneRxn.bindRadius);
    }
    const double createdRadius { molTemplateList[createdMol.molTypeIndex].radius };

    // Every molecule, by reference.  This used to return at the first molecule
    // it did not overlap, so it checked little beyond moleculeList[0].
    for (const auto& partMol : moleculeList) {
        // An emptied slot has no complex (myComIndex is -1) and no interfaces,
        // and an implicit lipid has no position to overlap; the time step's
        // overlap checks skip both too.
        if (partMol.index == createdMol.index || partMol.isEmpty || partMol.isImplicitLipid)
            continue;

        // Each interface lies within its template's radius of the molecule's
        // center, so from this far apart no pair can be within any bindRadius.
        Vec3D comVec { createdMol.comCoord - partMol.comCoord };
        if (comVec.length() >= createdRadius + molTemplateList[partMol.molTypeIndex].radius + maxBindRadius)
            continue;

        for (const auto& oneRxn : forwardRxns) {
            if (oneRxn.reactantListNew.size() != 2)
                continue;
            for (const auto& createdIface : createdMol.interfaceList) {
                for (const auto& partIface : partMol.interfaceList) {
                    bool canReact { (isReactant(createdIface, createdMol, oneRxn.reactantListNew[0])
                                        && isReactant(partIface, partMol, oneRxn.reactantListNew[1]))
                        || (isReactant(createdIface, createdMol, oneRxn.reactantListNew[1])
                            && isReactant(partIface, partMol, oneRxn.reactantListNew[0])) };
                    if (canReact && Vec3D { createdIface.coord - partIface.coord }.length() < oneRxn.bindRadius)
                        return true;
                }
            }
        }
    }

    return false;
}