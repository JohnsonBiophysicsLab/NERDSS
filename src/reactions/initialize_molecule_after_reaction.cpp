#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"
#include <numeric>

namespace {

/*!
 * \brief Puts the product's interface states on, as the reaction names them.
 *
 * Both draws below leave every interface in its template's first state, so this
 * has to follow each one.
 */
void set_product_iface_states(Molecule& mol, const MolTemplate& molTemplate, const CreateDestructRxn& currRxn)
{
    // TODO: Need to finish this, state should be set to reaction's product's state
    if (currRxn.productMolList.back().molTypeIndex != molTemplate.molTypeIndex)
        return; // TODO: do something?

    for (const auto& rxnIface : currRxn.productMolList.back().interfaceList) {
        mol.interfaceList[rxnIface.relIfaceIndex].stateIden = rxnIface.requiresState;
        // set the stateIndex by finding the state in the MolTemplate Interface's stateList
        for (unsigned stateItr { 0 }; stateItr < molTemplate.interfaceList[rxnIface.relIfaceIndex].stateList.size();
             ++stateItr) {
            if (molTemplate.interfaceList[rxnIface.relIfaceIndex].stateList[stateItr].iden == rxnIface.requiresState) {
                mol.interfaceList[rxnIface.relIfaceIndex].stateIndex = stateItr;
                break;
            }
        }
    }
}

} // namespace

void draw_coords_after_zeroth_reaction(
    Molecule& mol, const MolTemplate& molTemplate, const CreateDestructRxn& currRxn, const Membrane& membraneObject)
{
    mol.create_random_coords(molTemplate, membraneObject);
    set_product_iface_states(mol, molTemplate, currRxn);
}

void draw_coords_after_uni_reaction(
    Molecule& mol, const Molecule& parentMol, const MolTemplate& molTemplate, const CreateDestructRxn& currRxn)
{
    // put the specie in a random position around the specie which created it
    double theta { rand_gsl() * 2 * M_PI };
    double phi { std::acos(rand_gsl() * 2 - 1) };
    double cosTheta { std::cos(theta) };
    double sinTheta { std::sin(theta) };
    double cosPhi { std::cos(phi) };
    double sinPhi { std::sin(phi) };
    Vec3D transVec { currRxn.creationRadius * cosTheta * sinPhi, currRxn.creationRadius * sinTheta * sinPhi,
        currRxn.creationRadius * cosPhi };

    // create the coordinates
    mol.comCoord = parentMol.comCoord + transVec;
    for (unsigned ifaceItr { 0 }; ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
        mol.interfaceList[ifaceItr].coord = molTemplate.interfaceList[ifaceItr].iCoord + mol.comCoord;
        mol.interfaceList[ifaceItr].index = molTemplate.interfaceList[ifaceItr].stateList[0].index;
        mol.interfaceList[ifaceItr].relIndex = ifaceItr;
        mol.interfaceList[ifaceItr].stateIden = molTemplate.interfaceList[ifaceItr].stateList[0].iden;
        mol.interfaceList[ifaceItr].stateIndex = 0;
        mol.interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
    }

    // set the interface states to the states defined in the reaction
    set_product_iface_states(mol, molTemplate, currRxn);
}

Molecule initialize_molecule_after_zeroth_reaction(
    int index, Parameters& params, MolTemplate& molTemplate, const CreateDestructRxn& currRxn, const Membrane& membraneObject)
{
    // TRACE();
    /*!
     * \brief Creates a molecule according to a creation from concentration reaction (CreateDestructRxn), assigns
     * coords.
     *
     * params[in] molTemplate MolTemplate of the molecule to be created
     */

    Molecule tmp {};
    tmp.molTypeIndex = molTemplate.molTypeIndex;
    tmp.mass = molTemplate.mass;
    tmp.isLipid = molTemplate.isLipid;
    tmp.isPromoter = molTemplate.isPromoter;
    // Set up interface state vectors
    tmp.freelist = std::vector<int>(molTemplate.interfaceList.size());
    std::iota(tmp.freelist.begin(), tmp.freelist.end(), 0);
    tmp.interfaceList = std::vector<Molecule::Iface>(molTemplate.interfaceList.size());

    // Create center of mass
    draw_coords_after_zeroth_reaction(tmp, molTemplate, currRxn, membraneObject);

    // clean up
    tmp.isEmpty = false;

    // iterate number of molecules in the system and set index
    tmp.index = index;
    ++Molecule::numberOfMolecules;
    params.numTotalUnits = params.numTotalUnits + molTemplate.interfaceList.size() + 1;

    tmp.id = Molecule::maxID++;

    // keep track of molecule types
    ++MolTemplate::numEachMolType[molTemplate.molTypeIndex];

    return tmp;
}

Molecule initialize_molecule_after_uni_reaction(int index, const Molecule& parentMol, Parameters& params,
    MolTemplate& molTemplate, const CreateDestructRxn& currRxn)
{
    // TRACE();
    Molecule tmp {};
    tmp.molTypeIndex = molTemplate.molTypeIndex;
    tmp.mass = molTemplate.mass;
    tmp.isLipid = molTemplate.isLipid;
    tmp.isPromoter = molTemplate.isPromoter;

    // Set up interface state vectors
    tmp.freelist = std::vector<int>(molTemplate.interfaceList.size());
    std::iota(tmp.freelist.begin(), tmp.freelist.end(), 0);
    tmp.interfaceList = std::vector<Molecule::Iface>(molTemplate.interfaceList.size());

    draw_coords_after_uni_reaction(tmp, parentMol, molTemplate, currRxn);

    // clean up
    tmp.isEmpty = false;

    // iterate number of molecules in the system and set index
    tmp.index = index;
    ++Molecule::numberOfMolecules;
    params.numTotalUnits = params.numTotalUnits + molTemplate.interfaceList.size() + 1;

    tmp.id = Molecule::maxID++;

    // keep track of molecule types
    ++MolTemplate::numEachMolType[molTemplate.molTypeIndex];

    return tmp;
}
