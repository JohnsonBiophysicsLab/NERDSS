/*! \file association_dispatch.cpp
 * \brief Picks the box or sphere implementation for each association-time reaction.
 *
 * Four files held one of these each, and they are the same twelve lines with
 * different names substituted in.  Together they make the pattern visible: the
 * geometry is chosen once, at the top, and the `_box` / `_sphere` pair below is
 * handed the arguments untouched.
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"

#include <cmath>
#include <iomanip>

void associate(long long int iter,
    int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, const Parameters& params,
    ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList,
    copyCounters& counterArrays, std::vector<Complex>& complexList,
    Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, std::ofstream& assocDissocFile)
{
    if (membraneObject.isSphere == true) {
        associate_sphere(iter, ifaceIndex1, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2, params,
            currRxn, moleculeList, molTemplateList, observablesList,
            counterArrays, complexList, membraneObject, forwardRxns, backRxns, assocDissocFile);
    } else {
        associate_box(iter, ifaceIndex1, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2, params,
            currRxn, moleculeList, molTemplateList, observablesList,
            counterArrays, complexList, membraneObject, forwardRxns, backRxns, assocDissocFile);
    }
}

void associate_implicitlipid(long long int iter,
    int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, const Parameters& params,
    ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList,
    copyCounters& counterArrays, std::vector<Complex>& complexList,
    Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, std::ofstream& assocDissocFile)
{
    if (membraneObject.isSphere == true) {
        associate_implicitlipid_sphere(iter, ifaceIndex1, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2, params,
            currRxn, moleculeList, molTemplateList, observablesList,
            counterArrays, complexList, membraneObject, forwardRxns, backRxns, assocDissocFile);
    } else {
        associate_implicitlipid_box(iter, ifaceIndex1, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2, params,
            currRxn, moleculeList, molTemplateList, observablesList,
            counterArrays, complexList, membraneObject, forwardRxns, backRxns, assocDissocFile);
    }
}

void perform_bimolecular_state_change(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject)
{
    if (membraneObject.isSphere == true) {
        perform_bimolecular_state_change_sphere(stateChangeIface, facilitatorIface, rxnItr,
            stateChangeMol, facilitatorMol, stateChangeCom, facilitatorCom,
            counterArrays, params, forwardRxns,
            backRxns, moleculeList, complexList,
            molTemplateList, observablesList, membraneObject);
    } else {
        perform_bimolecular_state_change_box(stateChangeIface, facilitatorIface, rxnItr,
            stateChangeMol, facilitatorMol, stateChangeCom, facilitatorCom,
            counterArrays, params, forwardRxns,
            backRxns, moleculeList, complexList,
            molTemplateList, observablesList, membraneObject);
    }
}

void perform_implicitlipid_state_change(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject)
{
    if (membraneObject.isSphere == true) {
        perform_implicitlipid_state_change_sphere(stateChangeIface, facilitatorIface, rxnItr,
            stateChangeMol, facilitatorMol, stateChangeCom, facilitatorCom,
            counterArrays, params, forwardRxns,
            backRxns, moleculeList, complexList,
            molTemplateList, observablesList, membraneObject);
    } else {
        perform_implicitlipid_state_change_box(stateChangeIface, facilitatorIface, rxnItr,
            stateChangeMol, facilitatorMol, stateChangeCom, facilitatorCom,
            counterArrays, params, forwardRxns,
            backRxns, moleculeList, complexList,
            molTemplateList, observablesList, membraneObject);
    }
}
