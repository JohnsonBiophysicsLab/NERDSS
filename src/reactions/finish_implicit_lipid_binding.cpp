/*! \file finish_implicit_lipid_binding.cpp
 * \brief The bookkeeping that turns an accepted protein-to-implicit-lipid
 * association into a bound state.
 *
 * associate_implicitlipid_box() carried this twice, once as "mol2 is the
 * implicit lipid" and once as "mol1 is implicit lipid", and
 * associate_implicitlipid_sphere() carried a third copy that only knew the
 * first role.  None of it is geometry: every step is a copy, a counter, a flag
 * or a list edit, so the two roles differ only in which of the pair is the
 * protein and which the lipid.  This is the one copy, written in that order.
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"

#include <algorithm>

void finish_implicit_lipid_binding(long long int iter, Molecule& protein, Molecule& lipid, Complex& proteinCom,
    Complex& lipidCom, int proteinIface, int lipidIface, const ForwardRxn& currRxn, double RS3D,
    std::vector<Molecule>& moleculeList, std::vector<MolTemplate>& molTemplateList,
    const Membrane& membraneObject, std::ofstream& assocDissocFile)
{
    // write temporary to real coords and clear temporary coordinates
    for (auto memMol : proteinCom.memberList) {
        moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
        for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
            moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
        moleculeList[memMol].clear_tmp_association_coords();
        moleculeList[memMol].trajStatus = TrajStatus::propagated;
        moleculeList[memMol].need_to_send = true;
    }
    // add link to surface for the molecule
    moleculeList[protein.index].linksToSurface++;
    // update complexes
    proteinCom.iLipidIndex = lipid.index; // index in molecule list of the implicit lipid.
    proteinCom.linksToSurface++;
    proteinCom.update_properties(moleculeList, molTemplateList); // recalculate the properties of the protein's complex
    proteinCom.OnSurface = true; // the complex is on membrane now.
    proteinCom.D.z = 0;
    /*Clear implicit lipid properties*/
    moleculeList[lipid.index].clear_tmp_association_coords();
    lipidCom.memberList.clear();
    lipidCom.memberList.push_back(lipid.index);
    lipidCom.update_properties(moleculeList, molTemplateList); // recalculate the properties of the lipid's complex

    // Enforce boundary conditions.  For the sphere system, many reflections may
    // be needed to move the complex back inside the sphere.
    reflect_complex_rad_rot(membraneObject, proteinCom, moleculeList, RS3D, false);

    // Binding takes a destroyable monomer out of its template's monomerList.
    remove_from_monomer_list(protein, molTemplateList);

    // update Molecule interface statuses
    protein.interfaceList[proteinIface].interaction.partnerIndex = lipid.index;
    protein.interfaceList[proteinIface].interaction.partnerIfaceIndex = lipidIface;
    if (currRxn.isReversible) {
        protein.interfaceList[proteinIface].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
    }
    protein.interfaceList[proteinIface].isBound = true;
    protein.interfaceList[proteinIface].index = currRxn.productListNew[0].absIfaceIndex;

    // The record is (protein, lipid), the order break_interaction_implicitlipid()
    // writes the matching BREAK line in.
    if (assocDissocFile.is_open()) {
        assocDissocFile << "ITR:" << iter << "," << "BOND,"
                        << molTemplateList[protein.molTypeIndex].molName << "," << protein.index << "," << proteinIface << ","
                        << molTemplateList[lipid.molTypeIndex].molName << "," << lipid.index << "," << lipidIface << std::endl;
    }

    // add to the list of bound interfaces and remove from the list of free interfaces
    protein.bndlist.push_back(proteinIface);
    protein.bndpartner.push_back(lipid.index);
    {
        size_t tmpItr { 0 };
        for (unsigned i { 0 }; i < protein.freelist.size(); ++i) {
            if (protein.freelist[i] == proteinIface)
                tmpItr = i;
        }
        protein.freelist[tmpItr] = protein.freelist.back();
        protein.freelist.pop_back();
    }

    // Set probability of this protein to zero in all reactions so it doesn't try to
    // react again but the partners still will avoid overlapping.
    zero_partner_probvec(protein, moleculeList, true);

    // Update the crossed molecule lists so that the current molecules won't avoid anything, but others will.
    proteinCom.ncross = -1;
    protein.crossings.clear();
}
