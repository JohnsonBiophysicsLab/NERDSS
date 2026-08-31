#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "tracing.hpp"
#include "reactions/shared_reaction_functions.hpp"

void break_interaction_implicitlipid(long long int iter, size_t relIface1, size_t relIface2, Molecule& reactMol1, Molecule& reactMol2,
    const BackRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList, std::ofstream& assocDissocFile)
{
    // TRACE();
    /*Find the new absIface for mol1*/
    int absIface1 { -1 };
    int absIface2 { -1 };
    if (currRxn.isSymmetric) {
        absIface1 = currRxn.productListNew[0].absIfaceIndex;
        absIface2 = currRxn.productListNew[1].absIfaceIndex;
    } else {
        /*Here, if both proteins are the same protein, same interface, but distinct states, need to correct for that.*/
        // std::cout << "State of mol1: " << reactMol1.interfaceList[relIface1].stateIden << " of mol2: " << reactMol2.interfaceList[relIface2].stateIden << "\n";
        // std::cout << " Product requires state: " << currRxn.productListNew[0].requiresState << " " << currRxn.productListNew[1].requiresState << std::endl;
        if (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex
            && relIface1 == currRxn.productListNew[0].relIfaceIndex && currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
            //matched protein, interface, and state of product[0] to mol1
            absIface1 = currRxn.productListNew[0].absIfaceIndex;
            absIface2 = currRxn.productListNew[1].absIfaceIndex;

        } else if (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex
            && relIface1 == currRxn.productListNew[1].relIfaceIndex && currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
            //matched protein, interface and state of product[1] to mol1
            absIface2 = currRxn.productListNew[0].absIfaceIndex;
            absIface1 = currRxn.productListNew[1].absIfaceIndex;
        } else {
            std::cout << " IN BREAK INTERACTION, DID NOT MATCH protein, interface, and state of a product to Mol1 " << reactMol1.index << std::endl;
            exit(1);
        }
    }

    reactMol1.interfaceList[relIface1].interaction.clear();
    reactMol1.interfaceList[relIface1].isBound = false;
    reactMol1.interfaceList[relIface1].index = absIface1;
    if (assocDissocFile.is_open()) {
        assocDissocFile << "ITR:" << iter << "," << "BREAK," 
        << molTemplateList[reactMol1.molTypeIndex].molName << "," << reactMol1.index << "," << relIface1 << "," 
        << molTemplateList[reactMol2.molTypeIndex].molName << "," << reactMol2.index << "," << relIface2 << std::endl;
    }
    //Add these protein into the bimolecular association list
    reactMol1.freelist.push_back(relIface1);
    // Was two independent find_if + erase calls, neither guarded: if the
    // interface was not bound, find_if returns end() and erase(end()) is
    // undefined.  They also selected by different keys, so they could remove
    // entries describing different bonds.  erase_bond() locates once through
    // bndlist and removes the matching pair, or does nothing.
    erase_bond(reactMol1, static_cast<int>(relIface1));

    //------------------------START UPDATE MONOMERLIST-------------------------
    // update oneTemp.monomerList when oneTemp.canDestroy is true and mol is monomer, add to monomerList if new monomer produced
    // reactMol1
    {
        Molecule& oneMol { reactMol1 };
        MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
        bool isMonomer { oneMol.bndpartner.empty() };
        bool canDestroy { oneTemp.canDestroy };
        if (isMonomer && canDestroy) {
            //add to monomerList
            oneTemp.monomerList.emplace_back(oneMol.index);
        }
    }
    //------------------------END UPDATE MONOMERLIST---------------------------
}
