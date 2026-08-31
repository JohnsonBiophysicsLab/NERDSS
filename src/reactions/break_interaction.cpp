#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

bool break_interaction(long long int iter, size_t relIface1, size_t relIface2, Molecule& reactMol1, Molecule& reactMol2,
    const BackRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList, int ILindexMol, 
    const ForwardRxn& conjForwardRxn, bool& breakLinkComplex, double timeStep, std::ofstream& assocDissocFile)
{
    // const std::vector<ForwardRxn> &forwardRxns?

    bool cancelDissociation = false; // this can be true during loop breaking due to correction term.
    breakLinkComplex = false; // default is full dissociation
    /* Identify which complex each protein will be in following dissociation.
           Simultaneously will be checking if they'll still be in the same complex after dissociation,
           if so, they will all remain in the same single complex they started in.
     */

    /*Before checking which complex, break the existing bond stored in the bndpartner list, which is used in determine_parent_complex_IL*/
    // std::cout << " reactMol1 index: " << reactMol1.index << " size of bndpartner: " << reactMol1.bndpartner.size() << " reactMol2 index: " << reactMol2.index << " size of bndpartner: " << reactMol2.bndpartner.size() << " first partner: " << reactMol2.bndpartner[0] << " first partner of 1: " << reactMol1.bndpartner[0] << std::endl;

    // The bond is located through bndlist, which is the only key that names it
    // uniquely: an interface has at most one partner, while bndpartner can hold
    // the same molecule twice when two interfaces bind the same neighbour.  The
    // previous code erased by partner index with std::remove, which drops every
    // matching entry rather than the one being broken.
    //
    // bndlist itself is not erased here -- that happens further down, only on
    // the branch that goes through with the dissociation -- so the slot stays
    // valid for both the erase below and the restore in the cancel path.
    const size_t bondSlot1 { find_bond_slot(reactMol1, static_cast<int>(relIface1)) };
    const size_t bondSlot2 { reactMol2.isImplicitLipid
            ? reactMol2.bndlist.size()
            : find_bond_slot(reactMol2, static_cast<int>(relIface2)) };

    if (bondSlot1 < reactMol1.bndpartner.size())
        reactMol1.bndpartner.erase(reactMol1.bndpartner.begin() + bondSlot1);
    if (reactMol2.isImplicitLipid == false && bondSlot2 < reactMol2.bndpartner.size())
        reactMol2.bndpartner.erase(reactMol2.bndpartner.begin() + bondSlot2);

    // std::cout << "within break interaction " << std::endl;
    // std::cout << "After erasing the bonds: reactMol1 index: " << reactMol1.index << " size of bndpartner: " << reactMol1.bndpartner.size() << " reactMol2 index: " << reactMol2.index << " size of bndpartner: " << reactMol2.bndpartner.size() << " first partner: " << reactMol2.bndpartner[0] << " first partner of 1: " << reactMol1.bndpartner[0] << std::endl;

    bool keepSameComplex;

    // construct the new complex that will be created
    unsigned newComIndex = complexList.size();
    // std::cout << "empty complexes: " << Complex::emptyComList.size() << "\nMembers:";
    // for (auto com : Complex::emptyComList)
    //     std::cout << ' ' << com;
    // std::cout << '\n';
    if (Complex::emptyComList.size() != 0 && complexList[Complex::emptyComList.back()].isEmpty) {
        // if there is an empty complex slot, make the new Complex in it
        newComIndex = Complex::emptyComList.back();
        Complex::emptyComList.pop_back(); // remove the index from the list
    } else {
        // if we're making a new Complex, create an empty one at the end (index complexList.size())
        complexList.emplace_back();
    }
    // std::cout << "New Com Index: " << newComIndex << '\n';

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // This is commented out when merging
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    /*assign each protein in original complex c1 to one of the two new
     complexes, if the complex forms a loop, they will be put back together in
     c1, and the individual interfaces that dissociated freed.
     */
    // find the new absolute interfaces
    // int absIface1 { -1 };
    // int absIface2 { -1 };
    // if (currRxn.isSymmetric) {
    //     absIface1 = currRxn.productListNew[0].absIfaceIndex;
    //     absIface2 = currRxn.productListNew[1].absIfaceIndex;
    // } else {
    //     /*Here, if both proteins are the same protein, same interface, but distinct states, need to correct for that.*/
    //     // std::cout << "State of mol1: " << reactMol1.interfaceList[relIface1].stateIden << " of mol2: " << reactMol2.interfaceList[relIface2].stateIden << "\n";
    //     // std::cout << " Product requires state: " << currRxn.productListNew[0].requiresState << " " << currRxn.productListNew[1].requiresState << std::endl;
    //     if (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex
    //         && relIface1 == currRxn.productListNew[0].relIfaceIndex && currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
    //         //matched protein, interface, and state of product[0] to mol1
    //         absIface1 = currRxn.productListNew[0].absIfaceIndex;
    //         absIface2 = currRxn.productListNew[1].absIfaceIndex;

    //     } else if (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex
    //         && relIface1 == currRxn.productListNew[1].relIfaceIndex && currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
    //         //matched protein, interface and state of product[1] to mol1
    //         absIface2 = currRxn.productListNew[0].absIfaceIndex;
    //         absIface1 = currRxn.productListNew[1].absIfaceIndex;
    //     } else {
    //         std::cout << " IN BREAK INTERACTION, DID NOT MATCH protein, interface, and state of a product to Mol1 " << reactMol1.index << std::endl;
    //         std::cout << "reactMol1.molTypeIndex: " << reactMol1.molTypeIndex << std::endl;
    //         std::cout << "currRxn.productListNew[0].molTypeIndex: " << currRxn.productListNew[0].molTypeIndex << std::endl;
    //         std::cout << "currRxn.productListNew[1].molTypeIndex: " << currRxn.productListNew[1].molTypeIndex << std::endl;
    //         std::cout << "reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex: " << (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex) << std::endl;
    //         std::cout << "reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex: " << (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex) << std::endl;
    //         std::cout << "relIface1: " << relIface1 << std::endl;
    //         std::cout << "currRxn.productListNew[0].relIfaceIndex: " << currRxn.productListNew[0].relIfaceIndex << std::endl;
    //         std::cout << "currRxn.productListNew[1].relIfaceIndex: " << currRxn.productListNew[1].relIfaceIndex << std::endl;
    //         std::cout << "relIface1 == currRxn.productListNew[0].relIfaceIndex: " << (relIface1 == currRxn.productListNew[0].relIfaceIndex) << std::endl;
    //         std::cout << "relIface1 == currRxn.productListNew[1].relIfaceIndex: " << (relIface1 == currRxn.productListNew[1].relIfaceIndex) << std::endl;
    //         std::cout << "currRxn.productListNew[0].requiresState: " << currRxn.productListNew[0].requiresState << std::endl;
    //         std::cout << "currRxn.productListNew[1].requiresState: " << currRxn.productListNew[1].requiresState << std::endl;
    //         std::cout << "reactMol1.interfaceList[relIface1].stateIden: " << reactMol1.interfaceList[relIface1].stateIden << std::endl;
    //         std::cout << "currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden: " << (currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden) << std::endl;
    //         std::cout << "currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden: " << (currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) << std::endl;
    //         exit(1);
    //     }
    // }

    // reactMol1.interfaceList[relIface1].index = absIface1;
    // reactMol1.interfaceList[relIface1].interaction.clear();
    // The commented-out block that stood here held a 2023-01-04 attempt at the
    // same fix -- locate through bndlist, erase all lists at that position --
    // which was disabled during a merge and never restored.  It is implemented
    // above, in erase_bond() and the bondSlot handling, so the dead copy is
    // gone.  Its own reactMol2 half called std::distance on an iterator erase
    // had already invalidated, so it could not have been restored as written.
    // // check if they'll still be in the same complex after dissociation,
    // // if not, move them slightly apart
    // // if so, change complex identities back and don't move
    // bool keepSameComplex;
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (reactMol2.isImplicitLipid == false)
        keepSameComplex = determine_parent_complex_IL(reactMol1.index, reactMol2.index, newComIndex, moleculeList, complexList, ILindexMol);
    else
        keepSameComplex = false;

    /*if they are within the same complex, use the loop correction rate to decide whether to go forward with the dissociation*/
    if (keepSameComplex) {
	  //std::cout << " BREAKING A LOOP, check COrrectino term " << std::endl;
        double rateForward = conjForwardRxn.rateList[0].rate;
        double largestRate = -1;
        for (auto oneForwardRate : conjForwardRxn.rateList) {
            if (oneForwardRate.rate > largestRate)
                largestRate = oneForwardRate.rate;
            if (std::abs(oneForwardRate.rate - rateForward) / std::abs(rateForward) > 1E-2) {
                /*There is more than one Forward reaction rate for the conjugate reaction, conditional on other states of molecule
                  TODO: We will need to figure out what the exact rate is based on the binding states of both reactants that are now dissociating.
                */
                std::cerr << " BREAKING A LOOP, BUT REBINDING RATE HAS MULTIPLE VALUES, USING LARGEST ONE TO CORRECT FOR LARGE k*dt! " << std::endl;
            }
        }
        // std::cout << " Performing Dissociation on a CLOSED LOOP, not creating a new complex ! " << std::endl;
        double c0_nm3 = 0.602; // standard state 1M in units of /nm^3
        double coop = conjForwardRxn.loopCoopFactor;
        double rateClose = largestRate * c0_nm3 * coop;
        // rateClose will be in units of /us (due to ka units), so no need for 1E-6 factor, since timeStep has
        // units of us

        double poisson = timeStep * rateClose;
        double correctionRatio{(1 - exp(-poisson)) / poisson};
		
	    //std::cout <<"Correction Ratio: "<<correctionRatio<<std::endl;

        if (1.0 * rand_gsl() > correctionRatio) {
            /*Cancel the dissociation!!*/
            cancelDissociation = true;
        }
    }

    if (cancelDissociation) {
        /*cancel the dissociation, move both interfaces back into the bndpartner list, which is the only change we've made so far.*/
        // Reinserted at the slot each came from, not appended.  push_back left
        // bndpartner permuted against bndlist whenever the cancelled bond was
        // not the molecule's last, which is measurably what happened on
        // homoTrimer; see docs/bond_bookkeeping_defects.md.
        reactMol1.bndpartner.insert(
            reactMol1.bndpartner.begin() + std::min(bondSlot1, reactMol1.bndpartner.size()), reactMol2.index);
        reactMol2.bndpartner.insert(
            reactMol2.bndpartner.begin() + std::min(bondSlot2, reactMol2.bndpartner.size()), reactMol1.index);

        // reset empty complexList
        if (newComIndex + 1 == complexList.size())
            complexList.pop_back();
        else
            Complex::emptyComList.push_back(newComIndex);

        return true;
    } else {
        /*Continue on with dissociation as before.*/
	  //std::cout << "continue with dissociation " << std::endl;
        // record the previous lastNumberUpdateItrEachMol & numEachMol
        std::vector<int> numEachMolPrevious {};
        std::vector<long long int> lastNumberUpdateItrEachMolPrevious {};
        for (unsigned index = 0; index < molTemplateList.size(); index++) {
            numEachMolPrevious.emplace_back(complexList[reactMol1.myComIndex].numEachMol[index]);
            lastNumberUpdateItrEachMolPrevious.emplace_back(complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index]);
        }

        /*assign each protein in original complex c1 to one of the two new complexes,
         if the complex forms a loop, they will be put back together in c1, and the
         individual interfaces that dissociated freed.
         */
        // A block disabled "when merging" stood here, calling correct_structure()
        // on fiber complexes of two members.  That function is deleted: nothing
        // called it, it read bndRxnList[0] which was never maintained, and it
        // displaced a local copy of the molecule that it then discarded.  See
        // docs/bond_bookkeeping_defects.md.
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
                // matched protein, interface, and state of product[0] to mol1
                absIface1 = currRxn.productListNew[0].absIfaceIndex;
                absIface2 = currRxn.productListNew[1].absIfaceIndex;

            } else if (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex
                && relIface1 == currRxn.productListNew[1].relIfaceIndex && currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
                // matched protein, interface and state of product[1] to mol1
                absIface2 = currRxn.productListNew[0].absIfaceIndex;
                absIface1 = currRxn.productListNew[1].absIfaceIndex;
            } else {
                std::cout << " IN BREAK INTERACTION, DID NOT MATCH protein, interface, and state of a product to Mol1 " << reactMol1.index << std::endl;
                std::cout << "reactMol1.molTypeIndex: " << reactMol1.molTypeIndex << std::endl;
                std::cout << "currRxn.productListNew[0].molTypeIndex: " << currRxn.productListNew[0].molTypeIndex << std::endl;
                std::cout << "currRxn.productListNew[1].molTypeIndex: " << currRxn.productListNew[1].molTypeIndex << std::endl;
                std::cout << "reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex: " << (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex) << std::endl;
                std::cout << "reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex: " << (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex) << std::endl;
                std::cout << "relIface1: " << relIface1 << std::endl;
                std::cout << "currRxn.productListNew[0].relIfaceIndex: " << currRxn.productListNew[0].relIfaceIndex << std::endl;
                std::cout << "currRxn.productListNew[1].relIfaceIndex: " << currRxn.productListNew[1].relIfaceIndex << std::endl;
                std::cout << "relIface1 == currRxn.productListNew[0].relIfaceIndex: " << (relIface1 == currRxn.productListNew[0].relIfaceIndex) << std::endl;
                std::cout << "relIface1 == currRxn.productListNew[1].relIfaceIndex: " << (relIface1 == currRxn.productListNew[1].relIfaceIndex) << std::endl;
                std::cout << "currRxn.productListNew[0].requiresState: " << currRxn.productListNew[0].requiresState << std::endl;
                std::cout << "currRxn.productListNew[1].requiresState: " << currRxn.productListNew[1].requiresState << std::endl;
                std::cout << "reactMol1.interfaceList[relIface1].stateIden: " << reactMol1.interfaceList[relIface1].stateIden << std::endl;
                std::cout << "currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden: " << (currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden) << std::endl;
                std::cout << "currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden: " << (currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) << std::endl;
                exit(1);
            }
        }

        reactMol1.interfaceList[relIface1].index = absIface1;
        reactMol1.interfaceList[relIface1].interaction.clear();
        reactMol1.interfaceList[relIface1].isBound = false;
        if (reactMol2.isImplicitLipid == false) {
            reactMol2.interfaceList[relIface2].index = absIface2;
            reactMol2.interfaceList[relIface2].interaction.clear();
            reactMol2.interfaceList[relIface2].isBound = false;
        }
        if (assocDissocFile.is_open()) {
            assocDissocFile << "ITR:" << iter << "," << "BREAK," 
            << molTemplateList[reactMol1.molTypeIndex].molName << "," << reactMol1.index << "," << relIface1 << "," 
            << molTemplateList[reactMol2.molTypeIndex].molName << "," << reactMol2.index << "," << relIface2 << std::endl;
        }
        // Add these protein into the bimolecular association list
        reactMol1.freelist.push_back(relIface1);
        // bndpartner lost its entry at bondSlot1 at the top of this function;
        // this removes the matching bndlist entry, keeping the two in step.
        if (bondSlot1 < reactMol1.bndlist.size())
            reactMol1.bndlist.erase(reactMol1.bndlist.begin() + bondSlot1);
        if (reactMol2.isImplicitLipid == false) {
            reactMol2.freelist.push_back(relIface2);
            if (bondSlot2 < reactMol2.bndlist.size())
                reactMol2.bndlist.erase(reactMol2.bndlist.begin() + bondSlot2);
        }

        if (!keepSameComplex) {
            /*continue on with the dissociation that creates two complexes*/
            complexList[reactMol1.myComIndex].update_properties(moleculeList, molTemplateList);

            // add a lifetime to the previous larger cluster
            // compare current cluster size with the previous ones
            for (unsigned index = 0; index < molTemplateList.size(); index++) {
                if (molTemplateList[index].countTransition == true && complexList[reactMol1.myComIndex].numEachMol[index] < numEachMolPrevious[index]) {
                    // shrink
                    molTemplateList[index].lifeTime[numEachMolPrevious[index] - 1].emplace_back((iter - lastNumberUpdateItrEachMolPrevious[index]) * Parameters::dt / 1E6);
                    complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index] = iter;
                }
            }

            double small = 1E-9;
            double dx;
            if (reactMol2.isImplicitLipid == false)
                dx = reactMol1.interfaceList[relIface1].coord.x - reactMol2.interfaceList[relIface2].coord.x;
            else
                dx = 0.01;
            dx = (dx > 0) ? small : -small; // fix for precision issues
            Vec3D chg1 { dx, 0, 0 };

            /*Update the positions of each protein. Then calculate the COMs of each
             of the complexes. Then calculated the radius (uses the Complexlist COM).
             update rotational diffusion, translational diffusion was already updated.
             */
            complexList[reactMol1.myComIndex].translate(chg1, moleculeList);
            complexList[reactMol1.myComIndex].update_properties(moleculeList, molTemplateList);
            if (reactMol2.isImplicitLipid == false) {
                complexList[reactMol2.myComIndex].update_properties(moleculeList, molTemplateList);

                // update lastNumberUpdateItrEachMol for new complex
                complexList[reactMol2.myComIndex].lastNumberUpdateItrEachMol.resize(molTemplateList.size());
                for (unsigned index = 0; index < molTemplateList.size(); index++) {
                    complexList[reactMol2.myComIndex].lastNumberUpdateItrEachMol[index] = iter;
                }

                complexList[reactMol2.myComIndex].isEmpty = false;
            }

            // update transition matrix
            // compare current cluster size with the previous ones
            for (unsigned index = 0; index < molTemplateList.size(); index++) {
                if (molTemplateList[index].countTransition == true && complexList[reactMol1.myComIndex].numEachMol[index] < numEachMolPrevious[index]) {
                    // shrink
                    if (numEachMolPrevious[index] - 1 >= 0)
                        molTemplateList[index].transitionMatrix[numEachMolPrevious[index] - 1][numEachMolPrevious[index] - 1] += iter - Parameters::lastUpdateTransition[index] - 1;
                    if (complexList[reactMol1.myComIndex].numEachMol[index] - 1 >= 0 && numEachMolPrevious[index] - 1 >= 0)
                        molTemplateList[index].transitionMatrix[complexList[reactMol1.myComIndex].numEachMol[index] - 1][numEachMolPrevious[index] - 1] += 1;
                    if (complexList[reactMol2.myComIndex].numEachMol[index] - 1 >= 0 && numEachMolPrevious[index] - 1 >= 0)
                        molTemplateList[index].transitionMatrix[complexList[reactMol2.myComIndex].numEachMol[index] - 1][numEachMolPrevious[index] - 1] += 1;

                    // update diagonal elements for unchanged complexes
                    for (unsigned indexCom = 0; indexCom < complexList.size(); indexCom++) {
                        if ((indexCom != reactMol1.myComIndex) && (indexCom != reactMol2.myComIndex)) {
                            if (complexList[indexCom].numEachMol[index] - 1 >= 0) {
                                molTemplateList[index].transitionMatrix[complexList[indexCom].numEachMol[index] - 1][complexList[indexCom].numEachMol[index] - 1] += iter - Parameters::lastUpdateTransition[index];
                            }
                        }
                    }

                    // update lastUpdateTransition
                    Parameters::lastUpdateTransition[index] = iter;
                }
            }

            ++Complex::numberOfComplexes;
            complexList[newComIndex].id = Complex::maxID++;
            for (auto i : complexList[newComIndex].memberList) {
                moleculeList[i].complexId = complexList[newComIndex].id;
            }
            // std::cout << "new number of complexes: " << Complex::numberOfComplexes << '\n';
            // std::cout << "New Complexes:\n";
            // std::cout << "\tComplex " << reactMol1.myComIndex << " of " << complexList[reactMol1.myComIndex].memberList.size() << " molecules.\n";
            // std::cout << "\t\tMembers:";
            // for (auto memMol : complexList[reactMol1.myComIndex].memberList)
            //     std::cout << ' ' << memMol;
            // std::cout << '\n';
            // std::cout << "\tComplex " << newComIndex << " of " << complexList[newComIndex].memberList.size() << " molecules.\n";
            // std::cout << "\t\tMembers:";
            // for (auto memMol : complexList[newComIndex].memberList)
            //     std::cout << ' ' << memMol;
            // std::cout << '\n';
        } else {
            /*All proteins remain in complex c1, dissociation
                  will break the product state of the two proteins that dissociated but here they
                  are linked in a closed loop so it will not create a new complex.
                  positions don't change
            */
            // reset empty complexList
            breakLinkComplex = true;
            if (newComIndex + 1 == complexList.size())
                complexList.pop_back();
            else
                Complex::emptyComList.push_back(newComIndex);
        }
        //------------------------START UPDATE MONOMERLIST-------------------------
        // update oneTemp.monomerList when oneTemp.canDestroy is true and mol is monomer, add to monomerList if new monomer produced
        // reactMol1
        {
            Molecule& oneMol { reactMol1 };
            MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
            bool isMonomer { oneMol.bndpartner.empty() };
            bool canDestroy { oneTemp.canDestroy };
            if (isMonomer && canDestroy) {
                // add to monomerList
                oneTemp.monomerList.emplace_back(oneMol.index);
            }

            // std::cout << "For mol " << oneMol.index << ": "
            //           << "canDestory is " << oneTemp.canDestroy << "\t"
            //           << "isMonomer is " << isMonomer << std::endl;
            // std::cout << "Now the monomerList is: ";
            // for (auto one : oneTemp.monomerList) {
            //     std::cout << one << "\t";
            // }
            // std::cout << std::endl;
        }
        // reactMol2
        {
            Molecule& oneMol { reactMol2 };
            MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
            bool isMonomer { oneMol.bndpartner.empty() };
            bool canDestroy { oneTemp.canDestroy };
            if (isMonomer && canDestroy) {
                // add to monomerList
                oneTemp.monomerList.emplace_back(oneMol.index);
            }
            // std::cout << "For mol " << oneMol.index << ": "
            //           << "canDestory is " << oneTemp.canDestroy << "\t"
            //           << "isMonomer is " << isMonomer << std::endl;
            // std::cout << "Now the monomerList is: ";
            // for (auto one : oneTemp.monomerList) {
            //     std::cout << one << "\t";
            // }
            // std::cout << std::endl;
        }
        //------------------------END UPDATE MONOMERLIST---------------------------
    } // end if to actually perform dissociation, so not cancelled
    return cancelDissociation;
}
