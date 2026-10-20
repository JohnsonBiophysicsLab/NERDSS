#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include <cmath>
#include <iomanip>

// Either 1 or 2 is the implicit lipid.  Only works for the sphere system; the
// protein's complex binds to the surface by 3D->2D or 2D->2D.
//
void associate_implicitlipid_sphere(
    long long int iter, int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, const Parameters& params,
    ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList,
    copyCounters& counterArrays, std::vector<Complex>& complexList,
    Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, std::ofstream& assocDissocFile)
{
    double RS3D { -1.0 };
    for (int RS3Di = 0; RS3Di < 100; RS3Di++) {
        if ((std::abs(membraneObject.RS3Dvect[RS3Di] - currRxn.bindRadius) < 1E-15) && (std::abs(membraneObject.RS3Dvect[RS3Di + 100] - currRxn.rateList[0].rate) < 1E-15) && std::abs(membraneObject.RS3Dvect[RS3Di + 200] - (1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.x + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.x) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.y + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.y) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.z + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.z))) < 1E-15) {
            RS3D = membraneObject.RS3Dvect[RS3Di + 300];
            break;
        }
    }

    // The reaction may list the implicit lipid first or second.  Everything
    // below runs on the pair as (protein, lipid), whichever way it was listed:
    // the protein's complex is the one that moves to sigma, sticks to the
    // sphere and is reflected, the lipid's is the one the membrane orientation
    // is read from, and the rotations pivot on whichever interface is passed
    // first, so running them in listed order gives a different floating-point
    // sequence per listing.  The angles belong to the listed order, so they
    // are sorted into roles along with the molecules.
    const bool lipidListedSecond { reactMol2.isImplicitLipid };
    Molecule& protein { lipidListedSecond ? reactMol1 : reactMol2 };
    Molecule& lipid { lipidListedSecond ? reactMol2 : reactMol1 };
    Complex& proteinCom { lipidListedSecond ? reactCom1 : reactCom2 };
    Complex& lipidCom { lipidListedSecond ? reactCom2 : reactCom1 };
    int proteinIface { lipidListedSecond ? ifaceIndex1 : ifaceIndex2 };
    int lipidIface { lipidListedSecond ? ifaceIndex2 : ifaceIndex1 };
    const double proteinTheta { lipidListedSecond ? currRxn.assocAngles.theta1 : currRxn.assocAngles.theta2 };
    const double lipidTheta { lipidListedSecond ? currRxn.assocAngles.theta2 : currRxn.assocAngles.theta1 };
    const double proteinPhi { lipidListedSecond ? currRxn.assocAngles.phi1 : currRxn.assocAngles.phi2 };
    const double lipidPhi { lipidListedSecond ? currRxn.assocAngles.phi2 : currRxn.assocAngles.phi1 };
    const Vec3D& proteinNorm { lipidListedSecond ? currRxn.norm1 : currRxn.norm2 };
    const Vec3D& lipidNorm { lipidListedSecond ? currRxn.norm2 : currRxn.norm1 };

    // set the implicit lipid's temporary position according to the protein.
    lipid.create_position_implicit_lipid(protein, lipidIface, currRxn.bindRadius, membraneObject);
    lipidCom.comCoord = lipid.comCoord;

    // set up temporary coordinates
    for (auto& memMol : proteinCom.memberList)
        moleculeList[memMol].set_tmp_association_coords();
    for (auto& mol : lipidCom.memberList)
        moleculeList[mol].set_tmp_association_coords();
    // std::cout << " In Associate binidng to sphere surface " << std::endl;

    // create references to reacting interfaces
    Vec3D& proteinIfaceCrd = protein.tmpICoords[proteinIface];
    Vec3D& lipidIfaceCrd = lipid.tmpICoords[lipidIface];

    // orientation corrections for membrane bound components
    bool isOnMembrane = false;
    bool transitionToSurface = false;

    /* MOVE PROTEIN TO SIGMA */
    {
        Vec3D sigma { proteinIfaceCrd - lipidIfaceCrd };
        Vec3D transVec1 { 0.0, 0.0, 0.0 };
        double displaceFrac {};
        double sigmaMag;
        // if both in 2D, ignore the z-component
        // if (proteinCom.D.z < 1E-14) {
        if (proteinCom.OnSurface) {
            isOnMembrane = true;
        } else { // note not on the membrane
            transitionToSurface = true;
            // std::cout << "TRANSITIONING FROM 3D->2D " << std::endl;
            sigmaMag = sqrt((sigma.x * sigma.x) + (sigma.y * sigma.y) + (sigma.z * sigma.z));
            displaceFrac = (sigmaMag - currRxn.bindRadius) / sigmaMag;
            transVec1.x = -sigma.x * displaceFrac;
            transVec1.y = -sigma.y * displaceFrac;
            transVec1.z = -sigma.z * displaceFrac;
        }
        // update the temporary coordinates
        for (auto& mp : proteinCom.memberList)
            moleculeList[mp].update_association_coords(transVec1);
        // std::cout << "Position after pushed to sigma: " << std::endl;
    } // Move protein to sigma
    // now rotate the protein's complex
    // std::cout << "ImplicitLipid model, need to rotate the complex, and make it stick on the sphere.  \n";
    // std::cout << "Pre-MEMBRANE ROT, COORDS: " << std::endl;
    if (molTemplateList[protein.molTypeIndex].isPoint) {
        // If the binding protein is a point, no orientations to specify
        // std::cout << " Move two point particles to contact along current separation vector, NO ORIENTATION \n";
    } else {
        // THETA
        // std::cout << std::setw(8) << std::setfill('-') << ' ' << std::endl
        //           << "THETA 1" << std::endl
        //           << std::setw(8) << ' ' << std::setfill(' ') << std::endl;
        theta_rotation(proteinIfaceCrd, lipidIfaceCrd, protein, lipid,
            proteinTheta, proteinCom, lipidCom, moleculeList,
            params.numerics.associationAngles);
        // std::cout << std::setw(30) << std::setfill('-') << ' '
        //           << std::setfill(' ') << std::endl;
        // std::cout << "THETA 2" << std::endl
        //           << std::setw(8) << std::setfill('-') << ' ' << std::setfill(' ')
        //           << std::endl;
        theta_rotation(lipidIfaceCrd, proteinIfaceCrd, lipid, protein,
            lipidTheta, lipidCom, proteinCom, moleculeList,
            params.numerics.associationAngles);

        // OMEGA
        // if protein has theta M_PI, uses protein norm instead of com_iface vector
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "OMEGA" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
        if (!std::isnan(currRxn.assocAngles.omega)) {
            omega_rotation(proteinIfaceCrd, lipidIfaceCrd, lipidIface, protein,
                lipid, proteinCom, lipidCom,
                currRxn.assocAngles.omega, currRxn, moleculeList,
                molTemplateList, params.numerics.associationAngles);
        } //else
        // std::cout << "P1 or P2 is a rod-type protein, no dihedral for associated complex." << std::endl;

        // PHI
        // PHI 1
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "PHI 1" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

        if (!std::isnan(proteinPhi)) {
            phi_rotation(proteinIfaceCrd, lipidIfaceCrd, lipidIface, protein, lipid, proteinCom, lipidCom,
                proteinNorm, proteinPhi, currRxn, moleculeList, molTemplateList,
                params.numerics.associationAngles);
        } //else
        // std::cout << "P1 has no valid phi angle." << std::endl;

        // PHI 2
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "PHI 2" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

        if (!std::isnan(lipidPhi)) {
            phi_rotation(lipidIfaceCrd, proteinIfaceCrd, proteinIface, lipid, protein, lipidCom, proteinCom,
                lipidNorm, lipidPhi, currRxn, moleculeList, molTemplateList,
                params.numerics.associationAngles);
        } //else
        // std::cout << "P2 has no valid phi angle." << std::endl;
    } // end of if points.

    // need to adjust Lipid's orientation, make it verticle to sphere
    {
        Molecule memProtein;
        Molecule Lipid;
        set_memProtein_sphere(lipidCom, memProtein, moleculeList, membraneObject);
        find_lipid_sphere(lipidCom, Lipid, moleculeList, membraneObject);
        Quat memRot = save_mem_orientation(memProtein, Lipid, molTemplateList[Lipid.molTypeIndex]);
        Vec3D pivot = Lipid.tmpComCoord;
        /*rotate the molecules and their complexes.*/
        rotate(pivot, memRot, proteinCom, moleculeList);
        rotate(pivot, memRot, lipidCom, moleculeList);
        //make sure the protein's interface is sticking on sphere.
        double lamda = (membraneObject.sphereR - proteinIfaceCrd.length()) / proteinIfaceCrd.length();
        Vec3D dtrans = Vec3D { lamda * proteinIfaceCrd };
        // std::cout << " In ImplicitLipid model, protein interface is off spherical membrane, shift up by: " << dtrans.x << " " << dtrans.y << " " << dtrans.z
        //           << std::endl;
        // update the temporary coordinates of the protein's complex
        for (auto& mp : proteinCom.memberList)
            moleculeList[mp].update_association_coords(dtrans);
    }

    // std::cout << "Pos-MEMBRANE ROT, COORDS: " << std::endl;

    // std::cout << "OVERLAP CHECK  AND REFLECT OFF SPHERE: " << std::endl;
    // Reflect off the sphere.
    std::array<double, 3> traj;
    for (int mm = 0; mm < 3; mm++)
        traj[mm] = 0;

    /*This needs to evaluate the traj update, based on it initially being zero.
    And here, it should be called based on the tmpCoords, not the full
    coordinates. also requires updating the COM of this temporary new position
  */
    update_complex_tmp_com_crds(proteinCom, moleculeList);
    proteinCom.tmpOnSurface = true;
    reflect_traj_tmp_crds(params, moleculeList, proteinCom, traj, membraneObject, RS3D, false); // uses tmpCoords to calculate traj.

    if (std::abs(traj[0] + traj[1] + traj[2]) > 1E-14) {
        // update the temporary coordinates of the protein's complex
        Vec3D vtraj { traj[0], traj[1], traj[2] };
        for (auto& mp : proteinCom.memberList)
            moleculeList[mp].update_association_coords(vtraj);

        // std::cout << "CRDS after reflecting off of the SPHERE by " << traj[0] << ' '
        //           << traj[1] << ' ' << traj[2] << std::endl;
    }
    /* CHECKS */
    bool cancelAssoc { false };

    check_if_spans_sphere(cancelAssoc, params, proteinCom, lipidCom, moleculeList, membraneObject, membraneObject.sphereR);
    if (cancelAssoc == true)
        counterArrays.nCancelSpanBox++;
    run_association_rejection_checks(cancelAssoc, proteinCom, lipidCom, moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns, isOnMembrane,
        transitionToSurface, counterArrays);

    if (cancelAssoc == true) {
        // std::cout << "Canceling association, returning complexes to original state.\n";
        proteinCom.tmpOnSurface = false;
        for (auto memMol : proteinCom.memberList)
            moleculeList[memMol].clear_tmp_association_coords();
        for (auto memMol : lipidCom.memberList)
            moleculeList[memMol].clear_tmp_association_coords();
        return;
        // EXIT FROM THIS ROUTINE.
    }
    counterArrays.nAssocSuccess++; //keep track of total number of successful association moves.
    /*Keep track of the sizes of complexes that associated*/
    track_association_events(proteinCom, lipidCom, transitionToSurface, isOnMembrane, counterArrays);

    finish_implicit_lipid_binding(iter, protein, lipid, proteinCom, lipidCom, proteinIface, lipidIface, currRxn,
        RS3D, moleculeList, molTemplateList, membraneObject, assocDissocFile);

    // update the number of bound species
    update_Nboundpairs(reactMol1.molTypeIndex, reactMol2.molTypeIndex, 1, params,
        counterArrays);

    //Update species copy numbers
    counterArrays.copyNumSpecies[currRxn.reactantListNew[0].absIfaceIndex] -= 1; // decrement ifaceIndex1
    counterArrays.copyNumSpecies[currRxn.reactantListNew[1].absIfaceIndex] -= 1; // decrement ifaceIndex2
    //counterArrays.copyNumSpecies[reactMol1.interfaceList[ifaceIndex1].index] -= 1; // decrement ifaceIndex1
    //counterArrays.copyNumSpecies[reactMol2.interfaceList[ifaceIndex2].index] -= 1; // decrement ifaceIndex2
    counterArrays.copyNumSpecies[currRxn.productListNew[0].absIfaceIndex] += 1; // increment product state

    //-------------------------UPDATE BINDPAIRLIST----------------------------
    // auto molIndex { reactMol1.index };
    // auto comIndex { reactMol1.myComIndex };
    // if (molTemplateList[reactMol1.molTypeIndex].isImplicitLipid == true) {
    //     molIndex = reactMol2.index;
    //     comIndex = reactMol2.myComIndex;
    // }

    // if (counterArrays.canDissociate[currRxn.productListNew[0].absIfaceIndex]) {
    //     if (complexList[comIndex].linksToSurface == 1) {
    //         counterArrays.bindPairListIL3D[currRxn.productListNew[0].absIfaceIndex].emplace_back(molIndex);
    //     } else {
    //         counterArrays.bindPairListIL2D[currRxn.productListNew[0].absIfaceIndex].emplace_back(molIndex);
    //     }
    // }

    // if (complexList[comIndex].linksToSurface == 2) {
    //     // one interface in the complex bind to surface changing the linksToSurface of complex from one to two will impact the bindPairListIL of the other linking interface
    //     // need to remove the impacted interface's mol index from the bindPairListIL3D, then add it to bindPairListIL2D
    //     // try to find out the other linking interface
    //     int theOtherLinkingIfaceAbsIndex { -1 };
    //     if (moleculeList[molIndex].linksToSurface == 2) { // the other interface is in the same molecule
    //         // loop over the interfaceList to find the other interface
    //         for (auto oneIface : moleculeList[molIndex].interfaceList) {
    //             if (oneIface.interaction.partnerIndex == 0 && oneIface.index != currRxn.productListNew[0].absIfaceIndex) {
    //                 theOtherLinkingIfaceAbsIndex = oneIface.index;
    //                 if (counterArrays.canDissociate[theOtherLinkingIfaceAbsIndex]) {
    //                     counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].erase(std::find_if(counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].begin(), counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].end(), [&](const size_t& mol) { return mol == molIndex; }));
    //                     counterArrays.bindPairListIL2D[theOtherLinkingIfaceAbsIndex].emplace_back(molIndex);
    //                 }
    //             }
    //         }
    //     } else { // the other interface is in another molecule
    //         for (auto oneMember : complexList[comIndex].memberList) {
    //             if ((oneMember != molIndex) && (moleculeList[oneMember].linksToSurface == 1)) {
    //                 // loop over the interfaceList to find the other interface
    //                 for (auto oneIface : moleculeList[oneMember].interfaceList) {
    //                     if (oneIface.interaction.partnerIndex == 0) {
    //                         theOtherLinkingIfaceAbsIndex = oneIface.index;
    //                         if (counterArrays.canDissociate[theOtherLinkingIfaceAbsIndex]) {
    //                             counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].erase(std::find_if(counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].begin(), counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].end(), [&](const size_t& mol) { return mol == oneMember; }));
    //                             counterArrays.bindPairListIL2D[theOtherLinkingIfaceAbsIndex].emplace_back(oneMember);
    //                         }
    //                     }
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    // }
    //-----------------------END UPDATE BINDPAIRLIST--------------------------

    //Update free number of lipids of IL for each state
    int relStateIndex { implicit_lipid_state_index(
        implicit_lipid_state_list(moleculeList, molTemplateList, membraneObject),
        implicit_lipid_iface(currRxn.reactantListNew, molTemplateList)) };
    membraneObject.numberOfFreeLipidsEachState[relStateIndex] -= 1;

    // TODO: Insert species tracking here
    if (currRxn.isObserved) {
        auto obsItr = observablesList.find(currRxn.observeLabel);
        if (obsItr != observablesList.end())
            ++obsItr->second;
    }
    //    reactCom1.display();
}
