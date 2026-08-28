#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <algorithm>
#include <vector>

void check_bimolecular_reactions(int pro1Index, int pro2Index, int simItr, double* tableIDs, unsigned& DDTableIndex,
    const Parameters& params, std::vector<gsl_matrix*>& normMatrices, std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject)
{
    // TRACE();
    //  int track1 = 143;
    //int track2= 182;
    int pro1MolType = moleculeList[pro1Index].molTypeIndex;
    /* if(pro1Index== track1)
       std::cout <<"In check bimolecular Reaction for protein track1 !"<<track1<<" to  "<<pro2Index<<"mytype: "<<pro1MolType<<" pro2type: "<<moleculeList[pro2Index].molTypeIndex<<std::endl;
     if(pro1Index== track2)
       std::cout <<"In check bimolecular Reaction for protein track2 !" <<track2<< " to "<<pro2Index<<" mytype: "<<pro1MolType<<" pro2type: "<<moleculeList[pro2Index].molTypeIndex<<std::endl;
    */
    // only consider when pro2 is NOT implicit-lipid.  This rules the pair out on
    // its own and is a single flag test, so it is checked before scanning pro1's
    // partner list rather than after.
    bool canInteract { !moleculeList[pro2Index].isImplicitLipid };
    if (canInteract) {
        canInteract = false;
        for (auto partner : molTemplateList[pro1MolType].rxnPartners) {
            if (partner == moleculeList[pro2Index].molTypeIndex) {
                canInteract = true;
                break;
            }
        }
    }

    bool canExclude { false };
    if ((moleculeList[pro1Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro1Index].molTypeIndex].excludeVolumeBound == true)
        || (moleculeList[pro2Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro2Index].molTypeIndex].excludeVolumeBound == true)) {
        canExclude = true;
    }

    // If this pair of proteins are already bound together, don't test for binding OR overlap
    if (canInteract) {
        /*If this pair of proteins are already bound together, don't test for binding OR overlap, set canInteract=0*/
        if (moleculeList[pro1Index].bndlist.size() > 0 && moleculeList[pro2Index].bndlist.size() > 0) {
            bool boundPro1 { false };
            for (auto& partner : moleculeList[pro1Index].bndpartner) {
                if (partner == pro2Index) {
                    boundPro1 = true;
                    break;
                }
            }
            // pro2's partner list only matters when pro1 already named pro2, so
            // it is not scanned otherwise.
            if (boundPro1) {
                bool boundPro2 { false };
                for (auto& partner : moleculeList[pro2Index].bndpartner) {
                    if (partner == pro1Index) {
                        boundPro2 = true;
                        break;
                    }
                }
                // TODO: add other case
                if (boundPro2)
                    canInteract = false;
            }
        }
    }

    // Molecules further apart than rMaxLimit cannot bind, so the interface
    // loops below have nothing to find.  rMaxLimit is the largest, over every
    // bimolecular reaction, of Rmax plus both interface-to-COM arms, and the
    // triangle inequality gives
    //
    //     |COM1 - COM2| <= |iface1 - iface2| + r1 + r2 < Rmax + r1 + r2
    //
    // for any pair close enough to react.  That is the same bound the cell list
    // already rests on when it sizes its sub-volumes by rMaxLimit, so this adds
    // no assumption; it only stops the pair one step earlier.  Worth doing
    // because the cell list is deliberately over-inclusive: of the pairs
    // reaching this point, 78% (clathrin) to 99% (compartment) are outside
    // rMaxLimit yet still walk the whole freelist x rxnPartners x freelist
    // search and find_which_reaction() for every combination.
    //
    // Two pairs are exempt.
    //
    // The volume-exclusion path below is not gated at all: it carries its own
    // radii, and its 2D branch reaches to RMax * 10, far outside rMaxLimit.
    //
    // Pairs with both complexes on the surface are not gated either, because
    // rMaxLimit is not a bound for them.  set_rMaxLimit() estimates a 2D
    // reaction's reach as 3 sqrt(6 Dtot dt), while
    // determine_2D_bimolecular_reaction_probability() uses 3.5 sqrt(4 Dtot dt)
    // over a Dtot that add_2D_rotational_diffusion() and discretize_2D_Dtot()
    // have both revised.  On rev_2D, whose rMaxLimit is set by a 2D reaction,
    // 1550 of 20220 reacting pairs sit beyond rMaxLimit, by up to 1.109x.  On
    // rev_3D, rev_3Dto2D, clathrin and mem_localization, whose rMaxLimit comes
    // from a 3D reaction, not one of 16.7 million does.
    if (canInteract) {
        const Vec3D comSep { moleculeList[pro1Index].comCoord - moleculeList[pro2Index].comCoord };
        if (comSep.x * comSep.x + comSep.y * comSep.y + comSep.z * comSep.z
            > params.rMaxLimit * params.rMaxLimit) {
            const bool bothOnSurface { complexList[moleculeList[pro1Index].myComIndex].OnSurface
                && complexList[moleculeList[pro2Index].myComIndex].OnSurface };
            if (!bothOnSurface)
                canInteract = false;
        }
    }

    if (canInteract) {
        /* CALCULATE ASSOCIATION PROBABILITIES */
        /* if(pro1Index== track1)
 	    std::cout <<" calculate Association prob to ! "<<pro2Index<<std::endl;
 	    if(pro1Index== track2)
 	    std::cout <<" calculate Association prob to ! "<<pro2Index<<std::endl;
      */
        for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro1Index].freelist.size(); ++relIface1Itr) {
            /*test all of i1's binding partners to see whether they are on protein pro2 */
            int relIface1 { moleculeList[pro1Index].freelist[relIface1Itr] };
            int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
            int stateIndex1 { moleculeList[pro1Index].interfaceList[relIface1].stateIndex };
            const Interface::State& state {
                molTemplateList[pro1MolType].interfaceList[relIface1].stateList[stateIndex1]
            };
            //if(pro1Index == track1 || pro1Index ==track2)
            //std::cout <<"reliface, abs, stateindex, npartners: "<<relIface1<<' '<<absIface1<<' '<<stateIndex1<<' '<<state.rxnPartners.size()<<std::endl;
            for (auto statePartner : state.rxnPartners) {
                //if(pro1Index == track1 || pro1Index ==track2)
                //  std::cout <<"statePertner: "<<statePartner<<std::endl;
                for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].freelist.size(); ++relIface2Idx) {
                    int relIface2 { moleculeList[pro2Index].freelist[relIface2Idx] };
                    int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };

                    if (absIface2 == statePartner) { // both binding interfaces are available!
                        /*if(pro1Index== track1)
                        std::cout <<" Both Ifaces available! "<<pro2Index<<" abs1: "<<absIface1<<" abs2: "<<absIface2<<std::endl;
                        if(pro1Index== track2)
                        std::cout <<" Both Ifaces available!! "<<pro2Index<<" abs1: "<<absIface1<<" abs2: "<<absIface2<<std::endl;
                        */
                        /*Here now we evaluate the probability of binding*/
                        /*Different interfaces can have different diffusion constants if they
                        also rotate. <theta^2>=6Drparams.timeStep.
                        In that case, <d>=sin(sqrt(6Drparams.timeStep)/2)*2R
                        so add in <d>^2=4R^2sin^2(sqrt(6Drparams.timeStep)/2)
                        for a single clathrin, R=arm length,
                        otherwise R will be from pivot point of rotation, COM,
                        calculate distance from interface to the complex COM
                        */
                        int rxnIndex { -1 };
                        int rateIndex { -1 };
                        bool isStateChangeBackRxn { false };

                        find_which_reaction(relIface1, relIface2, rxnIndex, rateIndex, isStateChangeBackRxn, state,
                            moleculeList[pro1Index], moleculeList[pro2Index], forwardRxns, backRxns, molTemplateList);
                        if (rxnIndex != -1 && rateIndex != -1) {
                            // Read each molecule's complex index once.  Nothing
                            // below changes it, and it was previously reloaded
                            // roughly twenty times per candidate pair, which the
                            // compiler cannot elide across the calls that take
                            // moleculeList by non-const reference.
                            const int com1Index { moleculeList[pro1Index].myComIndex };
                            const int com2Index { moleculeList[pro2Index].myComIndex };

                            if (com1Index == com2Index) {
                                evaluate_binding_within_complex(pro1Index, pro2Index, relIface1, relIface2, rxnIndex,
                                    rateIndex, isStateChangeBackRxn, params, moleculeList,
                                    complexList, molTemplateList,
                                    forwardRxns[rxnIndex], backRxns, membraneObject, counterArrays);
                            } else {
                                const Complex& com1 = complexList[com1Index];
                                const Complex& com2 = complexList[com2Index];
                                Vec3D ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                    - com1.comCoord };
                                Vec3D ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                    - com2.comCoord };
                                double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                    + ifaceVec.z * ifaceVec.z };
                                double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                    + ifaceVec2.z * ifaceVec2.z };

                                // binding with explicit-lipid model.
                                //write_rng_state();
                                if (com1.onFiber && com2.onFiber) {
                                    // both complexes are on the fiber
                                    // Assume diffusion only occurs in x direction (1 D)
                                    double Dtot = com1.D.x + com2.D.x;

                                    BiMolData biMolData { pro1Index, pro2Index, com1Index, com2Index, relIface1, relIface2,
                                        absIface1, absIface2, Dtot, magMol1, magMol2 };
                                    // std::cout << "Evaluate 1-D binding " << pro1Index << ", " << pro2Index << " Dtot " << Dtot << std::endl;
                                    determine_1D_bimolecular_reaction_probability(
                                        simItr, rxnIndex, rateIndex,
                                        isStateChangeBackRxn, biMolData, params,
                                        moleculeList, complexList, forwardRxns,
                                        backRxns);
                                }
                                // else if (std::abs(com1.D.z) < 1E-16 && std::abs(com2.D.z) < 1E-16) {
                                else if (com1.OnSurface && com2.OnSurface) {
                                    // both Complexes are on the membrane, evaluate as 2D reaction
                                    double Dtot = 1.0 / 2.0 * (com1.D.x + com2.D.x)
                                        + 1.0 / 2.0 * (com1.D.y + com2.D.y);

                                    BiMolData biMolData { pro1Index, pro2Index, com1Index, com2Index, relIface1, relIface2,
                                        absIface1, absIface2, Dtot, magMol1, magMol2 };
                                    // std::cout << "Evaluate 2-D binding " << pro1Index << ", " << pro2Index << " Dtot " << Dtot << std::endl;
                                    determine_2D_bimolecular_reaction_probability(
                                        simItr, rxnIndex, rateIndex,
                                        isStateChangeBackRxn, DDTableIndex,
                                        tableIDs, biMolData, params,
                                        moleculeList, complexList, forwardRxns,
                                        backRxns, membraneObject, normMatrices,
                                        survMatrices, pirMatrices);
                                }
                                else {
                                    //3D reaction
                                    double Dtot = 1.0 / 3.0 * (com1.D.x + com2.D.x)
                                        + 1.0 / 3.0 * (com1.D.y + com2.D.y)
                                        + 1.0 / 3.0 * (com1.D.z + com2.D.z);

                                    BiMolData biMolData { pro1Index, pro2Index, com1Index, com2Index, relIface1, relIface2,
                                        absIface1, absIface2, Dtot, magMol1, magMol2 };
                                    // std::cout << "Evaluate 3-D binding " << pro1Index << ", " << pro2Index << " Dtot " << Dtot << std::endl;
                                    determine_3D_bimolecular_reaction_probability(
                                        simItr, rxnIndex, rateIndex,
                                        isStateChangeBackRxn, biMolData, params,
                                        moleculeList, complexList, forwardRxns,
                                        backRxns);
                                } //end else 3D
                                //read_rng_state();
                            }
                        }
                    }
                }
            }
        }
    } // These protein partners interact

    if (canExclude == true) {
        // pro1 is bound
        if (moleculeList[pro1Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro1Index].molTypeIndex].excludeVolumeBound == true) {
            // loop the bindlist to find the exclude partner of each interface
            for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro1Index].bndlist.size(); ++relIface1Itr) {
                int relIface1 { moleculeList[pro1Index].bndlist[relIface1Itr] };
                int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
                if (moleculeList[pro1Index].interfaceList[relIface1].isBound && molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList.empty() == false) { // make sure it's actually bound and need excludeVolumeBound
                    // loop the interfaceList of pro2
                    for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].interfaceList.size(); ++relIface2Idx) {
                        int relIface2 { relIface2Idx };
                        int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };
                        int molTypeIndex2 { moleculeList[pro2Index].molTypeIndex };
                        for (int indexItr { 0 }; indexItr < molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList.size(); ++indexItr) {
                            if (molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList[indexItr] == molTypeIndex2) {
                                double bindRadius { molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeRadiusList[indexItr] };
                                int rxnIndex { molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundReactList[indexItr] };
                                if (molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundIfaceList[indexItr] == relIface2) {
                                    // calculate the distance between the two infaces
                                    Vec3D ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                        - complexList[moleculeList[pro1Index].myComIndex].comCoord };
                                    Vec3D ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                        - complexList[moleculeList[pro2Index].myComIndex].comCoord };
                                    double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                        + ifaceVec.z * ifaceVec.z };
                                    double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                        + ifaceVec2.z * ifaceVec2.z };
                                    if (complexList[moleculeList[pro1Index].myComIndex].onFiber == true && complexList[moleculeList[pro2Index].myComIndex].onFiber == true){
                                        // this is on 1D (a fiber)
                                        double Dtot = complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x;
                                        double RMax { 4.0 * sqrt(2.0 * Dtot * params.timeStep) + bindRadius };
                                        if (moleculeList[pro1Index].isPromoter && !moleculeList[pro2Index].isPromoter){
                                            RMax = 4.0 * sqrt(2.0 * Dtot * params.timeStep);
                                        } else if (!moleculeList[pro1Index].isPromoter && moleculeList[pro2Index].isPromoter){
                                            RMax = 4.0 * sqrt(2.0 * Dtot * params.timeStep);
                                        }
                                        double R1 = abs(moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x);
                                        if (R1 < RMax) {
                                            record_crossing_pair(pro1Index, pro2Index, relIface1, relIface2,
                                                std::array<int, 3> { rxnIndex, 0, false }, moleculeList, complexList);
                                        }
                                    }
                                    // else if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-16 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z) < 1E-16) {
                                    else if (complexList[moleculeList[pro1Index].myComIndex].OnSurface && 
                                        complexList[moleculeList[pro2Index].myComIndex].OnSurface) {
                                        // both Complexes are on the membrane, evaluate as 2D reaction
                                        double Dtot = 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };
                                        add_2D_rotational_diffusion(biMolData, complexList, params);
                                        discretize_2D_Dtot(biMolData);

                                        double RMax { 3.5 * sqrt(4.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        if (membraneObject.isSphere == true) {
                                            Vec3D iface11 = moleculeList[pro1Index].interfaceList[relIface1].coord;
                                            Vec3D iface22 = moleculeList[pro2Index].interfaceList[relIface2].coord;
                                            double r1 = iface11.length();
                                            double r2 = iface22.length();
                                            double r = (r1 + r2) / 2.0; //membraneObject.sphereR; //
                                            double theta = acos((iface11.x * iface22.x + iface11.y * iface22.y + iface11.z * iface22.z) / r1 / r2);
                                            R1 = r * theta;
                                        } else {
                                            double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                            double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                            double dz = 0;
                                            R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        }
                                        if (R1 < RMax * 10.0) {
                                            record_crossing_pair(pro1Index, pro2Index, relIface1, relIface2,
                                                std::array<int, 3> { rxnIndex, 0, false }, moleculeList, complexList);
                                        }
                                    } else {
                                        //3D reaction
                                        double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.z + complexList[moleculeList[pro2Index].myComIndex].D.z);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };
                                        add_3D_rotational_diffusion(biMolData, complexList, params,
                                            params.numerics.classification.explicitLipidFlatDiffusion);

                                        double RMax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                        double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                        // double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
                                        //                 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
                                        double dz { (complexList[moleculeList[pro1Index].myComIndex].OnSurface &&
                                                    complexList[moleculeList[pro2Index].myComIndex].OnSurface)
                                                ? 0
                                                : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
                                        R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        if (R1 < RMax) {
                                            record_crossing_pair(pro1Index, pro2Index, relIface1, relIface2,
                                                std::array<int, 3> { rxnIndex, 0, false }, moleculeList, complexList);
                                        }
                                    } //end else 3D
                                }
                            }
                        }
                    }
                }
            }
        }

        // pro2 is bound
        if (moleculeList[pro2Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro2Index].molTypeIndex].excludeVolumeBound == true) {
            // loop the bindlist to find the exclude partner of each interface
            for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro2Index].bndlist.size(); ++relIface1Itr) {
                int relIface2 { moleculeList[pro2Index].bndlist[relIface1Itr] };
                int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };
                int molTypeIndex2 { moleculeList[pro2Index].molTypeIndex };
                if (moleculeList[pro2Index].interfaceList[relIface2].isBound && molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundList.empty() == false) { // make sure it's actually bound and need excludeVolumeBound
                    // loop the interfaceList of pro2
                    for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro1Index].interfaceList.size(); ++relIface2Idx) {
                        int relIface1 { relIface2Idx };
                        int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
                        int molTypeIndex1 { moleculeList[pro1Index].molTypeIndex };
                        for (int indexItr { 0 }; indexItr < molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundList.size(); ++indexItr) {
                            if (molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundList[indexItr] == molTypeIndex1) {
                                double bindRadius { molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeRadiusList[indexItr] };
                                int rxnIndex { molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundReactList[indexItr] };
                                if (molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundIfaceList[indexItr] == relIface1) {
                                    // calculate the distance between the two infaces
                                    Vec3D ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                        - complexList[moleculeList[pro1Index].myComIndex].comCoord };
                                    Vec3D ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                        - complexList[moleculeList[pro2Index].myComIndex].comCoord };
                                    double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                        + ifaceVec.z * ifaceVec.z };
                                    double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                        + ifaceVec2.z * ifaceVec2.z };
                                    // if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-16 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z) < 1E-16) {
                                    if (complexList[moleculeList[pro1Index].myComIndex].OnSurface && complexList[moleculeList[pro2Index].myComIndex].OnSurface) {
                                        // both Complexes are on the membrane, evaluate as 2D reaction
                                        double Dtot = 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };
                                        add_2D_rotational_diffusion(biMolData, complexList, params);
                                        discretize_2D_Dtot(biMolData);

                                        double RMax { 3.5 * sqrt(4.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        if (membraneObject.isSphere == true) {
                                            Vec3D iface11 = moleculeList[pro1Index].interfaceList[relIface1].coord;
                                            Vec3D iface22 = moleculeList[pro2Index].interfaceList[relIface2].coord;
                                            double r1 = iface11.length();
                                            double r2 = iface22.length();
                                            double r = (r1 + r2) / 2.0; //membraneObject.sphereR; //
                                            double theta = acos((iface11.x * iface22.x + iface11.y * iface22.y + iface11.z * iface22.z) / r1 / r2);
                                            R1 = r * theta;
                                        } else {
                                            double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                            double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                            // double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
                                            //                 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
                                            double dz { (complexList[moleculeList[pro1Index].myComIndex].OnSurface
                                                            && complexList[moleculeList[pro2Index].myComIndex].OnSurface)
                                                    ? 0
                                                    : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
                                            R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        }
                                        if (R1 < RMax * 10.0) {
                                            record_crossing_pair(pro1Index, pro2Index, relIface1, relIface2,
                                                std::array<int, 3> { rxnIndex, 0, false }, moleculeList, complexList);
                                        }
                                    } else {
                                        //3D reaction
                                        double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.z + complexList[moleculeList[pro2Index].myComIndex].D.z);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };
                                        add_3D_rotational_diffusion(biMolData, complexList, params,
                                            params.numerics.classification.explicitLipidFlatDiffusion);

                                        double RMax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                        double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                        // double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
                                        //                 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
                                        double dz { (complexList[moleculeList[pro1Index].myComIndex].OnSurface
                                                        && complexList[moleculeList[pro2Index].myComIndex].OnSurface)
                                                ? 0
                                                : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
                                        R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        if (R1 < RMax) {
                                            record_crossing_pair(pro1Index, pro2Index, relIface1, relIface2,
                                                std::array<int, 3> { rxnIndex, 0, false }, moleculeList, complexList);
                                        }
                                    } //end else 3D
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
