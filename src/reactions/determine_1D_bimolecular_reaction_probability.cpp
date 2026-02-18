#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "tracing.hpp"
#include <iostream>
#include <sstream>

void determine_1D_bimolecular_reaction_probability(
    int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    BiMolData &biMolData, const Parameters &params,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    const std::vector<ForwardRxn> &forwardRxns,
    const std::vector<BackRxn> &backRxns) {

  // rotation
  // double Dr1{};
  // double cf{cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep))};
  // Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
  // double Dr2{};
  // cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
  // Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
  // biMolData.Dtot += (Dr1 + Dr2) / (4.0 * params.timeStep); // add in contributions from rotation

  
  // the maximum distance along x axis for reaction to occur
  // for site + site or protein + protein, RMax = 4 * sqrt(2 * Dtot * dt) + sigma
  double RMax{4.0 * sqrt(2.0 * biMolData.Dtot * params.timeStep) + forwardRxns[rxnIndex].bindRadius};
  // Use 3D diffusion to identify if reactans are on the same fiber
  // Assume D is isotropic in 3D so that Dtot3D = Dtot1D
  double RMax3D{3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + forwardRxns[rxnIndex].bindRadius};
  // for protein binding to DNA sites, consider Sigma_1D = 0 when calculating RMax and reaction rate
  // protein + DNA site 
  // TODO: use xor in future
  // proteins always bind on top of DNA sites. Therefore, we do not need to add the binding radius here.
  // Namely, we assume sigma_x = 0
  bool bindfrom3Dto1D {false};
  if (moleculeList[biMolData.pro1Index].isPromoter && !moleculeList[biMolData.pro2Index].isPromoter) {
    RMax = 4.0 * sqrt(2.0 * biMolData.Dtot * params.timeStep);
    bindfrom3Dto1D = true; 
  } else if (!moleculeList[biMolData.pro1Index].isPromoter && moleculeList[biMolData.pro2Index].isPromoter) {
    RMax = 4.0 * sqrt(2.0 * biMolData.Dtot * params.timeStep);
    bindfrom3Dto1D = true;
  }
  // double sep{};
  // double R1{};
  // bool withinRmax{get_distance(
  //     biMolData.pro1Index, biMolData.pro2Index, biMolData.relIface1,
  //     biMolData.relIface2, rxnIndex, rateIndex, isStateChangeBackRxn, sep, R1,
  //     RMax, complexList, forwardRxns[rxnIndex], moleculeList, false)};

  // We do not use get_distance here since we need two Rmax values but get_distance 
  // only takes one Rmax as input

  // calculate the 3D distance to determine if the reaction is on the same fiber
  // calculate the x-axis distance and sep
  double R3D{};
  double R1{};
  double dx = moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord.x -
              moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord.x;
  double dy = moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord.y -
              moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord.y;
  double dz = moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord.z -
              moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord.z;
  R3D = sqrt((dx * dx) + (dy * dy) + (dz * dz));
  R1 = abs(dx);

  bool withinRmax{R1 <= RMax && R3D <= RMax3D};

  // if (R1 <= RMax && R3D > RMax3D) {
  //   std::cout << "Molecules not on the same fiber, but dX < Rmax1D." << std::endl;
  //   std::cout << " R1: " << R1 << " RMax: " << RMax << " R3D: " << R3D << " RMax3D: " << RMax3D << std::endl;
  //   std::cout << " Mol1Index: " << biMolData.pro1Index << " Mol2Index: " << biMolData.pro2Index << std::endl; 
  // }
  
  // Make sure R1 is at least bindRadius
  if (bindfrom3Dto1D) {
    // when binding from 3D to 1D, set any R1 is acceptable (sigma_x = 0)
  } else {
    if  (R1 < forwardRxns[rxnIndex].bindRadius){
      R1 = forwardRxns[rxnIndex].bindRadius;
    }
  }
    

  if (withinRmax) {
    // std::cout << "RMax3D " << RMax3D << " Mol1Index: " << biMolData.pro1Index << " Mol2Index: " << biMolData.pro2Index << std::endl;
    /*in this case we evaluate the probability of this reaction*/
    moleculeList[biMolData.pro1Index].crossbase.push_back(biMolData.pro2Index);
    moleculeList[biMolData.pro2Index].crossbase.push_back(biMolData.pro1Index);
    moleculeList[biMolData.pro1Index].mycrossint.push_back(biMolData.relIface1);
    moleculeList[biMolData.pro2Index].mycrossint.push_back(biMolData.relIface2);
    moleculeList[biMolData.pro1Index].crossrxn.push_back(
        std::array<int, 3> { rxnIndex, rateIndex, isStateChangeBackRxn });
    moleculeList[biMolData.pro2Index].crossrxn.push_back(
        std::array<int, 3> { rxnIndex, rateIndex, isStateChangeBackRxn });
    ++complexList[moleculeList[biMolData.pro1Index].myComIndex].ncross;
    ++complexList[moleculeList[biMolData.pro2Index].myComIndex].ncross;
    // in case they dissociated
    // moleculeList[biMolData.pro1Index].display_all();
    // moleculeList[biMolData.pro2Index].display_all();
    moleculeList[biMolData.pro1Index].probvec.push_back(0);
    moleculeList[biMolData.pro2Index].probvec.push_back(0);
  } 

  // if (moleculeList[biMolData.pro1Index].isDissociated == true ||
  //     moleculeList[biMolData.pro2Index].isDissociated == true) {
  //   std::cout << "One of the molecules is dissociated. Skip reaction probability calculation." << std::endl;
  //   std::cout << " Mol1Index: " << biMolData.pro1Index << " Mol2Index: " << biMolData.pro2Index << std::endl; 
  // }
  
  if (moleculeList[biMolData.pro1Index].isDissociated != true &&
      moleculeList[biMolData.pro2Index].isDissociated != true) {
    /*This movestat check is if you allow just dissociated proteins to avoid
     * overlap*/
    if (withinRmax && forwardRxns[rxnIndex].rateList[rateIndex].rate > 0) {
      // int probMatrixIndex{0};
      /*Evaluate probability of reaction, with reweighting*/
      // double ratio = forwardRxns[rxnIndex].bindRadius / R1;

      // declare intrinsic binding rate of 1D->1D case.
      double kact{forwardRxns[rxnIndex].rateList[rateIndex].rate / forwardRxns[rxnIndex].area3Dto1D};
      if (forwardRxns[rxnIndex].isSymmetric == false)
        kact /= 2.0; // for A(a)+B(b)->A(a!).B(b!) case
        // This is different from 3D case since molecules can only approach from one side in 1D
        // This is irrelevant of bypassing is allowed or not.

      double currnorm{1.0};
      double p0_ratio{1.0};

      /*protein pro1 is molTypeIndex wprot and pro2 is wprot2*/
      int proA = biMolData.pro1Index;
      int ifaceA = biMolData.relIface1;
      int proB = biMolData.pro2Index;
      int ifaceB = biMolData.relIface2;

      if (biMolData.pro1Index > biMolData.pro2Index) {
        proA = biMolData.pro2Index;
        ifaceA = biMolData.relIface2;
        proB = biMolData.pro1Index;
        ifaceB = biMolData.relIface1;
      }

      double rxnProb{};
      for (int s{0}; s < moleculeList[proA].prevlist.size(); ++s) {
        if (moleculeList[proA].prevlist[s] == proB &&
            moleculeList[proA].prevmyface[s] == ifaceA &&
            moleculeList[proA].prevpface[s] == ifaceB) {
          if (moleculeList[proA].prevsep[s] >= RMax) {
            // BEcause previous reweighting was for 3D, now
            p0_ratio = 1.0;
            // restart reweighting in 2D.
            currnorm = 1.0;
          } else {
            if (bindfrom3Dto1D) {
              p0_ratio = pirr_pfree_ratio_psF_1D(
                R1, moleculeList[proA].prevsep[s], params.timeStep,
                biMolData.Dtot, 0.0, kact,
                moleculeList[proA].ps_prev[s]);
            } else {
              p0_ratio = pirr_pfree_ratio_psF_1D(
                R1, moleculeList[proA].prevsep[s], params.timeStep,
                biMolData.Dtot, forwardRxns[rxnIndex].bindRadius, kact,
                moleculeList[proA].ps_prev[s]);
            }
            currnorm = moleculeList[proA].prevnorm[s] * p0_ratio;
          }
          break;
        }
      }
      if (bindfrom3Dto1D) {
        rxnProb = passocF_1D(R1, params.timeStep, biMolData.Dtot, 0.0, kact);
      } else {
        rxnProb = passocF_1D(R1, params.timeStep, biMolData.Dtot, forwardRxns[rxnIndex].bindRadius, kact);
      }
      // debug info
      // if (bindfrom3Dto1D) {
      //   if (rxnProb * currnorm > 1e-3){
      //     std::cout << "bindfrom3Dto1D is true." << std::endl;
      //     std::cout << "In determine_1D_prob proB " << proB << " proA " << proA
      //           << " R1 " << R1 << " kact " << kact << " rxnProb " << rxnProb
      //           << " p0_ratio " << p0_ratio << " currnorm " << currnorm << std::endl;
      //   } 
      // } 
      // END

      moleculeList[biMolData.pro1Index].probvec.back() = rxnProb * currnorm;
      moleculeList[biMolData.pro2Index].probvec.back() = rxnProb * currnorm;
      if (rxnProb > 1.000001) {
        std::cerr << "Error: prob of reaction is: " << rxnProb
                  << " > 1. Avoid this using a smaller time step." << std::endl;
        // exit(1);
      }
      if (rxnProb > 0.5) {
        std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction "
                     "for a bimolecular binding with multiple binding sites, "
                     "please use a smaller time step."
                  << std::endl;
      }
      moleculeList[proA].currprevsep.push_back(R1);
      moleculeList[proA].currlist.push_back(proB);
      moleculeList[proA].currmyface.push_back(ifaceA);
      moleculeList[proA].currpface.push_back(ifaceB);
      moleculeList[proA].currprevnorm.push_back(currnorm);
      moleculeList[proA].currps_prev.push_back(1.0 - rxnProb * currnorm);
    } // Within reaction zone
  }
}