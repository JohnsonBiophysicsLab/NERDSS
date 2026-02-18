#include <algorithm>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void sweep_separation_complex_rot_fiber_box(
    int simItr, int selfIndex, Parameters &params,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    const std::vector<ForwardRxn> &forwardRxns,
    const std::vector<MolTemplate> &molTemplateList,
    const Membrane &membraneObject) {
  /*
    The molecule "selfIndex" is on Fiber
  */

  // std::cout << "Iter num: " << simItr << " Molecule Index: " << selfIndex << std::endl;
  // TRACE();
  int comIndex1{moleculeList[selfIndex].myComIndex};
  int startProIndex = selfIndex;
  int com1Size = complexList[comIndex1].memberList.size();

  int maxRows{1};
  for (auto memMol : complexList[comIndex1].memberList) {
    if (moleculeList[memMol].crossbase.size() > maxRows)
      maxRows = moleculeList[memMol].crossbase.size();
  }

  std::vector<int> ifaceList(maxRows * com1Size, -1);
  std::vector<int> overlapList(maxRows * com1Size, -1);
  std::vector<int> memCheckList(maxRows * com1Size, -1);
  // int ifaceList[maxRows * com1Size];
  // int overlapList[maxRows * com1Size];
  // int memCheckList[maxRows * com1Size];

  /*The sampled displacement for p1 is stored in traj. the position from the
   previous step is still stored in bases[p1].xcom, etc, and will be updated
   at the end of this routine*/

  /*figure out i2*/
  for (int c{0}; c < com1Size; ++c) {
    int proIndex = complexList[comIndex1].memberList[c];
    for (int i{0}; i < moleculeList[proIndex].crossbase.size(); ++i) {
      int p2{moleculeList[proIndex].crossbase[i]};
      int k2{moleculeList[p2].myComIndex};
      // if (complexList[k2].D.z < 1E-15) {
      if (complexList[k2].onFiber) {
        memCheckList[maxRows * c + i] = 1;
      } else
        memCheckList[maxRows * c + i] = 0;
      int i1{moleculeList[proIndex].mycrossint[i]};
      std::array<int, 3> rxnItr = moleculeList[proIndex].crossrxn[i];

      // get the partner interface
      ifaceList[maxRows * c + i] =
          (forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex == i1)
              ? forwardRxns[rxnItr[0]].reactantListNew[1].relIfaceIndex
              : forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex;
    }
  }

  // determine RS3Dinput
  double RS3Dinput{0.0};
  Complex targCom{complexList[comIndex1]};
  if (membraneObject.implicitLipid) {
    for (auto &molIndex : targCom.memberList) {
      for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
        if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] -
                     moleculeList[molIndex].molTypeIndex) < 1E-2) {
          RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
          break;
        }
      }
    }
  }

  // int tsave = 0; // This is not actually used

  int itr{0};
  int maxItr{50};
  int saveit{0};
  while (itr < maxItr) {
    int numOverlap{0};
    bool hasOverlap{false};
    int comIndex2{};
    double dr2{};
    for (unsigned memMolItr{0}; memMolItr < complexList[comIndex1].memberList.size(); ++memMolItr) 
    {
      int pro1Index = complexList[comIndex1].memberList[memMolItr];
      // if (simItr==3226 && selfIndex==52){
      //   for (int crossMemItr{0}; crossMemItr < moleculeList[pro1Index].crossbase.size();
      //      ++crossMemItr) {
      //       std::cout << "Mol 2: " << moleculeList[pro1Index].crossbase[crossMemItr] << std::endl;
      //      }
      // }
      for (int crossMemItr{0}; crossMemItr < moleculeList[pro1Index].crossbase.size();
           ++crossMemItr) {
        int pro2Index{moleculeList[pro1Index].crossbase[crossMemItr]};
        if (moleculeList[pro2Index].isImplicitLipid) continue;

        // For co-localized proteins, they do not affect other 1D objects
        if (moleculeList[pro1Index].isPromoter && !moleculeList[pro2Index].isPromoter) {
          continue; 
        } else if (!moleculeList[pro1Index].isPromoter && moleculeList[pro2Index].isPromoter) {
          continue; 
        }

        comIndex2 = moleculeList[pro2Index].myComIndex;
        /*Do not sweep for overlap if proteins are in the same complex, they
         * cannot move relative to one another!
         */
        if (comIndex1 != comIndex2) {
          int relIface1{moleculeList[pro1Index].mycrossint[crossMemItr]};
          int rxnItr{moleculeList[pro1Index].crossrxn[crossMemItr][0]};
          int relIface2{ifaceList[maxRows * memMolItr + crossMemItr]};

          Vector iface1Vec{
              moleculeList[pro1Index].interfaceList[relIface1].coord -
              complexList[comIndex1].comCoord};
          std::array<double, 9> M = create_euler_rotation_matrix(complexList[comIndex1].trajRot);
          iface1Vec = matrix_rotate(iface1Vec, M);

          double dx1{complexList[comIndex1].comCoord.x + iface1Vec.x +
                     complexList[comIndex1].trajTrans.x};
          double dy1{complexList[comIndex1].comCoord.y + iface1Vec.y +
                     complexList[comIndex1].trajTrans.y};
          double dz1{complexList[comIndex1].comCoord.z + iface1Vec.z +
                     complexList[comIndex1].trajTrans.z};

          /*Now complex 2*/
          // if (reflectList[comIndex2] == 0) {
          //     reflect_traj_complex_rad_rot(params, moleculeList,
          //     complexList[comIndex2], membraneObject, RS3Dinput);
          //     reflectList[comIndex2] = 1;
          // }

          Vector iface2Vec{
              moleculeList[pro2Index].interfaceList[relIface2].coord -
              complexList[comIndex2].comCoord};
          std::array<double, 9> M2 =
              create_euler_rotation_matrix(complexList[comIndex2].trajRot);
          iface2Vec = matrix_rotate(iface2Vec, M2);
          double dx2{complexList[comIndex2].comCoord.x + iface2Vec.x +
                     complexList[comIndex2].trajTrans.x};
          double dy2{complexList[comIndex2].comCoord.y + iface2Vec.y +
                     complexList[comIndex2].trajTrans.y};
          double dz2{complexList[comIndex2].comCoord.z + iface2Vec.z +
                     complexList[comIndex2].trajTrans.z};

          /*separation*/
          double dfx{dx1 - dx2};
          double dfy{dy1 - dy2};
          double dfz{dz1 - dz2};

          double sepLimit_sq = forwardRxns[rxnItr].bindRadius * forwardRxns[rxnItr].bindRadius;

          dr2 = (dfx * dfx);
          if (memCheckList[maxRows * memMolItr + crossMemItr] != 1) {
            dr2 += (dfy * dfy) + (dfz * dfz);
            if (dr2 < sepLimit_sq) {
              /*reselect positions for protein p2*/
              overlapList[numOverlap] = pro2Index;
              numOverlap++;
              hasOverlap = true;
            }
          } else {
            double dfx_origin = complexList[comIndex1].comCoord.x + iface1Vec.x -
                                complexList[comIndex2].comCoord.x - iface2Vec.x;
            if (!moleculeList[pro1Index].isPromoter && !moleculeList[pro2Index].isPromoter) {
              // Both molecules are not DNA sites
              // no hop and no overlap
              if (dfx_origin * dfx < 0 || dr2 < sepLimit_sq) {
                /*reselect positions for protein p2*/
                overlapList[numOverlap] = pro2Index;
                numOverlap++;
                hasOverlap = true;
              }
            } else {
              // At least one molecule is a DNA site
              // no overlap
              if (dr2 < sepLimit_sq) {
                /*reselect positions for protein p2*/
                overlapList[numOverlap] = pro2Index;
                numOverlap++;
                hasOverlap = true;
              }
            }
          }
          //   std::cout << "dr^2:" << dr2 << std::endl;
          //   std::cout << "sigma:" << forwardRxns[rxnItr].bindRadius <<
          //   std::endl; std::cout << "overlapSep:" << params.overlapSepLimit
          //   << std::endl;
        }  // ignore proteins within the same complex.
      }
    }
    /*Now resample positions of p1 and overlapList, if t>0, otherwise no
     overlap, so break from loop*/
    if (hasOverlap) {
      // std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<
      // std::endl; std::cout << complexList[comIndex1].comCoord.x +
      // complexList[comIndex1].trajTrans.x << ","
      //           << complexList[comIndex2].comCoord.x +
      //           complexList[comIndex2].trajTrans.x << std::endl;
      ++itr;
      complexList[comIndex1].trajTrans.x =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].D.x) * GaussV();
      complexList[comIndex1].trajTrans.y =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].D.y) * GaussV();
      complexList[comIndex1].trajTrans.z =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].D.z) * GaussV();
      complexList[comIndex1].trajRot.x =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.x) * GaussV();
      complexList[comIndex1].trajRot.y =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.y) * GaussV();
      complexList[comIndex1].trajRot.z =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.z) * GaussV();

      // reflectList[comIndex1] = 0;
      reflect_traj_complex_rad_rot(params, moleculeList, complexList[comIndex1],
                                   membraneObject, RS3Dinput,
                                   false);  // set isInsideCompartment = false.
      // reflectList[comIndex1] = 1;

      std::vector<int> resampleList(complexList.size(), 0);
      // int resampleList[complexList.size()];  // if this is 0, we need resample
      // for (auto &i : resampleList) {         // initialize
      //   i = 0;
      // }

      for (int checkMolItr{0}; checkMolItr < numOverlap; checkMolItr++) {
        int p2{overlapList[checkMolItr]};
        comIndex2 = moleculeList[p2].myComIndex;

        // avoid repeated resample
        if (resampleList[comIndex2] == 0) {
          // if (p2 > startProIndex && (moleculeList[p2].trajStatus ==
          // TrajStatus::none || moleculeList[p2].trajStatus ==
          // TrajStatus::canBeResampled)) {
          if ((moleculeList[p2].trajStatus == TrajStatus::none ||
               moleculeList[p2].trajStatus == TrajStatus::canBeResampled)) {
                /*
                We loop over proteins sequentially, so earlier proteins have already
                moved and avoided their neighbors and should not be moved again. These
                new positions selected for proteins not yet moved will be stored and
                then used when they test for overlap themselves.
                */
                // if (simItr==3226 && selfIndex==52){
                //   std::cout << "resample the overlaped moelcule: " << p2  << std::endl;
                // }

                /*If p2 just dissociated, also don't try to move again*/
                complexList[comIndex2].trajTrans.x =
                    sqrt(2.0 * params.timeStep * complexList[comIndex2].D.x) * GaussV();
                complexList[comIndex2].trajTrans.y =
                    sqrt(2.0 * params.timeStep * complexList[comIndex2].D.y) * GaussV();
                complexList[comIndex2].trajTrans.z =
                    sqrt(2.0 * params.timeStep * complexList[comIndex2].D.z) * GaussV();
                complexList[comIndex2].trajRot.x =
                    sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.x) * GaussV();
                complexList[comIndex2].trajRot.y =
                    sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.y) * GaussV();
                complexList[comIndex2].trajRot.z =
                    sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.z) * GaussV();

                // reflectList[comIndex2] = 0;
                reflect_traj_complex_rad_rot(
                    params, moleculeList, complexList[comIndex2], membraneObject,
                    RS3Dinput, false);  // set isInsideCompartment = false.
                // reflectList[comIndex2] = 1;
                resampleList[comIndex2] = 1;
          }
        }
      }
      // tsave = numOverlap;
    } else {
      saveit = itr + 1;
      itr = maxItr;  // break from loop
    }
  }  // end maximum iterations
  
  if (saveit == 0) {
    // std::cout << "Overlap Check Reached Iteration Limitation" << "(" << maxItr << ") at itr " << 
    // simItr << " for molecule " << selfIndex << std::endl;
    // throw std::runtime_error("Debug error: Overlap not fixed!");

    // Choice 1 (default): accept any movement. Moving is better than stuck. 
    // Although excluded volumes can prohibit most of unphysical overlaps,
    // in crowded systems, some overlaps may still happen.
    // So we accept any propagation to avoid deadlock.
    // DEBUG: Warning is given to inform users that overlaps happen.
    // std::cout << "Warning: Overlap Check Reached Iteration Limitation" << "(" << maxItr
    //           << ") at itr " << simItr << " for molecule " << selfIndex << std::endl;

    // Choice 2: cancel propagation if no available place found
    // complexList[comIndex1].trajTrans.zero_crds();
    // complexList[comIndex1].trajRot.zero_crds();

    /* Explain why there was no propagation cancel
    previously, in case to form a fake loop on a straight 1D fiber,
    (like GAGA factor which forms a hexamer while all subunits bind with the same DNA)
    proteins are allowed to form a "3D" structure while they are still considered as 
    1D objects. For 1D objects, the distance (both for overlap and reactions) is dx = |x2 - x1|,
    and the 3D structure may bring two proteins too close such that their 1D distance is 
    smaller than the overlap limit, sigma. To avoid proteins being trapped after dissociation
    (proteins are not reset to a reasonable place after dissociation), we forced proteins to 
    be propagated even the propagation may result in a place with overlap issue. However, there 
    are two problems:
    1. The no cancel method messed up kinetics.
    2. For normal crowded proteins, the no cancel method creates huge overlap issue.
    So, now we returned to the "cancel propagation" method to make sure crowded proteins have 
    correct movements.
    3. PRICE: have to apply excluded volumes for all 1D interactions, unless for co-localized proteins.
    */ 

  }

  // accept any propagation
  complexList[comIndex1].propagate(moleculeList, membraneObject, molTemplateList);

  // Reset displacements to zero so distance is measured to your current
  // updated position that won't change again this turn
  complexList[comIndex1].trajTrans.zero_crds();
  complexList[comIndex1].trajRot.zero_crds();
}