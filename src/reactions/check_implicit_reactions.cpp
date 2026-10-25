#include "classes/class_Membrane.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include <algorithm>

// pro2Index is a lipid's index.

#ifdef NERDSS_CHECK_LIPID_LOOKUP
#include <cstdlib>
#include <iostream>
namespace {
// How many lookups the both-ways build compared, reported at exit.
struct LookupTally {
  unsigned long long checked{0};
  ~LookupTally() {
    std::cerr << "NERDSS_CHECK_LIPID_LOOKUP: " << checked
              << " lookups agree\n";
  }
} lookupTally;

// Both-ways build: repeats the lookup the old way -- the lipid's interface
// index moved by the state offset, find_which_reaction(), moved back -- and
// stops at the first difference, including in the totMatches count.
void check_lipid_lookup(int relIface1, int relIface2, int lipidState,
                        unsigned long newCount, int rxnIndex, int rateIndex,
                        bool isStateChangeBackRxn,
                        const Interface::State &state, const Molecule &pro1,
                        const Molecule &lipid,
                        const std::vector<ForwardRxn> &forwardRxns,
                        const std::vector<BackRxn> &backRxns) {
  int oldRxn{-1};
  int oldRate{-1};
  bool oldBack{false};
  const unsigned long before{totMatches};
  Molecule &bumped{const_cast<Molecule &>(lipid)};
  bumped.interfaceList[relIface2].index += lipidState;
  find_which_reaction(relIface1, relIface2, oldRxn, oldRate, oldBack, state,
                      pro1, lipid, forwardRxns, backRxns);
  bumped.interfaceList[relIface2].index -= lipidState;
  const unsigned long oldCount{totMatches - before};
  totMatches = before;
  if (oldRxn != rxnIndex || oldRate != rateIndex ||
      oldBack != isStateChangeBackRxn || oldCount != newCount) {
    std::cerr << "NERDSS_CHECK_LIPID_LOOKUP: mismatch after "
              << lookupTally.checked
              << " lookups: old (" << oldRxn << ", " << oldRate << ", "
              << oldBack << ", " << oldCount << "), new (" << rxnIndex << ", "
              << rateIndex << ", " << isStateChangeBackRxn << ", " << newCount
              << ")\n";
    std::abort();
  }
  ++lookupTally.checked;
}
} // namespace
#endif

void check_implicit_reactions(
    int pro1Index, int pro2Index, int simItr, const Parameters &params,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    const std::vector<MolTemplate> &molTemplateList,
    const std::vector<ForwardRxn> &forwardRxns,
    const std::vector<BackRxn> &backRxns, copyCounters &counterArrays,
    Membrane &membraneObject, std::vector<double> &IL2DbindingVec,
    std::vector<double> &IL2DUnbindingVec, std::vector<double> &ILTableIDs) {
  // only consider when pro2 IS implicit-lipid
  if (moleculeList[pro2Index].isImplicitLipid == false)
    return;

  int pro1MolType = moleculeList[pro1Index].molTypeIndex;

  bool canInteract{true};
  int index = moleculeList[pro2Index]
                  .interfaceList[0]
                  .index;   // protein2 must be implicit-lipid, and has only one
                            // interface.
  int number_of_lipids = 0; // sum of all states of IL
  for (int i = 0; i < membraneObject.numberOfFreeLipidsEachState.size(); i++) {
    number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
  }
  if (number_of_lipids <= 0) {
    canInteract = false;
  }

  // allow multiple links to the surface.
  /*      if (canInteract) {
      // If this protein is already bound to surface, don't test for binding OR
 overlap, set canInteract= false/ if (moleculeList[pro1Index].bndlist.size() > 0
 ) { bool boundPro1 { false }; for (auto& partner :
 moleculeList[pro1Index].bndpartner) { if (partner == pro2Index) { boundPro1 =
 true; break;
              }
          }
          if (boundPro1)
              canInteract = false;
      }
 }
 */
  if (canInteract) {
    /* CALCULATE ASSOCIATION PROBABILITIES */
    for (int relIface1Itr{0};
         relIface1Itr < moleculeList[pro1Index].freelist.size();
         ++relIface1Itr) {
      // test all of i1's binding partners to see whether they are on protein
      // pro2
      int relIface1{moleculeList[pro1Index].freelist[relIface1Itr]};
      int absIface1{moleculeList[pro1Index].interfaceList[relIface1].index};
      int stateIndex1{
          moleculeList[pro1Index].interfaceList[relIface1].stateIndex};
      const Interface::State &state{molTemplateList[pro1MolType]
                                        .interfaceList[relIface1]
                                        .stateList[stateIndex1]};

      for (std::size_t partnerItr{0}; partnerItr < state.rxnPartners.size();
           ++partnerItr) {
        const unsigned statePartner{state.rxnPartners[partnerItr]};
        for (int relIface2Idx = 0;
             relIface2Idx < moleculeList[pro2Index].freelist.size();
             ++relIface2Idx) {
          int relIface2{moleculeList[pro2Index].freelist[relIface2Idx]};
          int absIface2{moleculeList[pro2Index].interfaceList[relIface2].index};
          // The lipid states used to be tried one by one, each passing when
          // absIface2 == statePartner - state in unsigned arithmetic.  Only
          // state == statePartner - absIface2 can, so that one is tested alone.
          const unsigned lipidState{statePartner -
                                    static_cast<unsigned>(absIface2)};
          if (membraneObject.nStates <= 0 ||
              lipidState >= static_cast<unsigned>(membraneObject.nStates))
            continue;
          const int tmpImplicitLipidStateIndex{static_cast<int>(lipidState)};
          if (membraneObject.numberOfFreeLipidsEachState
                  [tmpImplicitLipidStateIndex] >=
              0) { // both binding interfaces are available!
            int rxnIndex{-1};
            int rateIndex{-1};
            bool isStateChangeBackRxn{false};

            if (!moleculeList[pro1Index].isImplicitLipid) {
              // What find_which_reaction() returns with the lipid's
              // interface index moved to statePartner.
#ifdef NERDSS_CHECK_LIPID_LOOKUP
              const unsigned long matchesBefore{totMatches};
#endif
              find_which_reaction_for_partner(
                  relIface1, partnerItr, rxnIndex, rateIndex,
                  isStateChangeBackRxn, state, moleculeList[pro1Index],
                  moleculeList[pro2Index], forwardRxns, backRxns);
#ifdef NERDSS_CHECK_LIPID_LOOKUP
              check_lipid_lookup(relIface1, relIface2,
                                 tmpImplicitLipidStateIndex,
                                 totMatches - matchesBefore, rxnIndex,
                                 rateIndex, isStateChangeBackRxn, state,
                                 moleculeList[pro1Index],
                                 moleculeList[pro2Index], forwardRxns,
                                 backRxns);
#endif
            } else {
              moleculeList[pro2Index].interfaceList[relIface2].index +=
                  tmpImplicitLipidStateIndex;
              find_which_reaction(relIface1, relIface2, rxnIndex, rateIndex,
                                  isStateChangeBackRxn, state,
                                  moleculeList[pro1Index],
                                  moleculeList[pro2Index], forwardRxns,
                                  backRxns);
              moleculeList[pro2Index].interfaceList[relIface2].index -=
                  tmpImplicitLipidStateIndex;
            }

            if (rxnIndex != -1 && rateIndex != -1) {
              // int com1Index { moleculeList[pro1Index].myComIndex };
              // int com2Index { moleculeList[pro2Index].myComIndex };

              {
                Vec3D ifaceVec{
                    moleculeList[pro1Index].interfaceList[relIface1].coord -
                    complexList[moleculeList[pro1Index].myComIndex].comCoord};
                Vec3D ifaceVec2{
                    moleculeList[pro2Index].interfaceList[relIface2].coord -
                    complexList[moleculeList[pro2Index].myComIndex].comCoord};
                double magMol1{ifaceVec.x * ifaceVec.x +
                               ifaceVec.y * ifaceVec.y +
                               ifaceVec.z * ifaceVec.z};
                double magMol2{ifaceVec2.x * ifaceVec2.x +
                               ifaceVec2.y * ifaceVec2.y +
                               ifaceVec2.z * ifaceVec2.z};

                // find the state of IL
                int forwardRxnIndex{rxnIndex};
                const ForwardRxn &currRxn = forwardRxns[forwardRxnIndex];
                int relStateIndex{implicit_lipid_state_index(
                    implicit_lipid_state_list(moleculeList, molTemplateList,
                                              membraneObject),
                    implicit_lipid_iface(currRxn.reactantListNew,
                                         molTemplateList))};

                // binding with implicit-lipid model
                if (std::abs(
                        complexList[moleculeList[pro1Index].myComIndex].D.z) <
                    1E-10) {
                  // both Complexes are on the membrane, evaluate as 2D
                  // reaction
                  double Dtot =
                      1.0 / 2.0 *
                          (complexList[moleculeList[pro1Index].myComIndex]
                               .D.x +
                           complexList[moleculeList[pro2Index].myComIndex]
                               .D.x) +
                      1.0 / 2.0 *
                          (complexList[moleculeList[pro1Index].myComIndex]
                               .D.y +
                           complexList[moleculeList[pro2Index].myComIndex]
                               .D.y);

                  BiMolData biMolData{pro1Index,
                                      pro2Index,
                                      moleculeList[pro1Index].myComIndex,
                                      moleculeList[pro2Index].myComIndex,
                                      relIface1,
                                      relIface2,
                                      absIface1,
                                      absIface2,
                                      Dtot,
                                      magMol1,
                                      magMol2};
                  determine_2D_implicitlipid_reaction_probability(
                      simItr, rxnIndex, rateIndex, isStateChangeBackRxn,
                      ILTableIDs, biMolData, params, moleculeList,
                      complexList, forwardRxns, backRxns, IL2DbindingVec,
                      IL2DUnbindingVec, membraneObject, relStateIndex);
                } else {
                  // 3D->2D reaction
                  double Dtot =
                      1.0 / 3.0 *
                          (complexList[moleculeList[pro1Index].myComIndex]
                               .D.x +
                           complexList[moleculeList[pro2Index].myComIndex]
                               .D.x) +
                      1.0 / 3.0 *
                          (complexList[moleculeList[pro1Index].myComIndex]
                               .D.y +
                           complexList[moleculeList[pro2Index].myComIndex]
                               .D.y) +
                      1.0 / 3.0 *
                          (complexList[moleculeList[pro1Index].myComIndex]
                               .D.z +
                           complexList[moleculeList[pro2Index].myComIndex]
                               .D.z);

                  BiMolData biMolData{pro1Index,
                                      pro2Index,
                                      moleculeList[pro1Index].myComIndex,
                                      moleculeList[pro2Index].myComIndex,
                                      relIface1,
                                      relIface2,
                                      absIface1,
                                      absIface2,
                                      Dtot,
                                      magMol1,
                                      magMol2};

                  determine_3D_implicitlipid_reaction_probability(
                      simItr, rxnIndex, rateIndex, isStateChangeBackRxn,
                      biMolData, params, moleculeList, complexList,
                      forwardRxns, backRxns, membraneObject, relStateIndex);
                }
              }
            }
          }
        }
      }
    } // These protein partners interact
  }
}
