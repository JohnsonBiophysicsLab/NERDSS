#include "math/constants.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

/* Per-molecule unimolecular creation, destruction and dissociation, called by
 * the MPI step only (check_perform_zeroth_first_order_reactions()); the serial
 * loop replaced its call with check_for_unimolecular_reactions_population(),
 * which draws per reaction from the templates' monomerLists, lists the MPI
 * exchange does not maintain.  Unimolecular state changes are not here: the
 * MPI step runs check_for_unimolstatechange_reactions() for those, the same
 * scan the serial loop runs.
 *
 * Unlike the serial routines it mirrors, this does not maintain
 * counterArrays.copyNumSpecies, and that is deliberate: in the MPI build the
 * array is a scratch buffer that src/io_mpi/write_all_species.cpp zeroes and
 * recounts from moleculeList at every write, and nothing reads it between
 * writes (the free-protein count in check_dissociation_implicitlipid.cpp
 * counts moleculeList itself under mpi_).  A tally kept here would be
 * overwritten before anyone read it. */
void check_for_unimolecular_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<double>& IL2DbindingVec,
    std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs, std::ofstream& assocDissocFile) {
  // Note: cannot use a for-range or vector iterator loop here, since we are
  // changing the moleculeList vector while iterating over it
  for (unsigned molItr{0}; molItr < moleculeList.size(); ++molItr) {
    // only do checks if the Molecule exists
    if (moleculeList[molItr].isEmpty || moleculeList[molItr].isGhosted)
      continue;
#ifdef mpi_
    if (moleculeList[molItr].firstOrderHeld) continue;
#endif

    // now check for unimolecular creation and destruction
    // TODO: this is currently constrained to only consider a Molecule and only
    // have one rate
    bool hasDestroyed{false};
    if (moleculeList[molItr].isImplicitLipid == true) {
      for (const auto& oneRxn : createDestructRxns) {
        for (int tmpImplicitLipidStateIndex = 0;
             tmpImplicitLipidStateIndex < membraneObject.nStates;
             tmpImplicitLipidStateIndex++) {
          moleculeList[molItr].interfaceList[0].index +=
              tmpImplicitLipidStateIndex;
          moleculeList[molItr].interfaceList[0].relIndex +=
              tmpImplicitLipidStateIndex;
          moleculeList[molItr].interfaceList[0].stateIndex +=
              tmpImplicitLipidStateIndex;
          if (oneRxn.rxnType == ReactionType::destruction) {
            if (isReactant(moleculeList[molItr],
                           complexList[moleculeList[molItr].myComIndex], oneRxn,
                           moleculeList)) {
              // desruction of IL, , IL only can have one inerface and must be
              // the first mol type
              int tmpOneStateLipid =
                  membraneObject
                      .numberOfFreeLipidsEachState[tmpImplicitLipidStateIndex];

              for (int tmpItr = 0; tmpItr < tmpOneStateLipid; tmpItr++) {
                double prob{first_order_probability(oneRxn.rateList[0].rate,
                                                   params.timeStep)};

                // make sure that the time step is resonable according to the
                // prob of reaction
                if (prob > 1.000001) {
                  std::cerr << "Error: prob of reaction > 1. Avoid this using "
                               "a smaller time step."
                            << std::endl;
                  // exit(1);
                }
                if (prob > 0.5) {
                  // std::cout << "WARNING: prob of reaction > 0.5. If this is a
                  // reaction for a bimolecular binding with multiple binding
                  // sites, please use a smaller time step." << std::endl;
                }

                if (prob > rand_gsl()) {
                  // double rNum2 { rNum + rand_gsl() * Constants::iRandMax };
                  // // to get higher resolution if (prob > rNum2) { std::cout
                  // << "Destroying molecule of type IL"
                  //           << " with reaction "
                  //           << oneRxn.absRxnIndex << " at iteration " <<
                  //           simItr << '\n';

                  --membraneObject.numberOfFreeLipidsEachState
                        [tmpImplicitLipidStateIndex];

                  if (oneRxn.isObserved) {
                    auto observeItr = observablesList.find(oneRxn.observeLabel);
                    if (observeItr == observablesList.end()) {
                      // std::cerr << "WARNING: Observable " <<
                      // oneRxn.observeLabel << " not defined.\n";
                    } else {
                      --observeItr->second;
                    }
                  }
                  break;
                  // }
                }
              }
            }
          }
          if (oneRxn.rxnType == ReactionType::uniMolCreation) {
            if (isReactant(moleculeList[molItr],
                           complexList[moleculeList[molItr].myComIndex], oneRxn,
                           moleculeList)) {
              // creation of IL, IL only can have one inerface and must be the
              // first mol type
              int tmpOneStateLipid =
                  membraneObject
                      .numberOfFreeLipidsEachState[tmpImplicitLipidStateIndex];
              for (int tmpItr = 0; tmpItr < tmpOneStateLipid; tmpItr++) {
                // double rNum { 1.0 * rand_gsl() };
                double prob{first_order_probability(oneRxn.rateList.at(0).rate,
                                                   params.timeStep)};

                // make sure that the time step is resonable according to the
                // prob of reaction
                if (prob > 1.000001) {
                  std::cerr << "Error: prob of reaction > 1. Avoid this using "
                               "a smaller time step."
                            << std::endl;
                  // exit(1);
                }
                if (prob > 0.5) {
                  // std::cout << "WARNING: prob of reaction > 0.5. If this is a
                  // reaction for a bimolecular binding with multiple binding
                  // sites, please use a smaller time step." << std::endl;
                }

                if (prob > rand_gsl()) {
                  // std::cout << "Creating molecule of type IL"
                  //           << " from reaction " << oneRxn.absRxnIndex <<
                  //           '\n';

                  ++membraneObject.numberOfFreeLipidsEachState
                        [tmpImplicitLipidStateIndex];

                  if (oneRxn.isObserved) {
                    auto observeItr = observablesList.find(oneRxn.observeLabel);
                    if (observeItr == observablesList.end()) {
                      // std::cerr << "WARNING: Observable " <<
                      // oneRxn.observeLabel << " not defined.\n";
                    } else {
                      ++observeItr->second;
                    }
                  }
                  break;
                }
              }
            }
          }
          moleculeList[molItr].interfaceList[0].index -=
              tmpImplicitLipidStateIndex;
          moleculeList[molItr].interfaceList[0].relIndex -=
              tmpImplicitLipidStateIndex;
          moleculeList[molItr].interfaceList[0].stateIndex -=
              tmpImplicitLipidStateIndex;
        }
      }
    } else {
      for (const auto& oneRxn : createDestructRxns) {
        if (oneRxn.rxnType == ReactionType::destruction &&
            moleculeList[molItr].trajStatus == TrajStatus::none &&
            isReactant(moleculeList[molItr],
                       complexList[moleculeList[molItr].myComIndex], oneRxn,
                       moleculeList)) {
          double prob{
              first_order_probability(oneRxn.rateList[0].rate, params.timeStep)};

          // make sure that the time step is resonable according to the prob of
          // reaction
          if (prob > 1.000001) {
            std::cerr << "Error: prob of reaction > 1. Avoid this using a "
                         "smaller time step."
                      << std::endl;
            //  exit(1);
          }
          if (prob > 0.5) {
            // std::cout << "WARNING: prob of reaction > 0.5. If this is a
            // reaction for a bimolecular binding with multiple binding sites,
            // please use a smaller time step." << std::endl;
          }

          if (prob > rand_gsl()) {
            // double rNum2 { rNum + rand_gsl() * Constants::iRandMax }; // to
            // get higher resolution if (prob > rNum2) { std::cout <<
            // "Destroying molecule of type "
            //           <<
            //           molTemplateList[moleculeList[molItr].molTypeIndex].molName
            //           << " with reaction "
            //           << oneRxn.absRxnIndex << " at iteration " << simItr <<
            //           '\n';

            complexList[moleculeList[molItr].myComIndex].destroy(moleculeList,
                                                                 complexList);
            // remove the molecule from the SimulVolume subsCellList
            // have this here to avoid circular header calls with SimulVolume
            // and Molecule_Complex
            for (auto itr =
                     simulVolume.subCellList[moleculeList[molItr].mySubVolIndex]
                         .memberMolList.begin();
                 itr !=
                 simulVolume.subCellList[moleculeList[molItr].mySubVolIndex]
                     .memberMolList.end();
                 ++itr) {
              if (*itr == moleculeList[molItr].index) {
                simulVolume.subCellList[moleculeList[molItr].mySubVolIndex]
                    .memberMolList.erase(itr);
                break;
              }
            }
            moleculeList[molItr].mySubVolIndex = -1;

            if (oneRxn.isObserved) {
              auto observeItr = observablesList.find(oneRxn.observeLabel);
              if (observeItr == observablesList.end()) {
                // std::cerr << "WARNING: Observable " << oneRxn.observeLabel <<
                // " not defined.\n";
              } else {
                --observeItr->second;
              }
            }
            hasDestroyed = true;
            break;  // don't do anything else with the molecule this timestep
                    // }
          }
        } else if (oneRxn.rxnType == ReactionType::uniMolCreation &&
                   moleculeList[molItr].trajStatus == TrajStatus::none &&
                   isReactant(moleculeList[molItr],
                              complexList[moleculeList[molItr].myComIndex],
                              oneRxn, moleculeList)) {
          // double rNum { 1.0 * rand_gsl() };
          double prob{first_order_probability(oneRxn.rateList.at(0).rate,
                                             params.timeStep)};

          // make sure that the time step is resonable according to the prob of
          // reaction
          if (prob > 1.000001) {
            std::cerr << "Error: prob of reaction > 1. Avoid this using a "
                         "smaller time step."
                      << std::endl;
            // exit(1);
          }
          if (prob > 0.5) {
            // std::cout << "WARNING: prob of reaction > 0.5. If this is a
            // reaction for a bimolecular binding with multiple binding sites,
            // please use a smaller time step." << std::endl;
          }

          if (prob > rand_gsl()) {
            int newMolIndex{0};
            int newComIndex{0};
            // std::cout << "Creating molecule of type " <<
            // oneRxn.productMolList.back().molName
            //           << " from reaction " << oneRxn.absRxnIndex << '\n';

            create_molecule_and_complex_from_rxn(
                molItr, newMolIndex, newComIndex, true,
                molTemplateList[oneRxn.productMolList.back().molTypeIndex],
                params, oneRxn, simulVolume, moleculeList, complexList,
                molTemplateList, forwardRxns, membraneObject);
            moleculeList[molItr].trajStatus = TrajStatus::propagated;
            complexList[moleculeList[molItr].myComIndex].trajStatus =
                TrajStatus::propagated;

            if (oneRxn.isObserved) {
              auto observeItr = observablesList.find(oneRxn.observeLabel);
              if (observeItr == observablesList.end()) {
                // std::cerr << "WARNING: Observable " << oneRxn.observeLabel <<
                // " not defined.\n";
              } else {
                ++observeItr->second;
              }
            }
            break;  // don't do anything else with the molecule this timestep
          }
        }
      }
    }

    if (hasDestroyed == true) continue;

    if (moleculeList[molItr].isImplicitLipid == true) continue;

    if (moleculeList[molItr].isDissociated == true) {
      continue;
    }

    // now check for dissociation
    check_dissociation(simItr, params, simulVolume, molTemplateList, observablesList, molItr, moleculeList, complexList,
                backRxns, forwardRxns, createDestructRxns, counterArrays, membraneObject,assocDissocFile);
  }
}
