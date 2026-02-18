/**
 * @file
 * nerdss_mpi.cpp
 *
 * @brief
 * Main function for simulation with parallel nerdss.
 *
 */

#ifndef mpi_
#define mpi_
#endif

#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "boundary_conditions/reflect_functions.hpp"
#include "debug/debug.hpp"
#include "error/error.hpp"
#include "io/io.hpp"
#include "macro.hpp"
#include "math/constants.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "mpi.h"
#include "mpi/mpi_function.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "split.cpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#ifdef ENABLE_PROFILING
#include <gperftools/profiler.h>
#endif

using namespace std;

/**
 * @def LOC(…)
 * @brief Macro to demark code sections for identification purposes. Comment out
 * to print LOC.
 */
#define LOC(...)

/**
 * @def PROC_ALL
 * @brief Macro for LOC argument to specify printing for every rank.
 */
#define PROC_ALL mpiContext.rank

/**
 * @def PROC_0
 * @brief Macro for LOC argument to specify printing for only rank 0.
 */
#define PROC_0 (mpiContext.rank == 0) ? -1 : -2

/**
 * @var gsl_rng* r
 * @brief Global random number generator.
 */
gsl_rng* r;

/**
 * @typedef MDTimer
 * @brief Alias for std::chrono::system_clock used for timing.
 */
using MDTimer = std::chrono::system_clock;

/**
 * @typedef timeDuration
 * @brief Alias for std::chrono::duration used for timing.
 */
using timeDuration = std::chrono::duration<double, std::chrono::seconds>;

/**
 * @var randNum
 * @brief Global variable to store a random number.
 */
long long randNum = 0;

/**
 * @var totMatches
 * @brief Global variable to store the total number of matches.
 */
unsigned long totMatches = 0;

/**
 * @brief Prints the elements of a vector and an associated message.
 *
 * @tparam T Type of the vector elements.
 * @param vec Vector to print.
 * @param s Message to print before the vector.
 */
template <typename T>
void printVec(vector<T>& vec, string s) {
  cout << s << ":" << endl;
  for (auto elem : vec) {
    cout << elem << "\t";
  }
  cout << endl;
}

/**
 * @brief Main function that serves as an entrypoint to the program.
 *
 * @param argc The number of command line arguments passed to the program.
 * @param argv An array of C-strings holding the command line arguments.
 * @return An integer indicating the exit status of the program.
 */
int main(int argc, char* argv[]) {
  // Initialize MPI
  // MPI::Init(argc, argv);
  MPI_Init(&argc, &argv);
  // Initialize the MPI context
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // MpiContext mpiContext(MPI::COMM_WORLD.Get_size(), NEIGHBOR_BUFFER_SIZE);
  MpiContext mpiContext(size, NEIGHBOR_BUFFER_SIZE);

  // Set the rank of the current MPI process
  // mpiContext.rank = MPI::COMM_WORLD.Get_rank();
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiContext.rank);

  double start_computation_time, end_computation_time, total_computation_time,
      start_communication_time, end_communication_time,
      total_communication_time;

#ifdef ENABLE_PROFILING
  std::string profileFileName{"profile_output_" +
                              std::to_string(mpiContext.rank) + ".prof"};
  ProfilerStart(profileFileName.c_str());
#endif

  mkdir("PDB", 0777);
  mkdir("RESTARTS", 0777);
  mkdir("DATA", 0777);
  mkdir("mergePDB", 0777);
  mkdir("mergeOUT", 0777);

  if (DEBUG) {
    LOC(1.0, " @ MAIN ", PROC_0);
  }

  // Initialize random_device object to generate seed for the generator
  std::random_device rd{};
  unsigned seed{rd()};
  // seed = 12345;

  // Start timer
  MDTimer::time_point totalTimeStart = MDTimer::now();

  // This will store the list of molecule templates
  std::vector<MolTemplate> molTemplateList{};

  // These will store the list of reactions
  std::vector<ForwardRxn> forwardRxns{};
  std::vector<BackRxn> backRxns{};
  std::vector<CreateDestructRxn> createDestructRxns{};
  std::vector<TransmissionRxn>
      transmissionRxns{}; // list of transmission reactions

  // Pre-allocate memory for the reactions vector to avoid reallocation during
  // runtime The `reserve()` function pre-allocates memory for the vector to
  // hold 10 elements.
  forwardRxns.reserve(10);
  backRxns.reserve(10);
  createDestructRxns.reserve(10);

  // This will store the system parameters
  Parameters params{};

  // Set filenames using string interpolation
  std::string observablesFileName{"DATA/observables_time_" +
                                  std::to_string(mpiContext.rank) + ".dat"};
  std::string trajFileName{"DATA/trajectory_" + std::to_string(mpiContext.rank) +
                           ".xyz"};
  std::string transitionFileName{"DATA/transition_matrix_time_" +
                                 std::to_string(mpiContext.rank) + ".dat"};
  std::string restartFileName{"DATA/restart_" + std::to_string(mpiContext.rank) +
                              ".dat"};
  std::string addFileNameInput{};
  std::string paramFile{};
  std::string coordinateFileName{};
  std::string restartFileNameInput;
  std::string debugFileName{"DATA/debug_output_" + std::to_string(mpiContext.rank) +
                            ".dat"};

  std::string outFileName{"DATA/output_" + std::to_string(mpiContext.rank) + ".txt"};
  if (freopen(outFileName.c_str(), "a", stdout) == nullptr) {
    std::cerr << "Failed to open the output file.\n";
    return 1;
  }

  // Create and write into debug output file
  std::ofstream debugFile{debugFileName};
  debugFile << "# This file is used for the debug output." << std::endl;
  debugFile.close();

  if (DEBUG) {
    LOC(2.0, " @ PARSE_COMMAND_LINE ", PROC_0);
  }

  // Command line flag parser
  parse_command(argc, argv, params, paramFile, restartFileNameInput,
                addFileNameInput, coordinateFileName, seed);

  if (DEBUG) {
    LOC(2.1, " @ INSTATIATE_FILE_SCOPE_VARIABLES ", PROC_0);
  }

  auto startTime = MDTimer::to_time_t(totalTimeStart);
  char charTime[24];
  std::cout << "\nStart date: ";
  if (0 <
      strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&startTime)))
    std::cout << charTime << '\n';
  std::cout << "RNG Seed: " << seed << std::endl;

  // Random number generator
  const gsl_rng_type* T;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);

  // 2D reaction probability tables, used to calculate reaction probability
  std::vector<gsl_matrix*> survMatrices;
  std::vector<gsl_matrix*> normMatrices;
  std::vector<gsl_matrix*> pirMatrices;
  double* tableIDs = new double[params.max2DRxns * 2];

  // list of the observables
  std::map<std::string, int> observablesList;

  // simulVolume object inlcuding the information of the division of the system
  // into subVolumes
  SimulVolume simulVolume{};

  // list of all molecules in the system
  std::vector<Molecule> moleculeList{};

  // list of all complexes in the system
  std::vector<Complex> complexList{};

  // class structure that contains boundary conditions, implicit lipid model.
  Membrane membraneObject;

  // contains arrays tracking bound molecule pairs, and species copy nums.
  copyCounters counterArrays;

  // implicit-lipid index, which is also stored in
  // membraneObject.implicitlipidIndex.
  int implicitlipidIndex{0};

  // simulation iteration count
  long long int simItr{0};

  int totalSpeciesNum{0};  ///< total species num after add

  // initialize event counters to zero. Restart will update any non-zero values.
  init_association_events(counterArrays);

  if (DEBUG) {
    LOC(2.2, " @ PARSE_NEW_SIMULATION_INPUT ", PROC_0);
  }

  std::cout << "\nParsing Input: " << std::endl;
  if (!params.fromRestart && paramFile != "") {
    parse_input_for_a_new_simulation(
        paramFile, params, observablesList, forwardRxns, backRxns,
        createDestructRxns, molTemplateList, membraneObject, moleculeList,
        complexList, simulVolume, observablesFileName, implicitlipidIndex,
        simItr, mpiContext, trajFileName, transitionFileName);
  } else if (params.fromRestart) {
    if (DEBUG) {
      LOC(2.3, " @ RESTART_INPUT ", PROC_0)
    };

    parse_input_for_a_restart_simulation(
        restartFileNameInput, addFileNameInput, params, observablesList,
        forwardRxns, backRxns, createDestructRxns, molTemplateList,
        membraneObject, moleculeList, complexList, simulVolume,
        observablesFileName, implicitlipidIndex, simItr, mpiContext,
        trajFileName, transitionFileName, seed, counterArrays);
  } else {
    error(
        "Please provide a parameter and/or restart file. Parameter file Syntax "
        "is : ./nerdss -f parameterfile.inp or ./nerdss -r restartfile.dat \n");
  }

  params.checkUnimoleculeReactionPopulation = false;
  mpiContext.checkUnimoleculeReactionPopulation =
      params.checkUnimoleculeReactionPopulation;
  params.hasRankCommunicationForLargeComplex = false;
  mpiContext.hasRankCommunicationForLargeComplex =
      params.hasRankCommunicationForLargeComplex;

  if (membraneObject.implicitLipid == true) params.implicitLipid = true;
  for (auto& oneReaction : forwardRxns) {
    if (oneReaction.rxnType == ReactionType::uniMolStateChange) {
      params.hasUniMolStateChange = true;
      break;
    }
  }
  if (createDestructRxns.empty() == false) {
    params.hasCreationDestruction = true;
    params.isNonEQ = true;

    // Set molTemp.canDestroy = ture
    for (auto& oneReaction : createDestructRxns) {
      if (oneReaction.rxnType == ReactionType::destruction) {
        molTemplateList[oneReaction.reactantMolList.at(0).molTypeIndex]
            .canDestroy = true;
      }
    }
  }

  // Setup output files
  char fnameProXYZ[100];

  std::string fileName{};
  fileName =
      "DATA/histogram_complexes_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream assemblyfile{fileName};
  fileName = "DATA/mono_dimer_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream dimerfile{fileName};
  fileName = "DATA/event_counters_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream eventFile{fileName};
  fileName = "DATA/bound_pair_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream pairOutfile{fileName};
  fileName = "DATA/copy_numbers_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream speciesFile1{fileName};
  fileName = "DATA/assoc_dissoc_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream assocDissocFile(fileName);
  if (params.assocDissocWrite == false) {
    assocDissocFile.close();
  }

  // Set the interval to write checkpoint and transition matrix file if it is
  // not specified
  if (params.checkPoint == -1) {
    params.checkPoint = params.nItr / 10;
  }
  if (params.transitionWrite == -1) {
    params.transitionWrite = params.nItr / 10;
  }

  set_exclude_volume_bound(forwardRxns, molTemplateList);

  for (auto& oneComplex : complexList) {
    oneComplex.update_properties(moleculeList, molTemplateList);
  }

  Parameters::dt = params.timeStep;
  Parameters::lastUpdateTransition.resize(molTemplateList.size());

  if (DEBUG) {
    LOC(3.0, " @ PARALLEL_PARTITION ", PROC_0);
  }

  if (mpiContext.nprocs) {  // if parallel execution
    if (!params.fromRestart) {
      prepare_data_structures_for_parallel_execution(
          moleculeList, simulVolume, membraneObject, molTemplateList, params,
          forwardRxns, backRxns, createDestructRxns, counterArrays, mpiContext,
          complexList, pairOutfile);
    }
  }

  int meanComplexSize{0};

  totalSpeciesNum = init_speciesFile(speciesFile1, counterArrays,
                                     molTemplateList, forwardRxns, params);
  init_counterCopyNums(counterArrays, moleculeList, complexList,
                       molTemplateList, membraneObject, totalSpeciesNum,
                       params);

  write_all_species((simItr - params.itrRestartFrom) * params.timeStep *
                            Constants::usToSeconds +
                        params.timeRestartFrom,
                    moleculeList, speciesFile1, counterArrays, membraneObject,
                    mpiContext, simulVolume);

  init_print_dimers(dimerfile, params, molTemplateList);
  init_NboundPairs(
      counterArrays, pairOutfile, params, molTemplateList,
      moleculeList);  // initializes to zero, re-calculated for a restart!!
  write_NboundPairs(counterArrays, pairOutfile, simItr, params, moleculeList);
  print_dimers(complexList, dimerfile, simItr, params, molTemplateList);
  print_association_events(counterArrays, eventFile, simItr, params);

  int number_of_lipids = 0;  ///< sum of all states of IL
  for (size_t i = 0; i < membraneObject.numberOfFreeLipidsEachState.size();
       i++) {
    number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
  }
  meanComplexSize = print_complex_hist(
      complexList, assemblyfile, simItr, params, molTemplateList, moleculeList,
      number_of_lipids, mpiContext, simulVolume);

  if (DEBUG) {
    LOC(4.0, " @ SIMULATION_PREP ", PROC_0);
  }

  // Print out system information
  std::cout << "\nSimulation Parameters\n";
  params.display();
  membraneObject.display();
  std::cout << "\nMolecule Information\n";
  display_all_MolTemplates(molTemplateList);
  std::cout << "\nReactions\n";
  display_all_reactions(forwardRxns, backRxns, createDestructRxns);

  std::cout << "*************** BEGIN SIMULATION **************** "
            << std::endl;

  //  begin the timer
  MDTimer::time_point simulTimeStart = MDTimer::now();
  int numSavedDurations{1000};
  std::vector<std::chrono::duration<double>> durationList(numSavedDurations);

  if (DEBUG) {
    LOC(4.1, " @ FILL_mList_GAPS ", PROC_0);
  }

  std::fill(durationList.begin(), durationList.end(),
            std::chrono::duration<double>(simulTimeStart - totalTimeStart));

  unsigned DDTableIndex{0};

  // Vectors to store binding probabilities for implicit lipids in 2D.
  std::vector<double> IL2DbindingVec{};
  std::vector<double> IL2DUnbindingVec{};
  std::vector<double> ILTableIDs{};

  if (params.fromRestart == true) {
    if (DEBUG) {
      LOC(4.1, " @ PRESERVE RNG STATE FROM RESTART", PROC_0);
    }

    try {
      read_rng_state(mpiContext.rank);  // read the current RNG state
    } catch (std::exception& e) {
      std::cout << "Read rng_state failed. Set the rng_state with the seed: "
                << seed << std::endl;
      gsl_rng_set(r, seed);
    }

    std::ofstream restartFile{restartFileName,
                              std::ios::out};  // to show different from append

    write_rng_state(mpiContext.rank);  // write the current RNG state
    write_restart(simItr, restartFile, params, simulVolume, moleculeList,
                  complexList, molTemplateList, forwardRxns, backRxns,
                  createDestructRxns, transmissionRxns, observablesList, membraneObject,
                  counterArrays);
    restartFile.close();
  }

  long long int startSimItr = simItr + 1;
  long long int stopSimItr = params.nItr - 1;

  if (DEBUG) {
    LOC(4.2, " @ START SIMULATION_LOOP ", PROC_ALL);
  }

  // write the pdb of the initial system
  write_pdb(simItr, simItr, params, moleculeList, molTemplateList,
            membraneObject, mpiContext.rank, mpiContext, simulVolume);

  for (simItr += 1; simItr <= params.nItr; ++simItr) {
    if (1) {
      // if (mpiContext.rank == 0) {
      start_computation_time = MPI_Wtime();
    }

    if (VERBOSE) {
      printf("simItr (%lld) start on rank (%d)...\n", simItr, mpiContext.rank);
    }

    mpiContext.simItr = simItr;

    for (auto& oneMol : moleculeList) {
      if (oneMol.isEmpty || oneMol.isImplicitLipid) continue;

      oneMol.trajStatus = TrajStatus::none;

      complexList[oneMol.myComIndex].ncross = 0;
      complexList[oneMol.myComIndex].trajStatus = TrajStatus::none;
    }
    if (VERBOSE) {
      printf(
          "Reset the TrajStatus of mol and TrajStatus and ncross of complex "
          "done.\n");
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("4.3 simItr BEGINS");
      DEBUG_FIND_COMPLEX("4.3 simItr BEGINS");
    }

    if (DEBUG) {
      LOC(4.3, " @ PARALLEL_RECEIVE_DATA ", PROC_ALL, simItr);
    }

    propCalled = 0;
    MDTimer::time_point startStep = MDTimer::now();
    // destruct, unimol create, and dissociation (explicit) based on population

    if (DEBUG) {
      LOC(5.1, " @ CX_UNIMOL_RXN_POP ", PROC_ALL, simItr);
    }

    if (VERBOSE) {
      printf("Check and perform 0th and 1st order reactions...\n");
    }
    check_perform_zeroth_first_order_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
        ILTableIDs, mpiContext);

    if (DEBUG) {
      DEBUG_FIND_MOL("4.5");
      DEBUG_FIND_COMPLEX("4.5");
    }

    // Measure separations between proteins in neighboring cells to identify all
    // possible reactions.
    if (DEBUG) {
      LOC(6.1, " @ MEASURE_SEPARATIONS-CX_IMPLICIT/BIMOL_RXNS ", PROC_ALL,
          simItr);
    }
    if (VERBOSE) {
      printf("Measure separations between two reactants...\n");
    }

    simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                      molTemplateList, membraneObject, simItr,
                                      mpiContext);

    measure_separations_to_identify_possible_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
        ILTableIDs, normMatrices, survMatrices, pirMatrices, implicitlipidIndex,
        tableIDs, DDTableIndex);

    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//(5)");
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("5 (End Measure separations)");
      DEBUG_FIND_COMPLEX("5 (End Measure separations)");
    }

    // Now that separations and reaction probabilities are calculated, decide
    // whether to perform reactions for each protein.
    if (DEBUG) {
      LOC(7.0, " @ REACTIONS ", PROC_ALL, simItr);
    }

    // Divide the subVolume into two parts and execute each part
    // synchronously with all processes
    if (VERBOSE) {
      printf("Assgin molecules to left and right parts...\n");
    }
    std::vector<int> left{};
    std::vector<int> right{};
    for (int i = 0; i < moleculeList.size(); i++) {
      // reset the need_to_send flag to false
      moleculeList[i].need_to_send = false;

      int x_bin = get_x_bin(mpiContext, moleculeList[i]);
      if (mpiContext.rank == 0) {
        if (x_bin < (simulVolume.numSubCells.x - 1) / 2)
          left.push_back(i);
        else
          right.push_back(i);
      } else if (mpiContext.rank == mpiContext.nprocs - 1) {
        if (x_bin < (simulVolume.numSubCells.x + 1) / 2)
          left.push_back(i);
        else
          right.push_back(i);
      } else {
        if (x_bin < (simulVolume.numSubCells.x) / 2)
          left.push_back(i);
        else
          right.push_back(i);
      }
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("5.5 (Before update molecules)");
      DEBUG_FIND_COMPLEX("5.5 (Before update molecules)");
    }
    // do the left part, send to left and recieve from right
    if (VERBOSE) {
      printf("Update the states of molecules in the left part...\n");
    }

    perform_bimolecular_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, left);

    if (DEBUG) {
      DEBUG_FIND_MOL("6 (End decide whether to perform reactions 1)");
      DEBUG_FIND_COMPLEX("6 (End decide whether to perform reactions 1)");
    }
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//(6)");
    }

    if (DEBUG) {
      LOC(7.1, " @ OVERLAP_SWEEP ", PROC_ALL, simItr);
    }

    if (VERBOSE) {
      printf("Checking overlap and propagate molecules in the left part...\n");
    }
    check_overlap(left, simItr, params, moleculeList, complexList, simulVolume,
                  forwardRxns, backRxns, createDestructRxns, molTemplateList,
                  observablesList, counterArrays, membraneObject, mpiContext);

    if (DEBUG) {
      LOC(8.0, " @ AFTER_PROPAGATE ", PROC_ALL, simItr);
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("7 (after propagate)");
      DEBUG_FIND_COMPLEX("7 (after propagate)");
    }
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//(7)");
    }

    if (DEBUG) {
      LOC(8.2, " @ CX_MOL_COORDS ", PROC_ALL, simItr);
    }

    if (VERBOSE) {
      printf(
          "Checking coords of molecules in the left part to determine if need "
          "communication...\n");
    }
    // check_molecule_coordinates(mpiContext, moleculeList, left, complexList,
    //                            molTemplateList, simulVolume, membraneObject,
    //                            counterArrays);
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//(10)");
    }

    if (DEBUG) {
      LOC(8.3, " @ COM_MGT1 ", PROC_ALL, simItr);
    }

    for (auto& oneMolIndex : left) {
      auto& oneMol{moleculeList[oneMolIndex]};
      if (oneMol.isEmpty || oneMol.isImplicitLipid) continue;

      clear_reweight_vecs(oneMol);
      // oneMol.trajStatus = TrajStatus::none;
      oneMol.crossbase.clear();
      oneMol.mycrossint.clear();
      oneMol.crossrxn.clear();
      oneMol.probvec.clear();
    }

    if (DEBUG) {
      LOC(8.4, " @ PARALLEL_SEND_DATA ", PROC_ALL, simItr);
    }

    // Sending data to neighboring ranks:
    if (DEBUG) {
      debug_firstEmptyIndex(mpiContext, "BEFORE SEND");
    }
    if (DEBUG) {
      DEBUG_FIND_MOL("BEFORE SEND");
      DEBUG_FIND_COMPLEX("BEFORE SEND");
    }

    simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                      molTemplateList, membraneObject, simItr,
                                      mpiContext);

    if (1) {
      // if (mpiContext.rank == 0) {
      end_computation_time = MPI_Wtime();
      total_computation_time += (end_computation_time - start_computation_time);
      start_communication_time = MPI_Wtime();
    }

    // send to left
    if (mpiContext.rank > 0) {
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d)...\n", mpiContext.rank,
               mpiContext.rank - 1);
      }
      send_data_to_left_neighboring_ranks(
          mpiContext, simItr, params, simulVolume, left, moleculeList,
          complexList, molTemplateList, membraneObject);
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d) done.\n",
               mpiContext.rank, mpiContext.rank - 1);
      }
    }

    if (DEBUG) {
      check_complex_index(mpiContext, moleculeList, complexList);
    }

    // recieve from right
    if (mpiContext.rank < mpiContext.nprocs - 1) {
      if (VERBOSE) {
        printf("Receiving data from rank (%d) to rank (%d)...\n",
               mpiContext.rank + 1, mpiContext.rank);
      }
      receive_right_neighborhood_zones(
          mpiContext, simItr, simulVolume, right, moleculeList, complexList,
          molTemplateList, membraneObject, counterArrays);
      if (VERBOSE) {
        printf("Receiving data from rank (%d) to rank (%d) done\n",
               mpiContext.rank + 1, mpiContext.rank);
      }
    }

    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//after receive right");
      check_complex_index(mpiContext, moleculeList, complexList);
    }

    if (1) {
      // if (mpiContext.rank == 0) {
      end_communication_time = MPI_Wtime();
      total_communication_time +=
          (end_communication_time - start_communication_time);
      start_computation_time = MPI_Wtime();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    // we need to reassign the molecules to left and right here because the
    // right region members may be changed after the communication
    left.clear();
    right.clear();
    for (int i = 0; i < moleculeList.size(); i++) {
      // reset the need_to_send flag to false
      moleculeList[i].need_to_send = false;

      int x_bin = get_x_bin(mpiContext, moleculeList[i]);
      if (mpiContext.rank == 0) {
        if (x_bin < (simulVolume.numSubCells.x - 1) / 2)
          left.push_back(i);
        else
          right.push_back(i);
      } else if (mpiContext.rank == mpiContext.nprocs - 1) {
        if (x_bin < (simulVolume.numSubCells.x + 1) / 2)
          left.push_back(i);
        else
          right.push_back(i);
      } else {
        if (x_bin < (simulVolume.numSubCells.x) / 2)
          left.push_back(i);
        else
          right.push_back(i);
      }
    }

    // do the right part, send to right and recieve from left
    if (VERBOSE) {
      printf("Update the states of molecules in the right part...\n");
    }
    perform_bimolecular_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, right);

    if (DEBUG) {
      DEBUG_FIND_MOL("6 (End decide whether to perform reactions 2)");
      DEBUG_FIND_COMPLEX("6 (End decide whether to perform reactions 2)");
    }
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//(6)");
    }

    if (DEBUG) {
      LOC(7.1, " @ OVERLAP_SWEEP ", PROC_ALL, simItr);
    }

    if (VERBOSE) {
      printf("Checking overlap and propagate molecules in the right part...\n");
    }
    check_overlap(right, simItr, params, moleculeList, complexList, simulVolume,
                  forwardRxns, backRxns, createDestructRxns, molTemplateList,
                  observablesList, counterArrays, membraneObject, mpiContext);

    if (DEBUG) {
      LOC(8.0, " @ AFTER_PROPAGATE ", PROC_ALL, simItr);
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("7 (after propagate)");
      DEBUG_FIND_COMPLEX("7 (after propagate)");
    }
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//(7)");
    }

    if (DEBUG) {
      LOC(8.2, " @ CX_MOL_COORDS ", PROC_ALL, simItr);
    }

    if (VERBOSE) {
      printf(
          "Checking coords of molecules in the right part to determine if need "
          "communication...\n");
    }
    // check_molecule_coordinates(mpiContext, moleculeList, right, complexList,
    //                            molTemplateList, simulVolume, membraneObject,
    //                            counterArrays);
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "//(10)");
    }

    // Clear lists used for reweighting and encounter tracking

    if (DEBUG) {
      LOC(8.3, " @ COM_MGT1 ", PROC_ALL, simItr);
    }

    for (auto& oneMolIndex : right) {
      auto& oneMol{moleculeList[oneMolIndex]};
      if (oneMol.isEmpty || oneMol.isImplicitLipid) continue;

      clear_reweight_vecs(oneMol);
      oneMol.trajStatus = TrajStatus::none;
      oneMol.crossbase.clear();
      oneMol.mycrossint.clear();
      oneMol.crossrxn.clear();
      oneMol.probvec.clear();
      oneMol.justBoundThisStep = false;
      oneMol.justUnboundThisStep = false;

      // update complexes. if the complex has more than one member, it'll be
      // updated more than once, but this is probably more efficient than
      // iterating over the complexes afterwards
      complexList[oneMol.myComIndex].ncross = 0;
      complexList[oneMol.myComIndex].trajStatus = TrajStatus::none;
    }

    for (auto& oneMolIndex : left) {
      auto& oneMol{moleculeList[oneMolIndex]};
      if (oneMol.isEmpty || oneMol.isImplicitLipid) continue;

      oneMol.trajStatus = TrajStatus::none;
      oneMol.justBoundThisStep = false;
      oneMol.justUnboundThisStep = false;

      complexList[oneMol.myComIndex].ncross = 0;
      complexList[oneMol.myComIndex].trajStatus = TrajStatus::none;
    }

    if (DEBUG) {
      LOC(8.4, " @ PARALLEL_SEND_DATA ", PROC_ALL, simItr);
    }

    // Sending data to neighboring ranks:
    if (DEBUG) {
      debug_firstEmptyIndex(mpiContext, "BEFORE SEND");
    }
    if (DEBUG) {
      DEBUG_FIND_MOL("BEFORE SEND");
      DEBUG_FIND_COMPLEX("BEFORE SEND");
    }

    simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                      molTemplateList, membraneObject, simItr,
                                      mpiContext);

    if (1) {
      // if (mpiContext.rank == 0) {
      end_computation_time = MPI_Wtime();
      total_computation_time += (end_computation_time - start_computation_time);
      start_communication_time = MPI_Wtime();
    }

    // send to right
    if (mpiContext.rank < mpiContext.nprocs - 1) {
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d)...\n", mpiContext.rank,
               mpiContext.rank + 1);
      }
      send_data_to_right_neighboring_ranks(
          mpiContext, simItr, params, simulVolume, right, moleculeList,
          complexList, molTemplateList, membraneObject);
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d) done.\n",
               mpiContext.rank, mpiContext.rank + 1);
      }
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("AFTER SEND");
      DEBUG_FIND_COMPLEX("AFTER SEND");
    }
    if (DEBUG) {
      debug_firstEmptyIndex(mpiContext, "AFTER SEND");
    }
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "// AFTER SEND");
    }

    // recieve from left
    if (mpiContext.rank > 0) {
      if (VERBOSE) {
        printf("Receiving data from rank (%d) to rank (%d)...\n",
               mpiContext.rank - 1, mpiContext.rank);
      }
      receive_left_neighborhood_zones(
          mpiContext, simItr, simulVolume, left, moleculeList, complexList,
          molTemplateList, membraneObject, counterArrays);
      if (VERBOSE) {
        printf("Receiving data from rank (%d) to rank (%d) done\n",
               mpiContext.rank - 1, mpiContext.rank);
      }
    }

    if (1) {
      //  if (mpiContext.rank == 0) {
      end_communication_time = MPI_Wtime();
      total_communication_time +=
          (end_communication_time - start_communication_time);
    }

    // LOC( 8.6, " @ COM_MGT2 ", PROC_ALL, simItr);

    if (DEBUG) {
      LOC(9.0, " @ IO_TRAJ_+_PDB_+_TMATRIX ", PROC_ALL, simItr);
    }

    // debug_bndpartner_interface(mpiContext, "34: (before remove empty
    // complexes)");

    if (DEBUG) {
      LOC(9.2, " @ RM_EMPTY_COMPLEXES ", PROC_ALL, simItr);
    }

    remove_empty_slots(simItr, params, moleculeList, complexList, simulVolume,
                       forwardRxns, backRxns, createDestructRxns,
                       molTemplateList, observablesList, counterArrays,
                       membraneObject, mpiContext);

    if (DEBUG) {
      DEBUG_FIND_MOL("7: before 4.10");
      DEBUG_FIND_COMPLEX("7: before 4.10");
    }
    if (DEBUG) {
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "// before 4.10");
    }

    //------------------------------------------------------------------------------------

    if (DEBUG) {
      LOC(9.4, " @ CLEAR_REWEIGHTING ", PROC_ALL, simItr);
    }

    if (DEBUG) {
      LOC(9.5, " @ IO_RESTART ", PROC_ALL, simItr);
    }

    write_output(simItr, params, trajFileName, moleculeList, molTemplateList,
                 membraneObject, mpiContext, transitionFileName, counterArrays,
                 speciesFile1, debugFile, debugFileName, restartFileName,
                 complexList, simulVolume, forwardRxns, backRxns,
                 createDestructRxns, durationList, startStep, totalTimeStart,
                 observablesList, assemblyfile);

    if (simItr % params.timeWrite == 0) {
      std::cout << "total computation time: " << total_computation_time
                << " seconds\n";
      std::cout << "total communication time: " << total_communication_time
                << " seconds\n";

      // MERGE THE OUTPUTS FROM ALL THE PROCESSORS
      // if (mpiContext.rank == 0) {
      //   merge_outputs(mpiContext.nprocs, molTemplateList.size());
      // }
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("8: end iterating over time steps");
      DEBUG_FIND_COMPLEX("8: end iterating over time steps");
    }
  }  // end iterating over time steps

  // Write out final result
  std::cout << llinebreak << "End simulation\n";
  auto endTime = MDTimer::now();
  auto endTimeFormat = MDTimer::to_time_t(endTime);
  std::cout << "End date: ";
  //<< std::put_time(std::localtime(&endTimeFormat), "%F %T") << '\n';
  if (0 < strftime(charTime, sizeof(charTime), "%F %T",
                   std::localtime(&endTimeFormat)))
    std::cout << charTime << '\n';
  std::chrono::duration<double> wallTime = endTime - totalTimeStart;
  std::cout << "\tWall Time: ";
  std::cout << wallTime.count() << " seconds\n";
  std::cout << "total computation time: " << total_computation_time
            << " seconds\n";
  std::cout << "total communication time: " << total_communication_time
            << " seconds\n";

  // MERGE THE OUTPUTS FROM ALL THE PROCESSORS
  if (mpiContext.rank == 0) {
    merge_outputs(mpiContext.nprocs, molTemplateList.size());
  }

  delete[] tableIDs;
  gsl_rng_free(r);

#ifdef ENABLE_PROFILING
  ProfilerStop();
#endif

  MPI_Finalize();
  return 0;
}  // end main
