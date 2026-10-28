#pragma once

#include <array>
#include <cstring>
#include <map>
#include <string>
#include <set>
#include <unordered_map>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#ifdef mpi_
#include "mpi.h"
#endif
#include "mpi_context.hpp"

using namespace std;

/*! \defgroup mpi_function
 * \brief Functions for mpi implementation.
 */

/*! \ingroup mpi_function
 * \brief move complex
 */
void move_complexes_based_on_propagation(MpiContext &mpiContext, vector<Molecule> &moleculeList, SimulVolume &simulVolume, Membrane &membraneObject);

/*! \ingroup mpi_function
 * \brief Search for the molecule with the same ID within moleculeList
 */
int find_molecule(vector<Molecule> &moleculeList, int id);

/*! \ingroup mpi_function
 * \brief Verify that every molecule is owned by exactly one rank.
 *
 * Reports ids owned by two ranks and ids held only as ghosts, and returns false
 * if either set is non-empty.  Collective: every rank must call it.
 */
bool trace_id(int id);
void trace_molecule(MpiContext &mpiContext, const char *where, const std::vector<Molecule> &moleculeList,
                    const std::vector<Complex> &complexList, const Molecule &mol);
bool check_ownership_invariant(MpiContext &mpiContext,
                               vector<Molecule> &moleculeList,
                               vector<Complex> &complexList, long long simItr,
                               const char *label);

/*! \ingroup mpi_function
 * \brief Set every molecule's isGhosted from its complex's ownerRank: a ghost
 * is a member of a complex another rank owns.  Called after each receive and
 * after each bimolecular pass, the two places where complex membership or
 * ownership changes.  With skipUnclaimed, members of a complex nobody has
 * claimed yet (ownerRank -1) keep their flag; see the definition.
 */
void derive_ghost_from_ownership(MpiContext &mpiContext,
                                 vector<Molecule> &moleculeList,
                                 vector<Complex> &complexList,
                                 bool skipUnclaimed);

/*! \ingroup mpi_function
 * \brief Search for the complex with the same ID within complexList
 */
int find_complex(vector<Complex> &complexList, int id);

/*! \ingroup mpi_function
 * \brief update the copycounters when complex enter ghosted zone
 */
void update_copyCounters_enter_ghosted_zone(MpiContext &mpiContext, Molecule &mol, vector<MolTemplate> &molTemplateList,
    copyCounters &counterArrays, vector<Molecule> &moleculeList, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the copycounters when complex exit ghosted zone
 */
void update_copyCounters_exit_ghosted_zone(MpiContext &mpiContext, Molecule &mol, vector<MolTemplate> &molTemplateList,
    copyCounters &counterArrays, vector<Molecule> &moleculeList, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the properties of the molecule when it enters the ghosted zone
 */
void enter_ghosted_zone(MpiContext &mpiContext, Molecule &mol, vector<MolTemplate> &molTemplateList, copyCounters &counterArrays,
    vector<Molecule> &moleculeList, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the properties of the molecule when it exits the ghosted zone
 */
void exit_ghosted_zone(MpiContext &mpiContext, Molecule &mol, vector<MolTemplate> &molTemplateList, copyCounters &counterArrays,
    vector<Molecule> &moleculeList, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the properties of the molecule when it moves across the shared
 * zones
 */
void check_molecule_coordinates(MpiContext &mpiContext, vector<Molecule> &moleculeList, vector<int> &molIndexList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, SimulVolume &simulVolume, Membrane &membraneObject, copyCounters &counterArrays);

/*! \ingroup mpi_function
 * \brief map local index to unique global id
 */
void indices_to_IDs(Molecule &mol, vector<Molecule> &moleculeList, vector<Complex> &complexList);

/*! \ingroup mpi_function
 * \brief map mySubVolIndex of one molecule
 */
void update_mySubVolIndex(MpiContext &mpiContext, Molecule &mol, Membrane &membraneObject, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief check if a recieved mol enter the ghosted zone
 */
bool check_enter_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew, vector<int> &indicesEnteringMolecules, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief check if a recieved mol exit the ghosted zone
 */
bool check_exit_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew, vector<int> &indicesEnteringMolecules, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief serialize molecules
 */
void serialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume, vector<Molecule> &moleculeList, vector<Complex> &complexList, unsigned char *arrayRank,
    int &nArrayRank, set<int> &complexesSet, set<int> &moleculesSet);

/*! \ingroup mpi_function
 * \brief serialize some molecules
 */
void serialize_molecules_partial(MpiContext &mpiContext, SimulVolume &simulVolume, vector<int> &indices, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, unsigned char *arrayRank, int &nArrayRank, set<int> &complexesSet, int cellX);

/*! \ingroup mpi_function
 * \brief What a molecule looked like at one point of the step: which complex
 * it was in and how many of its interfaces were bound.
 *
 * The second bimolecular pass is the only thing between the two exchanges of
 * a step that changes a bond or moves a molecule between complexes, so a
 * molecule whose entry differs from a snapshot taken before that pass was
 * changed by it.  trajStatus cannot tell: the left region keeps "propagated"
 * from the first half, and every flag is reset before the second send.
 */
struct MolSnapshot {
  int comId{-1};   //!< id of the molecule's complex, -1 for none
  int nBound{0};   //!< interfaces whose partnerIndex is set
};

/*! \ingroup mpi_function
 * \brief Record every molecule's MolSnapshot, indexed like moleculeList
 */
void snapshot_molecules(const vector<Molecule> &moleculeList, const vector<Complex> &complexList, vector<MolSnapshot> &snapshot);

/*! \ingroup mpi_function
 * \brief deserialize molecules
 */
void deserialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume, vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject, copyCounters &counterArrays, vector<int> &indices, vector<int> &indicesEnteringMolecules,
    vector<int> &indicesExitingMolecules, unsigned char *arrayRank, int &nArrayRank, const vector<MolSnapshot> &sentLeft,
    vector<char> &keptLocalCom, vector<char> &keptLocalMol);

void deserialize_molecules_right(MpiContext &mpiContext, SimulVolume &simulVolume, vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject, copyCounters &counterArrays, vector<int> &indices,
    vector<int> &indicesEnteringMolecules, vector<int> &indicesExitingMolecules, unsigned char *arrayRank, int &nArrayRank,
    const unordered_map<int, MolSnapshot> &sentRight, vector<char> &keptLocalCom, vector<char> &keptLocalMol);

/*! \ingroup mpi_function
 * \brief Serialize complexes from complexList with indices within complexesSet
 */
void serialize_complexes(MpiContext &mpiContext, SimulVolume &simulVolume, Membrane &membraneObject, set<int> &complexesSet,
    vector<Complex> &complexList, vector<Molecule> &moleculeList, unsigned char *arrayRank, int &nArrayRank);

/*! \ingroup mpi_function
 * \brief Deserialize complexes into complexList if they don't exist
 */
void deserialize_complexes(MpiContext &mpiContext, vector<Molecule> &moleculeList, vector<Complex> &complexList, unsigned char *arrayRank, int &nArrayRank,
    const vector<char> &keptLocalCom, const vector<char> &keptLocalMol);

void deserialize_complexes_right(MpiContext &mpiContext, vector<Molecule> &moleculeList, vector<Complex> &complexList, unsigned char *arrayRank, int &nArrayRank,
    const vector<char> &keptLocalCom, const vector<char> &keptLocalMol);

/*! \ingroup mpi_function
 * \brief Map global unique id to local index
 */
void IDs_to_indices(MpiContext &mpiContext, vector<Molecule> &moleculeList, vector<int> &indices);

/*! \ingroup mpi_function
 * \brief verify molecule/complex/partner links agree in both directions
 */
bool check_structure_invariant(MpiContext &mpiContext, std::vector<Molecule> &moleculeList,
                               std::vector<Complex> &complexList, long long simItr,
                               const char *label, int mode = 0);

/*! \ingroup mpi_function
 * \brief give an owner to every complex created since the last start-up seed
 */
void claim_unowned_complexes(MpiContext &mpiContext, vector<Molecule> &moleculeList,
                             vector<Complex> &complexList, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief hold the complexes the right neighbour can bind into back from the
 * top-of-step first-order pass (marks their members ghosts; see
 * defer_edge_reactions.cpp)
 */
void hold_back_right_edge_first_order(MpiContext &mpiContext,
                                      vector<Molecule> &moleculeList,
                                      vector<Complex> &complexList,
                                      vector<int> &heldBack);

/*! \ingroup mpi_function
 * \brief undo hold_back_right_edge_first_order()'s ghost marks
 */
void release_held_back(vector<Molecule> &moleculeList,
                       const vector<int> &heldBack);

/*! \ingroup mpi_function
 * \brief the first-order pass for the held-back complexes, after the first
 * exchange
 */
void perform_held_back_first_order(
    long long simItr, Parameters &params, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, SimulVolume &simulVolume,
    vector<ForwardRxn> &forwardRxns, vector<BackRxn> &backRxns,
    vector<CreateDestructRxn> &createDestructRxns,
    vector<MolTemplate> &molTemplateList,
    std::map<std::string, int> &observablesList, copyCounters &counterArrays,
    Membrane &membraneObject, vector<double> &IL2DbindingVec,
    vector<double> &IL2DUnbindingVec, vector<double> &ILTableIDs,
    MpiContext &mpiContext, const vector<int> &heldBackIds,
    std::ofstream &assocDissocFile);

/*! \ingroup mpi_function
 * \brief write the whole system's state as one restart file from rank 0, in
 * the serial format (src/io_mpi/write_gathered_restart.cpp); collective.
 * rngWrite: 0 no RNG state, 1 DATA/rng_state, 2 RESTARTS/rng_state<simItr>
 */
void write_gathered_restart(
    long long int simItr, const std::string &fileName, int rngWrite,
    Parameters &params, std::vector<Molecule> &moleculeList,
    std::vector<Complex> &complexList, std::vector<MolTemplate> &molTemplateList,
    std::vector<ForwardRxn> &forwardRxns, std::vector<BackRxn> &backRxns,
    std::vector<CreateDestructRxn> &createDestructRxns,
    const std::map<std::string, int> &observablesList,
    const std::map<std::string, int> &observablesAtStart,
    Membrane &membraneObject, copyCounters &counterArrays,
    const Membrane &membraneWhole, const SimulVolume &simulVolumeWhole,
    MpiContext &mpiContext);

/*! \ingroup mpi_function
 * \brief write the bonded complexes of the whole system as the serial build's
 * DATA/COMPLEXES/<step>.json, from rank 0 (src/io_mpi/write_gathered_restart.cpp);
 * collective
 */
void write_gathered_bonded_complexes(
    long long int simItr, const std::string &fileName,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    const std::vector<MolTemplate> &molTemplateList,
    const std::vector<ForwardRxn> &forwardRxns,
    const std::vector<BackRxn> &backRxns, const Membrane &membraneObject,
    const Membrane &membraneWhole, const SimulVolume &simulVolumeWhole,
    MpiContext &mpiContext);

/*! \ingroup mpi_function
 * \brief drop every stored reference to a molecule that has left this rank
 */
void drop_references_to_absent_molecules(vector<Molecule> &moleculeList,
                                         vector<Complex> &complexList);

/*! \ingroup mpi_function
 * \brief disconnect molecule partners
 */
void disconnect_molecule_partners(unsigned &targMolIndex, Molecule &mol, vector<Molecule> &moleculeList);

/*! \ingroup mpi_function
 * \brief The x bin at which this rank's columns split into the left region
 * (bins below it, acted on in the first half of the step) and the right
 * region (the rest, acted on in the second half).  One rule for the driver's
 * two assignment loops.
 */
int region_mid_index(MpiContext &mpiContext);

/*! \ingroup mpi_function
 * \brief Drop every candidate reaction between a molecule and a partner in
 * the right ghost column.  A pair across a column boundary is decided by the
 * rank to the right of it, in its first half; see region_split.cpp.
 */
void drop_crossings_decided_elsewhere(MpiContext &mpiContext, vector<Molecule> &moleculeList, vector<Complex> &complexList);

/*! \ingroup mpi_function
 * \brief candidate pairs this rank can act on, summed over ranks on rank 0
 */
std::array<long long, 2> count_candidate_pairs(MpiContext &mpiContext, vector<Molecule> &moleculeList);

/*! \ingroup mpi_function
 * \brief reweighting-history hit counts, edge columns vs interior, summed
 * over ranks on rank 0
 */
std::array<long long, 4> reweight_history_hits(MpiContext &mpiContext, vector<Molecule> &moleculeList);

/*! \ingroup mpi_function
 * \brief Delete all molecules that were in the shared zone before sending to
 * the neighbor rank, but were not received back:
 */
void delete_disappeared_molecules(MpiContext &mpiContext, SimulVolume &simulVolume, Membrane &membraneObject,
    vector<Molecule> &moleculeList, vector<Complex> &complexList, vector<int> &region, vector<MolTemplate> &molTemplateList,
    bool fromLeft);

/*! \ingroup mpi_function
 * \brief Delete all complexes that were in the shared zone before sending to
 * the neighbor rank, but were not received back:
 */
void delete_disappeared_complexes_partial(MpiContext &mpiContext, vector<Molecule> &moleculeList, vector<Complex> &complexList);

/*! \ingroup mpi_function
 * \brief Delete all complexes that were in the shared zone before sending to
 * the neighbor rank, but were not received back:
 */
void delete_disappeared_complexes(MpiContext &mpiContext, vector<Molecule> &moleculeList, vector<Complex> &complexList);

/*! \ingroup mpi_function
 * \brief Send data to left neighbor rank
 */
void send_data_to_left_neighboring_ranks(MpiContext &mpiContext, long long int simItr, Parameters &params, SimulVolume &simulVolume, 
    vector<int> &left, vector<Molecule> &moleculeList, vector<Complex> &complexList, vector<MolTemplate> &molTemplateList, Membrane &membraneObject);

/*! \ingroup mpi_function
 * \brief Send data to right neighbor rank
 */
void send_data_to_right_neighboring_ranks(MpiContext &mpiContext, long long int simItr, Parameters &params, SimulVolume &simulVolume, vector<int> &right,
    vector<Molecule> &moleculeList, vector<Complex> &complexList, vector<MolTemplate> &molTemplateList, Membrane &membraneObject,
    unordered_map<int, MolSnapshot> &sentRight);

/*! \ingroup mpi_function
 * \brief Receive shared zones from right neighboring rank, deserialize received
 * molecules, update interfaces to match molecules on this rank, and delete
 * those molecules that were sent, but not received back
 */
void receive_right_neighborhood_zones(MpiContext &mpiContext, long long int simItr, SimulVolume &simulVolume, vector<int> &right, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList, Membrane &membraneObject, copyCounters &counterArrays,
    const unordered_map<int, MolSnapshot> &sentRight);

/*! \ingroup mpi_function
 * \brief Receive shared zones from left neighboring rank, deserialize received
 * molecules, update interfaces to match molecules on this rank, and delete
 * those molecules that were sent, but not received back
 */
void receive_left_neighborhood_zones(MpiContext &mpiContext, long long int simItr, SimulVolume &simulVolume, vector<int> &left, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList, Membrane &membraneObject, copyCounters &counterArrays,
    const vector<MolSnapshot> &sentLeft);

/*! \ingroup mpi_function
 * \brief Prepare rank data
 */
int prepare_rank_data(int tempRank, vector<Molecule> &moleculeList, SimulVolume &simulVolume, Membrane &membraneObject, vector<MolTemplate> &molTemplateList, Parameters &params,
    vector<ForwardRxn> &forwardRxns, vector<BackRxn> &backRxns, vector<CreateDestructRxn> &createDestructRxns, copyCounters &counterArrays, MpiContext &mpiContext,
    vector<Complex> &complexList, ofstream &pairOutfile, unsigned char *arrayRank);

/*! \ingroup mpi_function
 * \brief Prepare data structures for parallel execution
 */
void prepare_data_structures_for_parallel_execution(vector<Molecule> &moleculeList, SimulVolume &simulVolume, Membrane &membraneObject, vector<MolTemplate> &molTemplateList,
    Parameters &params, vector<ForwardRxn> &forwardRxns, vector<BackRxn> &backRxns, vector<CreateDestructRxn> &createDestructRxns, copyCounters &counterArrays, MpiContext &mpiContext,
    vector<Complex> &complexList, ofstream &pairOutfile);

/*! \ingroup mpi_function
 * \brief Set a rank up to continue from a restart file (one rank only)
 */
void prepare_data_structures_for_parallel_restart(vector<Molecule> &moleculeList, SimulVolume &simulVolume, Membrane &membraneObject,
    vector<MolTemplate> &molTemplateList, Parameters &params, vector<ForwardRxn> &forwardRxns, vector<BackRxn> &backRxns,
    vector<CreateDestructRxn> &createDestructRxns, copyCounters &counterArrays, MpiContext &mpiContext, vector<Complex> &complexList,
    ofstream &pairOutfile);

/*! \ingroup mpi_function
 * \brief Used for the writing of output files
 */
bool is_ghosted(Molecule &mol, MpiContext &mpiContext, SimulVolume &simulVolume);

bool is_owned_by_processor(Molecule &mol, MpiContext &mpiContext, SimulVolume &simulVolume);