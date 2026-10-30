#include <fstream>
#include <unordered_set>

#include "debug/debug.hpp"
#include "mpi/mpi_function.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

// The first-order reactions at the top of a step -- dissociation, state
// change, destruction, unimolecular creation -- run on every rank at once,
// each on the complexes it owns.  The right neighbour's first half then binds
// into this rank's last owned column (bin n-2, its bin 0) through the ghost
// copy it was sent at the end of the previous step, which knows nothing of
// those reactions.  When both touch one complex -- this rank dissociates a
// dimer, the neighbour binds one of its monomers as if the dimer were whole --
// the first exchange can keep only one version (deserialize_molecules_right()
// keeps the neighbour's), so one reaction is lost and the bond the other made
// is left pointing at a partner that does not bond back.  A serial run has
// one order for the two: the first-order reaction, then the binding, which
// sees its result.
//
// So the complexes the neighbour can reach wait.  At the top of the step,
// every complex this rank owns with a member in bin n-2 or beyond has its
// members held back from the first-order pass, and the pass runs for them
// alone after the first exchange, when this rank holds whatever the
// neighbour did to them.  Nothing else touches them in between: this rank's
// own first half stops two columns short of bin n-2 (MpiContext::
// column_share(), region_mid_index()).  A complex the neighbour did change is
// left out of the late pass for this step, as a reactant of an association is
// in a serial step, where it cannot also dissociate or change state after it
// has bound; that costs the complex's older bonds one draw.
//
// The rank with no right neighbour has nothing to wait for.

// Members of the complexes to hold back, by index; marks them
// Molecule::firstOrderHeld, which every first-order loop skips.  Not
// isGhosted: the implicit lipid's unbinding counts this rank's free proteins
// by that flag, and held-back molecules are still this rank's.
void hold_back_right_edge_first_order(MpiContext &mpiContext,
                                      vector<Molecule> &moleculeList,
                                      vector<Complex> &complexList,
                                      vector<int> &heldBack) {
  heldBack.clear();
  if (mpiContext.rank >= mpiContext.nprocs - 1) return;
  const int edge = mpiContext.simulVolume->numSubCells.x - 2;
  for (auto &com : complexList) {
    if (com.isEmpty || com.ownerRank != mpiContext.rank) continue;
    bool reachable = false;
    for (int m : com.memberList) {
      if (m < 0 || m >= static_cast<int>(moleculeList.size())) continue;
      const Molecule &mol = moleculeList[m];
      if (mol.isEmpty || mol.isImplicitLipid) continue;
      if (get_x_bin(mpiContext, moleculeList[m]) >= edge) {
        reachable = true;
        break;
      }
    }
    if (!reachable) continue;
    for (int m : com.memberList) {
      if (m < 0 || m >= static_cast<int>(moleculeList.size())) continue;
      Molecule &mol = moleculeList[m];
      if (mol.isEmpty || mol.isImplicitLipid || mol.isGhosted) continue;
      mol.firstOrderHeld = true;
      heldBack.push_back(m);
    }
  }
}

void release_held_back(vector<Molecule> &moleculeList,
                       const vector<int> &heldBack) {
  for (int m : heldBack)
    if (m >= 0 && m < static_cast<int>(moleculeList.size()))
      moleculeList[m].firstOrderHeld = false;
}

// The late pass.  heldBackIds are the molecule ids recorded at the top of the
// step: the exchange in between may have appended molecules but leaves the
// ones that stayed at their index, and ids are what survive it regardless.
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
    std::ofstream &assocDissocFile) {
  if (heldBackIds.empty()) return;
  std::unordered_set<int> ids(heldBackIds.begin(), heldBackIds.end());

  // Everything outside the held-back set sits this pass out, held for the
  // length of the pass.  The pass can append molecules (unimolecular
  // creation) but never moves one, so the marks still line up by index.
  const size_t nBefore = moleculeList.size();
  vector<char> late(nBefore, 0);
  for (size_t i = 0; i < nBefore; i++) {
    Molecule &mol = moleculeList[i];
    if (mol.isEmpty) continue;
    bool isLate = !mol.isImplicitLipid && !mol.isGhosted && ids.count(mol.id);
    // The neighbour bound this complex in its first half: like a reactant of
    // an association in a serial step, it takes no first-order reaction after
    // that.  Any member carrying the mark rules out the whole complex.
    if (isLate && mol.myComIndex >= 0 &&
        mol.myComIndex < static_cast<int>(complexList.size())) {
      const Complex &com = complexList[mol.myComIndex];
      if (com.trajStatus == TrajStatus::propagated) isLate = false;
      for (int m : com.memberList)
        if (m >= 0 && m < static_cast<int>(moleculeList.size()) &&
            moleculeList[m].trajStatus == TrajStatus::propagated)
          isLate = false;
    }
    late[i] = isLate;
    if (!isLate) mol.firstOrderHeld = true;
  }

  check_for_unimolecular_reactions(
      simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
      backRxns, createDestructRxns, molTemplateList, observablesList,
      counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
      ILTableIDs, assocDissocFile);
  if (params.hasUniMolStateChange == true)
    check_for_unimolstatechange_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject);
  // The same loop as the top-of-step pass in
  // check_perform_zeroth_first_order_reactions().
  if (params.implicitLipid == true) {
    for (unsigned molItr{0}; molItr < moleculeList.size(); ++molItr) {
      if (moleculeList[molItr].isEmpty ||
          moleculeList[molItr].isImplicitLipid == true ||
          complexList[moleculeList[molItr].myComIndex].OnSurface == false ||
          moleculeList[molItr].isGhosted == true ||
          moleculeList[molItr].firstOrderHeld)
        continue;
      check_dissociation_implicitlipid(
          simItr, params, simulVolume, molTemplateList, observablesList, molItr,
          moleculeList, complexList, backRxns, forwardRxns, createDestructRxns,
          counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
          ILTableIDs, assocDissocFile);
    }
  }

  for (size_t i = 0; i < nBefore; i++) moleculeList[i].firstOrderHeld = false;

  // This step's candidate pairs were listed before the late pass, from the
  // state it has just changed.  A serial step draws first-order reactions
  // before listing, and lists a molecule that dissociated or changed state
  // with probability zero, for overlap avoidance only
  // (determine_3D_bimolecular_reaction_probability() and the 2D/1D arms skip
  // isDissociated); do the same to what the late pass reacted.  A molecule it
  // destroyed has to leave its partners' lists altogether, as a transmitted
  // one does, because the overlap sweeps look up every listed partner's
  // complex and a destroyed molecule's myComIndex is -1.
  for (size_t i = 0; i < nBefore; i++) {
    if (!late[i]) continue;
    Molecule &mol = moleculeList[i];
    if (mol.isEmpty) {
      remove_partner_crossings(mol, moleculeList, complexList);
      mol.crossings.clear();
      continue;
    }
    if (mol.trajStatus != TrajStatus::propagated) continue;
    zero_partner_probvec(mol, moleculeList);
    for (auto &cross : mol.crossings) cross.prob = 0;
  }
}

// The overlap sweep of a half-step acts on the molecules of its region, and
// the plain sweep moves the initiator's complex alone.  The cluster sweeps
// (clusterOverlapCheck) move every complex of the cluster, and a cluster
// reaches across the region boundary: in the first half it moved and marked
// owned molecules of the right region, whose columns the right neighbour
// binds into during the same half through the copy it holds.  The first
// exchange then kept this rank's moved copies (deserialize_molecules_right()
// keeps what this rank propagated) and the neighbour's bonds into them were
// lost -- on mem_localization with the cluster sweep, a third of the
// membrane binding.  Hold everything outside the region for the length of
// the sweep; the cluster sweeps leave a held molecule where it is, as they
// leave a ghost.
void hold_outside_region(vector<Molecule> &moleculeList,
                         const vector<int> &region) {
  for (auto &mol : moleculeList) mol.sweepHeld = !mol.isEmpty;
  for (int m : region)
    if (m >= 0 && m < static_cast<int>(moleculeList.size()))
      moleculeList[m].sweepHeld = false;
}

void release_sweep_hold(vector<Molecule> &moleculeList) {
  for (auto &mol : moleculeList) mol.sweepHeld = false;
}
