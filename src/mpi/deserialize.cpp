#include <algorithm>
#include <cstdio>
#include <iostream>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

// A received complex replaces the local memberList wholesale, so a molecule that
// was a member here and is not in the list that arrived keeps myComIndex
// pointing at the complex while the complex no longer lists it back.  Usually
// the molecule left that complex for another one that arrives in the same
// message and re-parents it; when it does not, the molecule is left claiming a
// complex that does not claim it, and delete_disappeared_complexes_partial() can
// then destroy the complex from under it.  Nothing downstream survives that:
// create_complex_propagation_vectors() walks the cleared memberList and
// check_bimolecular_reactions() and the ncross reset index complexList with the
// -1 it ends up holding.
//
// Put the molecule back in the list it still claims.  The complex is then stale
// rather than inconsistent, and the next exchange replaces its memberList again.
static void reconcile_complex_membership(vector<Molecule> &moleculeList,
                                         vector<Complex> &complexList) {
  for (auto &mol : moleculeList) {
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    if (mol.myComIndex < 0 ||
        mol.myComIndex >= static_cast<int>(complexList.size()))
      continue;
    Complex &com = complexList[mol.myComIndex];
    if (com.isEmpty) continue;
    if (std::find(com.memberList.begin(), com.memberList.end(), mol.index) ==
        com.memberList.end())
      com.memberList.push_back(mol.index);
  }
}

// Both ranks changed the same complex in the same half-step: this rank by
// the pass or the reactions it ran since it sent its copy, the neighbour by
// what it did to the copy it was sent.  The exchange keeps one of the two
// versions whole, so the other rank's reaction is lost -- and the partner it
// bound on that rank keeps a bond to a molecule that no longer has it.  A
// report, not a repair: the line is what run_mpi_suite.sh counts.
//
// keptLocalCom says which complexes this rank changed; neighbourChanged which
// ones the neighbour's copy says it changed.  On a rank that owns one column
// the first test includes a plain move, which the neighbour's bond replaces
// just as it would a reaction.
static void report_exchange_conflicts(MpiContext &mpiContext, const char *side,
                                      const vector<Molecule> &moleculeList,
                                      const vector<Complex> &complexList,
                                      const vector<char> &keptLocalCom,
                                      const vector<char> &neighbourChanged) {
  for (size_t ci = 0; ci < complexList.size() && ci < neighbourChanged.size();
       ci++) {
    if (!neighbourChanged[ci] || ci >= keptLocalCom.size() || !keptLocalCom[ci])
      continue;
    const Complex &com = complexList[ci];
    if (com.isEmpty) continue;
    fprintf(stderr,
            "EXCHANGE CONFLICT itr=%d rank=%d from=%s complex id=%d owner=%d "
            "members=%zu first member id=%d\n",
            mpiContext.simItr, mpiContext.rank, side, com.id, com.ownerRank,
            com.memberList.size(),
            com.memberList.empty() ? -1
                                   : moleculeList[com.memberList[0]].id);
  }
}

void snapshot_molecules(const vector<Molecule> &moleculeList,
                    const vector<Complex> &complexList,
                    vector<MolSnapshot> &snapshot) {
  snapshot.assign(moleculeList.size(), MolSnapshot{});
  for (size_t i = 0; i < moleculeList.size(); i++) {
    const Molecule &mol = moleculeList[i];
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    MolSnapshot &s = snapshot[i];
    if (mol.myComIndex >= 0 &&
        mol.myComIndex < static_cast<int>(complexList.size()))
      s.comId = complexList[mol.myComIndex].id;
    for (const auto &iface : mol.interfaceList)
      if (iface.interaction.partnerIndex != -1) s.nBound++;
  }
}

// The second exchange of a step carries the left neighbour's copy of this
// rank's left edge, and that copy is stale by whatever happened here after
// the copy was sent: the second bimolecular pass, whose initiator in the
// right region can bind a partner in the left region, rewiring the partner's
// interface and moving it into the initiator's complex; and, on a rank that
// owns one column, the receive from the right, which lands in the same column
// the left neighbour mirrors.  The neighbour knows nothing of either.
// trajStatus cannot tell the two apart -- the left region keeps "propagated"
// from the first half, so that a complex spanning both regions is not moved
// twice, and every flag is reset before the send -- so the test is against a
// snapshot taken as the copy was sent: a molecule whose complex id or
// bound-interface count differs was changed since, and so was one that did
// not exist then.  As at the first exchange, the decision is per complex: it
// is one object on the wire and one object here.
static bool changed_since(const Molecule &mol,
                          const vector<Complex> &complexList,
                          const vector<MolSnapshot> &before) {
  if (mol.index < 0 || mol.index >= static_cast<int>(before.size()))
    return true;
  const MolSnapshot &s = before[mol.index];
  int comId = -1;
  if (mol.myComIndex >= 0 &&
      mol.myComIndex < static_cast<int>(complexList.size()))
    comId = complexList[mol.myComIndex].id;
  if (comId != s.comId) return true;
  int nBound = 0;
  for (const auto &iface : mol.interfaceList)
    if (iface.interaction.partnerIndex != -1) nBound++;
  return nBound != s.nBound;
}

static bool changed_in_second_half(const Complex &com,
                                   const vector<Molecule> &moleculeList,
                                   const vector<Complex> &complexList,
                                   const vector<MolSnapshot> &before) {
  for (int m : com.memberList) {
    if (m < 0 || m >= static_cast<int>(moleculeList.size())) continue;
    const Molecule &mol = moleculeList[m];
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    if (changed_since(mol, complexList, before)) return true;
  }
  return false;
}

void deserialize_complexes(MpiContext &mpiContext,
                           vector<Molecule> &moleculeList,
                           vector<Complex> &complexList,
                           unsigned char *arrayRank, int &nArrayRank,
                           const vector<char> &keptLocalCom,
                           const vector<char> &keptLocalMol) {
  if (VERBOSE) cout << "deserialize_complexes begins" << endl;
  // Store number of complexes to be deserialized:
  size_t nComplexes;
  POP(nComplexes);
  if (VERBOSE) cout << "nComplexes = " << nComplexes << endl;
  // Deserialize received complexes:
  for (size_t i = 0; i < nComplexes; i++) {
    if (VERBOSE) cout << "i = " << i << endl;
    // Deserialize i-th received complex:
    Complex c;
    c.deserialize(arrayRank, nArrayRank);

    c.receivedFromNeighborRank = true;
    // cout << "deserialized c.id = " << c.id << endl;

    int complexIndex = find_complex(complexList, c.id);

    // Remap received IDs into moleculeList indices on this rank:
    vector<int> memberList{};
    if (VERBOSE) cout << "members: ";
    for (auto &it : c.memberList) {  // looping over IDs
      if (VERBOSE) cout << it;
      int molIndex = find_molecule(moleculeList, it);
      // A member this rank cannot resolve must not enter memberList: -1 is
      // indexed unguarded by Complex::propagate(), Complex::update_properties()
      // and the serializer, so it corrupts the heap rather than failing.  Drop
      // the member and say so; the sender is supposed to ship a complex whole,
      // so this means that invariant broke.
      //
      // An emptied slot is just as unusable and is reachable by find_molecule(),
      // which matches on id alone while Molecule::MPI_remove_from_one_rank()
      // and Complex::destroy() leave id set on the molecule they empty.  Such a
      // molecule carries molTypeIndex == -1, which update_properties() feeds
      // straight to numEachMol[].  All the incoming molecules are already
      // installed by the time complexes are deserialized, so a slot still empty
      // here is dead rather than about to be filled.
      if (molIndex == -1 || moleculeList[molIndex].isEmpty) {
        fprintf(stderr,
                "rank %d: complex id=%d arrived from the left listing member "
                "id=%d, which is %s on this rank; dropping the member\n",
                mpiContext.rank, c.id, it,
                molIndex == -1 ? "not" : "deleted");
        if (DEBUG) error(mpiContext, "5: complex member mol not found");
        continue;
      }
      // A member whose local copy was kept stays in the complex it is in
      // here; re-parenting it to the received complex would have it listed
      // by two complexes and naming the wrong one.
      if (molIndex < static_cast<int>(keptLocalMol.size()) &&
          keptLocalMol[molIndex])
        continue;
      memberList.push_back(molIndex);
      if (complexIndex == -1) {
        if (moleculeList[molIndex].justBoundThisStep) {
          complexIndex = moleculeList[molIndex].myComIndex;
        }
      }
      if (VERBOSE) cout << "(" << molIndex << "); ";
    }
    if (VERBOSE) cout << endl;
    // An empty memberList is worse than a stale one: delete_disappeared_
    // complexes_partial() reads memberList[0] of every complex that is not
    // flagged empty.  Leave this rank's own copy untouched instead.
    if (memberList.empty()) continue;
    c.memberList = memberList;

    // Extract unique complex identifier in the system,
    // and find index to local complex with the same ID:

    if (VERBOSE)
      cout << "complexIndex = " << complexIndex << "complexId = " << c.id
           << endl;

    if (complexIndex == -1) {  // new complex on this rank
      if (VERBOSE) printf("This is a new complex (id=%d) here\n", c.id);

      c.index = complexList.size();
      complexIndex = c.index;
      complexList.push_back(c);
    } else {
      if (VERBOSE)
        printf("This is a old complex (id=%d) here\n",
               complexList[complexIndex].id);
      // The decision deserialize_molecules() took for the members: every
      // member of a complex changed in the second half kept its local copy,
      // so the complex keeps its local copy too, or the members would name a
      // complex that no longer lists them.
      if (complexIndex < static_cast<int>(keptLocalCom.size()) &&
          keptLocalCom[complexIndex]) {
        if (VERBOSE)
          printf("Skip complex (id=%d) because this rank changed it.\n",
                 complexList[complexIndex].id);
        complexList[complexIndex].receivedFromNeighborRank = true;
        continue;
      }

      auto vecOld = complexList[complexIndex].memberList;
      if (VERBOSE)
        cout << "Old receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      // Only the rank that currently owns a complex may say where ownership
      // goes.  ownerRank travels inside the complex, so a rank holding a ghost
      // copy could reassign it just by sending that copy back -- and ranks do
      // send ghost copies back, because that is how a complex is handed over.
      // A middle rank would hand a complex to one neighbour and have the other
      // neighbour's stale copy hand it straight back, so both ended up owning
      // it, propagating it independently from then on.
      const int senderRank = mpiContext.rank - 1;
      const int currentOwner = complexList[complexIndex].ownerRank;
      const bool senderMayReassign =
          (currentOwner == senderRank) || (currentOwner == -1);
      complexList[complexIndex] = c;  // copy all deserialized properties
      if (!senderMayReassign)
        complexList[complexIndex].ownerRank = currentOwner;
      if (VERBOSE)
        cout << "New receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      complexList[complexIndex].index = complexIndex;
      auto &vecNew = complexList[complexIndex].memberList;

      // Update member molecules to match also old molecules:
      // for (auto &it : vecOld) {
      //   // cout << "add old..." << endl;
      //   //  If old memberList element isn't found in vecNew, add it:
      //   if (std::find(vecNew.begin(), vecNew.end(), it) == vecNew.end()) {
      //     // Push only if a member molecule isn't deleted:
      //     if ((it < moleculeList.size()) &&
      //         (moleculeList[it].myComIndex != -1)) {
      //       // Push "it" only if not in a shared zone:
      //       int xBin = get_x_bin(mpiContext, moleculeList[it]);
      //       if (((!mpiContext.rank) || (xBin > 1)) &&
      //           ((mpiContext.rank == mpiContext.nprocs - 1) ||
      //            (xBin < (*(mpiContext.simulVolume)).numSubCells.x - 2))) {
      //         vecNew.push_back(it);
      //       }
      //     }
      //   }
      // }
    }

    // Update myComIndex of member molecules to complexIndex:
    complexList[complexIndex].index = complexIndex;
    for (auto &it : complexList[complexIndex].memberList) {
      // cout << "updating member " << it << endl;
      if (DEBUG && ((it >= moleculeList.size()) || (it == -1)))
        error(mpiContext,
              "6: complex member mol not found: it=" + to_string(it));
      Molecule &mol = moleculeList[it];
      mol.myComIndex = complexIndex;
    }
  }
  reconcile_complex_membership(moleculeList, complexList);
  //    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
  //    "//end deserialize_complexes()");
  if (VERBOSE) cout << "deserialize_complexes ends" << endl;
}

// The first exchange of a step carries the right neighbour's copy of this
// rank's edge molecules, and that copy is stale by whatever this rank did to
// them at the top of the step: the unimolecular state changes and
// dissociations run on every rank at once, before anything is sent.  Those
// routines mark what they touched TrajStatus::propagated, so a complex this
// rank changed is one whose trajStatus, or any member's, says so.  The
// decision is per complex, because a complex is one object on the wire and
// one object here: taking some members and not others leaves the two halves
// naming different complexes.
static bool changed_here_this_step(const Complex &com,
                                   const vector<Molecule> &moleculeList) {
  if (com.trajStatus == TrajStatus::propagated) return true;
  for (int m : com.memberList) {
    if (m < 0 || m >= static_cast<int>(moleculeList.size())) continue;
    if (moleculeList[m].trajStatus == TrajStatus::propagated) return true;
  }
  return false;
}

void deserialize_complexes_right(MpiContext &mpiContext,
                                 vector<Molecule> &moleculeList,
                                 vector<Complex> &complexList,
                                 unsigned char *arrayRank, int &nArrayRank,
                                 const vector<char> &keptLocalCom,
                                 const vector<char> &keptLocalMol) {
  if (VERBOSE) cout << "deserialize_complexes begins..." << endl;
  // Store number of complexes to be deserialized:
  size_t nComplexes;
  POP(nComplexes);
  if (VERBOSE) cout << "nComplexes = " << nComplexes << endl;
  // Deserialize received complexes:
  for (size_t i = 0; i < nComplexes; i++) {
    if (VERBOSE) cout << "i = " << i << endl;
    // Deserialize i-th received complex:
    Complex c;
    c.deserialize(arrayRank, nArrayRank);

    c.receivedFromNeighborRank = true;
    int complexIndex = find_complex(complexList, c.id);

    // Remap received IDs into moleculeList indices on this rank:
    vector<int> memberList{};
    if (VERBOSE) cout << "members: ";
    for (auto &it : c.memberList) {  // looping over IDs
      if (VERBOSE) cout << it;
      int molIndex = find_molecule(moleculeList, it);
      // Same as the left-hand side.
      if (molIndex == -1 || moleculeList[molIndex].isEmpty) {
        fprintf(stderr,
                "rank %d: complex id=%d arrived from the right listing member "
                "id=%d, which is %s on this rank; dropping the member\n",
                mpiContext.rank, c.id, it,
                molIndex == -1 ? "not" : "deleted");
        if (DEBUG) error(mpiContext, "5: complex member mol not found");
        continue;
      }
      // A member whose local copy was kept stays in the complex it is in
      // here.  Re-parenting it to the received complex would have it listed
      // by two complexes and naming the wrong one, which is exactly the state
      // that the per-complex decision exists to prevent.
      if (molIndex < static_cast<int>(keptLocalMol.size()) &&
          keptLocalMol[molIndex])
        continue;
      memberList.push_back(molIndex);
      if (VERBOSE) cout << "(" << molIndex << "); ";
    }
    if (VERBOSE) cout << endl;
    if (memberList.empty()) continue;
    c.memberList = memberList;

    if (VERBOSE)
      cout << "complexIndex = " << complexIndex << "complexId = " << c.id
           << endl;

    if (complexIndex == -1) {  // new complex on this rank
      if (VERBOSE) printf("This is a new complex (id=%d) here\n", c.id);

      c.index = complexList.size();
      complexIndex = c.index;
      complexList.push_back(c);
    } else {
      if (VERBOSE)
        printf("This is a old complex (id=%d) here\n",
               complexList[complexIndex].id);
      // The decision deserialize_molecules_right() took for the members, read
      // from its snapshot rather than retaken here: by now the members that
      // were replaced carry the sender's trajStatus, which says nothing about
      // what this rank did.  Every member of a complex changed here kept its
      // local copy, so the complex keeps its local copy too, or the members
      // would name a complex that no longer lists them.
      if (complexIndex < static_cast<int>(keptLocalCom.size()) &&
          keptLocalCom[complexIndex]) {
        if (VERBOSE)
          printf("Skip complex (id=%d) because it is already propagated.\n",
                 complexList[complexIndex].id);
        complexList[complexIndex].receivedFromNeighborRank = true;
        continue;
      }

      auto vecOld = complexList[complexIndex].memberList;
      if (VERBOSE)
        cout << "Old receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      // Only the rank that currently owns a complex may say where ownership
      // goes.  ownerRank travels inside the complex, so a rank holding a ghost
      // copy could reassign it just by sending that copy back -- and ranks do
      // send ghost copies back, because that is how a complex is handed over.
      // A middle rank would hand a complex to one neighbour and have the other
      // neighbour's stale copy hand it straight back, so both ended up owning
      // it, propagating it independently from then on.
      const int senderRank = mpiContext.rank + 1;
      const int currentOwner = complexList[complexIndex].ownerRank;
      const bool senderMayReassign =
          (currentOwner == senderRank) || (currentOwner == -1);
      // ncross counts the candidate reactions of the members, which stay
      // this rank's (deserialize_molecules_right()); the sender's count is
      // of its own lists.
      const int localNcross = complexList[complexIndex].ncross;
      complexList[complexIndex] = c;  // copy all deserialized properties
      complexList[complexIndex].ncross = localNcross;
      if (!senderMayReassign)
        complexList[complexIndex].ownerRank = currentOwner;
      if (VERBOSE)
        cout << "New receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      complexList[complexIndex].index = complexIndex;
      auto &vecNew = complexList[complexIndex].memberList;

      // Update member molecules to match also old molecules:
      // for (auto &it : vecOld) {
      //   // cout << "add old..." << endl;
      //   //  If old memberList element isn't found in vecNew, add it:
      //   if (std::find(vecNew.begin(), vecNew.end(), it) == vecNew.end()) {
      //     // Push only if a member molecule isn't deleted:
      //     if ((it < moleculeList.size()) &&
      //         (moleculeList[it].myComIndex != -1)) {
      //       // Push "it" only if not in a shared zone:
      //       int xBin = get_x_bin(mpiContext, moleculeList[it]);
      //       if (((!mpiContext.rank) || (xBin > 1)) &&
      //           ((mpiContext.rank == mpiContext.nprocs - 1) ||
      //            (xBin < (*(mpiContext.simulVolume)).numSubCells.x - 2))) {
      //         vecNew.push_back(it);
      //       }
      //     }
      //   }
      // }
    }

    // Update myComIndex of member molecules to complexIndex:
    complexList[complexIndex].index = complexIndex;
    for (auto &it : complexList[complexIndex].memberList) {
      // cout << "updating member " << it << endl;
      if (DEBUG && ((it >= moleculeList.size()) || (it == -1)))
        error(mpiContext,
              "6: complex member mol not found: it=" + to_string(it));
      Molecule &mol = moleculeList[it];
      mol.myComIndex = complexIndex;
    }
  }
  reconcile_complex_membership(moleculeList, complexList);
  //    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
  //    "//end deserialize_complexes()");
  if (VERBOSE) cout << "deserialize_complexes ends" << endl;
}

// Deserializes stripe of molecules with cellX x-bin coordinate,
// and updates moleculeList.
// Example stripes x-bin coordinates:
// rank0 |rank1|rank2| rank3  <- ownership
// 0 1 2 | 3                  <- rank0 visibility (includes ghosted 3)
//     2 | 3 4 | 5            <- rank1 visibility (includes ghosted 2 and 5)
//           4 | 5 6 | 7      <- rank2 visibility (includes ghosted 4 and 7)
//                 6 | 7 8 9  <- rank3 visibility (includes ghosted 6)
void deserialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume,
                           vector<Molecule> &moleculeList,
                           vector<Complex> &complexList,
                           vector<MolTemplate> &molTemplateList,
                           Membrane &membraneObject,
                           copyCounters &counterArrays, vector<int> &indices,
                           vector<int> &indicesEnteringMolecules,
                           vector<int> &indicesExitingMolecules,
                           unsigned char *arrayRank, int &nArrayRank,
                           const vector<MolSnapshot> &sentLeft,
                           vector<char> &keptLocalCom,
                           vector<char> &keptLocalMol) {
  // Decide which complexes keep their local copy before any molecule is
  // overwritten: a received molecule carries the sender's bonds and complex
  // id, so once one member has been replaced the comparison against the
  // snapshot would read the neighbour's state as this rank's.
  keptLocalCom.assign(complexList.size(), 0);
  for (size_t ci = 0; ci < complexList.size(); ci++)
    if (!complexList[ci].isEmpty &&
        changed_in_second_half(complexList[ci], moleculeList, complexList,
                               sentLeft))
      keptLocalCom[ci] = 1;
  // A complex the right neighbour owns is not the left neighbour's to
  // replace: its copy of it, if it has one at all, is a mirror of a mirror.
  // Those meet only on a rank that owns one column, where the left
  // neighbour's echo of that column would put the right neighbour's complex
  // back to where it was at the first exchange.
  if (mpiContext.rank < mpiContext.nprocs - 1)
    for (size_t ci = 0; ci < complexList.size(); ci++)
      if (!complexList[ci].isEmpty &&
          complexList[ci].ownerRank == mpiContext.rank + 1)
        keptLocalCom[ci] = 1;
  keptLocalMol.assign(moleculeList.size(), 0);

  int nMols;   // number of molecules to deserialize
  POP(nMols);  // deserialize from arrayRank
  if (VERBOSE)
    printf("deserializing (%d) molecules from the left rank...\n", nMols);
  if (CHECK_STRUCTURE) {
    // The molecules are read one at a time below, so peek at the message
    // here for the conflict report and leave the read position alone.
    int peek = nArrayRank;
    vector<char> neighbourChanged(complexList.size(), 0);
    for (int i = 0; i < nMols; i++) {
      Molecule mol;
      mol.deserialize(arrayRank, peek);
      const int molIndex = find_molecule(moleculeList, mol.id);
      if (molIndex == -1 || molIndex >= static_cast<int>(sentLeft.size()))
        continue;
      int nBound = 0;
      for (const auto &iface : mol.interfaceList)
        if (iface.interaction.partnerIndex != -1) nBound++;
      // On the wire myComIndex carries the complex id.
      if (mol.myComIndex == sentLeft[molIndex].comId &&
          nBound == sentLeft[molIndex].nBound)
        continue;
      const int ci = moleculeList[molIndex].myComIndex;
      if (ci >= 0 && ci < static_cast<int>(complexList.size()))
        neighbourChanged[ci] = 1;
    }
    // keptLocalCom also holds the right neighbour's complexes by now, which
    // this rank did not change; ask the second-half test alone.
    vector<char> changedHere(complexList.size(), 0);
    for (size_t ci = 0; ci < complexList.size(); ci++)
      changedHere[ci] =
          !complexList[ci].isEmpty &&
          changed_in_second_half(complexList[ci], moleculeList, complexList,
                                 sentLeft);
    report_exchange_conflicts(mpiContext, "left", moleculeList, complexList,
                              changedHere, neighbourChanged);
  }
  // Deserialize received molecules and update structures:
  for (int iDeserialize = 0; iDeserialize < nMols; iDeserialize++) {
    // Deserialize each received molecule:
    Molecule mol;

    mol.deserialize(arrayRank, nArrayRank);

    if (VERBOSE) printf("deserializing mol (id=%d) ...\n", mol.id);
    if (CHECK_STRUCTURE) trace_molecule(mpiContext, "wire copy", moleculeList, complexList, mol);

    update_mySubVolIndex(mpiContext, mol, membraneObject, simulVolume);

    mol.receivedFromNeighborRank = true;

    // Ghost status is NOT decided here.  It used to be recomputed from
    // get_x_bin against this rank's own xOffset, while the sender had already
    // decided by column index -- two criteria on two ranks that did not
    // partition the boundary, so a molecule could end up owned by neither.
    // derive_ghost_from_ownership() sets it once, after the complexes carrying
    // the ownership stamp have arrived.

    int molIndex = find_molecule(moleculeList, mol.id);

    // Check whether the molecule was already in ghosted zone,
    // or it has just arrived (index == -1):
    if (molIndex != -1) {  // if found
      if (VERBOSE) printf("mol (id=%d) is already here.\n", mol.id);

      if (DEBUG && (moleculeList[molIndex].id != mol.id))
        error(mpiContext, "7 // deserialize_molecules()");
      mol.index = molIndex;

      // This rank's copy is the newer one when the second bimolecular pass
      // changed the molecule's complex here: the neighbour's copy dates from
      // the first exchange and would undo the bond on this side only, leaving
      // the initiator bound to a partner that is not bound back.
      const int localCom = moleculeList[molIndex].myComIndex;
      const bool keepLocal =
          changed_since(moleculeList[molIndex], complexList,
                        sentLeft) ||
          (localCom >= 0 && localCom < static_cast<int>(keptLocalCom.size()) &&
           keptLocalCom[localCom]);
      if (keepLocal) {
        if (VERBOSE)
          printf("mol (id=%d) was changed here this half-step, so skip the "
                 "deserialization.\n",
                 mol.id);
        moleculeList[molIndex].receivedFromNeighborRank = true;
        keptLocalMol[molIndex] = 1;
        if (CHECK_STRUCTURE) trace_molecule(mpiContext, "kept local", moleculeList, complexList, moleculeList[molIndex]);
        continue;
      }

      // Update mySubVolIndex if needed
      // (received molecule might have moved by another rank
      // since this rank serialized it)
      auto &oldCellIndex = moleculeList[molIndex].mySubVolIndex;

      // check if the mol cross the ghosted boundary
      // int xBinOld = get_x_bin(mpiContext, moleculeList[molIndex]);
      // int xBinNew = get_x_bin(mpiContext, mol);
      // if (check_enter_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                              indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) enter the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesEnteringMolecules.push_back(molIndex);
      // }
      // if (check_exit_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                             indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) exit the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesExitingMolecules.push_back(molIndex);
      // }

      // The reweighting history is not on the wire (Molecule::serialize()),
      // and it is this rank's: the record of the pairs it evaluated, filed
      // under its own indices.  Losing it at every exchange gave every pair
      // with an edge molecule as its lower-index member a fresh history each
      // step, i.e. no reweighting at all.
      mol.prevReweight = std::move(moleculeList[molIndex].prevReweight);
      mol.currReweight = std::move(moleculeList[molIndex].currReweight);
      moleculeList[molIndex] = mol;
      if (CHECK_STRUCTURE) trace_molecule(mpiContext, "overwritten", moleculeList, complexList, moleculeList[molIndex]);

      auto &newCellIndex = moleculeList[molIndex].mySubVolIndex;
      // Check if cell index has changed:
      if (newCellIndex != oldCellIndex) {
        // Update simulVolume cells memberMolList of both mySubVolIndex:
        if (oldCellIndex < simulVolume.subCellList.size()) {
          auto &vec = simulVolume.subCellList[oldCellIndex].memberMolList;
          vec.erase(remove(vec.begin(), vec.end(), molIndex), vec.end());
        }
        // Deserialized molecule will be at the same index as existing one:
        if (newCellIndex < simulVolume.subCellList.size()) {
          simulVolume.subCellList[newCellIndex].memberMolList.push_back(
              molIndex);
        }
      }
    } else {  // if not found
      if (VERBOSE) printf("mol (id=%d) is new here.\n", mol.id);
      //  The molecule doesn't exist in moleculeList. Add it:
      mol.index = moleculeList.size();
      molIndex = mol.index;
      moleculeList.push_back(mol);
      if (CHECK_STRUCTURE) trace_molecule(mpiContext, "appended", moleculeList, complexList, moleculeList.back());
    }
    indices.push_back(molIndex);
    // if(DEBUG) DEBUG_MOL("1DESERIALIZed A MOLECULE");
  }

  if (VERBOSE)
    cout << "Total deserialized number of molecules: " << nMols << endl;
}

void deserialize_molecules_right(
    MpiContext &mpiContext, SimulVolume &simulVolume,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject,
    copyCounters &counterArrays, vector<int> &indices,
    vector<int> &indicesEnteringMolecules, vector<int> &indicesExitingMolecules,
    unsigned char *arrayRank, int &nArrayRank,
    const unordered_map<int, MolSnapshot> &sentRight,
    vector<char> &keptLocalCom, vector<char> &keptLocalMol) {
  int nMols;   // number of molecules to deserialize
  POP(nMols);  // deserialize from arrayRank
  if (VERBOSE)
    printf("deserializing (%d) molecules from the right rank...\n", nMols);
  // The whole message is read before anything is decided: the decision per
  // complex below reads the received copies as well as the local ones.
  vector<Molecule> received(nMols);
  for (int i = 0; i < nMols; i++) received[i].deserialize(arrayRank, nArrayRank);

  // Decide which complexes keep their local copy before any molecule is
  // overwritten.  A received molecule carries the sender's trajStatus, so once
  // one member of a complex has been replaced the test would read the
  // neighbour's flag as this rank's and keep the rest of the complex stale.
  keptLocalCom.assign(complexList.size(), 0);
  for (size_t ci = 0; ci < complexList.size(); ci++)
    if (!complexList[ci].isEmpty &&
        changed_here_this_step(complexList[ci], moleculeList))
      keptLocalCom[ci] = 1;
  // ...unless the neighbour changed the complex too.  Its copy of a molecule
  // that no longer matches what this rank sent it is a bond the neighbour's
  // first half made into this rank's edge, or a change to a complex the
  // neighbour owns; either way the neighbour has acted on it and its copy is
  // the one both ranks have to end up with.  On a rank that owns two columns
  // or more this rank's edge is not moved before this exchange, so the two
  // tests never both hold; on a rank that owns one column its whole column
  // was moved in the first half, and keeping it would drop every bond the
  // right neighbour made into it, leaving the neighbour's initiator bound to
  // a partner that is not bound back.  A molecule this rank also bound in
  // the same half-step is lost either way; that is the one case neither
  // choice can make consistent.
  vector<char> neighbourChanged(complexList.size(), 0);
  for (const Molecule &mol : received) {
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    auto sent = sentRight.find(mol.id);
    if (sent == sentRight.end()) continue;
    int nBound = 0;
    for (const auto &iface : mol.interfaceList)
      if (iface.interaction.partnerIndex != -1) nBound++;
    // On the wire myComIndex carries the complex id.
    if (mol.myComIndex == sent->second.comId && nBound == sent->second.nBound)
      continue;
    const int molIndex = find_molecule(moleculeList, mol.id);
    if (molIndex == -1) continue;
    const int ci = moleculeList[molIndex].myComIndex;
    if (ci >= 0 && ci < static_cast<int>(complexList.size()))
      neighbourChanged[ci] = 1;
  }
  if (CHECK_STRUCTURE)
    report_exchange_conflicts(mpiContext, "right", moleculeList, complexList,
                              keptLocalCom, neighbourChanged);
  for (size_t ci = 0; ci < complexList.size(); ci++)
    if (neighbourChanged[ci]) keptLocalCom[ci] = 0;
  // A complex the left neighbour owns is not the right neighbour's to
  // replace; see deserialize_molecules() for the mirror of this.
  if (mpiContext.rank > 0)
    for (size_t ci = 0; ci < complexList.size(); ci++)
      if (!complexList[ci].isEmpty &&
          complexList[ci].ownerRank == mpiContext.rank - 1)
        keptLocalCom[ci] = 1;
  keptLocalMol.assign(moleculeList.size(), 0);

  // Apply received molecules and update structures:
  for (Molecule &mol : received) {
    if (VERBOSE) printf("deserializing mol (id=%d) ...\n", mol.id);
    if (CHECK_STRUCTURE) trace_molecule(mpiContext, "wire copy", moleculeList, complexList, mol);

    update_mySubVolIndex(mpiContext, mol, membraneObject, simulVolume);

    mol.receivedFromNeighborRank = true;

    // Ghost status is NOT decided here.  It used to be recomputed from
    // get_x_bin against this rank's own xOffset, while the sender had already
    // decided by column index -- two criteria on two ranks that did not
    // partition the boundary, so a molecule could end up owned by neither.
    // derive_ghost_from_ownership() sets it once, after the complexes carrying
    // the ownership stamp have arrived.

    int molIndex = find_molecule(moleculeList, mol.id);

    // Check whether the molecule was already in ghosted zone,
    // or it has just arrived (index == -1):
    if (molIndex != -1) {  // if found
      if (VERBOSE) printf("mol (id=%d) is already here.\n", mol.id);

      if (DEBUG && (moleculeList[molIndex].id != mol.id))
        error(mpiContext, "7 // deserialize_molecules()");
      mol.index = molIndex;

      // This rank's copy is the newer one when this rank changed the
      // molecule's complex at the top of the step.  The test used to be the
      // molecule's own trajStatus, but a state change marks only the molecule
      // it changed and the complex, so the other member of a dimer was
      // replaced with the stale copy while the complex itself was then kept:
      // the replaced member was left holding the complex ID the wire carries
      // in myComIndex, with nothing to map it back to an index.
      // The decision is the complex's: a propagated member of a complex the
      // neighbour changed is replaced with the rest of it.
      const int localCom = moleculeList[molIndex].myComIndex;
      const bool keepLocal =
          (localCom >= 0 && localCom < static_cast<int>(keptLocalCom.size()) &&
           keptLocalCom[localCom]);
      if (keepLocal) {
        if (VERBOSE)
          printf(
              "mol (id=%d) has been propagated, so skip the deserialization.\n",
              mol.id);
        moleculeList[molIndex].receivedFromNeighborRank = true;
        keptLocalMol[molIndex] = 1;
        if (CHECK_STRUCTURE) trace_molecule(mpiContext, "kept local", moleculeList, complexList, moleculeList[molIndex]);
        continue;
      }

      // Update mySubVolIndex if needed
      // (received molecule might have moved by another rank
      // since this rank serialized it)
      auto oldCellIndex = moleculeList[molIndex].mySubVolIndex;

      // check if the mol cross the ghosted boundary
      // int xBinOld = get_x_bin(mpiContext, moleculeList[molIndex]);
      // int xBinNew = get_x_bin(mpiContext, mol);
      // if (check_enter_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                              indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) enter the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesEnteringMolecules.push_back(molIndex);
      // }
      // if (check_exit_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                             indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) exit the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesExitingMolecules.push_back(molIndex);
      // }

      // A received copy carries no candidate reactions: the sender's refer
      // to its own indices and are not on the wire.  The candidates this
      // rank listed at the top of the step are still the ones its second
      // half draws from, so they stay -- unless the sender reacted or moved
      // the molecule, which it says by trajStatus, and then the molecule is
      // out of this step's draws the same way a reacted molecule is in a
      // serial run, where associate() clears both reactants' lists.
      if (mol.trajStatus != TrajStatus::propagated)
        mol.crossings = std::move(moleculeList[molIndex].crossings);
      // Nor is the reweighting history; see deserialize_molecules().
      mol.prevReweight = std::move(moleculeList[molIndex].prevReweight);
      mol.currReweight = std::move(moleculeList[molIndex].currReweight);
      moleculeList[molIndex] = mol;
      if (CHECK_STRUCTURE) trace_molecule(mpiContext, "overwritten", moleculeList, complexList, moleculeList[molIndex]);

      auto newCellIndex = moleculeList[molIndex].mySubVolIndex;
      // Check if cell index has changed:
      if (newCellIndex != oldCellIndex) {
        // Update simulVolume cells memberMolList of both mySubVolIndex:
        if (oldCellIndex < simulVolume.subCellList.size()) {
          auto &vec = simulVolume.subCellList[oldCellIndex].memberMolList;
          vec.erase(remove(vec.begin(), vec.end(), molIndex), vec.end());
        }
        // Deserialized molecule will be at the same index as existing one:
        if (newCellIndex < simulVolume.subCellList.size()) {
          simulVolume.subCellList[newCellIndex].memberMolList.push_back(
              molIndex);
        }
      }
    } else {  // if not found
      if (VERBOSE) printf("mol (id=%d) is new here.\n", mol.id);
      // The molecule doesn't exist in moleculeList. Add it:
      mol.index = moleculeList.size();
      molIndex = mol.index;
      moleculeList.push_back(mol);
      if (CHECK_STRUCTURE) trace_molecule(mpiContext, "appended", moleculeList, complexList, moleculeList.back());
    }
    indices.push_back(molIndex);
    // if(DEBUG) DEBUG_MOL("1DESERIALIZed A MOLECULE");
  }

  if (VERBOSE)
    cout << "Total deserialized number of molecules: " << nMols << endl;
}

void update_mySubVolIndex(MpiContext &mpiContext, Molecule &mol,
                          Membrane &membraneObject, SimulVolume &simulVolume) {
  if (VERBOSE) cout << "update_mySubVolIndex begins" << endl;
  // The molecule that a rank sends might return with new coordinates/cell
  // index; therefore, searching should be done at least around current cell:
  // ...for(auto &molIndex :
  // simulVolume.subCellList[mol.mySubVolIndex].memberList){... Another possible
  // optimization: receive action requrement to move the molecule to certain
  // cell Another possible optimization: for each molecule, keep left index and
  // right index for finding it in the moleculeList of neighbor rank
  // TODO:ANSW: (based on answer on another question)
  //         int currBin = mol.mySubVolIndex - mpiContext.xOffset;
  int xItr{int((mol.comCoord.x + membraneObject.waterBox.x / 2) /
               simulVolume.subCellSize.x) -
           mpiContext.xOffset};
  int yItr{int((mol.comCoord.y + membraneObject.waterBox.y / 2) /
               simulVolume.subCellSize.y)};
  int zItr;
  if (membraneObject.waterBox.z > 0)
    zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) /
               simulVolume.subCellSize.z);
  else
    zItr = 0;

  int currBin = xItr + (yItr * simulVolume.numSubCells.x) +
                (zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y);
  // Based on coordinates, assign index of cell:
  mol.mySubVolIndex = currBin;
  if (VERBOSE) cout << "update_mySubVolIndex ends" << endl;
}

bool check_enter_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew,
                              vector<int> &indicesEnteringMolecules,
                              SimulVolume &simulVolume) {
  if (mpiContext.rank && xBinOld == 1 && xBinNew == 0) {  // if not first rank
    return true;
  }
  if (mpiContext.rank < mpiContext.nprocs - 1 &&
      xBinOld == simulVolume.numSubCells.x - 2 &&
      xBinNew == simulVolume.numSubCells.x - 1) {  // if not last rank
    return true;
  }
  return false;
}

bool check_exit_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew,
                             vector<int> &indicesEnteringMolecules,
                             SimulVolume &simulVolume) {
  if (mpiContext.rank && xBinOld == 0 && xBinNew == 1) {  // if not first rank
    return true;
  }
  if (mpiContext.rank < mpiContext.nprocs - 1 &&
      xBinOld == simulVolume.numSubCells.x - 1 &&
      xBinNew == simulVolume.numSubCells.x - 2) {  // if not last rank
    return true;
  }
  return false;
}
