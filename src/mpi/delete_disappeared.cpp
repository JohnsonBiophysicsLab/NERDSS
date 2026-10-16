#include <iostream>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

// A complex created after start-up -- by an association or a dissociation --
// carries Complex::ownerRank's default of -1: prepare() seeds ownership from the
// initial decomposition and nothing assigns it again, because the reaction code
// is shared with the serial build and has no rank to record.  ownerRank is the
// authority for who integrates and who tallies a complex, so an unclaimed
// complex belongs to nobody.
//
// Claim the ones this rank made.  Two rules, either suffices, and both are
// only ever applied to a complex no rank has stamped yet, so neither can
// contradict a handover.
//
// By position: a complex is this rank's if it has a member this rank owns by
// position.  This is the rule prepare() seeds with, and the position test
// partitions the bins cleanly, so two ranks cannot both claim.
//
// By descent: a complex is this rank's if none of its members is a ghost.
// isGhosted is derived from the owning complex's stamp at each exchange, and
// the only reaction code that writes it is a creation, for the molecule it
// makes; so a member that is not a ghost was, at the last exchange, in a
// complex this rank owned, and an unstamped complex made of such members was
// split off one of this rank's complexes by a dissociation at the top of this
// step, or created here.  What this rank dissociates is this rank's, wherever
// the product's members sit.
//
// The position rule alone left one product unowned for half a step: a dimer
// this rank owns straddling its last owned column and the right neighbour's
// first, dissociated at the top of the step with the fresh complex going to
// the member in the neighbour's column.  Nothing here owns that column by
// position, so the product kept ownerRank -1 through the send to the left;
// the right neighbour's copy of it, still bound in the parent, is that
// rank's ghost of a complex owned here, and derive_ghost_from_ownership()
// at the receive from the right then marked this rank's copy a ghost too.
// Owned by nobody, it stayed that way through the second half until the
// send to the right stamped it over (the bin n-1 loop in
// send_data_to_right_*), one exchange late.  The mirror image on the left
// never showed: a product in bin 0 keeps its parent's non-ghost flag and is
// stamped over by the send to the left before any exchange re-derives the
// flag.  Claiming by descent before that first exchange makes the right side
// behave the same way: the product is this rank's through the second half,
// in which the right ghost column lies, and the next send hands it to the
// neighbour.  Nothing else about the step changes: the dissociation has
// already displaced the product and marked it propagated, and a pair across
// the right boundary is the neighbour's to decide
// (drop_right_ghost_crossings()), so an owned product does nothing more in
// that half-step than an unowned one did.  What the claim restores is the
// invariant everything else reads -- one owner for every complex at every
// checkpoint -- for the tallies, for the far-neighbour rules keyed on
// ownerRank, and for check_ownership_invariant().
//
// Two ranks cannot both claim by descent either: a member can be a non-ghost
// on one rank only, the rank whose complex it was in.  Nor can descent pick
// up a complex that is a neighbour's: an association's survivor is the
// initiator's complex and keeps its stamp, so the only unstamped complexes
// are this step's dissociation and creation products.
void claim_unowned_complexes(MpiContext &mpiContext,
                             vector<Molecule> &moleculeList,
                             vector<Complex> &complexList,
                             SimulVolume &simulVolume) {
  for (auto &com : complexList) {
    if (com.isEmpty || com.ownerRank != -1) continue;
    bool byPosition = false;
    bool byDescent = true;
    int counted = 0;
    for (int m : com.memberList) {
      if (m < 0 || m >= (int)moleculeList.size()) continue;
      if (moleculeList[m].isEmpty || moleculeList[m].isImplicitLipid) continue;
      counted++;
      if (moleculeList[m].isGhosted) byDescent = false;
      if (!byPosition &&
          is_owned_by_processor(moleculeList[m], mpiContext, simulVolume))
        byPosition = true;
    }
    // The implicit lipid is exempt from the ghost rule and is skipped above,
    // so a complex with no counted member says nothing about descent.
    if (counted == 0) byDescent = false;
    if (!byPosition && !byDescent) continue;
    com.ownerRank = mpiContext.rank;
    // isGhosted is derived from ownerRank, so it has to follow: write_all_species()
    // reads the flag, and check_ownership_invariant() classifies by it.
    for (int m : com.memberList) {
      if (m < 0 || m >= (int)moleculeList.size()) continue;
      if (moleculeList[m].isEmpty) continue;
      moleculeList[m].isGhosted = false;
    }
  }
}

void delete_disappeared_complexes(MpiContext &mpiContext,
                                  vector<Molecule> &moleculeList,
                                  vector<Complex> &complexList) {
  if (VERBOSE) cout << "delete_disappeared_complexes begins" << endl;
  int count{0};
  for (auto &com : complexList) {
    if (!com.receivedFromNeighborRank) {
      if (!com.deleteIfNotReceivedBack) {
        // TODO:ASK: (complex management):
        //  update myComIndex of com.memberList molecules to match new molecule.

      } else {
        int zone_index = get_x_bin(mpiContext, moleculeList[com.memberList[0]]);
        if (!com.isEmpty) {
          count++;
          if (VERBOSE)
            cout << "### Deleting complex id=" << com.id
                 << ", complex index=" << com.index << endl
                 << "        zone index=" << zone_index;
          com.memberList.clear();
          com.destroy(moleculeList, complexList);
        }
      }
    }
  }

  if (VERBOSE)
    cout << "delete_disappeared_complexes ends" << endl
         << count << " complexes are deleted" << endl;
}

void delete_disappeared_complexes_partial(MpiContext &mpiContext,
                                          vector<Molecule> &moleculeList,
                                          vector<Complex> &complexList) {
  if (VERBOSE) cout << "delete_disappeared_complexes_partial begins" << endl;
  int count{0};
  // How many molecules still on this rank claim each complex through
  // myComIndex.  Normally a complex that was not received back has none left:
  // its molecules were not received back either and
  // delete_disappeared_molecules() has already gone through them.  When some
  // did come back, destroying the complex would leave them claiming a complex
  // that no longer exists -- and memberList.clear() below means destroy() does
  // not take them with it.  Such a molecule reaches complexList through
  // myComIndex in create_complex_propagation_vectors(),
  // check_bimolecular_reactions() and the ncross reset at the top of the step.
  vector<int> claims(complexList.size(), 0);
  for (auto &mol : moleculeList) {
    if (mol.isEmpty) continue;
    if (mol.myComIndex >= 0 &&
        mol.myComIndex < static_cast<int>(complexList.size()))
      claims[mol.myComIndex]++;
  }
  // A complex the neighbour merged into another, or that moved beyond the
  // neighbour's edge, does not come back, and its molecules either did not
  // come back either -- delete_disappeared_molecules() has had them -- or
  // came back naming the complex they are in now.  Whether this rank's copy
  // goes is decided by whether any molecule here still names it.  One that
  // is still claimed is stale rather than gone, and the next exchange
  // replaces it.  One that nobody names is garbage wherever its last member
  // now sits, whoever owned it, and whatever receivedFromNeighborRank says:
  // that flag is one bit for both directions, and on a rank that owns one
  // column the left neighbour's echo of the column sets it after the send
  // to the right has cleared it, so a complex the right neighbour merged
  // away read as received.  The region gate this pass used to apply had the
  // same hole from the other side: a member the neighbour had just moved
  // across the regions' boundary put its old complex beyond both passes.
  for (auto &com : complexList) {
    if (com.isEmpty) continue;
    if (com.index >= 0 && com.index < static_cast<int>(claims.size()) &&
        claims[com.index] > 0)
      continue;
    count++;
    if (VERBOSE)
      cout << "### Deleting complex id=" << com.id
           << ", complex index=" << com.index << endl;
    com.memberList.clear();
    com.destroy(moleculeList, complexList);
  }

  if (VERBOSE)
    cout << "delete_disappeared_complexes_partial ends" << endl
         << count << " complexes are deleted" << endl;
}

void delete_disappeared_molecules(
    MpiContext &mpiContext, SimulVolume &simulVolume, Membrane &membraneObject,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<int> &region, vector<MolTemplate> &molTemplateList,
    bool fromLeft) {
  if (VERBOSE) cout << "delete_disappeared_molecules begins" << endl;
  int count{0};
  // The neighbour that just sent answers only for the complexes it owns.  A
  // ghost of the other neighbour's complex that sits in this region -- on a
  // rank with one owned column, that column is mirrored by both neighbours,
  // and a complex the right neighbour owns can have a member in it -- is not
  // in this message because the sender never had it, not because it is gone.
  const int farNeighbour = fromLeft ? mpiContext.rank + 1 : mpiContext.rank - 1;
  for (auto &molIndex : region) {
    Molecule &mol = moleculeList[molIndex];
    if (mol.myComIndex == -1 || mol.isImplicitLipid == true || mol.isEmpty == true)
      continue;  // no checking for deleted molecules, implicit lipid
    if (mol.myComIndex >= 0 &&
        mol.myComIndex < static_cast<int>(complexList.size()) &&
        complexList[mol.myComIndex].ownerRank == farNeighbour)
      continue;

    // Only a ghost may be dropped here.  isShared marks a molecule this rank
    // owns and a neighbour mirrors, and receivedFromNeighborRank is a single
    // flag for both directions: sending to the left clears it, and this pass --
    // which is checking what the right neighbour returned -- then read that as
    // "the right neighbour dropped it" and deleted a molecule this rank owns.
    //
    // A middle rank is where the two overlap.  With one owned cell its left edge
    // stripe and its right edge stripe are the same cell, and that cell is in
    // the left region (see region_mid_index()), so every molecule it owns is
    // offered to this test at the receive from the left every step.  A
    // molecule genuinely handed over is marked isGhosted by the ghost-stripe
    // loop before it is sent, so handovers are unaffected.
    if (mol.isGhosted) {
      if (mol.receivedFromNeighborRank == false) {
        // Every crossing this molecule holds was counted in its complex's
        // ncross by record_crossing_pair(), so uncount them while myComIndex
        // still says which complex that was.  The entries pointing the other
        // way are dropped by drop_references_to_absent_molecules(), which has
        // to wait until IDs_to_indices() has run.
        if (mol.myComIndex >= 0 &&
            mol.myComIndex < static_cast<int>(complexList.size()))
          complexList[mol.myComIndex].ncross -=
              static_cast<int>(mol.crossings.size());
        mol.crossings.clear();
        // Remove the molecule from this rank
        // and remove it from complex as well
        mol.MPI_remove_from_one_rank(moleculeList, complexList);
        count++;
      }
    }
  }
  if (VERBOSE)
    cout << "delete_disappeared_molecules ends" << endl
         << count << " molecules are deleted" << endl;
  // cout << "delete_disappeared_molecules ends" << endl << count << " molecules
  // are deleted" << endl;
}

// Nothing may still name a molecule that has left this rank.  Molecules go in
// delete_disappeared_molecules(), which runs before IDs_to_indices(), so it is
// in no position to do this itself: at that point the molecules just taken off
// the wire still carry the *sender's* partnerIndex values, and bndpartner is
// whatever the previous exchange left.  Only once IDs_to_indices() has mapped
// the received IDs back to local indices does every list mean what it says, and
// the ones that survived the exchange untouched are exactly the ones nothing
// has revisited -- which is how a molecule ends up still bound to a slot that
// MPI_remove_from_one_rank() emptied.
//
// A departed molecule's slot keeps its id and its size, so these references are
// not caught by a bounds check: they read a molecule with a cleared
// interfaceList and myComIndex == -1, which is complexList[-1] in
// sweep_separation_box() and in check_dissociation().
void drop_references_to_absent_molecules(vector<Molecule> &moleculeList,
                                         vector<Complex> &complexList) {
  auto absent = [&moleculeList](int idx) {
    return idx < 0 || idx >= static_cast<int>(moleculeList.size()) ||
           moleculeList[idx].isEmpty;
  };
  for (auto &mol : moleculeList) {
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    // -1 is what the rest of the code reads as "this interface's partner is not
    // visible on this rank".  partnerId is deliberately left alone, so the bond
    // itself survives and IDs_to_indices() can restore the index if the partner
    // comes back.
    for (auto &iface : mol.interfaceList) {
      if (iface.interaction.partnerIndex != -1 &&
          absent(iface.interaction.partnerIndex))
        iface.interaction.partnerIndex = -1;
    }
    // bndpartner and bndlist are parallel, so an entry leaves both.
    for (int i = static_cast<int>(mol.bndpartner.size()) - 1; i >= 0; --i) {
      if (mol.bndpartner[i] != -1 && !absent(mol.bndpartner[i])) continue;
      mol.bndpartner.erase(mol.bndpartner.begin() + i);
      if (i < static_cast<int>(mol.bndlist.size()))
        mol.bndlist.erase(mol.bndlist.begin() + i);
    }
    for (int i = static_cast<int>(mol.crossings.size()) - 1; i >= 0; --i) {
      if (!absent(mol.crossings[i].partner)) continue;
      mol.crossings.erase(mol.crossings.begin() + i);
      if (mol.myComIndex >= 0 &&
          mol.myComIndex < static_cast<int>(complexList.size()))
        --complexList[mol.myComIndex].ncross;
    }
  }
}

void disconnect_molecule_partners(unsigned &targMolIndex, Molecule &mol,
                                  vector<Molecule> &moleculeList) {
  if (VERBOSE) cout << "disconnect_molecule_partners begins" << endl;
  // Remove partnerIndex from mol.bndpartner
  for (auto &partnerIndex : mol.bndpartner) {
    if (partnerIndex != -1) {
      Molecule &partner = moleculeList[partnerIndex];

      // Erasing while iterating forwards without stepping back skipped the
      // second of two consecutive entries.  Walking backwards removes every
      // match and leaves the indices below the cursor untouched.  bndlist is
      // erased at the same slot, which is the invariant break_interaction.cpp
      // now maintains; the size guard keeps this safe if some other path has
      // not.
      for (int i = static_cast<int>(partner.bndpartner.size()) - 1; i >= 0; --i) {
        if (partner.bndpartner[i] == targMolIndex) {
          partner.bndpartner.erase(partner.bndpartner.begin() + i);
          if (i < static_cast<int>(partner.bndlist.size()))
            partner.bndlist.erase(partner.bndlist.begin() + i);
        }
      }

      // Set partnerIndex to -1 in mol.iFace[].interaction.partnerIndex
      for (auto &iFace : partner.interfaceList) {
        if (iFace.interaction.partnerIndex == targMolIndex) {
          iFace.interaction.partnerIndex = -1;
        }
      }
    }
  }
  if (VERBOSE) cout << "disconnect_molecule_partners ends" << endl;
}
