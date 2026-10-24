#include "debug/debug.hpp"
#include "mpi/mpi_function.hpp"

// Where this rank's visible columns split into the left region, acted on in
// the first half of the step, and the right region, acted on in the second.
// Column (x bin) indices below the returned value are in the left region.
//
// The split is tied to the two exchanges of the step, not to load.  After the
// first half a rank sends its bins 0 and 1 to the left neighbour, and that
// neighbour's copy of bin 1 comes back after the second half as the copy of
// record: positions and all, with only the bonds the second half made here
// kept local.  So bin 1 has to be finished by the end of the first half, or
// it is put back to where it was.  Symmetrically, the right neighbour binds
// into bin n-2 during ITS first half and its copy arrives at this rank's
// first exchange, where this rank keeps whatever it moved in the first half;
// so bin n-2 must not be moved before that exchange.  A rank that owns two
// columns or more satisfies both: bin 1 goes left, bin n-2 goes right, and
// the formulas below put the first floor(owned/2) owned columns on the left.
//
// A rank that owns one column has bin 1 == bin n-2, and the original
// formulas put it on the right: the column was moved in the second half,
// after the left neighbour had been sent its copy, and that copy came back
// and put every molecule back where it started.  Of the two requirements the
// first is the one no guard can make up for, so a single owned column goes
// left: it is moved before the left neighbour sees it, and the left
// neighbour's returned copy is current.  That used to leave the right
// neighbour's first-half bonds into the column meeting a column this rank
// was acting on in the same half-step, and the exchange could keep only one
// of the two; now only the last rank is ever given a single column
// (MpiContext::column_share() gives a middle rank at least three), and it has
// no right neighbour.
//
// Rank 0 never receives from the left, so its only constraint is the right
// one, and its formula already puts a single owned column (bin 0) on the
// right.  With one rank there is no exchange at all and the split is only an
// order of evaluation, which the rank-0 formula keeps as it was.
//
// Rank 0 with two owned columns is the exception the other way: (n - 1) / 2
// put bin 0 on the left, and a bin-0 initiator pairs with bin 1, which is
// bin n-2 -- the column the right neighbour binds into during the same first
// half.  Two ranks then bound one molecule at once.  With nothing on its left
// to satisfy, rank 0 can leave the first half empty, as it already did with
// one column.  Middle ranks cannot (their bin 1 must go first), so they are
// given enough columns not to meet the problem; see
// MpiContext::column_share().
int region_mid_index(MpiContext &mpiContext) {
  const int n = mpiContext.simulVolume->numSubCells.x;
  int mid;
  if (mpiContext.rank == 0)
    mid = (mpiContext.endCell - mpiContext.startCell + 1 <= 2) ? 0
                                                              : (n - 1) / 2;
  else if (mpiContext.rank == mpiContext.nprocs - 1)
    mid = (n + 1) / 2;
  else
    mid = n / 2;
  if (mpiContext.rank > 0 && mpiContext.endCell == mpiContext.startCell)
    mid = 2;  // the left ghost column and the one owned column
  return mid;
}

// A pair that straddles a column boundary is seen by both ranks: the right
// rank lists it for its bin-1 initiator against a bin-0 ghost, the left rank
// for its bin n-2 initiator against a bin n-1 ghost.  Both would draw for it,
// so the pair is the right rank's to decide, in its first half, and the left
// rank drops its candidates here, right after measure_separations() listed
// them.  That is the single evaluation a serial run gives every pair, and it
// keeps the left rank from binding a molecule that the right rank, in the
// same half-step, may already have bound to someone else.
//
// Until now this happened by accident: the right neighbour's copy of bins
// n-2 and n-1 replaced this rank's at the first exchange, and a received
// molecule carries no candidate list, so the bin n-2 initiators had nothing
// left to draw from in the second half.  A rank that owns one column does
// its first half before that exchange, so its candidates against the right
// ghost column were live, and the rule has to be explicit.
void drop_right_ghost_crossings(MpiContext &mpiContext,
                                vector<Molecule> &moleculeList,
                                vector<Complex> &complexList) {
  if (mpiContext.rank == mpiContext.nprocs - 1) return;
  const int rightGhost = mpiContext.simulVolume->numSubCells.x - 1;
  vector<char> inRightGhost(moleculeList.size(), 0);
  for (size_t i = 0; i < moleculeList.size(); i++) {
    Molecule &mol = moleculeList[i];
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    inRightGhost[i] = get_x_bin(mpiContext, mol) == rightGhost;
  }
  for (size_t i = 0; i < moleculeList.size(); i++) {
    Molecule &mol = moleculeList[i];
    if (mol.isEmpty || mol.isImplicitLipid || mol.crossings.empty()) continue;
    int dropped = 0;
    for (int k = static_cast<int>(mol.crossings.size()) - 1; k >= 0; --k) {
      const int p = mol.crossings[k].partner;
      const bool partnerInGhost = p >= 0 &&
                                  p < static_cast<int>(moleculeList.size()) &&
                                  inRightGhost[p];
      if (!inRightGhost[i] && !partnerInGhost) continue;
      mol.crossings.erase(mol.crossings.begin() + k);
      dropped++;
    }
    // Every entry was counted in the complex's ncross by
    // record_crossing_pair(); uncount what was dropped, as
    // remove_partner_crossings() does.
    if (dropped && mol.myComIndex >= 0 &&
        mol.myComIndex < static_cast<int>(complexList.size()))
      complexList[mol.myComIndex].ncross -= dropped;
  }
}
