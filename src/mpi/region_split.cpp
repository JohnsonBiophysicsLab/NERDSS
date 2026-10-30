#include <array>

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
//
// The left ghost column has the mirror rule.  A pair with both members in it
// lies inside the left neighbour's last owned column and is the left
// neighbour's, which lists it as its own and keeps it.  measure_separations()
// lists such a pair here too whenever one member is owned here -- a member of
// a complex straddling the boundary -- and since it pairs ghosts with owned
// molecules (it used to skip ghosts altogether, and with them most pairs
// across the boundary) it can list more; drop them, or both ranks draw.
void drop_crossings_decided_elsewhere(MpiContext &mpiContext,
                                      vector<Molecule> &moleculeList,
                                      vector<Complex> &complexList) {
  const bool hasRight = mpiContext.rank < mpiContext.nprocs - 1;
  const bool hasLeft = mpiContext.rank > 0;
  if (!hasRight && !hasLeft) return;
  const int rightGhost = mpiContext.simulVolume->numSubCells.x - 1;
  // 1: in the right ghost column, 2: in the left ghost column
  vector<char> where(moleculeList.size(), 0);
  for (size_t i = 0; i < moleculeList.size(); i++) {
    Molecule &mol = moleculeList[i];
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    const int bin = get_x_bin(mpiContext, mol);
    if (hasRight && bin == rightGhost)
      where[i] = 1;
    else if (hasLeft && bin == 0)
      where[i] = 2;
  }
  for (size_t i = 0; i < moleculeList.size(); i++) {
    Molecule &mol = moleculeList[i];
    if (mol.isEmpty || mol.isImplicitLipid || mol.crossings.empty()) continue;
    int dropped = 0;
    for (int k = static_cast<int>(mol.crossings.size()) - 1; k >= 0; --k) {
      const int p = mol.crossings[k].partner;
      const int pw = (p >= 0 && p < static_cast<int>(moleculeList.size()))
                         ? where[p]
                         : 0;
      const bool rightsPair = where[i] == 1 || pw == 1;
      const bool leftsPair = where[i] == 2 && pw == 2;
      if (!rightsPair && !leftsPair) continue;
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

// How many candidate interface pairs this rank holds that it can act on --
// listed after drop_crossings_decided_elsewhere() -- summed over the ranks
// on rank 0, and how many of them pair two ghosts (members of a
// neighbour's straddling complexes, which this rank decides since
// measure_separations() stopped skipping them).  At step 1 nothing has moved
// and every decomposition starts from rank 0's coordinates, so the total has
// to be the np=1 count: each pair listed by exactly one rank; and no ghost
// pairs with a ghost there, since nothing straddles yet.  A diagnostic for
// CHECK_STRUCTURE builds; collective.
std::array<long long, 2> count_candidate_pairs(MpiContext &mpiContext,
                                               vector<Molecule> &moleculeList) {
  std::array<long long, 2> local{0, 0};
  for (const auto &mol : moleculeList) {
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    for (const auto &cross : mol.crossings) {
      const int p = cross.partner;
      if (p < 0 || p >= static_cast<int>(moleculeList.size())) continue;
      const Molecule &partner = moleculeList[p];
      if (partner.isImplicitLipid) continue;
      // each pair sits in both lists; count it from the lower id's
      if (mol.id < partner.id) {
        local[0]++;
        if (mol.isGhosted && partner.isGhosted) local[1]++;
      }
    }
  }
  std::array<long long, 2> total{0, 0};
  MPI_Reduce(local.data(), total.data(), 2, MPI_LONG_LONG, MPI_SUM, 0,
             MPI_COMM_WORLD);
  return total;
}

// How often a pair evaluated this step found last step's reweighting history
// (determine_*_bimolecular_reaction_probability() looks it up in the lower
// index's prevReweight by partner and faces), split by where the molecule
// holding it sits: an edge column (0, 1, n-2, n-1), whose copies the
// exchanges overwrite, or the interior.  Without per-rank bias the two hit
// rates agree.  Accumulated over the run; returns the four sums over ranks
// on rank 0 as {edge hits, edge pairs, interior hits, interior pairs}.  A
// diagnostic for CHECK_STRUCTURE builds; collective.
std::array<long long, 4> reweight_history_hits(MpiContext &mpiContext,
                                               vector<Molecule> &moleculeList) {
  static std::array<long long, 4> local{0, 0, 0, 0};
  const int n = mpiContext.simulVolume->numSubCells.x;
  for (const auto &mol : moleculeList) {
    if (mol.isEmpty || mol.isImplicitLipid || mol.isGhosted) continue;
    const int bin = get_x_bin(mpiContext, const_cast<Molecule &>(mol));
    const bool edge = (mpiContext.rank > 0 && (bin == 0 || bin == 1)) ||
                      (mpiContext.rank < mpiContext.nprocs - 1 &&
                       (bin == n - 2 || bin == n - 1));
    for (const auto &curr : mol.currReweight) {
      bool hit = false;
      for (const auto &prev : mol.prevReweight)
        if (prev.partner == curr.partner && prev.myFace == curr.myFace &&
            prev.partnerFace == curr.partnerFace) {
          hit = true;
          break;
        }
      local[edge ? 0 : 2] += hit;
      local[edge ? 1 : 3] += 1;
    }
  }
  std::array<long long, 4> total{0, 0, 0, 0};
  MPI_Reduce(local.data(), total.data(), 4, MPI_LONG_LONG, MPI_SUM, 0,
             MPI_COMM_WORLD);
  return total;
}
