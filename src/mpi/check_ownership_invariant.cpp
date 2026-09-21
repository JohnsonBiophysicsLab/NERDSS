#include <algorithm>
#include <map>
#include <vector>

#include "io/io.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

/*!
 * The parallel scheme rests on one invariant: every molecule instance is owned
 * by exactly one rank, and the copies other ranks keep are ghosts.  Ownership
 * loss is silent -- an unowned molecule is integrated by nobody and skipped by
 * every tally that filters ghosts, so it shows up much later as a small drift in
 * copy numbers, or not at all.
 *
 * This checks the invariant directly rather than inferring it from totals.  A
 * reference count would be wrong for any system with creation or destruction
 * reactions, where the total legitimately changes; "owned by exactly one rank"
 * holds regardless.
 *
 * Every rank sends the ids it owns and the ids it holds as ghosts to rank 0,
 * which reports ids owned twice and ids held but owned nowhere.  That is O(N)
 * communication, so this is meant to run every few hundred steps under
 * CHECK_OWNERSHIP, not on the hot path.
 *
 * Set CHECK_OWNERSHIP_EVERY to 1 when hunting a leak.  A coarse period will
 * miss the step on which a molecule is dropped, and reports only tell you it
 * was already lost some time earlier.
 *
 * Returns true when the invariant holds.  The implicit lipid is exempt: it is
 * deliberately present on every rank.
 */
bool check_ownership_invariant(MpiContext &mpiContext,
                               std::vector<Molecule> &moleculeList,
                               std::vector<Complex> &complexList,
                               long long simItr, const char *label) {
  std::vector<int> ownedLocal, ghostLocal;
  ownedLocal.reserve(moleculeList.size());
  for (const auto &mol : moleculeList) {
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    if (mol.isGhosted)
      ghostLocal.push_back(mol.id);
    else
      ownedLocal.push_back(mol.id);
  }

  auto gather = [&](const std::vector<int> &local, std::vector<int> &all) {
    int n = static_cast<int>(local.size());
    std::vector<int> counts(mpiContext.nprocs, 0);
    MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> displs(mpiContext.nprocs, 0);
    int total = 0;
    if (!mpiContext.rank) {
      for (int r = 0; r < mpiContext.nprocs; r++) {
        displs[r] = total;
        total += counts[r];
      }
      all.assign(total, 0);
    }
    MPI_Gatherv(local.data(), n, MPI_INT, mpiContext.rank ? nullptr : all.data(),
                counts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);
    return counts;
  };

  std::vector<int> ownedAll, ghostAll;
  std::vector<int> ownedCounts = gather(ownedLocal, ownedAll);
  gather(ghostLocal, ghostAll);

  int bad = 0;
  std::vector<int> firstUnowned;
  if (!mpiContext.rank) {
    std::map<int, int> ownedTimes;
    for (int id : ownedAll) ++ownedTimes[id];

    std::vector<int> doubleOwned, unowned;
    for (const auto &kv : ownedTimes)
      if (kv.second > 1) doubleOwned.push_back(kv.first);
    for (int id : ghostAll)
      if (!ownedTimes.count(id)) unowned.push_back(id);
    std::sort(unowned.begin(), unowned.end());
    unowned.erase(std::unique(unowned.begin(), unowned.end()), unowned.end());

    if (!doubleOwned.empty() || !unowned.empty()) {
      bad = 1;
      firstUnowned = unowned.empty() ? doubleOwned : unowned;
      fprintf(stderr,
          "OWNERSHIP VIOLATION itr=%lld at=%s : %zu id(s) owned by more than "
          "one rank, %zu id(s) held only as ghosts and owned by nobody\n",
          simItr, label, doubleOwned.size(), unowned.size());
      size_t show = 8;
      if (!doubleOwned.empty()) {
        fprintf(stderr, "  owned twice:");
        for (size_t i = 0; i < doubleOwned.size() && i < show; i++)
          fprintf(stderr, " %d", doubleOwned[i]);
        if (doubleOwned.size() > show) fprintf(stderr, " ...");
        fprintf(stderr, "\n");
      }
      if (!unowned.empty()) {
        fprintf(stderr, "  owned by nobody:");
        for (size_t i = 0; i < unowned.size() && i < show; i++)
          fprintf(stderr, " %d", unowned[i]);
        if (unowned.size() > show) fprintf(stderr, " ...");
        fprintf(stderr, "\n");
      }
      fprintf(stderr, "  owned per rank:");
      for (int r = 0; r < mpiContext.nprocs; r++)
        fprintf(stderr, " %d", ownedCounts[r]);
      fprintf(stderr, "\n");
    }
  }
  // every rank leaves with the same verdict, so callers may branch on it
  MPI_Bcast(&bad, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // When the invariant fails, every rank reports what it knows about the first
  // offending id.  Ownership bugs are only diagnosable from both sides at once,
  // and reconstructing that from two separate logs afterwards is painful.
  if (bad) {
    int firstBad = -1;
    if (!mpiContext.rank && !firstUnowned.empty()) firstBad = firstUnowned.front();
    MPI_Bcast(&firstBad, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (firstBad >= 0) {
      for (int r = 0; r < mpiContext.nprocs; r++) {
        if (r == mpiContext.rank) {
          bool found = false;
          for (const auto &mol : moleculeList) {
            if (mol.id != firstBad) continue;
            found = true;
            int ci = mol.myComIndex;
            int owner = (ci >= 0 && ci < (int)complexList.size())
                            ? complexList[ci].ownerRank
                            : -999;
            int comId = (ci >= 0 && ci < (int)complexList.size())
                            ? complexList[ci].id
                            : -999;
            fprintf(stderr,
                    "    id=%d on rank %d: ghost=%d empty=%d myComIndex=%d "
                    "complexId=%d complexOwnerRank=%d nCom=%zu\n",
                    firstBad, mpiContext.rank, (int)mol.isGhosted,
                    (int)mol.isEmpty, ci, comId, owner, complexList.size());
          }
          if (!found)
            fprintf(stderr, "    id=%d on rank %d: not present\n", firstBad,
                    mpiContext.rank);
          fflush(stderr);
        }
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
  return bad == 0;
}
