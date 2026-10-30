#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "debug/debug.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

/*!
 * Every structural link between molecules and complexes on this rank has to
 * agree with its reverse: a molecule's myComIndex names a complex that lists
 * it, a complex's memberList names live molecules that name it back, and a
 * bound interface names a partner whose interface names it back and that sits
 * in the same complex.  The reaction code indexes through those links without
 * checking them, so a broken one is a heap read somewhere downstream, usually
 * much later than the step that broke it.
 *
 * Prints every violation (up to a cap) to stderr with the step and the label
 * of the phase just completed, then aborts the run, so the first report names
 * the phase that broke the structure.  Costs O(N) per call; runs only under
 * CHECK_STRUCTURE.
 */
// mode 0 checks everything.  Mid-exchange the structure is legitimately in
// transit: after deserialize_molecules() a received molecule's myComIndex holds
// the sender's complex ID until deserialize_complexes() maps it, and its
// partnerIndex values are the sender's until IDs_to_indices() runs, and a
// complex that the neighbour merged away still lists its member until
// delete_disappeared_complexes_partial() takes it.  mode 1 therefore checks
// only that memberLists name live slots; mode 2 adds the molecule -> complex
// direction and leaves the partner links and the complex -> molecule
// direction to mode 0.
// Molecule ids to trace through the exchanges, from NERDSS_TRACE_IDS (a
// comma-separated list; give it to every rank with mpirun -x).  Diagnostic
// for CHECK_STRUCTURE builds.
bool trace_id(int id) {
  static std::vector<int> ids;
  static bool read = false;
  if (!read) {
    read = true;
    const char *env = std::getenv("NERDSS_TRACE_IDS");
    if (env) {
      int v = 0;
      bool have = false;
      for (const char *c = env;; c++) {
        if (*c >= '0' && *c <= '9') {
          v = v * 10 + (*c - '0');
          have = true;
        } else {
          if (have) ids.push_back(v);
          v = 0;
          have = false;
          if (!*c) break;
        }
      }
    }
  }
  return std::find(ids.begin(), ids.end(), id) != ids.end();
}

void trace_molecule(MpiContext &mpiContext, const char *where,
                    const std::vector<Molecule> &moleculeList,
                    const std::vector<Complex> &complexList,
                    const Molecule &mol) {
  if (!trace_id(mol.id)) return;
  const int ci = mol.myComIndex;
  fprintf(stderr,
          "TRACE rank %d itr=%d %s: id=%d idx=%d x=%.2f ghost=%d recv=%d "
          "traj=%d com(idx=%d id=%d owner=%d) bonds:",
          mpiContext.rank, mpiContext.simItr, where, mol.id, mol.index,
          mol.comCoord.x, mol.isGhosted ? 1 : 0,
          mol.receivedFromNeighborRank ? 1 : 0, static_cast<int>(mol.trajStatus),
          ci,
          (ci >= 0 && ci < static_cast<int>(complexList.size())) ? complexList[ci].id : -999,
          (ci >= 0 && ci < static_cast<int>(complexList.size())) ? complexList[ci].ownerRank : -999);
  for (const auto &iface : mol.interfaceList)
    fprintf(stderr, " [idx=%d id=%d]", iface.interaction.partnerIndex,
            iface.interaction.partnerId);
  fprintf(stderr, "\n");
}

bool check_structure_invariant(MpiContext &mpiContext,
                               std::vector<Molecule> &moleculeList,
                               std::vector<Complex> &complexList,
                               long long simItr, const char *label, int mode) {
  int bad = 0;
  const int cap = 30;
  const int nMol = static_cast<int>(moleculeList.size());
  const int nCom = static_cast<int>(complexList.size());
  auto report = [&](const char *what, int molIdx, int comIdx) {
    if (bad++ >= cap) return;
    int molId = (molIdx >= 0 && molIdx < nMol) ? moleculeList[molIdx].id : -999;
    int comId = (comIdx >= 0 && comIdx < nCom) ? complexList[comIdx].id : -999;
    fprintf(stderr,
            "STRUCTURE VIOLATION rank %d itr=%lld at=%s : %s (mol idx=%d id=%d, "
            "complex idx=%d id=%d)\n",
            mpiContext.rank, simItr, label, what, molIdx, molId, comIdx, comId);
  };

  for (int i = 0; mode != 1 && i < nMol; i++) {
    const Molecule &mol = moleculeList[i];
    if (mol.isEmpty || mol.isImplicitLipid) continue;
    if (mol.index != i) report("mol.index != position", i, mol.myComIndex);
    int ci = mol.myComIndex;
    if (ci < 0 || ci >= nCom) {
      report("myComIndex out of range", i, ci);
      continue;
    }
    const Complex &com = complexList[ci];
    if (com.isEmpty) report("myComIndex names an empty complex", i, ci);
    else if (std::find(com.memberList.begin(), com.memberList.end(), i) ==
             com.memberList.end())
      report("complex does not list the molecule back", i, ci);
    for (size_t f = 0; mode == 0 && f < mol.interfaceList.size(); f++) {
      const auto &inter = mol.interfaceList[f].interaction;
      int p = inter.partnerIndex;
      if (p == -1) continue;
      if (p < 0 || p >= nMol) {
        report("partnerIndex out of range", i, ci);
        continue;
      }
      const Molecule &partner = moleculeList[p];
      if (partner.isImplicitLipid) continue;
      if (partner.isEmpty) {
        report("partner is an emptied slot", i, ci);
        continue;
      }
      int pf = inter.partnerIfaceIndex;
      if (pf < 0 || pf >= static_cast<int>(partner.interfaceList.size())) {
        report("partnerIfaceIndex out of range", i, ci);
        continue;
      }
      if (partner.interfaceList[pf].interaction.partnerIndex != i ||
          partner.myComIndex != ci) {
        // Both sides of the broken link, for the record.
        const int back = partner.interfaceList[pf].interaction.partnerIndex;
        const int pci = partner.myComIndex;
        fprintf(stderr,
                "  mol id=%d (ghost=%d, complex id=%d owner=%d) iface %zu -> "
                "partner idx=%d id=%d (ghost=%d, complex id=%d owner=%d) iface "
                "%d, whose partnerIndex=%d (id=%d)\n",
                mol.id, mol.isGhosted ? 1 : 0, com.id, com.ownerRank, f, p,
                partner.id, partner.isGhosted ? 1 : 0,
                (pci >= 0 && pci < nCom) ? complexList[pci].id : -999,
                (pci >= 0 && pci < nCom) ? complexList[pci].ownerRank : -999,
                pf, back,
                (back >= 0 && back < nMol) ? moleculeList[back].id : -999);
      }
      if (partner.interfaceList[pf].interaction.partnerIndex != i)
        report("partner interface does not point back", i, ci);
      if (partner.myComIndex != ci)
        report("bound partner is in another complex", i, ci);
    }
  }
  for (int c = 0; c < nCom; c++) {
    const Complex &com = complexList[c];
    if (com.isEmpty) continue;
    if (com.index != c) report("complex.index != position", -1, c);
    if (com.memberList.empty()) report("live complex with empty memberList", -1, c);
    for (size_t k = 0; k < com.memberList.size(); k++) {
      int m = com.memberList[k];
      if (m < 0 || m >= nMol) {
        report("member index out of range", m, c);
        continue;
      }
      const Molecule &mol = moleculeList[m];
      if (mol.isEmpty) report("member is an emptied slot", m, c);
      else if (mol.molTypeIndex < 0) report("member has molTypeIndex < 0", m, c);
      else if (mode == 0 && mol.myComIndex != c) report("member names another complex", m, c);
      if (std::find(com.memberList.begin() + k + 1, com.memberList.end(), m) !=
          com.memberList.end())
        report("member listed twice", m, c);
    }
  }
  if (bad) {
    fprintf(stderr, "STRUCTURE VIOLATION rank %d itr=%lld at=%s : %d violation(s); aborting\n",
            mpiContext.rank, simItr, label, bad);
    fflush(stderr);
    exit(3);
  }
  return true;
}

/*!
 * A complex is moved by one rank, its owner, and the two exchanges of every
 * step keep the neighbours' ghost copies current.  Ownership moves with the
 * stripe: send_data_to_*_neighboring_ranks() hands every complex that has a
 * member in the rank's ghost column to the neighbour on that side, whole.  A
 * complex that reaches both ghost columns of a middle rank is therefore
 * handed to the left neighbour at one exchange and to the right neighbour at
 * the next, every step; the two neighbours never hear from each other, each
 * stamp overrides the other on whichever rank receives it last, and the
 * complex ends up owned by nobody or by two ranks (CHECK_OWNERSHIP reports
 * it as ids held only as ghosts).  It is the width the stripe cannot carry:
 * more than a middle rank's owned columns, and nothing stops a chain or a
 * lattice from growing that wide.  (A member two ranks away from the owner,
 * the case the audit named, is the same thing: bonds are shorter than a
 * column, so a complex is contiguous in columns, and one reaching from rank
 * r - 1's columns into rank r + 1's covers both of rank r's ghost columns.)
 *
 * This reports it.  Once per step, after the step's last exchange, every
 * complex on the rank is measured member by member -- its lowest and highest
 * column, the x bin as get_x_bin() computes it made absolute with xOffset, so
 * members held outside the rank's own columns count too -- and compared with
 * every middle rank's owned columns (MpiContext::column_share(), the same
 * table on every rank): a complex whose columns include the one before and
 * the one after a middle rank's owned columns is reported, with the stamp it
 * carries here, whatever that stamp says.  At np <= 2 there is no middle
 * rank and nothing to check.
 *
 * Cost in the -O3 build: a complex of k members spans at most (k - 1) times
 * the longest bound-pair separation (twice the farthest interface from a
 * centre of mass, plus the largest binding radius), and reaching past both
 * ends of a middle rank's owned columns takes more than their width, so a
 * complex that cannot be that wide is skipped on its member count alone --
 * one integer compare per complex per step, and the member loop never runs
 * for trimers and hexamers.  Checker builds (CHECK_STRUCTURE) walk every
 * member of every complex.
 *
 * The -O3 build prints one warning line, on the first step a rank sees such
 * a complex, and report_complex_extent_summary() a count at the end of the
 * run; nothing is printed when it never happens.  Checker builds print every
 * offending complex (up to a cap per step) under the STRUCTURE VIOLATION
 * label the suite counts, and go on: the structure is intact, the physics is
 * what is wrong.
 */
namespace {
long long extentSteps = 0;       // steps on which this rank saw one
long long extentFirstStep = -1;
long long extentComplexes = 0;   // complex copies reported, all steps
}  // namespace

int check_complex_extent(MpiContext &mpiContext,
                         std::vector<Molecule> &moleculeList,
                         std::vector<Complex> &complexList,
                         const std::vector<MolTemplate> &molTemplateList,
                         const std::vector<ForwardRxn> &forwardRxns,
                         int totalColumns, long long simItr) {
  if (mpiContext.nprocs < 3) return 0;

  // Every middle rank's owned columns, and the shortest x-extent that can
  // reach past both ends of the narrowest of them, both fixed for the run.
  static std::vector<int> firstOwned, lastOwned;  // per middle rank
  static double minReach = 0.0;    // nm; a narrower complex cannot violate
  static double maxLink = 0.0;     // nm; longest bound-pair COM separation
  if (firstOwned.empty()) {
    int narrowest = -1;
    for (int r = 1; r < mpiContext.nprocs - 1; r++) {
      int first = 0, count = 0;
      mpiContext.column_share(totalColumns, r, first, count);
      firstOwned.push_back(first);
      lastOwned.push_back(first + count - 1);
      if (narrowest < 0 || count < narrowest) narrowest = count;
    }
    minReach = narrowest * mpiContext.simulVolume->subCellSize.x;
    double reach = 0.0, sigma = 0.0;
    for (const auto &oneTemp : molTemplateList)
      for (const auto &iface : oneTemp.interfaceList)
        reach = std::max(reach, iface.iCoord.length());
    for (const auto &rxn : forwardRxns)
      sigma = std::max(sigma, std::max(rxn.bindRadius, rxn.bindRadius2D));
    maxLink = 2.0 * reach + sigma;
  }

  const int nMol = static_cast<int>(moleculeList.size());
  const int cap = 20;
  int bad = 0;
  for (const auto &com : complexList) {
    if (com.isEmpty) continue;
    const int k = static_cast<int>(com.memberList.size());
    if (!CHECK_STRUCTURE && (k - 1) * maxLink <= minReach) continue;
    int lo = -1, hi = -1;
    for (int m : com.memberList) {
      if (m < 0 || m >= nMol) continue;
      const Molecule &mol = moleculeList[m];
      if (mol.isEmpty || mol.isImplicitLipid) continue;
      int col = get_x_bin(mpiContext, moleculeList[m]) + mpiContext.xOffset;
      if (col < 0) col = 0;
      if (col >= totalColumns) col = totalColumns - 1;
      if (lo < 0 || col < lo) lo = col;
      if (hi < 0 || col > hi) hi = col;
    }
    if (lo < 0) continue;
    int middle = -1;
    for (size_t i = 0; i < firstOwned.size(); i++)
      if (lo < firstOwned[i] && hi > lastOwned[i]) {
        middle = static_cast<int>(i) + 1;
        break;
      }
    if (middle < 0) continue;
    bad++;
    if (CHECK_STRUCTURE) {
      if (bad <= cap)
        fprintf(stderr,
                "STRUCTURE VIOLATION rank %d itr=%lld at=complex extent : "
                "complex id=%d (%d members, owner stamp %d) spans columns "
                "%d..%d, past both ends of rank %d's columns %d..%d; it is "
                "handed to both of that rank's neighbours in turn\n",
                mpiContext.rank, simItr, com.id, k, com.ownerRank, lo, hi,
                middle, firstOwned[middle - 1], lastOwned[middle - 1]);
    } else if (extentSteps == 0 && bad == 1) {
      fprintf(stderr,
              "WARNING rank %d itr=%lld: complex id=%d (%d members) spans "
              "columns %d..%d, past both ends of rank %d's columns %d..%d; "
              "the stripe hands it to both of that rank's neighbours in turn, "
              "so this complex and what binds it are unreliable from here on "
              "(one warning per rank; a count follows at the end of the run)\n",
              mpiContext.rank, simItr, com.id, k, lo, hi, middle,
              firstOwned[middle - 1], lastOwned[middle - 1]);
    }
  }
  if (bad) {
    if (extentSteps == 0) extentFirstStep = simItr;
    extentSteps++;
    extentComplexes += bad;
  }
  return bad;
}

void report_complex_extent_summary(MpiContext &mpiContext) {
  if (extentSteps == 0) return;
  fprintf(stderr,
          "WARNING rank %d: a complex wider than the stripe carries on %lld "
          "step(s), first at step %lld, %lld complex copies in all\n",
          mpiContext.rank, extentSteps, extentFirstStep, extentComplexes);
}
