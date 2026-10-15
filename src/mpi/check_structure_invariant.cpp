#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

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
