/*! \file rank_tallies.cpp
 * \brief The per-rank versions of the tallies merge_outputs() sums: a rank
 *        reports what it owns, so that the ranks' files add up to the
 *        serial build's.
 */
#include "io/io.hpp"
#include "mpi_context.hpp"

// The serial print_dimers() counts every non-empty complex in the list, which
// on a rank includes the ghost copies of the neighbours' straddling
// complexes: summed over the ranks, every such complex is counted twice.
// Count a complex on the rank that owns it, as print_complex_hist() does; the
// columns and the format are the serial build's.  Same classification as the
// serial function: a one-member complex is a monomer of its type, a
// two-member complex is a dimer for each type it holds.
void print_dimers(std::vector<Complex>& complexList, std::ofstream& outfile,
                  int it, Parameters params,
                  std::vector<MolTemplate>& molTemplateList,
                  MpiContext& mpiContext)
{
    const int nTypes = params.numMolTypes;
    std::vector<int> monomers(nTypes, 0), dimers(nTypes, 0);
    for (const auto& com : complexList) {
        if (com.isEmpty || com.ownerRank != mpiContext.rank) continue;
        if (static_cast<int>(com.numEachMol.size()) < nTypes) continue;
        if (com.memberList.size() == 1) {
            for (int j = 0; j < nTypes; j++)
                if (com.numEachMol[j] == 1) monomers[j]++;
        } else if (com.memberList.size() == 2) {
            for (int j = 0; j < nTypes; j++)
                if (com.numEachMol[j] == 1 || com.numEachMol[j] == 2) dimers[j]++;
        }
    }
    outfile << (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom << '\t';
    for (int j = 0; j < nTypes; j++)
        outfile << monomers[j] << '\t' << dimers[j] << '\t';
    outfile << std::endl;
}

// init_NboundPairs() counts a restart's bonds from every molecule on the
// rank, ghost copies included, so a bond of a complex straddling the boundary
// was counted on both of its ranks.  Count it again from the molecules this
// rank owns (both ends of a bond are in one complex, so both are owned here
// or neither); the rule is init_NboundPairs()'s.
void recount_owned_bound_pairs(copyCounters& counterArrays, const Parameters& params,
                               const std::vector<Molecule>& moleculeList)
{
    for (auto& n : counterArrays.nBoundPairs) n = 0;
    for (const auto& mol : moleculeList) {
        if (mol.isEmpty || mol.isImplicitLipid || mol.isGhosted) continue;
        for (int partner : mol.bndpartner) {
            if (mol.index < partner) continue;
            if (partner < 0 || partner >= static_cast<int>(moleculeList.size())) continue;
            int p1 = mol.molTypeIndex;
            int p2 = moleculeList[partner].molTypeIndex;
            if (p2 < p1) std::swap(p1, p2);
            const int index = p1 * params.numMolTypes + p2;
            if (index >= 0 && index < static_cast<int>(counterArrays.nBoundPairs.size()))
                counterArrays.nBoundPairs[index]++;
        }
    }
}
