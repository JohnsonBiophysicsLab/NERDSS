#include <iostream>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

void serialize_complexes(MpiContext &mpiContext, SimulVolume &simulVolume,
                         Membrane &membraneObject, set<int> &complexesSet,
                         vector<Complex> &complexList,
                         vector<Molecule> &moleculeList,
                         unsigned char *arrayRank, int &nArrayRank) {
  if (VERBOSE) cout << "serialize_complexes begins" << endl;
  // The count has to be the number of complexes actually written, not the size
  // of the set, because the loop below skips empty ones.  Reserve the slot now
  // and fill it in once the count is known, the way serialize_molecules() does:
  // the receiver reads exactly this many complexes off the stream, so counting
  // one too many makes it read into the bytes that follow and desynchronizes
  // the rest of the message.
  size_t nComplexes = 0;
  int writeLengthAt = nArrayRank;
  nArrayRank += sizeof(nComplexes);

  // TODO: AssumING 2 molecule complexes:
  for (auto &com : complexList) {
    com.deleteIfNotReceivedBack = true;
  }

  // The members go on the wire as ids, which the receiver maps back to its
  // own indices.  Written straight from the complex (Complex::serialize_as())
  // rather than from a copy with its memberList rewritten; the one id vector
  // is reused across the loop.
  static vector<int> memberIds;
  // Loop over set elements:
  for (auto itr = complexesSet.begin(); itr != complexesSet.end(); itr++) {
    Complex &com = complexList[*itr];
    if (com.isEmpty) continue;
    memberIds.clear();
    for (auto &it : com.memberList) memberIds.push_back(moleculeList[it].id);

    // Serialize the complex *itr into the arrayRank:
    com.serialize_as(arrayRank, nArrayRank, memberIds);

    com.receivedFromNeighborRank = false;
    nComplexes++;
  }

  // Serialize the number of complexes before the serialized complexes:
  *((size_t *)&(arrayRank[writeLengthAt])) = nComplexes;
  // debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
  // "//serialize_complexes()");
  if (VERBOSE) cout << "serialize_complexes ends" << endl;
}

// TODO: nenadko
//  update interfaces (partnerIndex) in serialize_molecules as molecules flow
//  into the rank. Update also on deserialize part
// TODO:ASK: DO WE NEED TO TRANSFER CELLS AS WELL?
void serialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume,
                         vector<Molecule> &moleculeList,
                         vector<Complex> &complexList, unsigned char *arrayRank,
                         int &nArrayRank, set<int> &complexesSet,
                          set<int> &moleculesSet) {
  if (VERBOSE) cout << "serialize_molecules begins" << endl;
  // Molecules are serialized after the number of molecules,
  // making it possible to deserialize.
  // Store index in arrayRank where number of molecules will be serialized once
  // known:
  int writeLengthAt = nArrayRank;
  // Start serializing molecules right after the position for the number of
  // molecules:
  nArrayRank += sizeof(int);
  // Count molecules:
  int nMols = 0;

  for (auto &molIndex : moleculesSet) {
    Molecule &mol = moleculeList[molIndex];
    if (mol.isImplicitLipid) continue;
    if (mol.isEmpty) continue;
    
    // The receiver finds partners and the complex by id: bring the partner
    // ids up to date and write the complex's id where myComIndex goes,
    // without touching the molecule's own myComIndex.
    refresh_partner_ids(mol, moleculeList);
    mol.serialize_as(arrayRank, nArrayRank, complexList[mol.myComIndex].id);
    mol.receivedFromNeighborRank = false;
    nMols++;
  }

  // Serialize number of molecules before serialized molecules:
  *((int *)&(arrayRank[writeLengthAt])) = nMols;
  if (VERBOSE) cout << "Total serialized number of molecules:" << nMols << endl;
}

void serialize_molecules_partial(MpiContext &mpiContext,
                                 SimulVolume &simulVolume, vector<int> &indices,
                                 vector<Molecule> &moleculeList,
                                 vector<Complex> &complexList,
                                 unsigned char *arrayRank, int &nArrayRank,
                                 set<int> &complexesSet, int cellX) {
  if (VERBOSE) cout << "serialize_molecules begins" << endl;
  // Molecules are serialized after the number of molecules,
  // making it possible to deserialize.
  // Store index in arrayRank where number of molecules will be serialized once
  // known:
  int writeLengthAt = nArrayRank;
  // Start serializing molecules right after the position for the number of
  // molecules:
  nArrayRank += sizeof(int);
  // Count molecules:
  int nMols = 0;
  // Iterate over all cells, and determine whether their molecules have to be
  // transfered
  // TODO: transfer also complexes which occupy these cells:
  for (unsigned cellItr{0}; cellItr < simulVolume.subCellList.size();
       ++cellItr) {  // for all cells
    auto &cell = simulVolume.subCellList[cellItr];
    // cout << "Testing whether cellItr = " << cellItr << " should be
    // serialized" << endl;
    if (cell.xIndex - mpiContext.xOffset == cellX) {
      // cout << "Looping over all molecules from cell: cellItr = " << cellItr
      // << endl;
      for (unsigned memItr{0}; memItr < cell.memberMolList.size(); ++memItr) {
        // Find molecule index:
        unsigned molIndex = cell.memberMolList[memItr];
        if (DEBUG && (molIndex >= moleculeList.size()))
          error("1: molIndex >= moleculeList.size()");
        // Find molecule:
        Molecule &mol = moleculeList[molIndex];
        if (moleculeList[molIndex].myComIndex == -1 ||
            moleculeList[molIndex].isImplicitLipid ==
                true)  // molecule to be deleted or implicit lipid
          continue;
        // cout << "Serializing molecule, index = " << molIndex << endl;
        DEBUG_MOL("SEND " + to_string(cellX));
        // Mark myComIndex complex to be serialized as well:
        complexesSet.insert(mol.myComIndex);

        // Based on interfaces' partnerIndex, assign partnerIDs, and write the
        // complex's id where myComIndex goes:
        refresh_partner_ids(mol, moleculeList);
        mol.serialize_as(arrayRank, nArrayRank, complexList[mol.myComIndex].id);

        // Before sending molecules to the neighbor rank,
        // mark each molecule as not received back from the other rank,
        // so that once molecules from the same zones are received back from
        // other rank, those molecules that were sent, but are not received
        // back, can be deleted:
        mol.receivedFromNeighborRank = false;
        nMols++;
        // TODO: it is a complex => determine whether the complex is broken, and
        // introduce COM rank for the complex.
      }
    }
  }
  // Serialize number of molecules before serialized molecules:
  *((int *)&(arrayRank[writeLengthAt])) = nMols;
  if (VERBOSE) cout << "Total serialized number of molecules:" << nMols << endl;
}
