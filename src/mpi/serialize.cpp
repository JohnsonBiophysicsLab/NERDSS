#include <iostream>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

void serialize_complexes(MpiContext &mpiContext, SimulVolume &simulVolume,
                         Membrane &membraneObject, const vector<int> &complexes,
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

  // This used to set deleteIfNotReceivedBack on every complex of the rank, a
  // pass over the whole list per exchange.  The flag's only reader was
  // delete_disappeared_complexes(), which nothing called and which has since
  // been removed, and nothing ever clears it, so the pass never changed a
  // value.

  // The members go on the wire as ids, which the receiver maps back to its
  // own indices.  Written straight from the complex (Complex::serialize_as())
  // rather than from a copy with its memberList rewritten; the one id vector
  // is reused across the loop.
  static vector<int> memberIds;
  for (int comIndex : complexes) {
    Complex &com = complexList[comIndex];
    if (com.isEmpty) continue;
    memberIds.clear();
    for (auto &it : com.memberList) memberIds.push_back(moleculeList[it].id);

    // Serialize the complex into the arrayRank:
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
                         int &nArrayRank, const vector<int> &molecules) {
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

  for (int molIndex : molecules) {
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
