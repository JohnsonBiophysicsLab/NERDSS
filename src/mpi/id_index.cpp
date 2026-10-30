#include <cstdio>
#include <iostream>
#include <unordered_map>

#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

// Id-to-index maps for the two lists a rank holds.
//
// find_molecule() and find_complex() used to be linear scans over
// moleculeList and complexList, called once per received molecule and once
// per bound partner of each (deserialize.cpp, IDs_to_indices()): about 7% of
// a rank's samples at np=2 on hetTrimer, and quadratic in the stripe
// population.  The maps below answer the same question in O(1).
//
// What the scans answered, exactly: the lowest index whose id field equals
// the id asked for, whether or not that slot is empty.  The distinction
// matters: Molecule::MPI_remove_from_one_rank() leaves id set on the slot it
// empties, and deserialize_molecules() relies on finding that slot again so a
// molecule that comes back in the same step is written into it rather than
// appended twice; deserialize_complexes() finds it too and then drops the
// member because it is empty.  Molecule::destroy() and Complex::destroy() set
// id to -1, so those slots never match.  A hit is therefore validated by
// `list[index].id == id` alone, and the callers keep their own isEmpty tests.
//
// How the maps are kept current:
//  * rebuild_id_index() after prepare() has built the rank's lists and after
//    remove_empty_slots() has compacted them (every MPI step, in the function
//    that moves the slots).
//  * note_molecule_installed() / note_complex_installed() from the deserialisers
//    when they append a slot.  In-place overwrites keep the id and the index.
//  * refresh_id_index() at the top of each receive, for what the reactions
//    created on this rank since the last rebuild.  Ids are handed out from
//    Molecule::maxID / Complex::maxID, which only ever grow, and prepare()
//    gives every rank its own id range, so the ids created here since the last
//    pass are exactly those in [watermark, maxID): one integer compare per
//    slot finds them, whatever slot the creation took.  Nothing else installs
//    an id below the watermark except the deserialisers, which report it.
// A stale entry -- a slot destroyed, or recycled under a new id -- fails the
// validation, which is what the scan would have reported too.  Under
// CHECK_STRUCTURE every lookup is also checked against the scan.

namespace {

struct IdIndex {
  std::unordered_map<int, int> index;
  int watermark{0};  // ids below this were in the map at the last pass
};

IdIndex molIds;
IdIndex comIds;

template <typename T>
void rebuild(IdIndex &ids, const vector<T> &list, int maxID) {
  ids.index.clear();
  ids.index.reserve(list.size() * 2);
  for (size_t i = 0; i < list.size(); i++) {
    if (list[i].id < 0) continue;
    ids.index.emplace(list[i].id, static_cast<int>(i));  // keeps the lowest
  }
  ids.watermark = maxID;
}

template <typename T>
void refresh(IdIndex &ids, const vector<T> &list, int maxID) {
  if (maxID > ids.watermark) {
    for (size_t i = 0; i < list.size(); i++) {
      const int id = list[i].id;
      if (id >= ids.watermark && id < maxID)
        ids.index.emplace(id, static_cast<int>(i));
    }
    ids.watermark = maxID;
  }
}

template <typename T>
int scan(const vector<T> &list, int id) {
  size_t i = 0;
  while ((i < list.size()) && (list[i].id != id)) i++;
  if (i == list.size()) return -1;
  return static_cast<int>(i);
}

template <typename T>
int lookup(const IdIndex &ids, const vector<T> &list, int id) {
  // Destroyed slots carry id -1 and are not in the map; a caller asking for
  // a negative id gets the scan's answer, as before.  Nothing does in a run
  // whose invariants hold.
  if (id < 0) return scan(list, id);
  auto it = ids.index.find(id);
  if (it == ids.index.end()) return -1;
  const int index = it->second;
  if (index < 0 || index >= static_cast<int>(list.size())) return -1;
  return list[index].id == id ? index : -1;
}

}  // namespace

void rebuild_id_index(const vector<Molecule> &moleculeList,
                      const vector<Complex> &complexList) {
  rebuild(molIds, moleculeList, Molecule::maxID);
  rebuild(comIds, complexList, Complex::maxID);
}

void refresh_id_index(const vector<Molecule> &moleculeList,
                      const vector<Complex> &complexList) {
  refresh(molIds, moleculeList, Molecule::maxID);
  refresh(comIds, complexList, Complex::maxID);
}

void note_molecule_installed(int id, int index) { molIds.index[id] = index; }

void note_complex_installed(int id, int index) { comIds.index[id] = index; }

int find_complex(vector<Complex> &complexList, int id) {
  const int index = lookup(comIds, complexList, id);
  if (CHECK_STRUCTURE) {
    const int scanned = scan(complexList, id);
    if (scanned != index)
      fprintf(stderr,
              "STRUCTURE VIOLATION: complex id index gives %d for id %d, the "
              "scan %d\n",
              index, id, scanned);
  }
  return index;
}

int find_molecule(vector<Molecule> &moleculeList, int id) {
  const int index = lookup(molIds, moleculeList, id);
  if (CHECK_STRUCTURE) {
    const int scanned = scan(moleculeList, id);
    if (scanned != index)
      fprintf(stderr,
              "STRUCTURE VIOLATION: molecule id index gives %d for id %d, the "
              "scan %d\n",
              index, id, scanned);
  }
  return index;
}

void IDs_to_indices(MpiContext &mpiContext, vector<Molecule> &moleculeList,
                    vector<int> &indices) {
  if (VERBOSE) cout << "IDs_to_indices begins" << endl;
  for (auto &index : indices) {
    Molecule &mol = moleculeList[index];
    //        DEBUG_MOL("4A");
    // Update interfaces to store also partner-molecule ID and
    // re-populate bndpartner and bndlist:
    mol.bndpartner.clear();
    mol.bndlist.clear();
    for (size_t i = 0; i < mol.interfaceList.size(); i++) {
      Molecule::Interaction &inter = mol.interfaceList[i].interaction;
      // Update only interfaces whose partner is visible on this rank:
      if (inter.partnerId != -1) {  // if it can be seen
        // cout << "4A0 + " << inter.partnerId << endl;
        if (inter.partnerId == -3) {
          // partner is implicit lipid
          inter.partnerIndex = 0;
          inter.partnerId = 0;
          mol.bndpartner.push_back(inter.partnerIndex);
          mol.bndlist.push_back(i);
          continue;
        }
        inter.partnerIndex = find_molecule(moleculeList, inter.partnerId);
        if (inter.partnerIndex != -1) {  // if found local partner
                                         // cout << "4AA" << endl;
          //  Both have to be updated:
          mol.bndpartner.push_back(inter.partnerIndex);
          mol.bndlist.push_back(i);

          // Update interfaces of a partner as well:
          Molecule &partnerMol = moleculeList[inter.partnerIndex];
          for (size_t iPartner = 0; iPartner < partnerMol.interfaceList.size();
               iPartner++) {
            Molecule::Interaction &interPartner =
                partnerMol.interfaceList[iPartner].interaction;
            if (interPartner.partnerId ==
                mol.id) {  // if it is connected to mol
              interPartner.partnerIndex = mol.index;
              //  Add this mol as a partner only if it doesn't exist:
              size_t iBndpartner = 0;
              while ((iBndpartner < partnerMol.bndlist.size()) &&
                     (partnerMol.bndlist[iBndpartner] != iPartner))
                iBndpartner++;
              if (iBndpartner == partnerMol.bndpartner.size()) {
                partnerMol.bndpartner.push_back(mol.index);
                partnerMol.bndlist.push_back(iPartner);
              } else {
                partnerMol.bndpartner[iBndpartner] = mol.index;
              }
            }
          }
        }
      }
    }

    // The reweighting history is not mapped here: it is not on the wire, and
    // a molecule received over a local copy keeps that copy's history, which
    // already names local indices (deserialize_molecules*()).
  }
  if (CHECK_STRUCTURE)
    for (auto &index : indices)
      trace_molecule(mpiContext, "after IDs_to_indices", moleculeList, *mpiContext.complexList, moleculeList[index]);
  if (VERBOSE) cout << "IDs_to_indices ends" << endl;
}

void indices_to_IDs(Molecule &mol, vector<Molecule> &moleculeList,
                    vector<Complex> &complexList) {
  if (VERBOSE) cout << "indices_to_IDs begins" << endl;
  // Update interfaces to store also partner-molecule ID:
  for (size_t i = 0; i < mol.interfaceList.size(); i++) {
    Molecule::Interaction &inter = mol.interfaceList[i].interaction;
    if (inter.partnerIndex != -1) {  // if it can be seen
      if (moleculeList[inter.partnerIndex].isImplicitLipid == false) {
        // Update only interfaces whose partner is visible on this rank:
        inter.partnerId = moleculeList[inter.partnerIndex].id;
        if (VERBOSE) cout << "Assigned: " << inter.partnerId << ". ";
      } else {
        inter.partnerId = -3;  // its partner is implicit lipid
      }
    }
  }
  // Set myComIndex to the ID of the complex:
  mol.myComIndex = complexList[mol.myComIndex].id;
  // Map the reweighting history's partner indices to IDs.  An entry can be -1
  // when the neighbour it refers to has moved to another rank; moleculeList[-1]
  // is a read off the front of the vector, so pass the sentinel through
  // untouched.
  for (auto &prevEntry : mol.prevReweight) {
    prevEntry.partner =
        (prevEntry.partner >= 0 &&
         prevEntry.partner < static_cast<int>(moleculeList.size()))
            ? moleculeList[prevEntry.partner].id
            : -1;
  }
  if (VERBOSE) cout << "indices_to_IDs ends" << endl;
}
