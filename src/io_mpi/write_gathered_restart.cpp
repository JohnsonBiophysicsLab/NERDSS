#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <unordered_map>
#include <vector>

#include "io/io.hpp"
#include "macro.hpp"
#include "math/rand_gsl.hpp"
#include "mpi/mpi_function.hpp"

/* A checkpoint of a run on any number of ranks, written by rank 0 as one
 * restart file in the format the serial build writes, so that `nerdss -r`
 * and `nerdss_mpi -np 1 -r` continue it.
 *
 * Each rank holds part of the system, and holds it under its own numbering:
 * molecule and complex indices are positions in that rank's lists, and every
 * reference between them -- a bond's partnerIndex, myComIndex, a complex's
 * memberList, the reweighting history's partners -- uses them.  The restart
 * format is index-based too, so the file cannot be any one rank's lists.
 * Every rank therefore sends rank 0 what it owns (molecules that are not
 * ghosts, and the complexes they belong to) with every reference written as
 * an id, which is unique system-wide, and rank 0 numbers the lot once and
 * maps each id back to its new index.
 *
 * What a rank keeps as a tally of its own share is summed: the implicit
 * lipid's free sites per state, the counters and event histograms, the
 * transition matrices, and each observable's change since the start (every
 * rank starts from the whole system's count); the lists of cluster lifetimes
 * are concatenated.  What the
 * decomposition rescaled -- the water box's x bounds and volume, the membrane
 * area and site count -- is taken from the whole system as it was parsed,
 * before the ranks split it.  What depends on the numbering is recomputed:
 * bndpartner from the interfaces (in bndlist's order), each template's monomerList,
 * the dissociable pairs of bindPairList, the count of each molecule type, and
 * mySubVolIndex on the whole system's grid (the reader rebuilds the
 * sub-volumes from coordinates anyway).
 *
 * Rank 0's own copy of the implicit lipid, of its complex, of the implicit
 * lipid's 2D binding table and of the reactions and parameters stands for all
 * of them: they are the same on every rank, or (the binding table) built from
 * a density every rank shares.
 *
 * Collective: every rank must call it, at the same step.
 */

namespace {

// The reweighting history is not part of Molecule::serialize(), whose wire
// format the exchanges use and which they do not need it in; a checkpoint
// does.
void push_reweight(const std::vector<Molecule::ReweightEntry> &history,
                   std::vector<unsigned char> &buffer) {
  const int n = static_cast<int>(history.size());
  const unsigned char *p = reinterpret_cast<const unsigned char *>(&n);
  buffer.insert(buffer.end(), p, p + sizeof(n));
  for (const auto &entry : history) {
    p = reinterpret_cast<const unsigned char *>(&entry);
    buffer.insert(buffer.end(), p, p + sizeof(entry));
  }
}

std::vector<Molecule::ReweightEntry> pop_reweight(const unsigned char *data,
                                                  int &at) {
  int n = 0;
  std::memcpy(&n, data + at, sizeof(n));
  at += sizeof(n);
  std::vector<Molecule::ReweightEntry> history(n);
  for (int i = 0; i < n; i++) {
    std::memcpy(&history[i], data + at, sizeof(history[i]));
    at += sizeof(history[i]);
  }
  return history;
}

// Serialize one object into the growing buffer.  The PUSH-based serializers
// write through a raw pointer, so give them room and trim to what they used.
template <typename T>
void push_object(T &object, std::vector<unsigned char> &buffer,
                 std::size_t room) {
  const std::size_t start = buffer.size();
  buffer.resize(start + room);
  int n = 0;
  object.serialize(buffer.data() + start, n);
  if (static_cast<std::size_t>(n) > room) {
    fprintf(stderr, "write_gathered_restart: an object needed %d bytes, %zu "
                    "were reserved\n",
            n, room);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  buffer.resize(start + n);
}

template <typename T>
void push_value(const T &value, std::vector<unsigned char> &buffer) {
  const unsigned char *p = reinterpret_cast<const unsigned char *>(&value);
  buffer.insert(buffer.end(), p, p + sizeof(value));
}

template <typename T>
T pop_value(const unsigned char *data, int &at) {
  T value;
  std::memcpy(&value, data + at, sizeof(value));
  at += sizeof(value);
  return value;
}

// Room for one serialized object: everything fixed-size is far below the
// base, and each vector costs its length times its element plus a length.
std::size_t molecule_room(const Molecule &mol) {
  return 4096 + mol.interfaceList.size() * (sizeof(Molecule::Iface) + 64) +
         (mol.freelist.size() + mol.bndlist.size() + mol.bndpartner.size()) *
             sizeof(int) * 2;
}

std::size_t complex_room(const Complex &com) {
  return 4096 + com.memberList.size() * sizeof(int) * 2 +
         com.numEachMol.size() * sizeof(int) * 2 +
         com.lastNumberUpdateItrEachMol.size() * sizeof(long long) * 2 +
         com.obs.size() * sizeof(int) * 2;
}

// Element-wise sum over ranks onto rank 0.
template <typename T>
void reduce_sum(std::vector<T> &values, MPI_Datatype type) {
  std::vector<T> total(values.size(), T{});
  MPI_Reduce(values.data(), total.data(), static_cast<int>(values.size()),
             type, MPI_SUM, 0, MPI_COMM_WORLD);
  values.swap(total);
}

// Every rank's values laid end to end on rank 0, in rank order.
std::vector<double> gather_concatenated(const std::vector<double> &mine,
                                        MpiContext &mpiContext) {
  const int n = static_cast<int>(mine.size());
  std::vector<int> counts(mpiContext.nprocs, 0), displs(mpiContext.nprocs, 0);
  MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<double> all;
  if (mpiContext.rank == 0) {
    int total = 0;
    for (int r = 0; r < mpiContext.nprocs; r++) {
      displs[r] = total;
      total += counts[r];
    }
    all.resize(total);
  }
  MPI_Gatherv(mine.data(), n, MPI_DOUBLE, all.data(), counts.data(),
              displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return all;
}

}  // namespace

void write_gathered_restart(
    long long int simItr, const std::string &fileName, int rngWrite,
    Parameters &params,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    std::vector<MolTemplate> &molTemplateList,
    std::vector<ForwardRxn> &forwardRxns, std::vector<BackRxn> &backRxns,
    std::vector<CreateDestructRxn> &createDestructRxns,
    const std::map<std::string, int> &observablesList,
    const std::map<std::string, int> &observablesAtStart,
    Membrane &membraneObject, copyCounters &counterArrays,
    const Membrane &membraneWhole, const SimulVolume &simulVolumeWhole,
    MpiContext &mpiContext) {
  const bool hasIL = membraneObject.implicitLipid;
  const int ilMol = hasIL ? membraneObject.implicitlipidIndex : -1;
  const int ilCom =
      (hasIL && ilMol >= 0 && ilMol < static_cast<int>(moleculeList.size()))
          ? moleculeList[ilMol].myComIndex
          : -1;

  // ---- this rank's share, with every reference written as an id
  std::vector<unsigned char> mine;
  std::vector<int> ownedMols, ownedComs;
  for (const auto &mol : moleculeList)
    if (!mol.isEmpty && !mol.isImplicitLipid && !mol.isGhosted)
      ownedMols.push_back(mol.index);
  for (const auto &com : complexList)
    if (!com.isEmpty && com.index != ilCom && com.ownerRank == mpiContext.rank)
      ownedComs.push_back(com.index);

  push_value(static_cast<int>(ownedMols.size()), mine);
  for (int m : ownedMols) {
    Molecule copy = moleculeList[m];
    indices_to_IDs(copy, moleculeList, complexList);
    push_object(copy, mine, molecule_room(copy));
    // not in Molecule::serialize()
    push_value(copy.isPromoter, mine);
    push_value(copy.enforceCompartmentBC, mine);
    push_reweight(copy.prevReweight, mine);
  }
  push_value(static_cast<int>(ownedComs.size()), mine);
  for (int c : ownedComs) {
    Complex copy = complexList[c];
    for (auto &member : copy.memberList)
      member = (member >= 0 && member < static_cast<int>(moleculeList.size()))
                   ? moleculeList[member].id
                   : -1;
    push_object(copy, mine, complex_room(copy));
    push_value(copy.onFiber, mine);
  }

  // ---- gather the shares onto rank 0
  const int myBytes = static_cast<int>(mine.size());
  std::vector<int> bytes(mpiContext.nprocs, 0), displs(mpiContext.nprocs, 0);
  MPI_Gather(&myBytes, 1, MPI_INT, bytes.data(), 1, MPI_INT, 0,
             MPI_COMM_WORLD);
  std::vector<unsigned char> all;
  if (mpiContext.rank == 0) {
    long long total = 0;
    for (int r = 0; r < mpiContext.nprocs; r++) {
      displs[r] = static_cast<int>(total);
      total += bytes[r];
    }
    all.resize(static_cast<std::size_t>(total));
  }
  MPI_Gatherv(mine.data(), myBytes, MPI_CHAR, all.data(), bytes.data(),
              displs.data(), MPI_CHAR, 0, MPI_COMM_WORLD);

  // ---- tallies kept per rank
  std::vector<int> freeLipids = membraneObject.numberOfFreeLipidsEachState;
  reduce_sum(freeLipids, MPI_INT);
  std::vector<long long> counts{
      counterArrays.nLoops,           counterArrays.nCancelOverlapPartner,
      counterArrays.nCancelOverlapSystem, counterArrays.nCancelDisplace2D,
      counterArrays.nCancelDisplace3D, counterArrays.nCancelDisplace3Dto2D,
      counterArrays.nCancelSpanBox,   counterArrays.nAssocSuccess};
  reduce_sum(counts, MPI_LONG_LONG);
  std::vector<int> events3D = counterArrays.events3D;
  std::vector<int> events2D = counterArrays.events2D;
  std::vector<int> events3Dto2D = counterArrays.events3Dto2D;
  reduce_sum(events3D, MPI_INT);
  reduce_sum(events2D, MPI_INT);
  reduce_sum(events3Dto2D, MPI_INT);
  std::vector<int> observed;
  for (const auto &observable : observablesList)
    observed.push_back(observable.second);
  reduce_sum(observed, MPI_INT);
  std::vector<long long> transitions;
  for (const auto &oneTemp : molTemplateList) {
    if (!oneTemp.countTransition) continue;
    for (const auto &row : oneTemp.transitionMatrix)
      transitions.insert(transitions.end(), row.begin(), row.end());
  }
  reduce_sum(transitions, MPI_LONG_LONG);
  // lifeTime is not a tally but a list per cluster size, appended to as
  // clusters change size, so each rank's is as long as its own history: the
  // system's lists are the ranks' laid end to end.
  std::vector<std::vector<std::vector<double>>> lifeTimes;
  for (const auto &oneTemp : molTemplateList) {
    if (!oneTemp.countTransition) continue;
    lifeTimes.emplace_back();
    for (const auto &row : oneTemp.lifeTime)
      lifeTimes.back().push_back(gather_concatenated(row, mpiContext));
  }

  if (mpiContext.rank == 0) {
    // Molecule::deserialize() and Complex::deserialize() pop the sender's
    // numberOfMolecules and numberOfComplexes into the statics; keep this
    // rank's to put back.
    const int rankMolecules = Molecule::numberOfMolecules;
    const int rankComplexes = Complex::numberOfComplexes;

    // ---- one numbering for the whole system
    std::vector<Molecule> mols;
    std::vector<Complex> coms;
    std::vector<char> molFiber;  // isPromoter per new molecule
    std::vector<char> comFiber;
    if (hasIL) {
      Molecule lipid = moleculeList[ilMol];
      lipid.index = 0;
      lipid.myComIndex = 0;
      mols.push_back(lipid);
      Complex lipidCom = complexList[ilCom];
      lipidCom.index = 0;
      lipidCom.memberList = {0};
      coms.push_back(lipidCom);
    }
    std::unordered_map<int, int> molIndexOfId, comIndexOfId;
    for (int r = 0; r < mpiContext.nprocs; r++) {
      const unsigned char *data = all.data() + displs[r];
      int at = 0;
      const int nMols = pop_value<int>(data, at);
      for (int i = 0; i < nMols; i++) {
        Molecule mol;
        mol.deserialize(const_cast<unsigned char *>(data), at);
        mol.isPromoter = pop_value<bool>(data, at);
        mol.enforceCompartmentBC = pop_value<bool>(data, at);
        mol.prevReweight = pop_reweight(data, at);
        mol.index = static_cast<int>(mols.size());
        if (!molIndexOfId.emplace(mol.id, mol.index).second)
          fprintf(stderr, "write_gathered_restart: molecule id %d is owned by "
                          "more than one rank\n",
                  mol.id);
        mols.push_back(mol);
      }
      const int nComs = pop_value<int>(data, at);
      for (int i = 0; i < nComs; i++) {
        Complex com;
        com.deserialize(const_cast<unsigned char *>(data), at);
        com.onFiber = pop_value<bool>(data, at);
        com.index = static_cast<int>(coms.size());
        if (!comIndexOfId.emplace(com.id, com.index).second)
          fprintf(stderr, "write_gathered_restart: complex id %d is owned by "
                          "more than one rank\n",
                  com.id);
        coms.push_back(com);
      }
    }

    // ---- map every id back to a new index
    int unresolved = 0;
    auto molIndex = [&](int id) {
      auto it = molIndexOfId.find(id);
      if (it == molIndexOfId.end()) {
        unresolved++;
        return -1;
      }
      return it->second;
    };
    const Vec3D cell = simulVolumeWhole.subCellSize;
    for (auto &mol : mols) {
      if (mol.isImplicitLipid) continue;
      // On the wire myComIndex carries the complex id.
      auto com = comIndexOfId.find(mol.myComIndex);
      if (com == comIndexOfId.end()) {
        unresolved++;
        mol.myComIndex = -1;
      } else {
        mol.myComIndex = com->second;
      }
      for (auto &iface : mol.interfaceList) {
        auto &inter = iface.interaction;
        inter.partnerIndex =
            !iface.isBound ? -1
                           : (inter.partnerId == -3 ? 0 : molIndex(inter.partnerId));
      }
      // bndlist keeps the order the bonds were made in: the step-top
      // dissociation pass draws a molecule's bonds in that order, so a file
      // that reordered them would continue differently.  bndpartner is
      // parallel to it and names the partner each interface now has.
      std::vector<int> order;
      for (int iface : mol.bndlist)
        if (iface >= 0 && iface < static_cast<int>(mol.interfaceList.size()) &&
            mol.interfaceList[iface].isBound)
          order.push_back(iface);
      mol.bndlist = order;
      mol.bndpartner.clear();
      for (int iface : mol.bndlist)
        mol.bndpartner.push_back(mol.interfaceList[iface].interaction.partnerIndex);
      std::vector<Molecule::ReweightEntry> history;
      for (auto entry : mol.prevReweight) {
        auto it = molIndexOfId.find(entry.partner);
        if (it == molIndexOfId.end()) continue;  // the partner left the system
        entry.partner = it->second;
        history.push_back(entry);
      }
      mol.prevReweight.swap(history);
      // The whole system's grid, as update_mySubVolIndex() computes it on a
      // rank with no x offset.
      const int xItr = int((mol.comCoord.x + membraneWhole.waterBox.x / 2) / cell.x);
      const int yItr = int((mol.comCoord.y + membraneWhole.waterBox.y / 2) / cell.y);
      const int zItr =
          (membraneWhole.waterBox.z > 0)
              ? int(-(mol.comCoord.z + 1E-6 - membraneWhole.waterBox.z / 2.0) / cell.z)
              : 0;
      mol.mySubVolIndex = xItr + yItr * simulVolumeWhole.numSubCells.x +
                          zItr * simulVolumeWhole.numSubCells.x *
                              simulVolumeWhole.numSubCells.y;
    }
    for (auto &com : coms) {
      if (hasIL && com.index == 0) continue;
      for (auto &member : com.memberList) member = molIndex(member);
    }
    // Every link has to be symmetric once renumbered; anything else means the
    // ranks disagreed about ownership, and the file would not restart.
    int broken = 0;
    for (const auto &mol : mols) {
      if (mol.isImplicitLipid) continue;
      if (mol.myComIndex < 0) continue;
      const auto &members = coms[mol.myComIndex].memberList;
      if (std::find(members.begin(), members.end(), mol.index) == members.end())
        broken++;
      for (const auto &iface : mol.interfaceList) {
        if (!iface.isBound || iface.interaction.partnerIndex < 0) continue;
        const Molecule &partner = mols[iface.interaction.partnerIndex];
        if (partner.isImplicitLipid) continue;
        const auto &back =
            partner.interfaceList[iface.interaction.partnerIfaceIndex].interaction;
        if (back.partnerIndex != mol.index) broken++;
      }
    }
    if (unresolved || broken)
      fprintf(stderr, "write_gathered_restart: itr=%lld %d unresolved and %d "
                      "asymmetric references in the gathered system\n",
              simItr, unresolved, broken);

    // ---- the counts that follow from the numbering
    std::vector<MolTemplate> templates = molTemplateList;
    std::vector<int> numEachMolType(templates.size(), 0);
    for (auto &oneTemp : templates) oneTemp.monomerList.clear();
    for (const auto &mol : mols) {
      if (mol.molTypeIndex < 0) continue;
      numEachMolType[mol.molTypeIndex]++;
      // add_to_monomer_list()'s rule
      MolTemplate &oneTemp = templates[mol.molTypeIndex];
      if (!mol.isImplicitLipid && mol.bndpartner.empty() && oneTemp.canDestroy)
        oneTemp.monomerList.push_back(mol.index);
    }
    std::size_t t = 0, l = 0;
    for (auto &oneTemp : templates) {
      if (!oneTemp.countTransition) continue;
      for (auto &row : oneTemp.transitionMatrix)
        for (auto &v : row) v = transitions[t++];
      oneTemp.lifeTime = lifeTimes[l++];
    }

    Membrane membrane = membraneWhole;
    membrane.numberOfFreeLipidsEachState = freeLipids;
    membrane.IL2DbindingVec = membraneObject.IL2DbindingVec;
    membrane.ILTableIDs = membraneObject.ILTableIDs;

    copyCounters counters = counterArrays;
    // The serial build's population-based dissociation draws its pairs from
    // bindPairList by molecule index, so it has to name the new indices; and
    // nerdss_mpi does not keep it at all.  associate_box()'s rule: a bond of
    // a reaction whose product can dissociate is listed once, under the
    // product's absIfaceIndex, by the molecule on the product's first side --
    // check_for_unimolecular_reactions_population() follows that molecule's
    // productListNew[0].relIfaceIndex interface to the partner.  (Both sides of
    // a bond can carry the product's index, so the side has to be read from
    // the reaction, not from which molecule comes first.)
    for (auto &list : counters.bindPairList) list.clear();
    std::unordered_map<int, std::pair<int, int>> firstSide;  // product -> (type, iface)
    for (const auto &rxn : forwardRxns)
      if (rxn.rxnType == ReactionType::bimolecular && !rxn.productListNew.empty())
        firstSide[rxn.productListNew[0].absIfaceIndex] = {
            rxn.productListNew[0].molTypeIndex, rxn.productListNew[0].relIfaceIndex};
    for (const auto &mol : mols) {
      if (mol.isImplicitLipid) continue;
      for (const auto &iface : mol.interfaceList) {
        if (!iface.isBound) continue;
        const int product = iface.index;
        if (product < 0 || product >= static_cast<int>(counters.canDissociate.size()) ||
            product >= static_cast<int>(counters.bindPairList.size()) ||
            !counters.canDissociate[product])
          continue;
        auto side = firstSide.find(product);
        if (side == firstSide.end() || side->second.first != mol.molTypeIndex ||
            side->second.second != iface.relIndex)
          continue;
        // A symmetric bond -- both sides on the first-side interface -- has
        // two molecules that qualify; list it once, by the lower index.  A
        // molecule may be listed for one bond and its partner for another
        // (a chain of A(c)+A(q) bonds), so this is decided per bond.
        const int partner = iface.interaction.partnerIndex;
        const bool symmetric =
            partner >= 0 && mols[partner].molTypeIndex == side->second.first &&
            iface.interaction.partnerIfaceIndex == side->second.second;
        if (symmetric && partner < mol.index) continue;
        counters.bindPairList[product].push_back(mol.index);
      }
    }

    std::map<std::string, int> observables;
    std::size_t o = 0;
    for (const auto &observable : observablesList) {
      auto start = observablesAtStart.find(observable.first);
      const int startValue =
          (start == observablesAtStart.end()) ? 0 : start->second;
      observables[observable.first] =
          observed[o++] - (mpiContext.nprocs - 1) * startValue;
    }

    // The legacy format writes the sub-volumes' member lists.
    SimulVolume volume = simulVolumeWhole;
    for (auto &subCell : volume.subCellList) subCell.memberMolList.clear();
    for (const auto &mol : mols)
      if (!mol.isImplicitLipid && mol.mySubVolIndex >= 0 &&
          mol.mySubVolIndex < static_cast<int>(volume.subCellList.size()))
        volume.subCellList[mol.mySubVolIndex].memberMolList.push_back(mol.index);

    // The writer reads these from the classes' statics, which hold this
    // rank's values for the rest of the run: swap the system's in for the
    // write and put this rank's back.
    std::swap(MolTemplate::numEachMolType, numEachMolType);
    std::vector<int> rankEmptyMols, rankEmptyComs;
    std::swap(Molecule::emptyMolList, rankEmptyMols);
    std::swap(Complex::emptyComList, rankEmptyComs);
    Molecule::numberOfMolecules = static_cast<int>(mols.size());
    Complex::numberOfComplexes = static_cast<int>(coms.size());

    std::vector<TransmissionRxn> transmissionRxns{};
    std::ofstream restartFile{fileName, std::ios::out};
    write_restart(simItr, restartFile, params, volume, mols, coms, templates,
                  forwardRxns, backRxns, createDestructRxns, transmissionRxns,
                  observables, membrane, counters);
    restartFile.close();

    std::swap(MolTemplate::numEachMolType, numEachMolType);
    std::swap(Molecule::emptyMolList, rankEmptyMols);
    std::swap(Complex::emptyComList, rankEmptyComs);
    Molecule::numberOfMolecules = rankMolecules;
    Complex::numberOfComplexes = rankComplexes;

    // Rank 0's generator, under the names the serial build writes and reads:
    // 1 for DATA/rng_state, 2 for the checkpoint's RESTARTS/rng_state<itr>.
    if (rngWrite == 1) write_rng_state();
    if (rngWrite == 2) write_rng_state_simItr(simItr);
  }
}
