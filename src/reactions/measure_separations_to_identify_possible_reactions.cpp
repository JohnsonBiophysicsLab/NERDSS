/*! \file measure_separations_to_identify_possible_reactions.cpp
 *
 * The pairwise search of one time step, shared by the serial build
 * (EXEs/nerdss.cpp) and the MPI build (EXEs/nerdss_mpi.cpp).  Until the two
 * were joined the serial build carried its own copy of this walk inline, and
 * the gates added to it since 1.3 -- the type mask of a (type, type) pair, the
 * type mask of a whole SubBox, the spherical shell index and the compartment
 * check -- never reached this file, so an MPI rank offered every pair its
 * cell stencil produced and spent about twice the serial build's time per
 * step here (mem_localization at np=1: 2.98x the serial time, with 99.6% of
 * the offered pairs lipid-lipid).
 *
 * Every gate skips only pairs check_bimolecular_reactions() would have
 * rejected without side effects, so the pairs that remain are offered in the
 * same order as without the gates and the trajectory does not change.  The
 * ghost rules of the MPI decomposition stay under mpi_, where the serial
 * build never sees them.
 */
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "boundary_conditions/reflect_functions.hpp"
#include "debug/debug.hpp"
#include "error/error.hpp"
#include "io/io.hpp"
#include "macro.hpp"
#include "math/constants.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#ifdef mpi_
#include "mpi.h"
#endif
#include "mpi/mpi_function.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "mpi_context.hpp"
#include "system_setup/system_setup.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// Bit t2 of interactionMask[t1] is set iff a (t1, t2) pair can do anything
// inside check_bimolecular_reactions().  The pairwise search offers every
// pair the cell stencil produces regardless of what the two molecules are,
// and on a membrane model that is nearly all of them: mem_localization holds
// 3755 lipids among 3955 molecules, offers about 165 000 pairs per step, and
// 699 of them -- 0.42% -- get past the molecule-type test inside the
// function.  The other 99.6% are lipid-lipid pairs that no reaction names.
//
// Skipping them in the walk rather than inside the function removes only
// pairs that the function would have rejected without side effects, so the
// pairs that remain are offered in the same order as before and the
// trajectory does not change.
//
// One bit per molecule type, so the table only exists for models with at
// most 64 of them; left empty, the search offers every pair as it did before.
// Upper bound on a SubBox's neighborList: a cubic SubBox looks at 13 forward
// neighbours (SimulVolume::maxNeighbors); 32 leaves room without a heap
// buffer per molecule.
static constexpr int kMaxNeighbors{32};

static std::vector<uint64_t> build_interaction_mask(
    const std::vector<MolTemplate>& molTemplateList) {
  std::vector<uint64_t> interactionMask;
  // Past 64 types the bits run out.  Every mask is then all ones, so the
  // test below passes every pair, as the walk did before the table existed:
  // the same answer as an empty table, without a flag to test per partner.
  if (molTemplateList.size() > 64) {
    interactionMask.assign(molTemplateList.size(), ~uint64_t(0));
    return interactionMask;
  }
  interactionMask.assign(molTemplateList.size(), 0);
  for (size_t typeItr1{0}; typeItr1 < molTemplateList.size(); ++typeItr1) {
    for (size_t typeItr2{0}; typeItr2 < molTemplateList.size(); ++typeItr2) {
      // check_bimolecular_reactions() looks for the partner's type in the
      // target's rxnPartners, in that order, so this table is not symmetric.
      bool pairMatters{false};
      for (int partner : molTemplateList[typeItr1].rxnPartners) {
        if (partner == static_cast<int>(typeItr2)) {
          pairMatters = true;
          break;
        }
      }
      // Its volume-exclusion path does not consult rxnPartners and turns on
      // an interface being bound, which is runtime state.  Testing the
      // template flag instead keeps every pair that path could ever want.
      if (molTemplateList[typeItr1].excludeVolumeBound ||
          molTemplateList[typeItr2].excludeVolumeBound)
        pairMatters = true;
      if (pairMatters) interactionMask[typeItr1] |= (uint64_t(1) << typeItr2);
    }
  }
  return interactionMask;
}

// Whether any two explicit molecules can react with each other at all.  The
// shell index stays inactive otherwise, since it would then hold only pairs
// the search discards anyway.
static bool any_explicit_pair_reaction(
    const std::vector<MolTemplate>& molTemplateList) {
  for (size_t typeItr1{0}; typeItr1 < molTemplateList.size(); ++typeItr1) {
    if (molTemplateList[typeItr1].isImplicitLipid) continue;
    for (int partner : molTemplateList[typeItr1].rxnPartners) {
      if (partner < 0 || size_t(partner) >= molTemplateList.size()) continue;
      if (molTemplateList[partner].isImplicitLipid) continue;
      return true;
    }
  }
  return false;
}

void PairSearchIndex::build(const std::vector<MolTemplate>& molTemplateList,
                            const Membrane& membraneObject,
                            const Parameters& params) {
  interactionMask = build_interaction_mask(molTemplateList);

  // A second neighbour index, over the membrane of a spherical system only.
  // SimulVolume's grid covers the (2R)^3 box that bounds the sphere, so for
  // molecules pinned to the shell it spends its resolution on interior volume
  // that is empty by construction, and it bins on Cartesian coordinates while
  // get_distance() measures those pairs by geodesic arc.  This index covers
  // surface-to-surface pairs; molecules stay in the Cartesian grid as well, so
  // surface-to-interior pairs are untouched.  It stays inactive for every
  // non-spherical system, which then behaves exactly as before.
  shellIndex.build(membraneObject.isSphere() ? membraneObject.sphereR : 0.0,
                   params.rMaxLimit, any_explicit_pair_reaction(molTemplateList));
  shellIndex.display();
}

void measure_separations_to_identify_possible_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<TransmissionRxn>& transmissionRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<double>& IL2DbindingVec,
    std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs,
    std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, int implicitlipidIndex,
    double* tableIDs, unsigned& DDTableIndex, PairSearchIndex& pairSearch) {
  // Only the SubBoxes that hold molecules are visited.  Walking all of
  // subCellList cost time proportional to the SubBox count however few
  // molecules the system held, and that count is typically far larger than
  // the number of occupied SubBoxes: rev_3D holds 2 000 molecules in 27 000
  // SubBoxes, of which about 1 800 are occupied, and clathrin fills 97 of
  // 2 744.  Sorted, occupiedSubCells visits exactly the same SubBoxes in
  // exactly the same order as the old sweep, so the output is unchanged.
  //
  // Both indices are refreshed here, after the dissociation and creation
  // blocks of the step, so that they describe the molecule set the search is
  // about to walk.  The caller's update_memberMolLists() has just rebuilt
  // every member list through add_member(), which is what keeps
  // SubBox::typeMask, memberTypeIndex and occupancyMask in step with the
  // lists; the MPI-only sites that push into memberMolList directly (the
  // deserializers, prepare.cpp) or renumber it (remove_empty_slots()) all run
  // after this walk and before the next rebuild, so none of them is ever
  // read here.
  simulVolume.refresh_occupied_cells();
  ShellIndex& shellIndex = pairSearch.shellIndex;
  shellIndex.rebin(moleculeList, complexList);
  const bool shellActive{shellIndex.active};
  const std::vector<uint64_t>& interactionMask{pairSearch.interactionMask};
  // Does the partner's type pass the target's mask.  The shift count is
  // taken modulo 64 so that a model with more than 64 types, whose masks
  // are all ones, is tested with the same code: the old walk tested a flag
  // for that case before every pair, and clang could keep the flag in a
  // register only while the walk was inline in main().
  auto typePasses = [](uint64_t targMask, int partType) {
    return (targMask >> (partType & 63)) & 1;
  };
  // Partner types for the type-mask tests below, without touching moleculeList.
  const std::vector<int>& memberTypeIndex{simulVolume.memberTypeIndex};
#ifdef mpi_
  const int columns{simulVolume.numSubCells.x};
#endif

  for (int cellItr : simulVolume.occupiedSubCells) {
    for (unsigned memItr{0};
         memItr < simulVolume.subCellList[cellItr].memberMolList.size();
         ++memItr) {
      int targMolIndex{simulVolume.subCellList[cellItr].memberMolList[memItr]};
      if (moleculeList[targMolIndex].isImplicitLipid) continue;
      // A ghost is still a partner.  The walk below is a half shell -- the
      // same cell's later members and 13 forward neighbour cells -- so each
      // pair is met from one side only, and for a pair of a ghost in the left
      // ghost column and an owned molecule in the first owned column that
      // side is usually the ghost's: (x+1, *, *) is forward, (x-1, *, *) only
      // with y+1 or z+1.  Skipping ghosts here therefore lost most pairs
      // across every rank boundary -- the left neighbour drops them as its
      // right ghost's (drop_crossings_decided_elsewhere()), so nobody drew
      // them -- and assembly near the boundaries ran slow.
      //
      // Two ghosts pair too, unless both sit in the left ghost column: that
      // column is the left neighbour's last owned one, where it lists and
      // keeps the pair.  Anywhere else a pair of ghosts is nobody's but this
      // rank's.  The ghosts in this rank's first owned columns are members
      // of the left neighbour's complexes straddling the boundary, and on
      // that neighbour they sit in its right ghost column, where it drops
      // every pair (drop_crossings_decided_elsewhere()); skipping every
      // ghost-ghost pair here, as this used to, left such a pair with no
      // initiator anywhere.  A pair of ghosts in the right ghost column is
      // listed and then dropped as the right neighbour's, like any other
      // pair there.  The serial build has no ghosts and never reaches the
      // test.
#ifdef mpi_
      const bool targGhost{moleculeList[targMolIndex].isGhosted};
      auto ghostsLeftNeighboursPair = [&](int partMolIndex) {
        return targGhost && moleculeList[partMolIndex].isGhosted &&
               moleculeList[targMolIndex].mySubVolIndex % columns == 0 &&
               moleculeList[partMolIndex].mySubVolIndex % columns == 0;
      };
#else
      constexpr bool targGhost{false};
      auto ghostsLeftNeighboursPair = [](int) { return false; };
#endif

      // Test bimolecular reactions, and binding to implicit-lipids
      if (moleculeList[targMolIndex].freelist.size() > 0 ||
          molTemplateList[moleculeList[targMolIndex].molTypeIndex]
                  .excludeVolumeBound == true) {
        // first, check for implicit-lipid binding
        int protype = moleculeList[targMolIndex].molTypeIndex;
        const uint64_t targMask{interactionMask[protype]};
        // Pairs the shell index owns are skipped here and offered there.
        const bool targOnShell{shellActive &&
                               shellIndex.isShellBinned[targMolIndex] != 0};
        if (!targGhost && molTemplateList[protype].bindToSurface == true) {
          check_implicit_reactions(
              targMolIndex, implicitlipidIndex, simItr, params, moleculeList,
              complexList, molTemplateList, forwardRxns, backRxns,
              counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
              ILTableIDs);
        }
        // Check for binding to the compartment in / or out. there must be a
        // reaction type that involves binding/unbinding from the compartment
        // surface. Implementation might be very similar to implicit lipid.
        // (nerdss_mpi refuses compartment models and its ranks never carry
        // hasCompartment, so this is the serial build's alone.)
        if (membraneObject.hasCompartment == true) {
          // only consider single molecule
          if (complexList[moleculeList[targMolIndex].myComIndex]
                  .memberList.size() == 1) {
            check_compartment_reaction(
                targMolIndex, -1, simItr, params, moleculeList, complexList,
                molTemplateList, transmissionRxns, backRxns, counterArrays,
                membraneObject, IL2DbindingVec, IL2DUnbindingVec, ILTableIDs);
          }
        }

        // secondly, loop over proteins in your same cell.
        for (unsigned memItr2{memItr + 1};
             memItr2 < simulVolume.subCellList[cellItr].memberMolList.size();
             ++memItr2) {
          int partMolIndex{
              simulVolume.subCellList[cellItr].memberMolList[memItr2]};
          if (targOnShell && shellIndex.isShellBinned[partMolIndex])
            continue;  // the shell index offers this pair
          if (!typePasses(targMask, memberTypeIndex[partMolIndex])) continue;
          if (ghostsLeftNeighboursPair(partMolIndex)) continue;

          check_bimolecular_reactions(
              targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex,
              params, normMatrices, survMatrices, pirMatrices, moleculeList,
              complexList, molTemplateList, forwardRxns, backRxns,
              counterArrays, membraneObject);
        }  // loop over protein partners in your same cell

        // thirdly, loop over all neighboring cells, and all proteins in those
        // cells. for PBC, all cells have maxnbor neighbor cells. For
        // reflecting, edge have fewer.
        // Which neighbouring SubBoxes can hold a partner of this molecule.
        // Scanned in a loop of its own, with no call in it: when the scan
        // shared a loop with the check_bimolecular_reactions() call, the
        // call's argument references kept the callee-saved registers and
        // clang spilled targMask, so each of the 13 mask tests per molecule
        // paid stack loads.  In the serial build that cost rev_3D 1.7% more
        // instructions than the old inline walk; this scan gives them back.
        // The accepted SubBoxes are visited in neighborList order, so the
        // pairs come out in the same order as before.
        const std::vector<int>& neighborList{
            simulVolume.subCellList[cellItr].neighborList};
        int liveCells[kMaxNeighbors];
        int numLiveCells{0};
        if (neighborList.size() > kMaxNeighbors)
          error("a SubBox has more neighbours than the pair walk allows");
        for (int neighCellItr : neighborList) {
          if (DEBUG && (neighCellItr > simulVolume.subCellList.size()))
            error("neighCellItr is higher than subCellList size!");
          if ((simulVolume.subCellList[neighCellItr].typeMask & targMask) != 0)
            liveCells[numLiveCells++] = neighCellItr;
        }
        for (int liveItr{0}; liveItr < numLiveCells; ++liveItr) {
          const int neighCellItr{liveCells[liveItr]};
          for (unsigned memItr2{0};
               memItr2 <
               simulVolume.subCellList[neighCellItr].memberMolList.size();
               ++memItr2) {
            int partMolIndex{
                simulVolume.subCellList[neighCellItr].memberMolList[memItr2]};
            if (DEBUG && (partMolIndex > moleculeList.size()))
              error("partMolIndex is higher than moleculeList size!");
            if (targOnShell && shellIndex.isShellBinned[partMolIndex])
              continue;  // the shell index offers this pair
            if (!typePasses(targMask, memberTypeIndex[partMolIndex])) continue;
            if (ghostsLeftNeighboursPair(partMolIndex)) continue;
            check_bimolecular_reactions(
                targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex,
                params, normMatrices, survMatrices, pirMatrices, moleculeList,
                complexList, molTemplateList, forwardRxns, backRxns,
                counterArrays, membraneObject);
          }  // loop over all proteins in this neighbor cell
        }    // loop over all neighbor cells
      }      // if protein i is free to bind
    }        // loop over all proteins in initial cell
  }          // End looping over all cells.

  // Second pass: the surface-to-surface pairs the Cartesian pass skipped.
  // Only check_bimolecular_reactions() runs here.  The implicit-lipid and
  // compartment checks are per-molecule, not per-pair, and already ran above.
  if (shellActive) {
    for (int shellCellItr : shellIndex.occupiedCells) {
      const std::vector<int>& shellMembers =
          shellIndex.memberMolList[shellCellItr];
      for (unsigned memItr{0}; memItr < shellMembers.size(); ++memItr) {
        const int targMolIndex{shellMembers[memItr]};
        if (moleculeList[targMolIndex].isImplicitLipid) continue;
        if (!(moleculeList[targMolIndex].freelist.size() > 0 ||
              molTemplateList[moleculeList[targMolIndex].molTypeIndex]
                      .excludeVolumeBound == true))
          continue;

        const int protype = moleculeList[targMolIndex].molTypeIndex;
        const uint64_t targMask{interactionMask[protype]};
#ifdef mpi_
        const bool targGhost{moleculeList[targMolIndex].isGhosted};
        auto ghostsLeftNeighboursPair = [&](int partMolIndex) {
          return targGhost && moleculeList[partMolIndex].isGhosted &&
                 moleculeList[targMolIndex].mySubVolIndex % columns == 0 &&
                 moleculeList[partMolIndex].mySubVolIndex % columns == 0;
        };
#else
        auto ghostsLeftNeighboursPair = [](int) { return false; };
#endif

        for (unsigned memItr2{memItr + 1}; memItr2 < shellMembers.size();
             ++memItr2) {
          const int partMolIndex{shellMembers[memItr2]};
          if (!typePasses(targMask, memberTypeIndex[partMolIndex])) continue;
          if (ghostsLeftNeighboursPair(partMolIndex)) continue;
          check_bimolecular_reactions(
              targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex,
              params, normMatrices, survMatrices, pirMatrices, moleculeList,
              complexList, molTemplateList, forwardRxns, backRxns,
              counterArrays, membraneObject);
        }

        for (int neighCellItr : shellIndex.neighborList[shellCellItr]) {
          if ((shellIndex.typeMask[neighCellItr] & targMask) == 0) continue;
          const std::vector<int>& neighMembers =
              shellIndex.memberMolList[neighCellItr];
          for (unsigned memItr2{0}; memItr2 < neighMembers.size(); ++memItr2) {
            const int partMolIndex{neighMembers[memItr2]};
            if (!typePasses(targMask, memberTypeIndex[partMolIndex])) continue;
            if (ghostsLeftNeighboursPair(partMolIndex)) continue;
            check_bimolecular_reactions(
                targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex,
                params, normMatrices, survMatrices, pirMatrices, moleculeList,
                complexList, molTemplateList, forwardRxns, backRxns,
                counterArrays, membraneObject);
          }
        }
      }
    }
  }
}
