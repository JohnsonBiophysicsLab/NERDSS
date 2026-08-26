/*! \file class_SimulVolume.cpp
 * \ingroup SimulClasses
 * ### Created on 10/19/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "classes/class_SimulVolume.hpp"
#include "io/io.hpp"
#include "error/error.hpp"

#include <chrono>
#include <classes/class_SimulVolume.hpp>
#include <iostream>

/* SIMULBOX::SUBBOX */
// Member Functions
void SimulVolume::SubVolume::display() {
  std::cout << "SubVolume " << absIndex << '\n';
  std::cout << "\tRel. Indices: [" << xIndex << ", " << yIndex << ", " << zIndex
            << "]\n";
  std::cout << "\tMolecule Members:";
  for (auto &mol : memberMolList)
    std::cout << ' ' << mol;
  std::cout << "\n\tNeighbors (abs. index):";
  for (auto &cell : neighborList)
    std::cout << ' ' << cell;
  std::cout << std::endl;
}

/* SIMULBOX::DIMENSIONS */
// Constructors
/*! \brief Largest number of SubBoxes that still leaves every edge >= cellLength.
 *
 * floor(L / cellLength) = n implies L / n >= cellLength, so this is exactly the
 * finest grid whose SubBoxes still cover the interaction range.  A box shorter
 * than cellLength gets one SubBox rather than none: one SubBox has no
 * neighbours to miss, so the whole dimension is searched.
 */
static int cells_along(double boxLength, double cellLength) {
  if (boxLength <= 0.0 || cellLength <= 0.0)
    return 1;
  return std::max(1, int(floor(boxLength / cellLength)));
}

SimulVolume::Dimensions::Dimensions(const Parameters &params,
                                    const Membrane &membraneObject) {
  // Every dimension is now derived the same way.  z used to carry a
  // std::max(4, ...) floor, which forced four SubBoxes into a box of any
  // thickness: below waterBox.z = 4 * rMaxLimit that made the SubBoxes thinner
  // than the interaction range, and the +/-1 neighbour stencil then stopped
  // covering it.  Pairs two SubBoxes apart in z were never offered to
  // check_bimolecular_reactions() at all -- 7.5% of the pairs within rMaxLimit
  // were lost in a 30 nm box against a 16.7 nm range, and 0.41% at 50 nm.
  const double cellLength{params.rMaxLimit};
  x = cells_along(membraneObject.waterBox.x, cellLength);
  y = cells_along(membraneObject.waterBox.y, cellLength);
  z = cells_along(membraneObject.waterBox.z, cellLength);

  tot = x * y * z;
}

// Member Functions
void SimulVolume::Dimensions::check_dimensions(const Parameters &params,
                                               const Membrane &membraneObject) {
  // The "is any dimension too small" repair that used to open this function is
  // gone: Dimensions() now derives every count from floor(L / rMaxLimit), so
  // no edge can come in below rMaxLimit.  The repair could not have fixed one
  // anyway -- it re-applied the same std::max(4, ...) floor that produced the
  // too-thin SubBoxes in the first place.
  //
  // What is left is the one cap on the total SubBox count, so params and
  // membraneObject are no longer read here.

  #ifdef mpi_
    // std::min because a cap must only ever remove SubBoxes.  Written without
    // it, this assignment could raise a count -- a 200 x 200 x 1 grid exceeds
    // 27000 and would come out 200 x 11 x 11 -- and raising a count shortens
    // that edge, which is exactly what create_simulation_volume() now refuses
    // to run with.  Untested here: this build has no MPI validation.
    if (x * y * z > 27000) {
      const int perAxis{static_cast<int>(std::sqrt(27000 / x))};
      y = std::min(y, perAxis);
      z = std::min(z, perAxis);
    }

    tot = x * y * z;
  #else
    // One cap on the total, applied to all three axes by the same factor.
    //
    // It replaces a flat 30 SubBoxes per dimension, which took no account of
    // the interaction range: on the 940 nm rev_3D box against a 16.70 nm range
    // the grid wants 56 SubBoxes per side, got 30, and so ran with SubBoxes
    // 1.88x too wide.  Candidate pairs grow with the cube of the SubBox edge,
    // and the measured over-inclusion doubled with it, from about 5x to 11.8x;
    // on the rev_2D slab, where the grid wants 172 per side, it reached 52x.
    // That cap was affordable only because the search walked every SubBox, and
    // it no longer does.
    //
    // Memory is the only thing left that justifies a cap, so the cap is
    // expressed in SubBoxes.  One costs about 145 bytes -- 64 for the struct
    // and the rest for the 13-entry neighbour list -- so this is roughly 72 MB
    // of grid.  No sample input comes close: the largest, rev_3Dto2D, wants
    // 59^3 = 205 379.
    //
    // Shrinking all three axes by the same cube-root factor keeps the SubBoxes
    // cubic, which the 13-neighbour stencil assumes, and can only lower a
    // count, so it cannot bring an edge back below rMaxLimit.
    const long long maxTotalCells{500000};
    const long long requested{1LL * x * y * z};
    long long current{requested};
    // Iterated, because one pass does not always land under the budget.  The
    // cube-root factor assumes all three axes absorb their share, and an axis
    // already down to one SubBox cannot: a 1 x 1 x 1197604 grid comes out of a
    // single pass at 1 x 1 x 895000, still nearly twice the cap.  Every pass
    // either lowers the product or leaves it untouched, and the loop stops on
    // the latter, so it terminates whatever the aspect ratio.
    while (current > maxTotalCells) {
      const double shrink{std::cbrt(double(maxTotalCells) / double(current))};
      x = std::max(1, int(x * shrink));
      y = std::max(1, int(y * shrink));
      z = std::max(1, int(z * shrink));
      const long long next{1LL * x * y * z};
      if (next == current)
        break; // every axis is down to one SubBox; nothing left to give
      current = next;
    }
    if (current != requested) {
      std::cout << "Sub-volume budget of " << maxTotalCells
                << " exceeded by the " << requested
                << " the interaction range asks for; using " << x << " x " << y
                << " x " << z << ".\n";
    }

    tot = x * y * z;
  #endif

  // The max(4000, 0.5 * N^2) SubBox budget that used to close this function is
  // gone.  It was there so that walking the SubBoxes could not cost more than
  // testing all N^2/2 pairs outright, and the previous commit removed that
  // reason: the pairwise search now visits only the occupied SubBoxes, so its
  // cost follows the molecule count, not the SubBox count.
  //
  // What the budget did cost was resolution, on the axes it happened to pick.
  // It shrank x and y by `scale` but z by `2 * scale`, so the SubBoxes it left
  // behind were not cubic -- which the 13-neighbour stencil assumes.  On the
  // compartment sample, 100 molecules put the budget at 5000 SubBoxes, the
  // 30-per-dimension grid of 27 000 exceeded it, and the loop settled on
  // 17 x 17 x 8: SubBoxes of 58.8 x 58.8 x 125 nm against a 28.35 nm
  // interaction range.  That grid offered 177 candidate pairs per step where
  // only 1.8 were within range, an over-inclusion of 98x, the worst of any
  // sample measured -- on the case whose own molecule count triggered the
  // guard that was meant to make it faster.
  //
  // N came from Molecule::numberOfMolecules read once here, at setup, and the
  // grid is never rebuilt, so a system that grows kept the resolution its
  // starting molecule count bought.

  tot = x * y * z;
}

/* SIMULBOX */
// Member Functions

void SimulVolume::display() {
  std::cout << "Simulation volume parameters:\n";
  std::cout << "Total sub-volumes: " << numSubCells.tot << '\n';
  std::cout << "\tDimensions: [" << numSubCells.x << ", " << numSubCells.y
            << ", " << numSubCells.z << "]\n";
  std::cout << "\tMaximum sub-volume neighbors: " << maxNeighbors << '\n';
  std::cout << "\tSub-volume size: [" << subCellSize.x << ", " << subCellSize.y
            << ", " << subCellSize.z << "]\n";
}

void SimulVolume::create_simulation_volume(const Parameters &params,
                                           const Membrane &membraneObject) {
  // Every SubBox count below is the box length divided by this, so it has to be
  // a real length.  set_rMaxLimit() leaves it at zero for a model that declares
  // forward reactions but no bimolecular one, and the division then produced a
  // negative SubBox count and the endless "CELL PAIR MAX EXCEEDED" report that
  // known_broken.tsv records against unimolecular_reverse.  Saying so is more
  // use than looping.
  if (!(params.rMaxLimit > 0.0)) {
    std::cerr << "ERROR: rMaxLimit is " << params.rMaxLimit
              << " nm, so the simulation volume cannot be divided into "
                 "sub-volumes.\n"
              << "       set_rMaxLimit() leaves it at zero when the model "
                 "declares no bimolecular reaction.\n";
    exit(1);
  }

  // Determine the number of boxes there will be in each dimension
  numSubCells = Dimensions(params, membraneObject);
  numSubCells.check_dimensions(params, membraneObject);

  // Calculate the cells' dimensions in nanometers
  if (membraneObject.waterBox.z > 0)
    subCellSize = Vec3D{membraneObject.waterBox.x / (numSubCells.x * 1.0),
                        membraneObject.waterBox.y / (numSubCells.y * 1.0),
                        membraneObject.waterBox.z / (numSubCells.z * 1.0)};
  else
    subCellSize = Vec3D{membraneObject.waterBox.x / (numSubCells.x * 1.0),
                        membraneObject.waterBox.y / (numSubCells.y * 1.0), 1};
  // The pairwise search offers same-SubBox pairs and +/-1-stencil pairs and
  // nothing else, so a SubBox edge shorter than rMaxLimit means molecules two
  // SubBoxes apart -- and therefore possibly within the interaction range --
  // are never compared.  A dimension holding a single SubBox is exempt: it has
  // no neighbours, so nothing can fall outside the stencil.  Nothing above can
  // reach this state any more; the check is here so that a future change to
  // the SubBox arithmetic cannot reintroduce silently missed reactions.
  const auto require_edge = [&](const char *axis, double edge, int count) {
    if (count > 1 && edge < params.rMaxLimit) {
      std::cerr << "ERROR: sub-volume edge along " << axis << " is " << edge
                << " nm across " << count << " sub-volumes, below the "
                << params.rMaxLimit << " nm interaction range.\n"
                << "       Reacting pairs would be missed. Exiting.\n";
      exit(1);
    }
  };
  require_edge("x", subCellSize.x, numSubCells.x);
  require_edge("y", subCellSize.y, numSubCells.y);
  if (membraneObject.waterBox.z > 0)
    require_edge("z", subCellSize.z, numSubCells.z);

  // Create cell neighborlists.
  subCellList = std::vector<SubVolume>(numSubCells.tot);
  occupancyMask.assign((numSubCells.tot + 63) / 64, 0);
  occupiedSubCells.clear();
  create_cell_neighbor_list_cubic();
}

void SimulVolume::create_cell_neighbor_list_cubic() {
  int cellNum{0};
  for (unsigned zItr{0}; zItr < numSubCells.z; ++zItr) {
    for (unsigned yItr{0}; yItr < numSubCells.y; ++yItr) {
      for (unsigned xItr{0}; xItr < numSubCells.x; ++xItr) {
        // set up SubVolume
        subCellList[cellNum].absIndex = cellNum;
        subCellList[cellNum].xIndex = xItr;
        subCellList[cellNum].yIndex = yItr;
        subCellList[cellNum].zIndex = zItr;

        /*For each cell figure out its 13 neighbors that are ~forward and up*/
        // This only works for cubic subvolumes
        if (xItr < (numSubCells.x - 1)) {
          // count all the plus x boxes
          subCellList[cellNum].neighborList.push_back(
              (xItr + 1) + yItr * numSubCells.x +
              zItr * (numSubCells.x * numSubCells.y));
          if (yItr < (numSubCells.y - 1)) {
            subCellList[cellNum].neighborList.push_back(
                (xItr + 1) + (yItr + 1) * numSubCells.x +
                zItr * (numSubCells.x * numSubCells.y));
            if (zItr < (numSubCells.z - 1)) {
              subCellList[cellNum].neighborList.push_back(
                  (xItr + 1) + (yItr + 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
          if (zItr < (numSubCells.z - 1)) {
            subCellList[cellNum].neighborList.push_back(
                (xItr + 1) + yItr * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
          }
        }
        if (yItr < (numSubCells.y - 1)) {
          subCellList[cellNum].neighborList.push_back(
              xItr + (yItr + 1) * numSubCells.x +
              zItr * (numSubCells.x * numSubCells.y));
          if (xItr > 0) {
            subCellList[cellNum].neighborList.push_back(
                (xItr - 1) + (yItr + 1) * numSubCells.x +
                zItr * (numSubCells.x * numSubCells.y));
            if (zItr < (numSubCells.z - 1)) {
              subCellList[cellNum].neighborList.push_back(
                  (xItr - 1) + (yItr + 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
          if (zItr < (numSubCells.z - 1)) {
            subCellList[cellNum].neighborList.push_back(
                xItr + (yItr + 1) * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
          }
        }
        if (zItr < (numSubCells.z - 1)) {
          subCellList[cellNum].neighborList.push_back(
              xItr + yItr * numSubCells.x +
              (zItr + 1) * (numSubCells.x * numSubCells.y));
          if (xItr > 0) {
            subCellList[cellNum].neighborList.push_back(
                xItr - 1 + yItr * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
            if (yItr > 0) {
              subCellList[cellNum].neighborList.push_back(
                  xItr - 1 + (yItr - 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
          if (yItr > 0) {
            subCellList[cellNum].neighborList.push_back(
                xItr + (yItr - 1) * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
            if (xItr < (numSubCells.x - 1)) {
              subCellList[cellNum].neighborList.push_back(
                  xItr + 1 + (yItr - 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
        }
        if (subCellList[cellNum].neighborList.size() > maxNeighbors) {
          std::cerr
              << "ERROR: Maximum number of neighbors exceeded for SubVolume "
              << cellNum << ". Exiting\n";
          exit(1);
        }
        ++cellNum;
      } // end looping over z cells
    }   // end looping over y cells
  }     // end looping over x cells
}

void SimulVolume::clear_member_lists() {
  for (size_t wordItr{0}; wordItr < occupancyMask.size(); ++wordItr) {
    uint64_t bits{occupancyMask[wordItr]};
    while (bits) {
      const int cellIndex{int(wordItr * 64) + lowest_set_bit(bits)};
      bits &= bits - 1;
      subCellList[cellIndex].memberMolList.clear();
      subCellList[cellIndex].typeMask = 0;
    }
    occupancyMask[wordItr] = 0;
  }
  occupiedSubCells.clear();
}

void SimulVolume::update_memberMolLists(
    const Parameters &params, std::vector<Molecule> &moleculeList,
    std::vector<Complex> &complexList,
    std::vector<MolTemplate> &molTemplateList, const Membrane &membraneObject,
    int simItr) {
  // make sure the list of member molecules is empty.  Sweeping all of
  // subCellList here would cost time proportional to the cell count on every
  // timestep no matter how few molecules there are, and the cell count is
  // typically far larger: a 494 nm box with a 33.7 nm interaction limit holds
  // 2744 cells, which a hundred-molecule system leaves almost entirely empty.
  clear_member_lists();

  int itrCheck =
      1000; // no need to check every step if it violates box boundaries.

  int itr{0};
  if (simItr % itrCheck != 0) {
    /*just assign bins, don't check bin limits/errors*/
    for (unsigned molItr{0}; molItr < moleculeList.size(); ++molItr) {
      Molecule &mol = moleculeList[molItr]; // just for legibility

      if (mol.isEmpty || mol.isImplicitLipid)
        continue;

      // get which box the Molecule belongs to
      int xItr{int((mol.comCoord.x + membraneObject.waterBox.x / 2) /
                   subCellSize.x)};
      int yItr{int((mol.comCoord.y + membraneObject.waterBox.y / 2) /
                   subCellSize.y)};
      int zItr;
      if (membraneObject.waterBox.z > 0)
        zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) /
                   subCellSize.z);
      else
        zItr = 0;

      // allow the modecule a bit out of the box
      if (xItr == -1)
        xItr = 0;
      if (xItr == numSubCells.x)
        xItr = numSubCells.x - 1;
      if (yItr == -1)
        yItr = 0;
      if (yItr == numSubCells.y)
        yItr = numSubCells.y - 1;
      if (zItr == -1)
        zItr = 0;
      if (zItr == numSubCells.z)
        zItr = numSubCells.z - 1;

      int currBin = xItr + (yItr * numSubCells.x) +
                    (zItr * numSubCells.x * numSubCells.y);

      mol.mySubVolIndex = currBin;
      if (currBin >= numSubCells.tot) {
        std::cerr << "Molecule " << mol.index
                  << " seems outside simulation volume, with center of mass "
                     "coordinates ["
                  << mol.comCoord << "].\n";
        exit(1);
      }
      add_member(currBin, mol.index, mol.molTypeIndex);
    }
  } else {
    /*make sure proteins are within bin limits, lipids are on membrane*/
    // Signed, so that a restart can leave the counter at -1 and have ++molItr
    // put it back to 0.  Written as `molItr = 0` in an ++molItr loop, each
    // restart below resumed at molecule 1 instead -- and since
    // clear_member_lists() had just emptied every member list, molecule 0 was
    // left out of the grid for that step, taking every pair it belongs to with
    // it.  None of the sample inputs reaches these branches, so no benchmark
    // case changes.
    for (int molItr{0}; molItr < int(moleculeList.size()); ++molItr) {
      Molecule &mol = moleculeList[molItr]; // just for legibility

      if (mol.isEmpty || mol.isImplicitLipid)
        continue;

      // get which box the Molecule belongs to
      int xItr{int((mol.comCoord.x + membraneObject.waterBox.x / 2) /
                   subCellSize.x)};
      int yItr{int((mol.comCoord.y + membraneObject.waterBox.y / 2) /
                   subCellSize.y)};
      int zItr;
      if (membraneObject.waterBox.z > 0)
        zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) /
                   subCellSize.z);
      else
        zItr = 0;

      if (xItr == -1)
        xItr = 0;
      if (xItr == numSubCells.x)
        xItr = numSubCells.x - 1;
      if (yItr == -1)
        yItr = 0;
      if (yItr == numSubCells.y)
        yItr = numSubCells.y - 1;
      if (zItr == -1)
        zItr = 0;
      if (zItr == numSubCells.z)
        zItr = numSubCells.z - 1;
      int currBin = xItr + (yItr * numSubCells.x) +
                    (zItr * numSubCells.x * numSubCells.y);

      // Make sure the Molecule is still on the membrane if its supposed to be
      if (molTemplateList[mol.molTypeIndex].isLipid) {
        // define RS3Dinput
        double RS3Dinput{0.0};

        if (membraneObject.implicitLipid == true) {
          for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
            if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] -
                         mol.molTypeIndex) < 1E-2) {
              RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
              //   std::cout << mol.molTypeIndex << "\t" <<
              //   membraneObject.RS3Dvect[RS3Dindex + 400] << "\t" << RS3Dindex
              //   << "\n";
              break;
            }
          }
        }

        if (mol.comCoord.z - 0.1 >
                -membraneObject.waterBox.z * 0.5 + RS3Dinput &&
            mol.isImplicitLipid == false) {
          //            && std::abs(mol.comCoord.z) -
          //            std::abs((membraneObject.waterBox.z / 2)) > 1E-6)
          std::cerr
              << "Molecule " << mol.index << " of type "
              << molTemplateList[mol.molTypeIndex].molName
              << " is off the membrane. Writing coordinates and exiting.\n";
          //  std::cout << mol.molTypeIndex << "\t" <<
          //  -membraneObject.waterBox.z * 0.5 + RS3Dinput << "\t" << RS3Dinput
          //  << "\n";
          write_xyz(std::string{"error_coord_dump.xyz"}, params, moleculeList,
                    molTemplateList);
          exit(1);
        }
      }

      // Now make sure the Molecule is still inside the box in all dimensions
      if (mol.comCoord.z > (membraneObject.waterBox.z / 2) ||
          mol.comCoord.z + 1E-6 < -(membraneObject.waterBox.z / 2)) {
        std::cout << "Molecule " << mol.index
                  << " is outside simulation volume in the z-dimension, with "
                     "center of mass coordinates ["
                  << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = -1; // ++molItr resumes at molecule 0
        clear_member_lists();
      } else if (mol.comCoord.y > (membraneObject.waterBox.y / 2) ||
                 mol.comCoord.y + 1E-6 < -(membraneObject.waterBox.y / 2)) {
        std::cout << "Molecule " << mol.index
                  << " is outside simulation volume in the y-dimension, with "
                     "center of mass coordinates ["
                  << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = -1; // ++molItr resumes at molecule 0
        clear_member_lists();
      } else if (mol.comCoord.x > (membraneObject.waterBox.x / 2) ||
                 mol.comCoord.x + 1E-6 < -(membraneObject.waterBox.x / 2)) {
        std::cout << "Molecule " << mol.index
                  << " is outside simulation volume in the x-dimension, with "
                     "center of mass coordinates ["
                  << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = -1; // ++molItr resumes at molecule 0
        clear_member_lists();
      } else if (currBin > (numSubCells.tot) || currBin < 0) {
        std::cout
            << "Molecule " << mol.index
            << " is outside simulation volume with center of mass coordinates ["
            << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = -1; // ++molItr resumes at molecule 0
        clear_member_lists();
      } else {
        // The Molecule is in the simulation volume, okay to proceed
        mol.mySubVolIndex = currBin;
        add_member(currBin, mol.index, mol.molTypeIndex);
      }
    } // loop over all molecules.
  }   // check all boundary limits are OK.
}

void SimulVolume::update_memberMolLists(
    const Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject,
    int simItr,
    MpiContext& mpiContext)  // xItr bin should be reduced by xOffset in order
                             // to calculate the box id in the current rank
{
  // make sure the list of member molecules is empty.  occupiedSubCells is
  // cleared alongside so that an MPI run, which reaches add_member() through
  // create_molecule_and_complex_from_rxn(), cannot accumulate stale entries.
  for (auto& subBox : subCellList) {
    subBox.memberMolList.clear();
    subBox.typeMask = 0;
  }
  std::fill(occupancyMask.begin(), occupancyMask.end(), 0);
  occupiedSubCells.clear();

  int itrCheck =
      1000;  // no need to check every step if it violates box boundaries.

  int itr{0};

  for (unsigned molItr{0}; molItr < moleculeList.size(); ++molItr) {
    Molecule& mol = moleculeList[molItr];  // just for legibility
    if (mol.isEmpty || mol.isImplicitLipid) continue;

    // Get which box the Molecule belongs to.
    // xOffset represents starting x coordinate of cell at current rank;
    // Namely, xItr was not customized for each rank, but left as it was.
    // Therefore, calculating currBin requires substracting this offset.
    int xItr{
        int((mol.comCoord.x + membraneObject.waterBox.x / 2) / subCellSize.x) -
        mpiContext.xOffset};
    int yItr{
        int((mol.comCoord.y + membraneObject.waterBox.y / 2) / subCellSize.y)};
    int zItr;
    if (membraneObject.waterBox.z > 0)
      zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) /
                 subCellSize.z);
    else
      zItr = 0;

    // allow the modecule a bit out of the box
    // TODO: Check whether this should be only for edge ranks? What happens if
    // xItr should belong to next or previous rank?
    //         if (xItr == -1)
    if ((xItr == -1) && (mpiContext.rank == 0)) xItr = 0;
    //        if (xItr == numSubCells.x)
    if ((xItr == numSubCells.x) &&
        (mpiContext.rank ==
         mpiContext.nprocs - 1))  // TODO (last ranks might be unoccupied)
      xItr = numSubCells.x - 1;
    if (yItr == -1) yItr = 0;
    if (yItr == numSubCells.y) yItr = numSubCells.y - 1;
    if (zItr == -1) zItr = 0;
    if (zItr == numSubCells.z) zItr = numSubCells.z - 1;

    // ignore the molecules that outside the rank
    if (mpiContext.rank > 0 && xItr < 0) {
      continue;
    }

    if (mpiContext.rank < mpiContext.nprocs - 1 && xItr >= numSubCells.x) {
      continue;
    }

    int currBin =
        xItr + (yItr * numSubCells.x) + (zItr * numSubCells.x * numSubCells.y);

    if (simItr % itrCheck != 0) {
      /*just assign bins, don't check bin limits/errors*/
      mol.mySubVolIndex = currBin;
      if (currBin >= numSubCells.tot) {
        std::cerr << "Molecule " << mol.index << " (ID=" << mol.id << ")"
                  << " seems outside simulation volume, with center of mass "
                     "coordinates ["
                  << mol.comCoord << "].\n";
        std::cerr << "numSubCells.x = " << numSubCells.x
                  << ", numSubCells.y = " << numSubCells.y
                  << ", numSubCells.z = " << numSubCells.z << "\n";
        std::cerr << "xItr = " << xItr << ", yItr = " << yItr
                  << ", zItr = " << zItr << "\n";
        std::cerr << "currBin = " << currBin
                  << ", numSubCells.tot = " << numSubCells.tot
                  << ", mpiContext.xOffset = " << mpiContext.xOffset << "\n";
        error("mol outside box.");
        exit(1);
      }
      add_member(currBin, molItr, mol.molTypeIndex);
    } else {
      /*make sure proteins are within bin limits, lipids are on membrane*/
      // Make sure the Molecule is still on the membrane if its supposed to be
      if (std::abs(molTemplateList[mol.molTypeIndex].D.z) < 1E-10) {
        // define RS3Dinput
        double RS3Dinput{0.0};

        if (membraneObject.implicitLipid == true) {
          for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
            if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] -
                         mol.molTypeIndex) < 1E-2) {
              RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
              //   std::cout << mol.molTypeIndex << "\t" <<
              //   membraneObject.RS3Dvect[RS3Dindex + 400] << "\t" << RS3Dindex
              //   << "\n";
              break;
            }
          }
        }

        if (mol.comCoord.z - 0.1 >
                -membraneObject.waterBox.z * 0.5 + RS3Dinput &&
            mol.isImplicitLipid == false) {
          //            && std::abs(mol.comCoord.z) -
          //            std::abs((membraneObject.waterBox.z / 2)) > 1E-6)
          std::cerr
              << "Molecule " << mol.index << " of type "
              << molTemplateList[mol.molTypeIndex].molName
              << " is off the membrane. Writing coordinates and exiting.\n";
          //  std::cout << mol.molTypeIndex << "\t" <<
          //  -membraneObject.waterBox.z * 0.5 + RS3Dinput << "\t" << RS3Dinput
          //  << "\n";
          write_xyz(std::string{"error_coord_dump.xyz"}, params, moleculeList,
                    molTemplateList);
          exit(1);
        }
      }

      // Now make sure the Molecule is still inside the box in all dimensions
      if (mol.comCoord.z > (membraneObject.waterBox.z / 2) ||
          mol.comCoord.z + 1E-6 < -(membraneObject.waterBox.z / 2)) {
        std::cout << "Molecule " << mol.index
                  << " is outside simulation volume in the z-dimension, with "
                     "center of mass coordinates ["
                  << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = 0;
        for (auto& subBox : subCellList) {
          subBox.memberMolList.clear();
          subBox.typeMask = 0;
        }
        std::fill(occupancyMask.begin(), occupancyMask.end(), 0);
        occupiedSubCells.clear();
      } else if (mol.comCoord.y > (membraneObject.waterBox.y / 2) ||
                 mol.comCoord.y + 1E-6 < -(membraneObject.waterBox.y / 2)) {
        std::cout << "Molecule " << mol.index
                  << " is outside simulation volume in the y-dimension, with "
                     "center of mass coordinates ["
                  << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = 0;
        for (auto& subBox : subCellList) {
          subBox.memberMolList.clear();
          subBox.typeMask = 0;
        }
        std::fill(occupancyMask.begin(), occupancyMask.end(), 0);
        occupiedSubCells.clear();
      } else if (mol.comCoord.x > (membraneObject.waterBox.x / 2) ||
                 mol.comCoord.x + 1E-6 < -(membraneObject.waterBox.x / 2)) {
        std::cout << "Molecule " << mol.index
                  << " is outside simulation volume in the x-dimension, with "
                     "center of mass coordinates ["
                  << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = 0;
        for (auto& subBox : subCellList) {
          subBox.memberMolList.clear();
          subBox.typeMask = 0;
        }
        std::fill(occupancyMask.begin(), occupancyMask.end(), 0);
        occupiedSubCells.clear();
      } else if (currBin > (numSubCells.tot) || currBin < 0) {
        std::cout
            << "Molecule " << mol.index
            << " is outside simulation volume with center of mass coordinates ["
            << mol.comCoord << "]. Attempting to fit back into box.\n";
        complexList[mol.myComIndex].put_back_into_SimulVolume(
            itr, mol, membraneObject, moleculeList, molTemplateList);
        // reset member search
        molItr = 0;
        for (auto& subBox : subCellList) {
          subBox.memberMolList.clear();
          subBox.typeMask = 0;
        }
        std::fill(occupancyMask.begin(), occupancyMask.end(), 0);
        occupiedSubCells.clear();
      } else {
        // The Molecule is in the simulation volume, okay to proceed
        mol.mySubVolIndex = currBin;
        add_member(currBin, mol.index, mol.molTypeIndex);
      }
    }  // check all boundary limits are OK.
  }    // loop over all molecules.
}
