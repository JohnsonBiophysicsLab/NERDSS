#include "classes/class_Vec3D.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "io/io.hpp"
#include "io/read_coordinate_file.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "system_setup/system_setup.hpp"
#include <cctype>
#include <cmath>
#include <iostream>

// Forward declarations of helper functions

int fixOverlappingMolecules(std::vector<Molecule> &moleculeList, const std::vector<MolTemplate> &molTemplateList, std::vector<Complex> &complexList, 
  const std::vector<ForwardRxn> &forwardRxns, const Membrane &membraneObject, const std::vector<int> &changedMoleculeIndices = {}, bool checkOnlyChanged = false);

std::vector<int> updateMoleculeCoordinates(std::vector<Molecule> &moleculeList,const std::vector<MolTemplate> &molTemplateList,const std::string &coordinateFileName);

bool checkAndResolveOverlap(const Molecule &mol1, Molecule &mol2, const std::vector<MolTemplate> &molTemplateList, const std::vector<ForwardRxn> &forwardRxns,
  const Membrane &membraneObject, std::vector<Complex> &complexList, bool mol2Changed,
  const std::vector<std::vector<double>> &maxBindRadius);

// Main function definition

/**
 * @brief Main function to generate coordinates for all molecules in the system
 */
void generate_coordinates(const Parameters &params,
                          std::vector<Molecule> &moleculeList,
                          std::vector<Complex> &complexList,
                          std::vector<MolTemplate> &molTemplateList,
                          const std::vector<ForwardRxn> &forwardRxns,
                          const Membrane &membraneObject,
                          std::string &coordinateFileName)
{
  // First create all the molecules and their corresponding complexes, with
  // random center of mass coordinates
  for (auto &oneTemp : molTemplateList)
  {
    if (oneTemp.isImplicitLipid == false)
    {
      for (unsigned itr{0}; itr < oneTemp.copies; ++itr)
      {
        moleculeList.emplace_back(initialize_molecule(Complex::numberOfComplexes, params, oneTemp, membraneObject));
        complexList.emplace_back(initialize_complex(moleculeList.back(), molTemplateList[moleculeList.back().molTypeIndex]));
        oneTemp.monomerList.emplace_back(moleculeList.back().index);
        moleculeList.back().complexId = complexList.back().id;
      }
    }
    else
    {
      oneTemp.isLipid = true;
      moleculeList.emplace_back(initialize_molecule(Complex::numberOfComplexes, params, oneTemp, membraneObject));
      complexList.emplace_back(initialize_complex(moleculeList.back(), molTemplateList[moleculeList.back().molTypeIndex]));
      moleculeList.back().complexId = complexList.back().id;
    }
  }
  std::cout << "NUMBER OF MOLECULES IN GEN COORDS: " << moleculeList.size()
            << std::endl;
  if (moleculeList.size() == 0 && params.numTotalUnits == 0)
  {
    std::cout << "No molecules present, skipping coordinate generation.\n";
    return;
  }
  // The interface states the starting copy numbers ask for ("A : 103
  // (ser~U), 5 (ser~P)"), so that the overlap passes test the pairs that can
  // react at the first step.  initialize_states() used to be the first to set
  // them, after placement, so every molecule was placed as if it started in
  // its interfaces' first states.  It still runs later, for the implicit
  // lipid's counts, and sets the same states again.
  for (const auto &oneTemp : molTemplateList)
    initialize_molecule_states(moleculeList, oneTemp);

  // Fix overlapping molecules in the initial configuration
  fixOverlappingMolecules(moleculeList, molTemplateList, complexList, forwardRxns, membraneObject);
  // Apply coordinates from the coordinate file if provided
  if (!coordinateFileName.empty())
  {
    std::cout << "\ncoordinate file: " << coordinateFileName << std::endl;

    // Update molecular coordinates from the coordinate file
    std::vector<int> changedMoleculeIndex = updateMoleculeCoordinates(moleculeList, molTemplateList, coordinateFileName);

    // Fix overlaps for molecules with new coordinates
    fixOverlappingMolecules(moleculeList, molTemplateList, complexList, forwardRxns, membraneObject, changedMoleculeIndex, true);
  }

  // Write final coordinates to XYZ file
  write_xyz("DATA/initial_crds.xyz", params, moleculeList, molTemplateList);
}

// Helper routines implementation

/**
 * @brief Check if two molecules overlap and resolve if needed
 * @return True if overlap was found, false otherwise
 */
bool checkAndResolveOverlap(
    const Molecule &mol1,
    Molecule &mol2,
    const std::vector<MolTemplate> &molTemplateList,
    const std::vector<ForwardRxn> &forwardRxns,
    const Membrane &membraneObject,
    std::vector<Complex> &complexList,
    bool mol2Changed,
    const std::vector<std::vector<double>> &maxBindRadius)
{
  // Two types with no bimolecular reaction between them cannot overlap, and
  // this runs for every ordered pair of molecules: most pairs of a large
  // model end here.
  const double pairBindRadius = maxBindRadius[mol1.molTypeIndex][mol2.molTypeIndex];
  if (pairBindRadius <= 0.0)
    return false;

  // Skip check if molecules are not in vicinity (optimization)
  if (!are_molecules_in_vicinity(mol1, mol2, molTemplateList, pairBindRadius))
    return false;

  // An implicit lipid has no position to overlap, and an emptied slot no
  // interfaces; the time step's overlap checks skip both too.  By the
  // template: initialize_paramters_for_implicitlipid_model() sets the
  // molecule's own flag only after placement.
  if (molTemplateList[mol1.molTypeIndex].isImplicitLipid || molTemplateList[mol2.molTypeIndex].isImplicitLipid
      || mol1.isEmpty || mol2.isEmpty)
    return false;

  // The test placement uses for created and added molecules too: an interface
  // pair overlaps when the two could react with each other and sit closer than
  // that reaction's bindRadius.  This used to take, for each interface of
  // mol1, the first reaction between the two molecule types, and hold it to
  // every interface of mol2 with that sigma, whichever interface or state the
  // reaction named: a mol2 interface that reacts with nothing pushed mol2 away,
  // and the answer depended on which molecule came first, which the -c pass
  // fixes as the file one.
  if (!interfacesOverlap(mol1, mol2, forwardRxns, pairBindRadius))
    return false;

  if (mol2Changed)
  {
    // Just warn about overlap if both molecules came from PDB file
    // not considered as overlap since this is defined by users
    std::cout << "|-WARNING!!! " << molTemplateList[mol1.molTypeIndex].molName
              << " at " << mol1.comCoord.x << " " << mol1.comCoord.y
              << " " << mol1.comCoord.z << " and " << molTemplateList[mol2.molTypeIndex].molName
              << " at " << mol2.comCoord.x << " " << mol2.comCoord.y
              << " " << mol2.comCoord.z << " overlaps. Please check your input!"
              << std::endl;
    return false;
  }

  // Relocate the molecule to resolve overlap.  create_random_coords() puts
  // every interface back in its first state; the molecule keeps the ones it
  // starts in, which the next pass tests it with.
  std::vector<Molecule::Iface> startingStates = mol2.interfaceList;
  mol2.create_random_coords(molTemplateList[mol2.molTypeIndex], membraneObject);
  for (unsigned int ifaceItr = 0; ifaceItr < mol2.interfaceList.size(); ++ifaceItr)
  {
    mol2.interfaceList[ifaceItr].index = startingStates[ifaceItr].index;
    mol2.interfaceList[ifaceItr].stateIden = startingStates[ifaceItr].stateIden;
    mol2.interfaceList[ifaceItr].stateIndex = startingStates[ifaceItr].stateIndex;
  }
  mol2.refresh_iface_state_bits();
  complexList[mol2.myComIndex].comCoord = mol2.comCoord;
  return true;
}

/**
 * @brief Updates the first N molecules of each type based on the coordinate file (PDB, XYZ or CSV)
 * @return Vec3D of indices of molecules that were moved
 */
std::vector<int> updateMoleculeCoordinates(
    std::vector<Molecule> &moleculeList,
    const std::vector<MolTemplate> &molTemplateList,
    const std::string &coordinateFileName)
{
  // Vec3D to store indices of modified molecules
  std::vector<int> changedMoleculeIndex;

  std::cout << "Parsing given coordinates..." << std::endl;

  // Key: molecule name, Value: COM positions from the file (PDB, XYZ or CSV;
  // see io/read_coordinate_file.hpp).  Names not in the system are skipped
  // below, and a name with more rows than copies loses the extra rows.
  std::map<std::string, std::vector<Vec3D>> moleculeCoordinates = read_coordinate_file(coordinateFileName);

  // Create index maps for efficient lookup of molecules by type
  std::map<std::string, std::vector<size_t>> moleculeIndicesByType;

  // First pass: build index of molecules by type (done only once)
  for (size_t i = 0; i < moleculeList.size(); ++i)
  {
    std::string molName = molTemplateList[moleculeList[i].molTypeIndex].molName;
    moleculeIndicesByType[molName].push_back(i);
  }

  // For each molecule type in the PDB file
  for (const auto &molCoordPair : moleculeCoordinates)
  {
    // Extract key and value from the map pair
    const std::string &molName = molCoordPair.first;
    const std::vector<Vec3D> &coordinates = molCoordPair.second;
    // Skip if there are no molecules of this type in the system
    if (moleculeIndicesByType.find(molName) == moleculeIndicesByType.end())
    {
      continue;
    }

    // Get the indices of molecules of this type
    const auto &indices = moleculeIndicesByType[molName];

    // Determine how many molecules to update (minimum of coordinates available and molecules present)
    size_t numToUpdate = std::min(coordinates.size(), indices.size());

    // Update only the first N molecules of this type
    for (size_t i = 0; i < numToUpdate; ++i)
    {
      size_t molIndex = indices[i];
      Molecule &movingMol = moleculeList[molIndex];

      // Get target coordinates
      double targetX = coordinates[i].x;
      double targetY = coordinates[i].y;
      double targetZ = coordinates[i].z;

      // Calculate translation vector
      Vec3D transVec{
          targetX - movingMol.comCoord.x,
          targetY - movingMol.comCoord.y,
          targetZ - movingMol.comCoord.z};

      // Update center of mass coordinate
      movingMol.comCoord = transVec + movingMol.comCoord;

      // Update all interface coordinates
      for (unsigned int j = 0; j < movingMol.interfaceList.size(); ++j)
      {
        movingMol.interfaceList[j].coord =
            transVec + movingMol.interfaceList[j].coord;
      }

      // Record this molecule as changed
      changedMoleculeIndex.emplace_back(molIndex);
    }
  }

  return changedMoleculeIndex;
}

/**
 * @brief Fix overlapping molecules in the system
 * @param checkOnlyChanged Flag to only check changed molecules (default: false)
 * @return Number of iterations performed
 */
int fixOverlappingMolecules(
    std::vector<Molecule> &moleculeList,
    const std::vector<MolTemplate> &molTemplateList,
    std::vector<Complex> &complexList,
    const std::vector<ForwardRxn> &forwardRxns,
    const Membrane &membraneObject,
    const std::vector<int> &changedMoleculeIndices,
    bool checkOnlyChanged)
{
  std::cout << (checkOnlyChanged ? "|-" : "")
            << "Finding and fixing overlapping proteins.\n";

  // How close two molecules' centers must be before their interfaces can overlap
  const std::vector<std::vector<double>> maxBindRadius = max_bind_radius_between_types(molTemplateList, forwardRxns);

  int currItr = 0;
  const int MAX_ITERATIONS = 50; // Maximum retry iterations

  // Molecules re-placed at random by the checkOnlyChanged pass that have not
  // yet been checked against every other molecule at their new position
  std::vector<int> replacedMoleculeIndices;

  while (currItr < MAX_ITERATIONS)
  {
    bool hasOverlap = false;
    ++currItr;
    int numOverlap = 0;

    if (checkOnlyChanged)
    {
      // Only check molecules that were changed from PDB
      for (int mol1Itr : changedMoleculeIndices)
      {
        auto &mol1 = moleculeList[mol1Itr];

        for (unsigned long mol2Itr = 0; mol2Itr < moleculeList.size(); ++mol2Itr)
        {
          if (mol1Itr == mol2Itr)
            continue; // Skip self comparison

          // Check if molecule 2 was also changed from PDB
          bool mol2Changed = (std::find(changedMoleculeIndices.begin(),
                                       changedMoleculeIndices.end(),
                                       mol2Itr) != changedMoleculeIndices.end());

          auto &mol2 = moleculeList[mol2Itr];

          // Check for overlap and resolve if needed
          if (checkAndResolveOverlap(mol1, mol2, molTemplateList, forwardRxns, membraneObject, complexList, mol2Changed, maxBindRadius))
          {
            hasOverlap = true;
            ++numOverlap;
            // mol2 was re-placed at random (overlaps between two file
            // molecules only warn), so it has to be checked again
            if (std::find(replacedMoleculeIndices.begin(), replacedMoleculeIndices.end(),
                          static_cast<int>(mol2Itr)) == replacedMoleculeIndices.end())
              replacedMoleculeIndices.push_back(mol2Itr);
          }
        }
      }

      // Check each re-placed molecule against all other molecules. It is
      // always the one moved again, so file molecules and molecules that did
      // not collide keep their positions
      std::vector<int> stillReplaced;
      for (int mol2Itr : replacedMoleculeIndices)
      {
        auto &mol2 = moleculeList[mol2Itr];

        for (unsigned long mol1Itr = 0; mol1Itr < moleculeList.size(); ++mol1Itr)
        {
          if (mol1Itr == static_cast<unsigned long>(mol2Itr))
            continue; // Skip self comparison

          if (checkAndResolveOverlap(moleculeList[mol1Itr], mol2, molTemplateList, forwardRxns, membraneObject, complexList, false, maxBindRadius))
          {
            hasOverlap = true;
            ++numOverlap;
            stillReplaced.push_back(mol2Itr);
            break; // Check the new position from the start next iteration
          }
        }
      }
      replacedMoleculeIndices.swap(stillReplaced);
    }
    else
    {
      // Check all molecules against each other
      unsigned long molListSize = moleculeList.size();

      for (unsigned long mol1Itr = 0; mol1Itr < molListSize; ++mol1Itr)
      {
        auto &mol1 = moleculeList[mol1Itr];

        for (unsigned long mol2Itr = 0; mol2Itr < molListSize; ++mol2Itr)
        {
          if (mol1Itr == mol2Itr)
            continue; // Skip self comparison

          auto &mol2 = moleculeList[mol2Itr];

          // Check for overlap and resolve if needed
          if (checkAndResolveOverlap(mol1, mol2, molTemplateList, forwardRxns, membraneObject, complexList, false, maxBindRadius))
          {
            hasOverlap = true;
            ++numOverlap;
          }
        }
      }
    }

    // Report progress
    std::cout << (checkOnlyChanged ? "|-" : "")
              << "Iteration: " << currItr << "\n"
              << (checkOnlyChanged ? "|-" : "")
              << "Overlapping Proteins (" << std::boolalpha << hasOverlap
              << "): " << numOverlap << '\n';

    // Exit if no overlaps found
    if (!hasOverlap)
    {
      std::cout << (checkOnlyChanged ? "|-" : "")
                << "No overlapping proteins found.\n";
      break;
    }
  }

  return currItr; // Return number of iterations performed
}