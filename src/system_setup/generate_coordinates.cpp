#include "classes/class_Coord.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vector.hpp"
#include "io/io.hpp"
#include "parser/parser_functions.hpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include <cctype>
#include <cmath>
#include <iostream>

// Forward declarations of helper functions

int fixOverlappingMolecules(std::vector<Molecule> &moleculeList, const std::vector<MolTemplate> &molTemplateList, std::vector<Complex> &complexList, 
  const std::vector<ForwardRxn> &forwardRxns, const Membrane &membraneObject, const std::vector<int> &changedMoleculeIndices = {}, bool checkOnlyChanged = false);

std::vector<int> updateMoleculeCoordinates(std::vector<Molecule> &moleculeList,const std::vector<MolTemplate> &molTemplateList,const std::string &coordinateFileName);

bool checkAndResolveOverlap(const Molecule &mol1, Molecule &mol2, const std::vector<MolTemplate> &molTemplateList, const std::vector<ForwardRxn> &forwardRxns,
  const Membrane &membraneObject, std::vector<Complex> &complexList, bool mol2Changed);

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
  // Fix overlapping molecules in the initial configuration
  fixOverlappingMolecules(moleculeList, molTemplateList, complexList, forwardRxns, membraneObject);
  // Apply coordinates from PDB file if provided
  if (!coordinateFileName.empty())
  {
    std::cout << "\ncoordinate file: " << coordinateFileName << std::endl;

    // Update molecular coordinates from PDB file
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
    bool mol2Changed)
{
  // Skip check if molecules are not in vicinity (optimization)
  if (!are_molecules_in_vicinity(mol1, mol2, molTemplateList))
    return false;
    
  bool overlapFound = false;
  const MolTemplate &mol2Temp = molTemplateList[mol2.molTypeIndex];

  // Check each interface pair for possible overlaps
  for (unsigned int iface1Itr = 0; iface1Itr < mol1.interfaceList.size(); ++iface1Itr)
  {
    const auto &iface1 = mol1.interfaceList[iface1Itr];
    
    for (unsigned int iface2Itr = 0; iface2Itr < mol2.interfaceList.size(); ++iface2Itr)
    {
      const auto &iface2 = mol2.interfaceList[iface2Itr];
      int rxnIndex = 0;
      bool theyInteract = false;
      
      // Check if molecules interact through a reaction
      for (auto rxnItr : molTemplateList[mol1.molTypeIndex].interfaceList[iface1Itr].stateList[0].myForwardRxns)
      {
        const ForwardRxn &rxn = forwardRxns[rxnItr];
        // Check if both molecules are reactants in the same reaction
        if ((rxn.reactantListNew[0].molTypeIndex == mol1.molTypeIndex &&
             rxn.reactantListNew[1].molTypeIndex == mol2.molTypeIndex) ||
            (rxn.reactantListNew[0].molTypeIndex == mol2.molTypeIndex &&
             rxn.reactantListNew[1].molTypeIndex == mol1.molTypeIndex))
        {
          theyInteract = true;
          rxnIndex = &rxn - &forwardRxns[0];
          break;
        }
      }
      
      if (theyInteract)
      {
        // Calculate distance between interfaces
        Vector tmpVec{iface1.coord - iface2.coord};
        double mag = tmpVec.x * tmpVec.x + tmpVec.y * tmpVec.y + tmpVec.z * tmpVec.z;
        
        // Check if interfaces are too close (overlap)
        if (mag < forwardRxns[rxnIndex].bindRadius * forwardRxns[rxnIndex].bindRadius)
        {
          if (mol2Changed)
          {
            // Just warn about overlap if both molecules came from PDB file
            // not considered as overlap since this is defined by users
            std::cout << "|-WARNING!!! " << molTemplateList[mol1.molTypeIndex].molName
                      << " at " << mol1.comCoord.x << " " << mol1.comCoord.y
                      << " " << mol1.comCoord.z << " and " << mol2Temp.molName 
                      << " at " << mol2.comCoord.x << " " << mol2.comCoord.y 
                      << " " << mol2.comCoord.z << " overlaps. Please check your input!" 
                      << std::endl;
          }
          else
          {
            overlapFound = true;
            // Relocate the molecule to resolve overlap
            mol2.create_random_coords(mol2Temp, membraneObject);
            complexList[mol2.myComIndex].comCoord = mol2.comCoord;
          }
          
          return overlapFound; // Exit early once overlap is found
        }
      }
    }
  }

  return overlapFound;
}

/**
 * @brief Updates the first N molecules of each type based on PDB file data
 * @return Vector of indices of molecules that were moved
 */
std::vector<int> updateMoleculeCoordinates(
    std::vector<Molecule> &moleculeList,
    const std::vector<MolTemplate> &molTemplateList,
    const std::string &coordinateFileName)
{
  // Vector to store indices of modified molecules
  std::vector<int> changedMoleculeIndex;

  // Open the PDB file for reading
  std::ifstream inputFile(coordinateFileName);
  if (!inputFile.is_open())
  {
    std::cerr << "Error: Could not open coordinate file: " << coordinateFileName << std::endl;
    return changedMoleculeIndex;
  }
  std::cout << "Parsing given coordinates..." << std::endl;
  
  // Create a map to store coordinates for each molecule type
  // Key: molecule name, Value: vector of coordinates from PDB
  std::map<std::string, std::vector<Vector>> moleculeCoordinates;

  // Parse the PDB file to extract COM coordinates
  std::string line;
  while (std::getline(inputFile, line))
  {
    // Skip lines that don't start with "ATOM"
    if (line.substr(0, 4) != "ATOM")
    {
      continue;
    }

    // Parse the line according to PDB format
    std::string atomType = line.substr(12, 4);
    // Trim whitespace from atomType
    atomType.erase(0, atomType.find_first_not_of(" "));
    atomType.erase(atomType.find_last_not_of(" ") + 1);

    // Only process COM atom type (center of mass)
    if (atomType != "COM")
    {
      continue;
    }

    std::string moleculeName = line.substr(17, 4);
    // Trim whitespace from moleculeName
    moleculeName.erase(0, moleculeName.find_first_not_of(" "));
    moleculeName.erase(moleculeName.find_last_not_of(" ") + 1);

    // Extract x, y, z coordinates
    double x = std::stod(line.substr(30, 8));
    double y = std::stod(line.substr(38, 8));
    double z = std::stod(line.substr(46, 8));

    // Store coordinates for this molecule type in a vector
    // Each entry in the vector represents one molecule of that type from the PDB
    moleculeCoordinates[moleculeName].push_back(Vector{x, y, z});
  }

  // Close the input file
  inputFile.close();

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
    const std::vector<Vector> &coordinates = molCoordPair.second;
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
      Vector transVec{
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

  int currItr = 0;
  const int MAX_ITERATIONS = 50; // Maximum retry iterations

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
          if (checkAndResolveOverlap(mol1, mol2, molTemplateList, forwardRxns, membraneObject, complexList, mol2Changed))
          {
            hasOverlap = true;
            ++numOverlap;
          }
        }
      }
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
          if (checkAndResolveOverlap(mol1, mol2, molTemplateList, forwardRxns, membraneObject, complexList, false))
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