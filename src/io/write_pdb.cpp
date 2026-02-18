#include "io/io.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_pdb(long long int simItr, 
    unsigned frameNum, 
    const Parameters& params, 
    const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, 
    const Membrane& membraneObject)
{
    std::string pdbFileName = "PDB/" + std::to_string(frameNum) + ".pdb";
    std::ofstream pdbFile(pdbFileName);
    if (!pdbFile.is_open()) {
        std::cerr << "Error: Unable to open PDB file for writing: " << pdbFileName << std::endl;
        return;
    }

    auto printTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    pdbFile << std::left << std::setw(6) << "TITLE" 
            << " PDB TIMESTEP " << simItr 
            << " CREATED " << std::ctime(&printTime);
    pdbFile << std::left << std::setw(6) << "CRYST1" 
            << std::setw(9) << membraneObject.waterBox.x 
            << std::setw(9) << membraneObject.waterBox.y 
            << std::setw(9) << membraneObject.waterBox.z 
            << std::setw(7) << 90 << std::setw(7) << 90 << std::setw(7) << 90 
            << " P 1" << std::endl;

    int complexIndex = 1;
    int moleculeIndex = 1;
    int atomIndex = 1;

    for (const auto& mol : moleculeList) {
        if (mol.isImplicitLipid || mol.isEmpty) {
            continue;
        }
        const MolTemplate& molTemplate = molTemplateList[mol.molTypeIndex];
        complexIndex = mol.myComIndex;

        // Write Center of Mass (COM)
        int comIndex = atomIndex;
        pdbFile << "ATOM  "
                << std::right << std::setw(5) << (atomIndex % 100000)
                << " " << std::left << std::setw(4) << "COM"
                << " " << std::right << std::setw(3) << molTemplate.molName.substr(0, 3)
                << " " << std::right << std::setw(4) << (moleculeIndex % 10000)
                << "    "
                << std::right << std::setw(8) << std::fixed << std::setprecision(1) 
                << (mol.comCoord.x + membraneObject.waterBox.x / 2)
                << std::right << std::setw(8) 
                << (mol.comCoord.y + membraneObject.waterBox.y / 2)
                << std::right << std::setw(8) 
                << (mol.comCoord.z + membraneObject.waterBox.z / 2)
                << std::right << std::setw(6) << 1.00  // Occupancy
                << std::right << std::setw(6) << complexIndex  // B-factor (Complex ID)
                << std::right << std::setw(12) << molTemplate.molName.substr(0, 4)  // Segment ID
                << std::endl;

        ++atomIndex;

        // Write Interfaces
        for (const auto& interface : mol.interfaceList) {
            int interfaceIndex = atomIndex;
            pdbFile << "ATOM  "
                    << std::right << std::setw(5) << (atomIndex % 100000)
                    << " " << std::left << std::setw(4) 
                    << molTemplate.interfaceList[interface.relIndex].name.substr(0, 4) 
                    << " " << std::right << std::setw(3) 
                    << molTemplate.molName.substr(0, 3)
                    << " " << std::right << std::setw(4) 
                    << (moleculeIndex % 10000)
                    << "    "
                    << std::right << std::setw(8) << std::fixed << std::setprecision(1) 
                    << (interface.coord.x + membraneObject.waterBox.x / 2)
                    << std::right << std::setw(8) 
                    << (interface.coord.y + membraneObject.waterBox.y / 2)
                    << std::right << std::setw(8) 
                    << (interface.coord.z + membraneObject.waterBox.z / 2)
                    << std::right << std::setw(6) << 1.00  // Occupancy
                    << std::right << std::setw(6) << complexIndex  // B-factor (Complex ID)
                    << std::right << std::setw(12) << molTemplate.molName.substr(0, 4)  // Segment ID
                    << std::endl;

            ++atomIndex;
        }

        ++moleculeIndex;
    }

    // Write END record
    pdbFile << "END" << std::endl;
    pdbFile.close();
}
