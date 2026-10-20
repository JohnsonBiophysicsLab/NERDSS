/*! \file class_moltemplate_functions

 * ### Created on 10/18/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */

#include "classes/class_MolTemplate.hpp"
#include "parser/parser_functions.hpp"

#include <algorithm>
#include <cmath>

#include <iostream>

/* INTERFACE */
// Static Variables
int Interface::State::totalNumOfStates = 0;
std::vector<int> MolTemplate::absToRelIface {};
unsigned MolTemplate::numMolTypes = 0;
std::vector<int> MolTemplate::numEachMolType {};

// Constructors
Interface::State::State(int index)
    : index(index)
{
}

Interface::Interface(std::string name, const Vec3D& iCoord)
    : iCoord(iCoord)
    , name(std::move(name))
{
}

Interface::State::State(const std::string& ifaceAndStateName, int index)
    : ifaceAndStateName(ifaceAndStateName)
    , index(index)
{
}

void Interface::build_excluded_ifaces()
{
    excludedIfaces.clear();
    for (std::size_t i { 0 }; i < excludeVolumeBoundList.size(); ++i)
        excludedIfaces.push_back(ExcludedIface { excludeVolumeBoundList[i], excludeVolumeBoundIfaceList[i],
            excludeRadiusList[i], excludeVolumeBoundReactList[i] });
    // Stable, so that entries naming the same interface keep their list order.
    std::stable_sort(excludedIfaces.begin(), excludedIfaces.end(),
        [](const ExcludedIface& a, const ExcludedIface& b) {
            return a.molTypeIndex != b.molTypeIndex ? a.molTypeIndex < b.molTypeIndex
                                                    : a.relIfaceIndex < b.relIfaceIndex;
        });
}

/*******************/

/* MOLTEMPLATE */

// Functions
void MolTemplate::display() const
{
    std::cout << "Molecule template " << molTypeIndex << '\n';
    std::cout << "Name: " << molName << '\n';
    std::cout << "Copy number:" << copies << '\n';
    std::cout << " Diffusion trans: " << D.x << ' ' << D.y << ' ' << D.z << '\n';
    std::cout << " Diffusion Rot: " << Dr.x << ' ' << Dr.y << ' ' << Dr.z << '\n';

    if (isLipid) {
        if (isImplicitLipid)
            std::cout << "Is a implicitLipid: " << std::boolalpha << isImplicitLipid << '\n';
        else
            std::cout << "Is a lipid: " << std::boolalpha << isLipid << '\n';
    }
    std::cout << "Is a promoter: " << std::boolalpha << isPromoter << '\n';
    std::cout << "Is a rod: " << std::boolalpha << isRod << '\n';
    std::cout << "Is a point: " << std::boolalpha << isPoint << '\n';
    std::cout << "Radius: " << radius << '\n';
    std::cout << "\nInterfaces:\n";
    for (auto& iface : interfaceList) {
        if (iface.stateList.size() == 1) {
            std::cout << "Name: " << iface.name << '\n';
            std::cout << "Relative index: " << iface.index << '\n';
            std::cout << "Absolute index: " << iface.stateList[0].index << '\n';
            if (!iface.stateList[0].myForwardRxns.empty()) {
                std::cout << "Forward Reactions: ";
                for (auto& rxn : iface.stateList[0].myForwardRxns)
                    std::cout << " (" << rxn << ")";
                std::cout << '\n';
            }
            if (!iface.stateList[0].myCreateDestructRxns.empty()) {
                std::cout << "Creation/Destruction Reactions: ";
                for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                    std::cout << " (" << rxn << ")";
                std::cout << '\n';
            }
            std::cout << '\n';
        } else {
            std::cout << "Name: " << iface.name << '\n';
            std::cout << "Relative index: " << iface.index << '\n';
            std::cout << "States:" << '\n';
            for (auto& state : iface.stateList) {
                std::cout << ""
                          << "(identity: " << state.iden << ", absolute index: " << state.index << ")\n";
                //                if (!state.myForwardRxns.empty()) {
                //                    std::cout << "Forward Reactions: ";
                //                    for (auto& rxn : state.myForwardRxns)
                //                        for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                //                    std::cout << '\n';
                //                }
                //                if (!state.stateChangeRxns.empty()) {
                //                    std::cout << "State Change Reactions: ";
                //                    for (auto& rxn : state.stateChangeRxns)
                //                        for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                //                    std::cout << '\n';
                //                }
                //                if (!state.myCreateDestructRxns.empty()) {
                //                    std::cout << "Creation/Destruction Reactions: ";
                //                    for (auto& rxn : state.myCreateDestructRxns)
                //                        for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                //                    std::cout << '\n';
                //                }
                std::cout << '\n';
            }
        }
    }
}

int MolTemplate::find_relIndex_from_absIndex(int targStateIndex) const
{
    /*!
     * \brief Finds the relative index of Interface given the absolute index of the State
     */

    // subroutine to find the Interface index from the State index
    int i { 0 };
    for (auto& iface : interfaceList) {
        for (auto& state : iface.stateList) {
            if (state.index == targStateIndex) {
                return i;
            }
        }
        ++i;
    }
    exit(1);
    //    not_in_cont_err(__PRETTY_FUNCTION__, __LINE__);
    return 0; // it'll terminate before it hits this, this is just so it
    // doesn't yell at me for not returning an int
}

void MolTemplate::cache_diffusion_derivatives()
{
    /*! \brief Precomputes the per-axis inverse cube root of Dr.
     *
     * Complex::update_properties() sums this quantity over the member molecules
     * of a complex on every timestep, and Dr is fixed for the whole run, so
     * evaluating pow() there repeats the same handful of results millions of
     * times.  The expression below is written exactly as update_properties()
     * used to write it so the cached double is bit-for-bit the value that used
     * to be computed inline.
     *
     * A zero Dr component is left as zero: update_properties() guards on
     * Dr != 0 and substitutes its own sentinel, so it never reads this value.
     */
    invCbrtDr.x = (Dr.x != 0) ? (1.0 / pow(Dr.x, 1.0 / 3.0)) : 0.0;
    invCbrtDr.y = (Dr.y != 0) ? (1.0 / pow(Dr.y, 1.0 / 3.0)) : 0.0;
    invCbrtDr.z = (Dr.z != 0) ? (1.0 / pow(Dr.z, 1.0 / 3.0)) : 0.0;
}

void MolTemplate::set_value(std::string& line, MolKeyword molKeyword)
{
    /*! \ingroup Parser
     * \brief sets the value of a MolTemplate's variable based on the keyword parsed from the parameter input file
     */

    auto key = static_cast<std::underlying_type<MolKeyword>::type>(molKeyword);
    switch (key) {
    case 0: {
        molName = line;
        break;
    }
    case 1: {
        copies = std::stoi(line);
        break;
    }
    case 2: {
        isRod = read_boolean(line);
        break;
    }
    case 3: {
        isLipid = read_boolean(line);
        std::cout << "Read in isLipid: " << std::boolalpha << isLipid << std::endl;
        break;
    }
    case 4: {
        D = Vec3D(parse_input_array(line));
        std::cout << "Read in D: [" << D.x << "um^2s^-1, " << D.y << "um^2s^-1, " << D.z << "um^2s^-1]" << std::endl;
        break;
    }
    case 5: {
        Dr = Vec3D(parse_input_array(line));
        cache_diffusion_derivatives();
        std::cout << "Read in Dr: [" << Dr.x << "rad^2s^-1, " << Dr.y << "rad^2s^-1, " << Dr.z << "rad^2s^-1]" << std::endl;
        break;
    }
    case 8: {
        mass = std::stod(line);
        std::cout << "Read in mass: " << std::boolalpha << mass << std::endl;
        break;
    }
    case 9: {
        checkOverlap = read_boolean(line);
        std::cout << "Read in checkOverlap: " << std::boolalpha << checkOverlap << std::endl;
        break;
    }
    case 11: {
        isPoint = read_boolean(line);
        break;
    }
    case 12: {
        isImplicitLipid = read_boolean(line);
        if (isImplicitLipid == true) {
            isLipid = true;
        }
        std::cout << "Read in isImplicitLipid: " << std::boolalpha << isImplicitLipid << std::endl;
        break;
    }
    case 13: {
        countTransition = read_boolean(line);
        std::cout << "Read in countTransition: " << std::boolalpha << countTransition << std::endl;
        break;
    }
    case 14: {
        transitionMatrixSize = std::stoi(line);
        break;
    }
    case 15: {
        outsideCompartment = read_boolean(line);
        if (outsideCompartment == true) {
            crossesCompartment = true;
        }
        break;
    }
    case 16: {
        insideCompartment = read_boolean(line);
        if (insideCompartment == true) {
            crossesCompartment = true;
        }
        break;
    }
    case 17: {
        isPromoter = read_boolean(line);
        std::cout << "Read in isPromoter: " << std::boolalpha << isPromoter << std::endl;
        break;
    }
    default: {
        std::cout << "Keyword [BLANK] is not a valid keyword, ignoring." << '\n';
    }
    }
}

/*******************/
