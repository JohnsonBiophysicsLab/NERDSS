#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <iostream>
#include <ostream>

bool determine_if_reaction_occurs(int& crossIndex1, int& crossIndex2, const double maxRandInt, Molecule& mol,
    std::vector<Molecule>& moleculeList, const std::vector<ForwardRxn>& forwardRxns)
{
    // TRACE();
    // bool willReact { false };
    for (unsigned crossMolItr { 0 }; crossMolItr < mol.crossings.size(); ++crossMolItr) {

        // make sure that the time step is resonable according to the prob of reaction
        // if (mol.crossings[crossMolItr].prob > 1.000001) {
        //     std::cerr << "Error: prob of reaction is: " << mol.crossings[crossMolItr].prob << " > 1. Avoid this using a smaller time step." << std::endl;
        //     exit(1);
        // }
        // if (mol.crossings[crossMolItr].prob > 0.5) {
        //     std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
        // }

        double rand1 { rand_gsl64() };
        // std::cout << mol.crossings[crossMolItr].prob << std::endl;
        // debug info
        // if (mol.crossings[crossMolItr].prob > 1e-3){
        //     std::cout << "Determining if reaction occurs for molecule index: " << mol.index
        //             << " on complex: " << mol.myComIndex << " crossMol: " << mol.crossings[crossMolItr].partner
        //             << " with reaction probability: " << mol.crossings[crossMolItr].prob
        //             << " random number: " << rand1 << std::endl;
        // }
        // END
        if (rand1 < mol.crossings[crossMolItr].prob) {

            crossIndex1 = crossMolItr;

            int mol2Index { mol.crossings[crossIndex1].partner };
            double pMatch { mol.crossings[crossIndex1].prob };
            /*Find the index on the partner protein's list of reactions that matches mol's.*/
            if (moleculeList[mol2Index].isImplicitLipid == true) {
                crossIndex2 = 0;
                return true;
                // std::cout << "DETERMINED REACTION TO OCCUR TO IL: " << mol.crossings[crossIndex1].rxn[0] << " prob: " << mol.crossings[crossMolItr].prob << " involving protein: " << mol.index << " on complex: " << mol.myComIndex << std::endl;
            } else {
                for (unsigned crossMolItr2 { 0 }; crossMolItr2 < moleculeList[mol2Index].crossings.size();
                     ++crossMolItr2) {
                    crossIndex2 = crossMolItr2;
                    bool isStateChangeBackRxn { mol.crossings[crossIndex1].rxn[2] == 1 };
                    int rxnItr { mol.crossings[crossIndex1].rxn[0] };
                    bool rxnMatches { false };
                    // check to make sure the reactant interfaces match the reaction. this is necessary if two reactions
                    // have the exact same probability (rare, but possible)
                    if (std::abs(moleculeList[mol2Index].crossings[crossMolItr2].prob - pMatch) < 1E-10
                        && mol.crossings[crossMolItr].rxn == moleculeList[mol2Index].crossings[crossMolItr2].rxn) {
                        if (mol.index == moleculeList[mol2Index].crossings[crossIndex2].partner) {
                            //proteins bind each other, using same reaction, and same probability. Will return crossIndex1 and

                            // corssIndex2 may be not desired if two interfaces have same reaction and same probability, for example cd1 + cd2

                            // determine mol is reactant1 or reactant2
                            if (mol.interfaceList[mol.crossings[crossIndex1].myIface].index == forwardRxns[rxnItr].reactantListNew[0].absIfaceIndex) {
                                //mol is reactant1
                                if (moleculeList[mol2Index].interfaceList[moleculeList[mol2Index].crossings[crossIndex2].myIface].index == forwardRxns[rxnItr].reactantListNew[1].absIfaceIndex)
                                    rxnMatches = true;
                            } else {
                                //mol is reactant2
                                if (moleculeList[mol2Index].interfaceList[moleculeList[mol2Index].crossings[crossIndex2].myIface].index == forwardRxns[rxnItr].reactantListNew[0].absIfaceIndex)
                                    rxnMatches = true;
                            }

                            if (rxnMatches == true) {
                                // std::cout << " DETERMINED WHICH REACTION OCCURED. "
                                //           << "Prbability: " << pMatch << " crossindex1: " << crossIndex1 << " crossindex2: " << crossIndex2 << " protein1: " << mol.index << " protein2: " << mol2Index << " check, protein2's partner: " << moleculeList[mol2Index].crossings[crossIndex2].partner << " isStateChangeBackRxn? 1 is true " << mol.crossings[crossIndex1].rxn[2] << std::endl;
                                return true;
                            }
                        }
                    }
                }
                return false;
            }
        }
    }
    return false;
}
