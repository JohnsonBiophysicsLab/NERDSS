#include "trajectory_functions/trajectory_functions.hpp"

void clear_reweight_vecs(Molecule& oneMol)
{
    /* Write reweight factor when exiting Rmax 
    We only care about the reweighting ratio when R > Rmax.
    To avoid counting when association has occured, 
    uncomment corresponding lines in nerdss.cpp
    by searching "Association happened, exit!"
    */
    // Kept in step with Molecule::ReweightEntry so that the instruction above --
    // uncomment this to get the reweighting factor -- still compiles.
    // for (const auto& prevEntry : oneMol.prevReweight) {
    //     // Search curr reweighting for the same item in prev reweighting
    //     bool found{false};
    //     for (const auto& currEntry : oneMol.currReweight){
    //         if (
    //             currEntry.partner == prevEntry.partner &&
    //             currEntry.myFace == prevEntry.myFace &&
    //             currEntry.partnerFace == prevEntry.partnerFace
    //         ){
    //             found = true;
    //             break;
    //         }
    //     }
    //     // If not found, write the previous reweighting factor
    //     if (!found){
    //         std::cout << "======= START REWEIGHTING INFO =======" << '\n';
    //         std::cout << prevEntry.partner << '\n';
    //         std::cout << prevEntry.myFace << '\n';
    //         std::cout << prevEntry.partnerFace << '\n';
    //         std::cout << prevEntry.norm << "\t # prev norm" << '\n';
    //         std::cout << prevEntry.survProb << "\t # prev SurvP" << '\n';
    //         std::cout << prevEntry.sep << '\n';
    //         std::cout << "======= END REWEIGHTING INFO =======" << std::endl;
    //         std::cout << "Exit for reweighting analysis with fixed initial separation." << std::endl;
    //         exit(0);
    //     }
    // }
    /* Curr values as prev values and clear curr values */
    // Six swaps and six clears became one of each when the twelve parallel
    // vectors were merged into Molecule::ReweightEntry.  It has to stay a swap
    // followed by clear(): clear() retains capacity, and that is what keeps the
    // sweep free of allocations in steady state (RESULTS.md section 16.2 put
    // the allocator at 0.94% of leaf samples).
    oneMol.prevReweight.swap(oneMol.currReweight);
    oneMol.currReweight.clear();
}
