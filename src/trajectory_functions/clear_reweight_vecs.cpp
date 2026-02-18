#include "trajectory_functions/trajectory_functions.hpp"

void clear_reweight_vecs(Molecule& oneMol)
{
    /* Write reweight factor when exiting Rmax 
    We only care about the reweighting ratio when R > Rmax.
    To avoid counting when association has occured, 
    uncomment corresponding lines in nerdss.cpp
    by searching "Association happened, exit!"
    */
    // for (int previd{0}; previd < oneMol.prevlist.size(); ++previd) {
    //     // Search curr reweighting for the same itme in prev reweighting
    //     bool found{false}; 
    //     for (int currid{0}; currid < oneMol.currlist.size(); ++currid){
    //         if (
    //             oneMol.currlist[currid] == oneMol.prevlist[previd] &&
    //             oneMol.currmyface[currid] == oneMol.prevmyface[previd] &&
    //             oneMol.currpface[currid] == oneMol.prevpface[previd]
    //         ){
    //             found = true;
    //             break;
    //         }   
    //     }
    //     // If not found, write the previous reweighting factor
    //     if (!found){
    //         std::cout << "======= START REWEIGHTING INFO =======" << '\n';
    //         std::cout << oneMol.prevlist[previd] << '\n';
    //         std::cout << oneMol.prevmyface[previd] << '\n';
    //         std::cout << oneMol.prevpface[previd] << '\n';
    //         std::cout << oneMol.prevnorm[previd] << "\t # prev norm" << '\n';
    //         std::cout << oneMol.ps_prev[previd] << "\t # prev SurvP" << '\n';
    //         std::cout << oneMol.prevsep[previd] << '\n';
    //         std::cout << "======= END REWEIGHTING INFO =======" << std::endl;
    //         std::cout << "Exit for reweighting analysis with fixed initial separation." << std::endl;
    //         exit(0);
    //     }
    // }
    /* Curr values as prev values and clear curr values */
    oneMol.prevlist = oneMol.currlist;
    oneMol.prevmyface = oneMol.currmyface;
    oneMol.prevpface = oneMol.currpface;
    oneMol.prevnorm = oneMol.currprevnorm;
    oneMol.ps_prev = oneMol.currps_prev;
    oneMol.prevsep = oneMol.currprevsep;
    oneMol.currlist.clear();
    oneMol.currmyface.clear();
    oneMol.currpface.clear();
    oneMol.currprevnorm.clear();
    oneMol.currps_prev.clear();
    oneMol.currprevsep.clear();
}
