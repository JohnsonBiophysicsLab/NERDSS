#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

/*
  In this _complex_ version, itr tests overlap not just for each protein, but for each complex, so all the proteins
 in a complex, before performing position updates.

  NEW 2018: IN THIS VERSION, IT DOES NOT ATTEMPT TO SOLVE OVERLAP FOR PROTEINS WITHIN THE SAME COMPLEX, SINCE THEY
 CANNOT DIFFUSE RELATIVE TO ONE ANOTHER!

  This routine checks whether protein p1 is overlapping any partners in its reaction
  zone at its new position that is given by its current position +traj. If itr does
 overlap, the displacement by traj is rejected and a new position for itself and any overlapping partners are
 selected. Once itr no longer overlaps anyone, this protein and its complex are moved and the partners retain their
 stored new displacements. If a protein has already updated its position (done in sequential order) then itr cannot
 resample a new position, the current protein must still continue to avoid overlapping, however.
*/
namespace {

/*! \brief Body shared by the solution and membrane box sweeps.
 *
 * `sweep_separation_complex_rot_box` and
 * `sweep_separation_complex_rot_memtest_box` were two files that differed in
 * exactly one expression: how the pair separation is measured.  In solution it
 * is the full 3D distance; on the membrane, when the *partner* complex is also
 * on the surface, the z term is dropped so that only the xy displacement is
 * judged - two complexes both stuck to the membrane cannot separate in z, so
 * counting z would report an overlap no resampling could fix.
 *
 * \param[in] ignoreZBetweenSurfacePartners the membrane behaviour above.  When
 * false the z term is always included, which is the same sum the solution
 * version wrote as one expression: `(a + b) + c` either way.
 */
void sweep_separation_box(int simItr, int pro1Index, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject,
    bool ignoreZBetweenSurfacePartners)
{
    // TRACE();
    int com1Index { moleculeList[pro1Index].myComIndex };
    size_t com1Size { complexList[com1Index].memberList.size() };

    int maxRows { 1 };
    for (auto memMol : complexList[com1Index].memberList) {
        if (moleculeList[memMol].crossbase.size() > maxRows)
            maxRows = moleculeList[memMol].crossbase.size();
    }

    int ifaceList[maxRows * com1Size];
    int overlapList[maxRows * com1Size];
    // Partner-is-on-the-surface flags; only read when ignoreZBetweenSurfacePartners.
    int memCheckList[maxRows * com1Size];

    /*The sampled displacement for p1 is stored in traj. the position from the
     previous step is still stored in bases[p1].xcom, etc, and will be updated
     at the end of this routine*/

    /*figure out i2*/
    for (int memMolItr { 0 }; memMolItr < com1Size; ++memMolItr) {
        pro1Index = complexList[com1Index].memberList[memMolItr];
        for (int i { 0 }; i < moleculeList[pro1Index].crossbase.size(); ++i) {
            if (ignoreZBetweenSurfacePartners) {
                int pro2Index { moleculeList[pro1Index].crossbase[i] };
                int com2Index { moleculeList[pro2Index].myComIndex };
                // if (complexList[com2Index].D.z < 1E-15) {
                memCheckList[maxRows * memMolItr + i] = complexList[com2Index].OnSurface ? 1 : 0;
            }

            int i1 { moleculeList[pro1Index].mycrossint[i] };
            std::array<int, 3> rxnItr = moleculeList[pro1Index].crossrxn[i];

            // get the partner interface
            ifaceList[maxRows * memMolItr + i] = (forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex == i1)
                ? forwardRxns[rxnItr[0]].reactantListNew[1].relIfaceIndex
                : forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex;
        }
    }

    // determine RS3Dinput
    double RS3Dinput { 0.0 };
    Complex targCom { complexList[com1Index] };
    if (membraneObject.implicitLipid) {
        for (auto& molIndex : targCom.memberList) {
            for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
                if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - moleculeList[molIndex].molTypeIndex) < 1E-2) {
                    RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
                    break;
                }
            }
        }
    }

    int index = moleculeList[pro1Index].molTypeIndex;
    bool isInsideCompartment = molTemplateList[index].insideCompartment;

    int itr { 0 };
    int maxItr { 10 };
    while (itr < maxItr) {
        int numOverlap { 0 };
        bool hasOverlap { false };
        int com2Index {};

        for (unsigned memMolItr { 0 }; memMolItr < com1Size; ++memMolItr) {
            pro1Index = complexList[com1Index].memberList[memMolItr];
            for (int crossMolItr { 0 }; crossMolItr < moleculeList[pro1Index].crossbase.size(); ++crossMolItr) {
                int pro2Index { moleculeList[pro1Index].crossbase[crossMolItr] };
                if (moleculeList[pro2Index].isImplicitLipid)
                    continue;

                com2Index = moleculeList[pro2Index].myComIndex;
                /*Do not sweep for overlap if proteins are in the same complex, they cannot move relative to one
                 * another!
                 */
                if (com1Index != com2Index) {
                    int i1 { moleculeList[pro1Index].mycrossint[crossMolItr] };
                    int rxnItr { moleculeList[pro1Index].crossrxn[crossMolItr][0] };
                    int i2 { ifaceList[maxRows * memMolItr + crossMolItr] };

                    Vec3D iface1Vec { moleculeList[pro1Index].interfaceList[i1].coord - complexList[com1Index].comCoord };
                    std::array<double, 9> M = create_euler_rotation_matrix(complexList[com1Index].trajRot);
                    iface1Vec = matrix_rotate(iface1Vec, M);

                    double dx1 { complexList[com1Index].comCoord.x + iface1Vec.x + complexList[com1Index].trajTrans.x };
                    double dy1 { complexList[com1Index].comCoord.y + iface1Vec.y + complexList[com1Index].trajTrans.y };
                    double dz1 { complexList[com1Index].comCoord.z + iface1Vec.z + complexList[com1Index].trajTrans.z };

                    /*Now complex 2*/
                    Vec3D iface2Vec { moleculeList[pro2Index].interfaceList[i2].coord - complexList[com2Index].comCoord };
                    std::array<double, 9> M2 = create_euler_rotation_matrix(complexList[com2Index].trajRot);
                    iface2Vec = matrix_rotate(iface2Vec, M2);
                    double dx2 { complexList[com2Index].comCoord.x + iface2Vec.x + complexList[com2Index].trajTrans.x };
                    double dy2 { complexList[com2Index].comCoord.y + iface2Vec.y + complexList[com2Index].trajTrans.y };
                    double dz2 { complexList[com2Index].comCoord.z + iface2Vec.z + complexList[com2Index].trajTrans.z };

                    /*separation*/
                    double df1 { dx1 - dx2 };
                    double df2 { dy1 - dy2 };
                    double df3 { dz1 - dz2 };

                    double dr2 = (df1 * df1) + (df2 * df2);
                    if (!(ignoreZBetweenSurfacePartners && memCheckList[maxRows * memMolItr + crossMolItr] == 1))
                        dr2 += (df3 * df3);

                    if (dr2 < forwardRxns[rxnItr].bindRadius * forwardRxns[rxnItr].bindRadius) {
                        /*reselect positions for protein pro2Index*/
                        overlapList[numOverlap] = pro2Index;
                        ++numOverlap;
                        hasOverlap = true;
                    }
                } // ignore proteins within the same complex.
            }
        }
        /*Now resample positions of p1 and overlapList, if numOverlap>0, otherwise no overlap, so
         break from loop*/
        if (!hasOverlap)
            break;

        ++itr;
        complexList[com1Index].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[com1Index].D.x) * GaussV();
        complexList[com1Index].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[com1Index].D.y) * GaussV();
        complexList[com1Index].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[com1Index].D.z) * GaussV();
        complexList[com1Index].trajRot.x = sqrt(2.0 * params.timeStep * complexList[com1Index].Dr.x) * GaussV();
        complexList[com1Index].trajRot.y = sqrt(2.0 * params.timeStep * complexList[com1Index].Dr.y) * GaussV();
        complexList[com1Index].trajRot.z = sqrt(2.0 * params.timeStep * complexList[com1Index].Dr.z) * GaussV();

        reflect_traj_complex_rad_rot(params, moleculeList, complexList[com1Index], membraneObject, RS3Dinput, isInsideCompartment);

        int resampleList[complexList.size()]; // if this is 0, we need resample
        for (auto& i : resampleList) { // initialize
            i = 0;
        }

        for (int checkMolItr { 0 }; checkMolItr < numOverlap; checkMolItr++) {
            int p2 { overlapList[checkMolItr] };
            com2Index = moleculeList[p2].myComIndex;

            // avoid repeated resample
            if (resampleList[com2Index] == 0) {
                if (moleculeList[p2].trajStatus == TrajStatus::none || moleculeList[p2].trajStatus == TrajStatus::canBeResampled) {
                    /*
                 We loop over proteins sequentially, so earlier proteins have already moved and avoided
                 their neighbors and should not be moved again.
                 These new positions selected for proteins not yet moved will be stored and
                 then used when they test for overlap themselves.
                 */

                    /*If p2 just dissociated, also don't try to move again*/
                    complexList[com2Index].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[com2Index].D.x) * GaussV();
                    complexList[com2Index].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[com2Index].D.y) * GaussV();
                    complexList[com2Index].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[com2Index].D.z) * GaussV();
                    complexList[com2Index].trajRot.x = sqrt(2.0 * params.timeStep * complexList[com2Index].Dr.x) * GaussV();
                    complexList[com2Index].trajRot.y = sqrt(2.0 * params.timeStep * complexList[com2Index].Dr.y) * GaussV();
                    complexList[com2Index].trajRot.z = sqrt(2.0 * params.timeStep * complexList[com2Index].Dr.z) * GaussV();

                    reflect_traj_complex_rad_rot(params, moleculeList, complexList[com2Index], membraneObject, RS3Dinput, isInsideCompartment);
                    resampleList[com2Index] = 1;
                }
            }
        }
    } // end maximum iterations

    // Both versions ended with a "can't solve overlap" report at itr == maxItr - 1
    // whose every print was commented out, leaving a walk over the cross list that
    // wrote nothing and read nothing the caller can see.

    complexList[com1Index].propagate(moleculeList, membraneObject, molTemplateList);

    // Reset displacements to zero so distance is measured to your current
    // updated position that won't change again this turn
    complexList[com1Index].trajTrans.zero();
    complexList[com1Index].trajRot.zero();
}

} // namespace

//! \brief Overlap sweep for a complex in solution: separations are measured in 3D.
void sweep_separation_complex_rot_box(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    sweep_separation_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList,
        membraneObject, false);
}

//! \brief Overlap sweep for a complex on the membrane: z is ignored against partners that are also on it.
void sweep_separation_complex_rot_memtest_box(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    sweep_separation_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList,
        membraneObject, true);
}
