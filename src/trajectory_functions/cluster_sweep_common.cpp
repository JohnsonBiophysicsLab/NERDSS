/*! \file cluster_sweep_common.cpp
 * \brief The two loops the cluster overlap sweeps repeat four times over.
 *
 * When a cluster is too big to resolve in one timestep, the sweep splits the
 * step into substeps.  Doing that means the complexes *around* the cluster have
 * to be kept in step with it: their partners are collected once, and their
 * trajectories are resampled at the new, shorter timestep - once when the
 * substepping is set up, and again before each substep after the first.
 *
 * Both loops were written out twice in each of
 * sweep_separation_complex_rot_memtest_cluster_box.cpp and
 * _sphere.cpp, so each existed four times.  They are the same loops, and the
 * copies had not drifted apart except in one dead spot noted below.
 *
 * The rest of those two files is *not* merged.  They differ in five live
 * places - the substep count, whether an interface on the membrane is placed by
 * rotation or by calculate_update_position_interface(), whether the
 * failed-to-converge warning prints, which pairwise sweep they fall back to,
 * and whether the trailing trajTrans/trajRot are zeroed - and a body carrying
 * five flags would be harder to read than the two bodies are.
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Cluster.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void collect_cluster_partners(int comIndex, const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList, std::vector<int>& partnerList,
    std::vector<TrajStatus>& movestatOrigPartnerList)
{
    for (unsigned memMolItr { 0 }; memMolItr < complexList[comIndex].memberList.size(); ++memMolItr) {
        int pro1Index = complexList[comIndex].memberList[memMolItr];
        for (int crossMemItr { 0 }; crossMemItr < moleculeList[pro1Index].crossbase.size(); ++crossMemItr) {
            int pro2Index { moleculeList[pro1Index].crossbase[crossMemItr] };
            if (moleculeList[pro2Index].isImplicitLipid)
                continue;
            partnerList.push_back(pro2Index);
            movestatOrigPartnerList.push_back(moleculeList[pro2Index].trajStatus);
        }
    }
}

void resample_partner_trajectories(const std::vector<int>& partnerList, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const Parameters& params, const Membrane& membraneObject, double RS3Dinput)
{
    // One resample per complex, not per molecule: partnerList holds molecules,
    // and several of them can belong to the same complex.
    std::vector<int> didMove;

    for (int i = 0; i < partnerList.size(); i++) {
        int p = partnerList[i];
        int k = moleculeList[p].myComIndex;
        if (moleculeList[p].trajStatus == TrajStatus::none || moleculeList[p].trajStatus == TrajStatus::canBeResampled) {
            int flag = 0;
            for (int d = 0; d < didMove.size(); d++) {
                if (didMove[d] == k)
                    flag = 1;
            }
            if (flag == 0) {

                // if (membraneObject.isSphere == true && complexList[k].D.z < 1E-15) { // complex on sphere surface
                if (membraneObject.isSphere == true && complexList[k].OnSurface) { // complex on sphere surface
                    // The box copies of this branch drew the translation from
                    // complexList[k1] - an index left over from an enclosing loop
                    // rather than the complex being resampled.  It cannot have
                    // mattered: the box sweep is only ever reached from
                    // sweep_separation_complex_rot_memtest_cluster() under
                    // !membraneObject.isSphere, so this branch never runs there.
                    Vec3D targTrans = create_complex_propagation_vectors_on_sphere(params, complexList[k]);
                    complexList[k].trajTrans.x = targTrans.x;
                    complexList[k].trajTrans.y = targTrans.y;
                    complexList[k].trajTrans.z = targTrans.z;
                    complexList[k].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k].Dr.x) * GaussV();
                    complexList[k].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k].Dr.y) * GaussV();
                    complexList[k].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k].Dr.z) * GaussV();
                } else {
                    complexList[k].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[k].D.x) * GaussV();
                    complexList[k].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[k].D.y) * GaussV();
                    complexList[k].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[k].D.z) * GaussV();
                    complexList[k].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k].Dr.x) * GaussV();
                    complexList[k].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k].Dr.y) * GaussV();
                    complexList[k].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k].Dr.z) * GaussV();
                }

                reflect_traj_complex_rad_rot(params, moleculeList, complexList[k], membraneObject, RS3Dinput, false);

                didMove.push_back(k);
            }
        }
    }
}
