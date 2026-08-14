/*! \file resample_complex_trajectory.cpp
 * \brief One Brownian step for one complex, redrawn.
 *
 * Six `GaussV()` draws - translation x, y, z then rotation x, y, z - each
 * scaled by the complex's own diffusion constant for that axis.  Every routine
 * that rejects a trial move and tries again writes this out, and the overlap
 * sweeps write it twice, once for each complex of the pair.
 *
 * The draw order is the whole of the observable behaviour: the six values come
 * off one shared random stream, so anything that reorders them, or that skips a
 * draw on some path, changes every later number in the simulation.  Having the
 * order in one place is the point of this file.
 */
#include "classes/class_Molecule_Complex.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void resample_complex_trajectory(Complex& targCom, const Parameters& params)
{
    targCom.trajTrans.x = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    targCom.trajTrans.y = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    targCom.trajTrans.z = sqrt(2.0 * params.timeStep * targCom.D.z) * GaussV();
    targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
    targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
    targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();
}
