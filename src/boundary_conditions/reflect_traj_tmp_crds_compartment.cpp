#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "tracing.hpp"

/*evaluates reflection into the compartment based on tmpCoords of all proteins. targCOM tmpCOM coords also must be consistent-updated.
 *Does not update molecule traj vectors, just updates the passed in vector traj, to collect for both complexes.
 *
*/
void reflect_traj_tmp_crds_compartment(
    const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    /*This routine updated April 2020 to test if a large complex that spans the sphere could extend out in both directions
    if so, it attempts to correct for this by resampling the complex's translational updates.
    */
    double RS3D { reflecting_surface_offset_tmp(targCom, RS3Dinput) };
    double sphereR = membraneObject.compartmentR + RS3D;

    const Vec3D base { targCom.tmpComCoord.x + traj[0], targCom.tmpComCoord.y + traj[1],
        targCom.tmpComCoord.z + traj[2] };

    // Identity: this routine performs no rotation, it only reflects the trial translation.
    std::array<double, 9> M {};
    M[0] = 1;
    M[4] = 1;
    M[8] = 1;

    /*This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.*/
    if (!((base.length() + targCom.radius) < sphereR))
        return;

    /*Now evaluate all interfaces distance from boundaries.*/
    // The compartment excludes rather than contains, so the point of interest is the
    // one that reaches deepest inside it.  Scored as a negated depth so the same
    // "largest score wins" accumulator applies; negation is exact.
    ExtremePoint deepest { Vec3D {}, 0.0 };
    for_each_complex_tmp_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D rot { matrix_rotate(point - targCom.tmpComCoord, M) };
        const Vec3D curr { base + rot };
        deepest.consider(curr, -(curr.length() - sphereR));
    });

    if (deepest.score > 0.0) {
        // Put back outside the compartment
        const double targR { deepest.point.length() };
        double lamda = -2.0 * (targR - sphereR) / targR;
        traj[0] = lamda * deepest.point.x;
        traj[1] = lamda * deepest.point.y;
        traj[2] = lamda * deepest.point.z;
    }
}
