#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "tracing.hpp"

/*evaluates reflection out of box based on tmpCoords of all proteins. targCOM tmpCOM coords also must be consistent-updated.
 *Does not update molecule traj vectors, just updates the passed in vector traj, to collect for both complexes.
 *
*/
void reflect_traj_tmp_crds_box(
    const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    /*This routine updated October 2019 to test if a large complex that spans the box could extend out in both directions
    if so, it attempts to correct for this by resampling the complex's translational and rotational updates.
    */
    double RS3D { reflecting_surface_offset_tmp(targCom, RS3Dinput) };

    // declare the six boundary sides of the system box;
    const Vec3D posSide { membraneObject.waterBox.x / 2.0, membraneObject.waterBox.y / 2.0,
        membraneObject.waterBox.z / 2.0 };
    const Vec3D negSide { -membraneObject.waterBox.x / 2.0, -membraneObject.waterBox.y / 2.0,
        -membraneObject.waterBox.z / 2.0 + RS3D };

    // Identity: this routine performs no rotation, it only reflects the trial translation.
    std::array<double, 9> M {};
    M[0] = 1;
    M[4] = 1;
    M[8] = 1;

    /*This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.*/
    for (int axis { 0 }; axis < 3; ++axis) {
        const double posWallSide { axis_value(posSide, axis) };
        const double negWallSide { axis_value(negSide, axis) };
        const double base { axis_value(targCom.tmpComCoord, axis) + traj[axis] };

        if (!((base + targCom.radius) > posWallSide || (base - targCom.radius) < negWallSide))
            continue;

        /*Now evaluate all interfaces distance from boundaries.*/
        const WallExtent extent { scan_axis_extent_tmp(
            targCom, moleculeList, M, axis, base, WallExtent { negWallSide, posWallSide }) };

        // check whether this complex is out of the box
        const bool outsidePos { extent.posWall > posWallSide };
        const bool outsideNeg { extent.negWall < negWallSide };

        // Put back inside the box.  A complex that reaches past both walls, or that
        // would be pushed out the far side by the correction, is left where it is:
        // reflect_traj_complex_rad_rot_box re-checks the span in that case, but this
        // tmp-coordinate routine never has.
        if (outsideNeg && !outsidePos)
            traj[axis] -= 2.0 * (extent.negWall - negWallSide);
        if (outsidePos && !outsideNeg)
            traj[axis] -= 2.0 * (extent.posWall - posWallSide);
    }
}
