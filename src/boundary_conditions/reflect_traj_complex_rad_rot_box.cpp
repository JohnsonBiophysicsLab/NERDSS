#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "tracing.hpp"

void reflect_traj_complex_rad_rot_box(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    // NOTE: it only works for a box system with the membrane surface located on the Z-bottom.

    // for implicit-lipid model, the boundary surface must consider the reflecting-surface RS3D
    // for explicit-lipid model, membraneObject.RS3D = 0. And, for those proteins that are bound on surface, they are not allowed to reflect along Z-axis
    double RS3D { reflecting_surface_offset(targCom, RS3Dinput) };

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    // declare the six boundary sides of the system box;
    const Vec3D posSide { membraneObject.waterBox.x / 2.0, membraneObject.waterBox.y / 2.0,
        membraneObject.waterBox.z / 2.0 };
    const Vec3D negSide { -membraneObject.waterBox.x / 2.0, -membraneObject.waterBox.y / 2.0,
        -membraneObject.waterBox.z / 2.0 + RS3D };

    // his routine updated March 2017 to test if a large complex that spans the box could extend out in both directions
    // if so, it attempts to correct for this by resampling the complex's translational and rotational updates.
    const Vec3D curr { targCom.comCoord + targCom.trajTrans };

    // This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.
    bool recheck { false };
    for (int axis { 0 }; axis < 3; ++axis) {
        const double posWallSide { axis_value(posSide, axis) };
        const double negWallSide { axis_value(negSide, axis) };
        const double base { axis_value(curr, axis) };

        if (!((base + targCom.radius) > posWallSide || (base - targCom.radius) < negWallSide))
            continue;

        // Now evaluate all interfaces distance from boundaries.
        const WallExtent extent { scan_axis_extent(
            targCom, moleculeList, M, axis, base, WallExtent { negWallSide, posWallSide }) };

        // Largest distance outside the positive side (marked +) and the negative side (marked -).
        const double posd { extent.overshoot(posWallSide) };
        const double negd { extent.undershoot(negWallSide) };
        const bool outsidePos { posd != 0.0 };
        const bool outsideNeg { negd != 0.0 };
        if (!outsidePos && !outsideNeg)
            continue;

        // Put back inside the box, extended out
        if (outsidePos)
            axis_ref(targCom.trajTrans, axis) -= 2.0 * posd;
        else if (outsideNeg)
            axis_ref(targCom.trajTrans, axis) -= 2.0 * negd;

        if (outsideNeg && outsidePos) {
            // For a large complex, test if it could be pushed back out the other side
            recheck = true;
        } else if (outsideNeg && extent.posWall - 2.0 * negd > posWallSide) {
            // Also need to check that update will not push you out the other side
            recheck = true;
        } else if (outsidePos && extent.negWall - 2.0 * posd < negWallSide) {
            recheck = true;
        }
    }

    if (recheck) {
        // Test that new coordinates have not pushed you out of the box for a very large complex,
        // if so, resample rotation matrix.
        reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, RS3Dinput);
    }
}
