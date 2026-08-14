#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "tracing.hpp"

void reflect_traj_complex_rad_rot_nocheck_box(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    double RS3D;
    if (targCom.OnSurface) {
        RS3D = 0;
    } else {
        RS3D = RS3Dinput;
    }

    std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);

    // declare the six boundary sides of the system box;
    const Vec3D posSide { membraneObject.waterBox.x / 2.0, membraneObject.waterBox.y / 2.0,
        membraneObject.waterBox.z / 2.0 };
    const Vec3D negSide { -membraneObject.waterBox.x / 2.0, -membraneObject.waterBox.y / 2.0,
        -membraneObject.waterBox.z / 2.0 + RS3D };

    const Vec3D curr { targCom.comCoord + targCom.trajTrans };

    /*Z is separate to allow the interfaces to approach to the membrane
     but don't need to test if the entire complex is far enough
     away from the boundary.
     */
    for (int axis { 0 }; axis < 3; ++axis) {
        const double posWallSide { axis_value(posSide, axis) };
        const double negWallSide { axis_value(negSide, axis) };
        const double base { axis_value(curr, axis) };

        if (!((base + targCom.radius) > posWallSide || (base - targCom.radius) < negWallSide))
            continue;

        const WallExtent extent { scan_axis_extent(
            targCom, moleculeList, M, axis, base, WallExtent { negWallSide, posWallSide }) };

        // Both can fire for a complex that reaches past either wall; this is not an else-if.
        if (extent.posWall > posWallSide)
            axis_ref(targCom.trajTrans, axis) -= 2.0 * (extent.posWall - posWallSide);
        if (extent.negWall < negWallSide)
            axis_ref(targCom.trajTrans, axis) -= 2.0 * (extent.negWall - negWallSide);
    }
}
