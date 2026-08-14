#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void reflect_traj_check_span_box(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    bool needsRecheck { true };
    int maxItr { 50 };
    int checkItr { 0 };

    double RS3D { reflecting_surface_offset(targCom, RS3Dinput) };

    // Built once, before the loop, and deliberately not rebuilt when the loop
    // resamples trajRot below: the resampled rotation is applied by
    // reflect_traj_complex_rad_rot_nocheck, which builds its own matrix.
    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    // declare the six boundary sides of the system box;
    const Vec3D posSide { membraneObject.waterBox.x / 2.0, membraneObject.waterBox.y / 2.0,
        membraneObject.waterBox.z / 2.0 };
    const Vec3D negSide { -membraneObject.waterBox.x / 2.0, -membraneObject.waterBox.y / 2.0,
        -membraneObject.waterBox.z / 2.0 + RS3D };

    while (checkItr < maxItr && needsRecheck) {
        needsRecheck = false; // without double span, this will stay 0

        bool moveFailed { false };
        bool outsideBox { false };
        bool outsidePos[3] { false, false, false };
        bool outsideNeg[3] { false, false, false };

        // these need to be what current positions
        // due to translation and rotation are
        BoxExtent extent;
        for (int axis { 0 }; axis < 3; ++axis)
            extent[axis] = WallExtent { axis_value(negSide, axis), axis_value(posSide, axis) };
        extent = scan_box_extent(targCom, moleculeList, M, targCom.comCoord + targCom.trajTrans, extent);

        for (int axis { 0 }; axis < 3; ++axis) {
            const double posWallSide { axis_value(posSide, axis) };
            const double negWallSide { axis_value(negSide, axis) };

            if (extent[axis].posWall > posWallSide) {
                outsidePos[axis] = true;
                outsideBox = true;
            }
            if (extent[axis].negWall < negWallSide) {
                outsideNeg[axis] = true;
                outsideBox = true;
            }

            // A complex that fills much of the box gets fewer attempts to find a
            // pose that fits.  Written per axis and in axis order, so the last
            // axis wide enough to trip a threshold is the one that sets the cap.
            const double span { std::abs(extent[axis].posWall - extent[axis].negWall) };
            const double boxSpan { std::abs(posWallSide - negWallSide) };
            if (span > 1.0 / 2.0 * boxSpan)
                maxItr = 20;
            if (span > 2.0 / 3.0 * boxSpan)
                maxItr = 10;
            if (span > 4.0 / 5.0 * boxSpan)
                maxItr = 5;
        }

        if (outsideBox) {
            for (int axis { 0 }; axis < 3; ++axis) {
                const double posWallSide { axis_value(posSide, axis) };
                const double negWallSide { axis_value(negSide, axis) };

                if (outsideNeg[axis] && outsidePos[axis]) {
                    // For a large complex, test if checkItr could be pushed back out the other side
                    moveFailed = true;
                }
                if (outsideNeg[axis] && !outsidePos[axis]) {
                    const double shift { 2.0 * (extent[axis].negWall - negWallSide) };
                    axis_ref(targCom.trajTrans, axis) -= shift;
                    // Also need to check that update will not push you out the other side
                    if (extent[axis].posWall - shift > posWallSide)
                        moveFailed = true;
                }
                if (!outsideNeg[axis] && outsidePos[axis]) {
                    const double shift { 2.0 * (extent[axis].posWall - posWallSide) };
                    axis_ref(targCom.trajTrans, axis) -= shift;
                    if (extent[axis].negWall - shift < negWallSide)
                        moveFailed = true;
                }
            }
        } // recheck span

        if (moveFailed == true) {
            // Resample, extends in x, y, and/or z
            resample_complex_trajectory(targCom, params);

            reflect_traj_complex_rad_rot_nocheck(params, targCom, moleculeList, membraneObject, RS3Dinput);
            ++checkItr;
            needsRecheck = true; // will need to recheck after resampling traj and trajR
        }
    } // loop over iterations and flag condition

    // A failure to converge used to walk the complex again to print the positions
    // it ended up with; every one of those prints is commented out, so the walk
    // did nothing.  The caller cannot see the outcome either way.
}
