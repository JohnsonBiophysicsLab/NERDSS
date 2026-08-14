/*! \file check_if_spans_box.cpp
 * ### Created on 2018-11-30 by Matthew Varga
 * ### TODO List
 * ***
 * can make this more efficient by just kicking out whenever cancelAssoc = true
 */
#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "tracing.hpp"

#include <cmath>
#include <iostream>

void check_if_spans_box(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Membrane& membraneObject)
{
    // TRACE();
    // Associating proteins have been moved to contact. Before assigning them to the complexsame complex,
    // test to see if the complex is too big to fit in the box.
    // declare the six boundary sides of the system box;
    // The reflecting surface is deliberately not considered here, because what is
    // being checked is whether the pair spans the box.
    const Vec3D posSide { membraneObject.waterBox.x / 2.0, membraneObject.waterBox.y / 2.0,
        membraneObject.waterBox.z / 2.0 };
    const Vec3D negSide { -membraneObject.waterBox.x / 2.0, -membraneObject.waterBox.y / 2.0,
        -membraneObject.waterBox.z / 2.0 };

    const double pairRadius { reactCom1.radius + reactCom2.radius };

    // Z first, then Y, then X: each axis shifts the tmp coordinates the next one reads.
    static const int axisOrder[3] = { 2, 1, 0 };

    for (int i { 0 }; i < 3; ++i) {
        const int axis { axisOrder[i] };
        const double posWallSide { axis_value(posSide, axis) };
        const double negWallSide { axis_value(negSide, axis) };

        // The approximate size of the pair puts it as outside; now test interface positions.
        if (!(pairRadius > posWallSide))
            continue;

        // Farthest tmp coordinate of either complex along this axis.
        WallExtent reach { -HUGE_VAL, HUGE_VAL };
        reach = scan_axis_extent_tmp_placed(reactCom1, moleculeList, axis, reach);
        reach = scan_axis_extent_tmp_placed(reactCom2, moleculeList, axis, reach);

        // How far the pair sticks out past each wall; seeded negative, so a pair
        // that is comfortably inside leaves both at the seed and trips nothing.
        double posWall { negWallSide };
        if (reach.posWall - posWallSide > posWall)
            posWall = reach.posWall - posWallSide;
        double negWall { negWallSide };
        if (negWallSide - reach.negWall > negWall)
            negWall = negWallSide - reach.negWall;

        // check whether to span the box
        const bool outsidePos { posWall > 0 };
        const bool outsideNeg { negWall > 0 };

        // translation or cancel
        if (outsideNeg && outsidePos)
            cancelAssoc = true;
        /*Also check if it sticks out far enough in one direction, that pushing back in will cause
          it to stick out the other side.
         */
        if (posWall + negWall > 0)
            cancelAssoc = true;

        // put back in the box. put at edge, rather than bouncing off.
        if (outsideNeg && !outsidePos) {
            translate_tmp_coords_along_axis(reactCom1, moleculeList, axis, -negWall);
            translate_tmp_coords_along_axis(reactCom2, moleculeList, axis, -negWall);
        }
        if (!outsideNeg && outsidePos) {
            translate_tmp_coords_along_axis(reactCom1, moleculeList, axis, posWall);
            translate_tmp_coords_along_axis(reactCom2, moleculeList, axis, posWall);
        }
    }
}
