/*! \file reflect_traj_complex_radial.cpp
 * \brief Reflects a complex's sampled translation back across a spherical
 *        boundary, on whichever side of it the complex must stay.
 *
 * This was `reflect_traj_complex_rad_rot_sphere` and
 * `reflect_traj_complex_compartment`.  As with the tmp-coordinate pair, the
 * shared body is the whole routine and the difference is one sign, carried here
 * by \ref RadialSide.
 *
 * Two things genuinely differ between the two, and both are arguments rather
 * than sign flips, so they stay visible at the call site:
 *
 *   * **The surface complexes.**  The sphere version skipped a complex with
 *     `OnSurface` set: such a complex lives *on* the membrane and is propagated
 *     by a different routine, so reflecting it off the membrane would be wrong.
 *     The compartment version had no such test, because a complex bound to the
 *     compartment surface is not tracked by `OnSurface`.  `skipOnSurface` keeps
 *     that split.
 *
 *   * **The span re-check.**  A reflection can push a large complex out through
 *     the far side of the boundary, so the sphere version re-samples by calling
 *     \ref reflect_traj_check_span_sphere.  The compartment version never did;
 *     its comment reads "assume the complex is not huge enough".  That
 *     assumption is unchecked, and this is the point in the plan where it would
 *     be revisited - but doing so changes results, so `recheckSpan` preserves
 *     it for now.  See docs/debranching_plan.md.
 *
 * ### Created on 02/25/2020 by Yiben Fu, as the two files this replaces.
 */
#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "tracing.hpp"

void reflect_traj_complex_radial(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom,
    const Membrane& membraneObject, double radius, RadialSide side, double RS3Dinput, bool skipOnSurface,
    bool recheckSpan)
{
    // TRACE();

    // for implicit-lipid model, the boundary surface must consider the reflecting-surface RS3D
    // for explicit-lipid model, membraneObject.RS3D = 0. And, for those proteins that are bound on surface, they are not allowed to reflect along Z-axis
    // for the complex on the surface, its propagation is different from those inside the sphere
    const double RS3D { reflecting_surface_offset(targCom, RS3Dinput) };
    const double boundaryR { radial_boundary(radius, RS3D, side) };

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    bool recheck = false;

    // if (targCom.D.z < 1E-14 || targCom.OnSurface) { // for the complex on the sphere surface
    if (!(skipOnSurface && targCom.OnSurface)) {
        /*calculate distance to the centre of the boundary. */
        /*first just test Complex COM+ radius, if it fits on the required side.
          if yes, then test if all interfaces do.
          then, if they do not, what is the max displacement past the boundary.
          move it along the radial direction back by that displacement.
         */
        /*assume the origin of the sphere is at zero. */
        const Vec3D base { targCom.comCoord + targCom.trajTrans };

        if (radial_may_escape(base, targCom.radius, boundaryR, side)) {
            /*Now evaluate all molecules and interfaces distance from boundaries.*/
            // Ranked by signed radius and seeded at the boundary, which is the
            // convention both of the routines this replaces were written with.
            // Not radial_escape(): see radial_signed_radius() for why the two
            // are not interchangeable.
            const double seed { radial_signed_boundary(boundaryR, side) };
            ExtremePoint escaped { Vec3D {}, seed };
            for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
                const Vec3D rot { matrix_rotate(point - targCom.comCoord, M) };
                const Vec3D curr { base + rot };
                escaped.consider(curr, radial_signed_radius(curr, side));
            });

            if (escaped.score > seed) {
                recheck = true;
                targCom.trajTrans = Vec3D(targCom.trajTrans + radial_reflection_shift(escaped.point, boundaryR));
            }
        }
    }

    if (recheck && recheckSpan) {
        // Test that new coordinates have not pushed the complex out through the
        // far side of the boundary; if so, resample the rotation matrix.
        reflect_traj_check_span_sphere(params, targCom, moleculeList, membraneObject, radius, RS3Dinput);
    }
}
