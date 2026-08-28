/*! \file reflect_traj_tmp_crds_radial.cpp
 * \brief Reflects a trial translation back across a spherical boundary, on
 *        whichever side of it the complex is required to stay.
 *
 * This was two files.  `reflect_traj_tmp_crds_sphere` kept a complex inside the
 * membrane and `reflect_traj_tmp_crds_compartment` kept it outside the
 * compartment, and apart from their names they differed in exactly three
 * expressions, each of which was the other with a sign flipped: the boundary
 * radius (`radius - RS3D` against `compartmentR + RS3D`), the bounding test
 * (`>` against `<`), and the score a point was ranked by.  Everything after
 * that - the identity rotation, the scan over tmp coordinates, the `score > 0`
 * test and the whole reflection tail - was identical, character for character.
 *
 * \ref RadialSide names the one thing that actually differs.  The reflection
 * tail needs no sign at all: see \ref radial_reflection_shift.
 *
 * ### Created on 02/25/2020 by Yiben Fu, as the two files this replaces.
 */
#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "tracing.hpp"

/*evaluates reflection across a spherical boundary based on tmpCoords of all proteins. targCOM tmpCOM coords also must be consistent-updated.
 *Does not update molecule traj vectors, just updates the passed in vector traj, to collect for both complexes.
 */
void reflect_traj_tmp_crds_radial(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom,
    std::array<double, 3>& traj, double radius, RadialSide side, double RS3Dinput)
{
    // TRACE();
    /*This routine updated April 2020 to test if a large complex that spans the sphere could extend out in both directions
    if so, it attempts to correct for this by resampling the complex's translational updates.
    */
    const double RS3D { reflecting_surface_offset_tmp(targCom, RS3Dinput) };
    const double boundaryR { radial_boundary(radius, RS3D, side) };

    const Vec3D base { targCom.tmpComCoord.x + traj[0], targCom.tmpComCoord.y + traj[1],
        targCom.tmpComCoord.z + traj[2] };

    // Identity: this routine performs no rotation, it only reflects the trial translation.
    std::array<double, 9> M {};
    M[0] = 1;
    M[4] = 1;
    M[8] = 1;

    /*This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.*/
    if (!radial_may_escape(base, targCom.radius, boundaryR, side))
        return;

    /*Now evaluate all interfaces distance from boundaries.*/
    // Scored by how far the point has escaped past the boundary, seeded at zero,
    // so a complex entirely on the correct side leaves the seed in place.
    ExtremePoint escaped { Vec3D {}, 0.0 };
    for_each_complex_tmp_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D rot { matrix_rotate(point - targCom.tmpComCoord, M) };
        const Vec3D curr { base + rot };
        escaped.consider(curr, radial_escape(curr, boundaryR, side));
    });

    if (escaped.score > 0.0) {
        // Put back on the required side of the boundary
        const Vec3D shift { radial_reflection_shift(escaped.point, boundaryR) };
        traj[0] = shift.x;
        traj[1] = shift.y;
        traj[2] = shift.z;
    }
}
