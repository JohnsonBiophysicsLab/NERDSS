/*! \file functions_for_spherical_system.hpp
 *
 * ### Created on 2/05/2020 by Yiben Fu
 *
 * \brief Geometry for a system whose membrane *is* a sphere.
 *
 * \section sphere_on_not_in On the sphere, not inside it
 *
 * Every function here treats the sphere as a two-dimensional world: the
 * molecules it moves live *on the spherical surface*, and their motion is
 * confined to it.  This is the curved counterpart of the flat membrane, and the
 * routines are the curved counterparts of the flat ones -- a step is an arc
 * along a great circle rather than a straight line, a binding radius is a
 * geodesic distance rather than a chord, and a rotation is about the local
 * outward normal rather than about a fixed axis.
 *
 * That is a different thing from a spherical *container*, which is what
 * "spherical system" reads as at first and what the reflecting boundary
 * conditions in reflect_functions.hpp actually implement: there, molecules
 * occupy the enclosed volume and the sphere only turns them back at the wall.
 * Nothing in this file will do anything sensible with an interior coordinate.
 * The radius the functions work with is not a parameter they are given, it is
 * `coordinate.length()` -- read off the coordinate itself, on the assumption
 * that it already sits on the surface.  Hand one an interior point and it will
 * silently answer for the smaller sphere that point happens to lie on.
 *
 * `Membrane::isSphere` is the flag that selects these routines; the callers are
 * the `_sphere` variants in reactions/ and trajectory_functions/.
 *
 * \section sphere_frames The two coordinate systems, and which is which
 *
 * Two different triples of doubles appear below and they are not interchangeable:
 *
 *   * \ref Vec3D is cartesian, always.  x, y, z are lengths in nm from the
 *     centre of the sphere, which is the origin.
 *   * \ref SphericalCoord is spherical.  Its members are named, so a polar angle
 *     cannot be read as an x coordinate.
 *
 * These functions used to pass spherical coordinates *as* `Vec3D`, packing
 * `(theta, phi, r)` into `(x, y, z)`.  Nothing marked which convention a given
 * `Vec3D` held; the two were told apart only by a trailing comment on the
 * declaration, and a cartesian vector handed to a function expecting the
 * spherical packing compiled and ran and produced coordinates.  Angles and
 * lengths now travel in distinct types, so that mistake is a compile error.
 *
 * \section sphere_frame_array The inner coordinate set
 *
 * `inner_coord_set()` and `inner_coord_set_new()` return an orthonormal frame
 * as nine doubles, three basis vectors laid out end to end:
 *
 *     crdSet[0..2]  i -- the outward radial direction at the complex's centre
 *     crdSet[3..5]  j -- tangent to the surface, along the direction of travel
 *     crdSet[6..8]  k -- tangent to the surface, i x j, completing the frame
 *
 * `i` is the local normal, so a rotation of the complex on the surface is a
 * rotation about `i`, and `j` and `k` span the local tangent plane in which
 * that rotation happens.
 *
 * \section sphere_params Why frames go by reference and vectors go by value
 *
 * The two are deliberately inconsistent, and the inconsistency is the
 * optimization.  A frame is nine doubles in a `std::array`, which the ABI hands
 * over in memory: by value the caller has to lay out and copy 72 bytes at every
 * call, and a frame is built once per complex per step and then read once per
 * molecule and once per interface, so those copies are all of the traffic.  It
 * is passed by reference.
 *
 * \ref Vec3D is three doubles in a trivially copyable aggregate, so it is a
 * homogeneous floating-point aggregate and travels in three FP registers -- see
 * the note on the type itself.  Passing one by reference is *slower*: it forces
 * the caller to spill a value it was holding in registers to the stack so that
 * there is an address to point at, and makes every read in the callee a load.
 * Measured on this file at -O3, `const Vec3D&` cost four extra instructions in
 * inner_coord_set() and five in translate_on_sphere().  \ref SphericalCoord is
 * the same shape, for the same reason.
 */
#pragma once

#include "attributes.hpp"
#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vec3D.hpp"

#include <array>
#include <cmath>

/*! \struct SphericalCoord
 * \brief A point on (or above) the sphere, in spherical coordinates.
 *
 * Distinct from \ref Vec3D so that the two conventions cannot be confused: see
 * \ref sphere_frames.  The physics convention, matching the one the ISO uses:
 * `theta` is measured down from +z, `phi` round from +x in the xy-plane.
 */
struct SphericalCoord {
    //! \brief Polar angle from +z, radians, in [0, pi].
    //!
    //! find_spherical_coords() reports the south pole as -pi rather than +pi;
    //! the comment on that function says why it is left that way.
    double theta { 0.0 };
    double phi { 0.0 }; //!< azimuthal angle from +x, radians, in [0, 2*pi)
    double r { 0.0 }; //!< distance from the centre of the sphere, nm

    SphericalCoord() = default;

    SphericalCoord(double theta, double phi, double r) noexcept
        : theta(theta)
        , phi(phi)
        , r(r)
    {
    }
};

/*!
 * \brief Cartesian to spherical.
 *
 * `r` is taken from the coordinate, so it is the radius of whatever sphere the
 * point happens to lie on -- see \ref sphere_on_not_in.  Both poles are special
 * cases, tested for exactly; `phi` is undefined there and reported as zero.
 */
NERDSS_NODISCARD SphericalCoord find_spherical_coords(Vec3D coord);

//! \brief Spherical to cartesian.  Inverse of find_spherical_coords().
NERDSS_NODISCARD Vec3D find_cartesian_coords(SphericalCoord coord);

/*!
 * \brief Where interface 1 ends up once the pair has closed along the surface.
 *
 * \param arc1        geodesic distance interface 1 travels, nm
 * \param iface1      the interface being moved, cartesian, on the surface
 * \param iface2      its partner, cartesian, on the surface
 * \param arcTotal    geodesic separation the pair starts from, nm
 * \param bindRadius  the target geodesic separation, nm; see calc_bindRadius2D()
 */
NERDSS_NODISCARD Vec3D find_position_after_association(
    double arc1, Vec3D iface1, Vec3D iface2, double arcTotal, double bindRadius);

//! \brief The frame at `com` before the step.  Layout: \ref sphere_frame_array.
NERDSS_NODISCARD std::array<double, 9> inner_coord_set(Vec3D com, Vec3D comNew);

//! \brief The frame at `comNew` after the step.  Layout: \ref sphere_frame_array.
NERDSS_NODISCARD std::array<double, 9> inner_coord_set_new(Vec3D com, Vec3D comNew);

//! \brief Carries `targ` from the old frame to the new one, both cartesian.
NERDSS_NODISCARD Vec3D translate_on_sphere(Vec3D targ, Vec3D com, Vec3D comNew,
    const std::array<double, 9>& crdSet, const std::array<double, 9>& crdSetNew);

//! \brief Rotates `targ` by `dAngle` about the local outward normal, both cartesian.
NERDSS_NODISCARD Vec3D rotate_on_sphere(
    Vec3D targ, Vec3D com, const std::array<double, 9>& crdSet, double dAngle);

/*!
 * \brief The chord binding radius re-expressed as a geodesic distance.
 *
 * A binding radius is quoted as a straight-line separation, but two molecules
 * confined to the surface can only approach each other along it, so the
 * distance the association actually has to close is the arc subtending that
 * chord: \f$ 2R \arcsin(\sigma / 2R) \f$.
 */
NERDSS_NODISCARD double calc_bindRadius2D(double bindRadius, Vec3D iFace);

/*!
 * \brief Builds the reference membrane protein the complex is re-oriented against.
 *
 * `memProtein` is filled with the complex's outermost lipid -- or, under the
 * implicit-lipid model where the complex has no lipid member, with the implicit
 * lipid one of its interfaces is bound to.  Its coordinates are then rewritten
 * so that its interface points at the centre of the sphere, which is what makes
 * it usable as a normal.  Exits if the complex has no lipid of either kind.
 */
void set_memProtein_sphere(const Complex& reactCom, Molecule& memProtein,
    const std::vector<Molecule>& moleculeList, const Membrane& membraneObject);

//! \brief As set_memProtein_sphere(), but keeps the lipid's own association coordinates.
void find_Lipid_sphere(const Complex& reactCom, Molecule& lipid,
    const std::vector<Molecule>& moleculeList, const Membrane& membraneObject);
