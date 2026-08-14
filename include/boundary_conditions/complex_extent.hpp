/*! \file complex_extent.hpp
 * \brief The scan every reflection routine in boundary_conditions was writing
 *        out by hand, once per Cartesian axis.
 *
 * Every routine in `src/boundary_conditions` answers the same question before
 * it does anything else: *how far does this complex reach, once its sampled
 * translation and rotation are applied?*  Each one answered it with the same
 * doubly nested loop - member molecule COM, then each of that molecule's
 * interfaces - written out three times, once for X, once for Y and once for Z,
 * with `x`/`currx`/`Posdx` search-and-replaced to `y` and `z`.  Six files held
 * eighteen copies of it between them.
 *
 * The loops are collected here.  The arithmetic is character for character what
 * the axis blocks computed, in the same order, so every value rounds and
 * contracts as it did before:
 *
 *   * the position of a point is `base + rot`, where `base` is the axis
 *     component of the complex COM plus its translation and `rot` is
 *     `row0*v.x + row1*v.y + row2*v.z` - two statements, as before, not one
 *     folded expression;
 *   * running extremes are updated with the same `>` / `<` comparisons in the
 *     same order, so ties and signed zeroes land the same way.
 *
 * \section extent_derived Quantities the axis blocks tracked twice
 *
 * Several callers also accumulated, inside the same loop, the largest distance
 * a point stuck out past a wall (`Posdx`, `Negdx`) and the flags saying whether
 * anything stuck out at all.  All four are functions of `posWall` and `negWall`
 * alone:
 *
 *   `outsidePos == posWall > posSide`, and when that holds `posWall` is the
 *   largest coordinate seen (the seed value `negSide` is below `posSide`), so
 *   `posWall - posSide` is the same subtraction on the same operands that the
 *   loop performed for that point.  The negative side likewise.
 *
 * So the callers derive them after the scan instead of carrying them through
 * it.  \ref WallExtent::overshoot and \ref WallExtent::undershoot spell that
 * out.
 */
#pragma once

#include "classes/class_Molecule_Complex.hpp"
#include "math/matrix.hpp"

#include <array>

/*! \ingroup BoundaryConditions
 * \brief Farthest reach of a complex along one axis: `[negWall, posWall]`.
 *
 * Seeded with the wall positions rather than with infinities, because that is
 * what the routines here have always done - a complex entirely inside the box
 * leaves the seed values in place.
 *
 * A plain aggregate with no default member initializers, so that
 * `WallExtent { negSide, posSide }` is brace initialization: this project is
 * built as C++11, where a default member initializer would disqualify it.
 */
struct WallExtent {
    double posWall;
    double negWall;

    //! \brief How far past `posSide` the complex reaches; 0 if it does not reach past it.
    double overshoot(double posSide) const { return posWall > posSide ? posWall - posSide : 0.0; }

    //! \brief How far past `negSide` the complex reaches, as a negative number; 0 if it does not.
    double undershoot(double negSide) const { return negWall < negSide ? negWall - negSide : 0.0; }
};

//! \brief Reach along all three axes, indexed 0/1/2 for x/y/z.
typedef std::array<WallExtent, 3> BoxExtent;

/*! \ingroup BoundaryConditions
 * \brief The point with the largest score seen so far, for the spherical
 *        boundaries, where "farthest out" is one number rather than six.
 *
 * The box routines all score a point the same way - its coordinate along one
 * axis - so \ref scan_axis_extent can do the whole scan.  The spherical ones do
 * not: one wants the largest radius, one the largest amount by which a point
 * pokes out past the membrane, and the compartment ones want the *smallest*
 * radius, each seeded differently and each with its own tie-breaking.  Those
 * differences decide real branches, so rather than hide them behind a flag the
 * callers keep their own one-line score expression and hand the point here.
 *
 * `score` and `point` are seeded by the caller, and `consider` uses a strict
 * `>` so that the first point to reach a given score wins the tie, as every one
 * of these loops has always done.
 */
struct ExtremePoint {
    Vec3D point;
    double score;

    void consider(const Vec3D& candidate, double candidateScore)
    {
        if (candidateScore > score) {
            point = candidate;
            score = candidateScore;
        }
    }
};

//! \brief Reads x, y or z by axis index, so an axis loop can index a Vec3D.
inline double axis_value(const Vec3D& vec, int axis)
{
    return (axis == 0) ? vec.x : ((axis == 1) ? vec.y : vec.z);
}

//! \brief Writable counterpart of \ref axis_value.
inline double& axis_ref(Vec3D& vec, int axis)
{
    return (axis == 0) ? vec.x : ((axis == 1) ? vec.y : vec.z);
}

/*! \ingroup BoundaryConditions
 * \brief Visit the current position of every member COM and every interface.
 *
 * The order - each member molecule's COM, then that molecule's interfaces - is
 * the order the reflection routines have always used, and it decides which
 * point wins a tie for the farthest position.
 */
template <typename F>
inline void for_each_complex_point(
    const Complex& targCom, const std::vector<Molecule>& moleculeList, F fn)
{
    for (auto& memMol : targCom.memberList) {
        fn(moleculeList[memMol].comCoord);
        for (const auto& iface : moleculeList[memMol].interfaceList)
            fn(iface.coord);
    }
}

/*! \ingroup BoundaryConditions
 * \brief `tmpCoords` twin of \ref for_each_complex_point, for the association paths.
 *
 * The interface loop walks `tmpICoords` directly.  `reflect_traj_tmp_crds_box`
 * used to bound it by `interfaceList.size()` while indexing `tmpICoords`, which
 * visits the same sequence: `Molecule::set_tmp_association_coords` and
 * `Molecule::update_association_coords` are the only things that fill
 * `tmpICoords`, and both append exactly one entry per interface.  Walking the
 * vector itself also cannot run off the end if the coordinates were never set.
 */
template <typename F>
inline void for_each_complex_tmp_point(
    const Complex& targCom, const std::vector<Molecule>& moleculeList, F fn)
{
    for (auto& memMol : targCom.memberList) {
        fn(moleculeList[memMol].tmpComCoord);
        for (const auto& iface : moleculeList[memMol].tmpICoords)
            fn(iface);
    }
}

/*! \ingroup BoundaryConditions
 * \brief Reach along one axis after the sampled translation and rotation.
 *
 * \param[in] M Euler matrix of the sampled rotation; row `axis` is the only one read.
 * \param[in] base Axis component of the complex COM plus its translation.
 * \param[in] seed Initial extremes, normally the two walls of that axis.
 */
inline WallExtent scan_axis_extent(const Complex& targCom, const std::vector<Molecule>& moleculeList,
    const std::array<double, 9>& M, int axis, double base, WallExtent seed)
{
    const double row0 { M[3 * axis] };
    const double row1 { M[3 * axis + 1] };
    const double row2 { M[3 * axis + 2] };
    const Vec3D com { targCom.comCoord };

    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D vec { point - com };
        const double rot { row0 * vec.x + row1 * vec.y + row2 * vec.z };
        const double curr { base + rot };

        if (curr > seed.posWall)
            seed.posWall = curr;
        if (curr < seed.negWall)
            seed.negWall = curr;
    });
    return seed;
}

/*! \ingroup BoundaryConditions
 * \brief `tmpCoords` twin of \ref scan_axis_extent.
 *
 * \param[in] base Axis component of the complex tmp COM plus the trial translation.
 */
inline WallExtent scan_axis_extent_tmp(const Complex& targCom, const std::vector<Molecule>& moleculeList,
    const std::array<double, 9>& M, int axis, double base, WallExtent seed)
{
    const double row0 { M[3 * axis] };
    const double row1 { M[3 * axis + 1] };
    const double row2 { M[3 * axis + 2] };
    const Vec3D com { targCom.tmpComCoord };

    for_each_complex_tmp_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D vec { point - com };
        const double rot { row0 * vec.x + row1 * vec.y + row2 * vec.z };
        const double curr { base + rot };

        if (curr > seed.posWall)
            seed.posWall = curr;
        if (curr < seed.negWall)
            seed.negWall = curr;
    });
    return seed;
}

/*! \ingroup BoundaryConditions
 * \brief Reach along all three axes in one pass over the points.
 *
 * `reflect_traj_check_span_box` rotates each point once with the full matrix
 * and tests all three axes, rather than making three passes with one matrix row
 * each; this keeps that.
 */
inline BoxExtent scan_box_extent(const Complex& targCom, const std::vector<Molecule>& moleculeList,
    const std::array<double, 9>& M, const Vec3D& base, BoxExtent seed)
{
    const Vec3D com { targCom.comCoord };

    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        const Vec3D vec { point - com };
        const Vec3D rot { matrix_rotate(vec, M) };

        for (int axis { 0 }; axis < 3; ++axis) {
            const double curr { axis_value(base, axis) + axis_value(rot, axis) };

            if (curr > seed[axis].posWall)
                seed[axis].posWall = curr;
            if (curr < seed[axis].negWall)
                seed[axis].negWall = curr;
        }
    });
    return seed;
}

/*! \ingroup BoundaryConditions
 * \brief Reach along one axis of the complex where it already stands.
 *
 * `reflect_complex_rad_rot_box` runs after the positions have been updated, so
 * it reads coordinates straight off the molecules with no rotation and no
 * translation to add.
 */
inline WallExtent scan_axis_extent_placed(
    const Complex& targCom, const std::vector<Molecule>& moleculeList, int axis, WallExtent seed)
{
    for_each_complex_point(targCom, moleculeList, [&](const Vec3D& point) {
        const double curr { axis_value(point, axis) };

        if (curr > seed.posWall)
            seed.posWall = curr;
        if (curr < seed.negWall)
            seed.negWall = curr;
    });
    return seed;
}

/*! \ingroup BoundaryConditions
 * \brief Reach of the *tmp* coordinates along one axis, where the complex stands.
 *
 * The association-time span check works in tmp coordinates that have already
 * been moved to contact, so like \ref scan_axis_extent_placed there is nothing
 * to rotate or translate.
 */
inline WallExtent scan_axis_extent_tmp_placed(
    const Complex& targCom, const std::vector<Molecule>& moleculeList, int axis, WallExtent seed)
{
    for_each_complex_tmp_point(targCom, moleculeList, [&](const Vec3D& point) {
        const double curr { axis_value(point, axis) };

        if (curr > seed.posWall)
            seed.posWall = curr;
        if (curr < seed.negWall)
            seed.negWall = curr;
    });
    return seed;
}

/*! \ingroup BoundaryConditions
 * \brief Shift a complex's tmp COM and tmp interface coordinates along one axis.
 */
inline void translate_tmp_coords_along_axis(
    Complex& targCom, std::vector<Molecule>& moleculeList, int axis, double shift)
{
    axis_ref(targCom.tmpComCoord, axis) -= shift;
    for (int memMol : targCom.memberList) {
        axis_ref(moleculeList[memMol].tmpComCoord, axis) -= shift;
        for (auto& iface : moleculeList[memMol].tmpICoords)
            axis_ref(iface, axis) -= shift;
    }
}

/*! \ingroup BoundaryConditions
 * \brief Shift a complex's COM and interface coordinates along one axis.
 */
inline void translate_coords_along_axis(
    Complex& targCom, std::vector<Molecule>& moleculeList, int axis, double shift)
{
    axis_ref(targCom.comCoord, axis) -= shift;
    for (auto memMol : targCom.memberList) {
        axis_ref(moleculeList[memMol].comCoord, axis) -= shift;
        for (auto& iface : moleculeList[memMol].interfaceList)
            axis_ref(iface.coord, axis) -= shift;
    }
}
