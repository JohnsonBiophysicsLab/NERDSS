/*! \file check_if_spans_sphere.cpp
 * ### Created on 2020-02-23 by Yiben Fu
 */
#include "boundary_conditions/complex_extent.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "tracing.hpp"

#include <iostream>

namespace {

/*! \brief Visit every tmp coordinate of `targCom` that is not an implicit lipid.
 *
 * The implicit lipid has no position of its own, so the span check has always
 * skipped it; \ref for_each_complex_tmp_point does not, which is why this one
 * lives here rather than in the header.
 */
template <typename F>
void for_each_real_tmp_point(const Complex& targCom, const std::vector<Molecule>& moleculeList, F fn)
{
    for (auto& memMol : targCom.memberList) {
        if (moleculeList[memMol].isImplicitLipid)
            continue;

        fn(moleculeList[memMol].tmpComCoord);
        for (const auto& iface : moleculeList[memMol].tmpICoords)
            fn(iface);
    }
}

//! \brief Shift a complex's tmp COM and tmp interface coordinates.
void translate_tmp_coords(Complex& targCom, std::vector<Molecule>& moleculeList, const Vec3D& trans)
{
    targCom.tmpComCoord += trans;
    for (int memMol : targCom.memberList) {
        moleculeList[memMol].tmpComCoord += trans;
        // update interface coords
        for (auto& iface : moleculeList[memMol].tmpICoords)
            iface += trans;
    }
}

} // namespace

void check_if_spans_sphere(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double radius)
{
    // TRACE();
    // Associating proteins have been moved to contact. Before assigning them to the complexsame complex,
    // test to see if the complex is too big to fit in the box.

    // declare the boundary side of the system;
    double sphereR = radius; // no considering the reflecting-surface, because here we are checking whether to span the box

    // find the new COM and new radius, to check whether the new radius is larger than the sphere radius
    Vec3D newCom;
    double newRadius = 0.0;
    com_of_two_tmp_complexes(reactCom1, reactCom2, newCom, moleculeList);

    auto growRadius = [&](const Vec3D& point) {
        Vec3D disVec { point - newCom };
        if (disVec.length() > newRadius)
            newRadius = disVec.length();
    };
    for_each_real_tmp_point(reactCom1, moleculeList, growRadius);
    for_each_real_tmp_point(reactCom2, moleculeList, growRadius);

    if (newRadius > sphereR) {
        // std::cout << "STICKS OUT THE SPHERE, CANCEL ASSOCIATION " << '\n';
        cancelAssoc = true;
        return;
    }

    // The approximate size of the complex (max size) puts it as outside, now test interface positions.
    // Scored by how far the point pokes out past the membrane, seeded at zero.
    ExtremePoint farthest { Vec3D {}, 0.0 };
    auto measure = [&](const Vec3D& point) { farthest.consider(point, point.length() - sphereR); };
    for_each_real_tmp_point(reactCom1, moleculeList, measure);
    for_each_real_tmp_point(reactCom2, moleculeList, measure);

    // put back in the box. put at edge, rather than bouncing off.
    if (farthest.score > 0.0) {
        double lamda = -farthest.score / farthest.point.length();
        Vec3D trans = lamda * farthest.point;
        translate_tmp_coords(reactCom1, moleculeList, trans);
        translate_tmp_coords(reactCom2, moleculeList, trans);
    }
}
