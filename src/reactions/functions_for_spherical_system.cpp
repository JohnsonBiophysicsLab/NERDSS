/*! \file functions_for_spherical_system.cpp
 *
 * ### Created on 2/05/2020 by Yiben Fu
 *
 * \brief Geometry for a system whose membrane *is* a sphere.
 *
 * The surface-not-volume reading of "spherical system", the cartesian/spherical
 * distinction these functions have to keep straight, and the layout of the
 * nine-double inner coordinate set are all documented once, on the header.
 * Start there.
 *
 * \section sphere_tolerances The degeneracy guards
 *
 * Almost every function here opens with a test for a degenerate input: a step
 * that did not move, a target that sits on top of the reference point, a
 * direction that lies along the axis it would have to be rotated about.  Each
 * one exists because the formula that follows divides by the quantity being
 * tested.  They are structural guards on the algebra, not physical cutoffs, and
 * they are not user-configurable for that reason -- unlike the tolerances in
 * numerics/numerical_settings.hpp, which decide when two *rates* or two
 * *angles* count as the same and are properties of the model rather than of the
 * arithmetic.  All three happen to be 1e-8; they are still written separately,
 * because they answer different questions and there is no reason a future
 * change to one should drag the others along.
 *
 * \section sphere_squared Comparing squared lengths, and where that is wrong
 *
 * Those guards used to take a square root and compare the result.  A square root
 * is monotonic, so `|v| < t` is `|v|^2 < t^2` and the root was never needed.
 * Each guard now compares against a `...Squared` constant formed from its
 * tolerance by multiplication, so the two cannot drift apart.
 *
 * That is only safe because the right-hand side is a *fixed tolerance*.  `sqrt`
 * is correctly rounded, so it can map two distinct squared values onto one
 * double; the squared comparison then separates a pair the rooted one called
 * equal.  Against 1e-8 that needs a displacement within about one part in 10^16
 * of the tolerance exactly, which no step in this simulation lands on.
 *
 * It is *not* safe for a comparison between two measured lengths, and this file
 * has three of those -- the nearer of two candidate positions in
 * find_position_after_association(), and the furthest-from-centre searches in
 * set_memProtein_sphere() and find_lipid_sphere().  Those keep their square
 * roots.  The searches are the reason: they run over molecules that all sit on
 * the same sphere, so the lengths being ranked are equal up to rounding by
 * construction, and the collapse is routine rather than rare.  Measured over
 * 2,000,000 pairs of points placed on one sphere of radius 100 by the same
 * construction the callers use:
 *
 *     |a| == |b| after sqrt                  34.42%
 *     (|a| > |b|) disagrees with (|a|^2 > |b|^2)   1.98%
 *
 * A disagreement there does not shift a coordinate by an ulp, it selects a
 * different lipid as the membrane reference, so it is a change in behaviour and
 * not in precision.  A differential test against the previous implementation
 * over 80,000 fixtures reported 2,628 such flips in each of the two searches;
 * none of the 18 cases in the validation and coverage suites shows them, because
 * each of their calls has only one candidate to rank.
 *
 * The saving given up is two square roots per association event and one per
 * lipid per event, against roughly 2,000 events in the sphere case -- far below
 * what the same run spends on the copies removed alongside them.
 */
#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vec3D.hpp"
#include "math/sincos.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {

/*!
 * \brief Below this, two points are the same point.  See \ref sphere_tolerances.
 *
 * Guards the frame construction against a step of zero length and the
 * inner-coordinate solve against a target that coincides with the reference
 * point: in both cases the direction the formula needs is undefined.
 */
const double sameCoordTolerance { 1E-8 };
const double sameCoordToleranceSquared { sameCoordTolerance * sameCoordTolerance };

/*!
 * \brief Below this, a coordinate is on the z axis.  See \ref sphere_tolerances.
 *
 * At either pole the fallback frame's seed direction `(0, 0, 1)` is parallel to
 * the radial direction, so their cross product is zero and cannot be normalized
 * into a basis vector.  `(-1, 0, 0)` is used instead.
 */
const double onAxisTolerance { 1E-8 };

/*!
 * \brief Below this, a projection counts as unit length.  See \ref sphere_tolerances.
 *
 * The second half of rotate_on_sphere()'s "nothing to rotate" test.  Compared
 * against 1 rather than 0, which is what the original wrote and is preserved
 * here: it is a length in nm being compared to unity, so it only ever fires for
 * a target one nanometre out along the rotation axis.
 */
const double unitProjectionTolerance { 1E-8 };

/*!
 * \brief Decomposes `target - com` onto the orthonormal frame `crdSet`.
 *
 * Returns the three coefficients (alpha, beta, gamma) with
 * `target - com == alpha*i + beta*j + gamma*k`, or all zeros when the target
 * coincides with `com`.  The branches ahead of the general solve pick off the
 * cases where the target is perpendicular to one or two of the basis vectors,
 * because the general formula divides by quantities that vanish there.
 *
 * File-local: translate_on_sphere() is the only caller and always has been.  It
 * was declared in the header, which meant it could not be inlined into that one
 * caller and had to take its frame by value.
 */
std::array<double, 3> calculate_inner_coord_coefficients(
    Vec3D target, Vec3D com, const std::array<double, 9>& crdSet)
{
    std::array<double, 3> coefficients {};

    const Vec3D offset { target - com };
    if (offset.length_squared() < sameCoordToleranceSquared) { // target is at com
        return coefficients;
    }

    // get the inner_coords_set
    const Vec3D i { crdSet[0], crdSet[1], crdSet[2] };
    const Vec3D j { crdSet[3], crdSet[4], crdSet[5] };
    const Vec3D k { crdSet[6], crdSet[7], crdSet[8] };

    double alpha, beta, gamma;
    // Whether the offset is perpendicular to a basis vector is asked of the
    // *direction*, so that the answer does not depend on how far away the
    // target is; the coefficients themselves are built from the unnormalized
    // offset.  The three projections are taken once here rather than in each
    // branch that reads them.
    const Vec3D offsetDir { offset.normalized() };
    const double dotI { offsetDir.dot(i) };
    const double dotJ { offsetDir.dot(j) };
    const double dotK { offsetDir.dot(k) };

    if (std::abs(dotI) < sameCoordTolerance) { // offset is perpendicular to i
        alpha = 0.0;
        if (std::abs(dotJ) < sameCoordTolerance) { // and to j, so it is along k
            beta = 0.0;
            gamma = offset.length();
            if (dotK < 0.0) {
                gamma = -gamma;
            }
        } else if (std::abs(dotK) < sameCoordTolerance) { // and to k, so it is along j
            gamma = 0.0;
            beta = offset.length();
            if (dotJ < 0.0) {
                beta = -beta;
            }
        } else {
            beta = (offset.x * k.y - offset.y * k.x) / (j.x * k.y - j.y * k.x);
            gamma = (offset.x * j.y - offset.y * j.x) / (k.x * j.y - k.y * j.x);
        }
    } else if (std::abs(dotJ) < sameCoordTolerance) { // offset is perpendicular to j
        beta = 0.0;
        if (std::abs(dotI) < sameCoordTolerance) { // and to i, so it is along k
            alpha = 0.0;
            gamma = offset.length();
            if (dotK < 0.0) {
                gamma = -gamma;
            }
        } else if (std::abs(dotK) < sameCoordTolerance) { // and to k, so it is along i
            gamma = 0.0;
            alpha = offset.length();
            if (dotI < 0.0) {
                alpha = -alpha;
            }
        } else {
            alpha = (offset.x * k.y - offset.y * k.x) / (i.x * k.y - i.y * k.x);
            gamma = (offset.x * i.y - offset.y * i.x) / (k.x * i.y - k.y * i.x);
        }
    } else if (std::abs(dotK) < sameCoordTolerance) { // offset is perpendicular to k
        gamma = 0.0;
        if (std::abs(dotI) < sameCoordTolerance) { // and to i, so it is along j
            alpha = 0.0;
            beta = offset.length();
            if (dotJ < 0.0) {
                beta = -beta;
            }
        } else if (std::abs(dotJ) < sameCoordTolerance) { // and to j, so it is along i
            beta = 0.0;
            alpha = offset.length();
            if (dotI < 0.0) {
                alpha = -alpha;
            }
        } else {
            alpha = (offset.x * j.y - offset.y * j.x) / (i.x * j.y - i.y * j.x);
            beta = (offset.x * i.y - offset.y * i.x) / (j.x * i.y - j.y * i.x);
        }
    } else {
        const double n1 { k.y * offset.x - k.x * offset.y };
        const double n2 { i.x * k.y - i.y * k.x };
        const double n3 { j.x * k.y - j.y * k.x };
        const double n4 { k.z * offset.x - k.x * offset.z };
        const double n5 { i.x * k.z - i.z * k.x };
        const double n6 { j.x * k.z - j.z * k.x };
        alpha = (n1 * n6 - n4 * n3) / (n2 * n6 - n5 * n3);
        beta = (n1 * n6 - n2 * n6 * alpha) / (n3 * n6);
        gamma = (offset.x - alpha * i.x - beta * j.x) / k.x;
    }

    coefficients[0] = alpha;
    coefficients[1] = beta;
    coefficients[2] = gamma;
    return coefficients;
}

} // namespace

SphericalCoord find_spherical_coords(Vec3D coord)
{
    const double r { coord.length() };
    if (coord.z == r) { // north pole: phi is arbitrary, so it is reported as zero
        return SphericalCoord { 0.0, 0.0, r };
    }
    if (coord.z == -r) {
        /* South pole.  theta is reported as -pi, not the +pi that acos() would
         * return and that \ref SphericalCoord documents as the range.  Either
         * maps back to the same cartesian point -- find_cartesian_coords()
         * differs between the two only in the sign of a sin(theta) that is
         * 1.2e-16 -- and this is what the original wrote, so it is left alone
         * rather than quietly moved.  The one caller that reads `theta` back
         * subtracts a small step from it, and does so identically either way.
         */
        return SphericalCoord { -M_PI, 0.0, r };
    }

    const double theta { std::acos(coord.z / r) }; // acos, gives angle range [0,pi]
    double phi { std::acos(coord.x / (r * std::sin(theta))) };
    if (std::isnan(phi)) {
        phi = 0.0;
    }
    if (coord.y < 0)
        phi = 2.0 * M_PI - phi;
    return SphericalCoord { theta, phi, r };
}

Vec3D find_cartesian_coords(SphericalCoord coord)
{
    // Both sines and both cosines come from two argument reductions rather than
    // four, and r*sin(theta) -- the cylindrical radius, shared by x and y -- is
    // formed once instead of twice.
    double sinTheta, cosTheta;
    sin_cos(coord.theta, sinTheta, cosTheta);
    double sinPhi, cosPhi;
    sin_cos(coord.phi, sinPhi, cosPhi);

    const double rSinTheta { coord.r * sinTheta };
    return Vec3D { rSinTheta * cosPhi, rSinTheta * sinPhi, coord.r * cosTheta };
}

Vec3D find_position_after_association(
    double arc1, Vec3D iface1, Vec3D iface2, double arcTotal, double bindRadius)
{
    const double x1 { iface1.x };
    const double y1 { iface1.y };
    const double z1 { iface1.z };
    const double x2 { iface2.x };
    const double y2 { iface2.y };
    const double z2 { iface2.z };

    const double sphereRadius { iface1.length() };
    // define the unit vector of the plane, Origin-Iface1-Iface2;
    double nx { 1.0 };
    double ny { (x2 * z1 - x1 * z2) / (y1 * z2 - y2 * z1) };
    double nz { (x2 * y1 - x1 * y2) / (z1 * y2 - z2 * y1) };
    const double mag { std::sqrt(nx * nx + ny * ny + nz * nz) };
    nx = nx / mag;
    ny = ny / mag;
    nz = nz / mag;
    /////////////////////////////////////////////////////////
    /* calculate the new position of Iface1.*/
    arc1 = std::abs(arc1);
    // The cosine of the subtended angle is the same for a11 and a22, so it is
    // taken once.
    const double cosArc1 { std::cos(arc1 / sphereRadius) };
    const double a1 { nz * x1 - nx * z1 };
    const double a11 { sphereRadius * sphereRadius * nz * cosArc1 };
    const double a2 { ny * x1 - nx * y1 };
    const double a22 { sphereRadius * sphereRadius * ny * cosArc1 };
    const double a3 { ny * z1 - nz * y1 };
    // function: quadA*x^2 + quadB*x + quadC = 0;
    const double quadA { a1 * a1 + a2 * a2 + a3 * a3 };
    const double quadB { -2.0 * (a1 * a11 + a2 * a22) };
    const double quadC { a11 * a11 + a22 * a22 - a3 * a3 * sphereRadius * sphereRadius };
    double delta { quadB * quadB - 4.0 * quadA * quadC };
    if (delta < 0.0)
        delta = 0.0;
    // The two roots differ only in the sign of this square root, and the
    // reciprocal is common to both.
    const double sqrtDelta { std::sqrt(delta) };
    const double invTwoA { 0.5 / quadA };
    // solution1
    const double root1 { invTwoA * (-quadB + sqrtDelta) };
    const Vec3D solution1 { root1, (a1 * root1 - a11) / a3, -(a2 * root1 - a22) / a3 };
    // solution2
    const double root2 { invTwoA * (-quadB - sqrtDelta) };
    const Vec3D solution2 { root2, (a1 * root2 - a11) / a3, -(a2 * root2 - a22) / a3 };
    // Which of the two the interface actually reaches: the nearer one to the
    // partner when the pair is closing, the further one when it is separating.
    // These keep their square roots, deliberately -- see \ref sphere_squared.
    const double distance1 { (solution1 - iface2).length() };
    const double distance2 { (solution2 - iface2).length() };
    Vec3D newPosition;
    if (bindRadius < arcTotal) {
        newPosition = (distance1 < distance2) ? solution1 : solution2;
    } else {
        newPosition = (distance1 > distance2) ? solution1 : solution2;
    }
    if (std::isnan(newPosition.x) || std::isnan(newPosition.y)) {
        std::cout << "WRONG: non position is generated in 'find_position_after_association'...EXIT! " << std::endl;
        exit(1);
    }
    return newPosition;
}

/*dtheta is the polar angle change, not the solid angle
  COM and targ has already been moved along the azimuth by dphi.
*/
std::array<double, 9> inner_coord_set(Vec3D com, Vec3D comNew)
{
    Vec3D i, j, k, v;
    if ((com - comNew).length_squared() < sameCoordToleranceSquared) {
        // The complex does not move, so there is no direction of travel to build
        // the frame's tangent from; seed it with an arbitrary vector instead.
        i = com;
        i.normalize();
        v = Vec3D { 0.0, 0.0, 1.0 };
        if (std::abs(std::abs(com.z) - com.length()) < onAxisTolerance) {
            v = Vec3D { -1.0, 0.0, 0.0 };
        }
        j = v.unit_cross(i);
        k = i.unit_cross(j);
    } else {
        i = com;
        i.normalize();
        v = comNew;
        v.normalize();
        k = i.unit_cross(v);
        j = k.unit_cross(i);
    }
    i.normalize();
    j.normalize();
    k.normalize();

    return std::array<double, 9> { { i.x, i.y, i.z, j.x, j.y, j.z, k.x, k.y, k.z } };
}

std::array<double, 9> inner_coord_set_new(Vec3D com, Vec3D comNew)
{
    Vec3D i, j, k, v;
    if ((com - comNew).length_squared() < sameCoordToleranceSquared) {
        // As in inner_coord_set(): no motion, so no direction of travel.
        i = comNew;
        i.normalize();
        v = Vec3D { 0.0, 0.0, 1.0 };
        if (std::abs(std::abs(comNew.z) - comNew.length()) < onAxisTolerance) {
            v = Vec3D { -1.0, 0.0, 0.0 };
        }
        j = v.unit_cross(i);
        k = i.unit_cross(j);
    } else {
        // The frame at the end of the step is seeded from a reference point one
        // further step along the same great circle, projected back onto the
        // sphere, so that its tangent points the way the complex was already
        // going rather than the way it came from.
        const Vec3D step { comNew - com };
        const double stepLength { step.length() };
        const double comRadius { com.length() };
        const double extendedLength { stepLength
            + stepLength * comRadius * comRadius / (comRadius * comRadius - stepLength * stepLength) };
        const Vec3D dRefNew { (extendedLength / stepLength) * step };
        Vec3D refNew { com + dRefNew };
        refNew = (comRadius / refNew.length()) * refNew;

        // Deliberately not normalized before the cross products, unlike
        // inner_coord_set(): unit_cross() normalizes its result, so i is the only
        // one that needs the pass below.
        i = comNew;
        v = refNew;
        k = i.unit_cross(v);
        j = k.unit_cross(i);
    }
    i.normalize();
    j.normalize();
    k.normalize();

    return std::array<double, 9> { { i.x, i.y, i.z, j.x, j.y, j.z, k.x, k.y, k.z } };
}

Vec3D translate_on_sphere(Vec3D targ, Vec3D com, Vec3D comNew,
    const std::array<double, 9>& crdSet, const std::array<double, 9>& crdSetNew)
{
    const Vec3D dCom { comNew - com };
    if (dCom.length_squared() < sameCoordToleranceSquared) { // no translation on sphere
        return targ;
    }
    // The point keeps its coefficients on the frame; it is the frame that moves.
    const std::array<double, 3> coefficients { calculate_inner_coord_coefficients(targ, com, crdSet) };
    const double alpha { coefficients[0] };
    const double beta { coefficients[1] };
    const double gamma { coefficients[2] };
    const Vec3D i { crdSetNew[0], crdSetNew[1], crdSetNew[2] };
    const Vec3D j { crdSetNew[3], crdSetNew[4], crdSetNew[5] };
    const Vec3D k { crdSetNew[6], crdSetNew[7], crdSetNew[8] };
    const Vec3D targNew { alpha * i + beta * j + gamma * k };
    return targNew + comNew;
}

Vec3D rotate_on_sphere(
    Vec3D targ, Vec3D com, const std::array<double, 9>& crdSet, double dAngle)
{
    Vec3D targNew;

    // targ will rotate along O-COM line, i.e. i axis. so its projection along i is not change
    const Vec3D i { crdSet[0], crdSet[1], crdSet[2] };
    const Vec3D j { crdSet[3], crdSet[4], crdSet[5] };
    const Vec3D k { crdSet[6], crdSet[7], crdSet[8] };
    const Vec3D offset { targ - com };

    const Vec3D offsetI { i * offset.dot(i) };
    const Vec3D offsetJK { offset - offsetI };
    // Squared, and first in the || so that it short-circuits: a target lying
    // straight out along the radius takes no square root at all on the way to
    // returning unchanged.  A lipid's interface is placed radially from the
    // centre of the sphere and a lone lipid's molecule COM is its complex's COM,
    // so that is the path most membrane molecules take.
    const double offsetJKLengthSquared { offsetJK.length_squared() };
    // if targ is on the line of O-COM, then no need to rotate
    if (offsetJKLengthSquared < sameCoordToleranceSquared
        || std::abs(offsetI.length() - 1.0) < unitProjectionTolerance) {
        targNew = targ;
    } else {
        const double offsetJKLength { std::sqrt(offsetJKLengthSquared) };
        double phi { std::acos(offsetJK.dot(j) / offsetJKLength) };
        if (offsetJK.dot(k) < 0.0) {
            phi = 2.0 * M_PI - phi;
        }
        phi = phi + dAngle;
        double sinPhi, cosPhi;
        sin_cos(phi, sinPhi, cosPhi);
        const Vec3D rotatedJK { j * (offsetJKLength * cosPhi) + k * (offsetJKLength * sinPhi) };
        const Vec3D rotated { offsetI + rotatedJK };
        targNew = Vec3D { rotated.x + com.x, rotated.y + com.y, rotated.z + com.z };
    }
    if (std::isnan(targNew.x)) {
        /* Only the rotating branch above can produce a NaN here: a NaN anywhere
         * in targ, com, crdSet or dAngle makes both tests on the way in false -
         * every comparison against a NaN is - so the non-rotating branch is
         * unreachable with a NaN result.
         *
         * The original re-tested the tangential component's cached magnitude at
         * this point, and the reassignment inside that branch had reset the cache
         * to zero, so this guard always fell through to returning targ unchanged.
         * The `exit(1)` arm it also had was dead code.
         */
        targNew = targ;
    }
    return targNew;
}

double calc_bindRadius2D(double bindRadius, Vec3D iFace)
{
    const double sphereRadius { iFace.length() };
    return sphereRadius * 2 * std::asin((0.5 * bindRadius) / sphereRadius);
}

void set_memProtein_sphere(const Complex& reactCom, Molecule& memProtein,
    const std::vector<Molecule>& moleculeList, const Membrane& membraneObject)
{
    // for explicit lipid model; lipid is a member of reactCom.  The outermost
    // lipid is the one whose tmp COM is furthest from the centre of the sphere.
    // The square roots stay: every candidate is on the same sphere, so this
    // compares near-equal lengths -- see \ref sphere_squared.
    double furthest { 0.0 };
    for (int mol : reactCom.memberList) {
        const Molecule& candidate = moleculeList[mol];
        if (!candidate.isLipid && !candidate.isImplicitLipid)
            continue;
        const double dist { candidate.tmpComCoord.length() };
        if (dist > furthest) {
            memProtein = candidate;
            furthest = dist;
        }
    }
    memProtein.comCoord = memProtein.tmpComCoord;
    // here memProtein is an Lipid, has only one interface
    const double comRadius { memProtein.comCoord.length() }; // same for every interface
    for (size_t i = 0; i < memProtein.interfaceList.size(); i++) {
        const Vec3D ifaceToCom { memProtein.tmpICoords[i] - memProtein.tmpComCoord };
        const double bond { ifaceToCom.length() };
        memProtein.interfaceList[i].coord = (comRadius - bond) / comRadius * memProtein.comCoord;
    }

    // for implicit lipid model, on 2D->2D case, reactCom has no implicitlipid member.
    if (memProtein.isLipid == false && memProtein.isImplicitLipid == false && membraneObject.implicitLipid == true) {
        Vec3D iface;
        Vec3D com;
        furthest = 0.0;
        for (int mol : reactCom.memberList) {
            const Molecule& member = moleculeList[mol];
            for (size_t i = 0; i < member.interfaceList.size(); i++) {
                if (member.interfaceList[i].isBound == true) {
                    const int index = member.interfaceList[i].interaction.partnerIndex;
                    // tmpICoords is indexed only once the partner is known to be
                    // the implicit lipid, as it was: the two vectors are not
                    // guaranteed to be the same length.
                    if (moleculeList[index].isImplicitLipid == true
                        && member.tmpICoords[i].length() > furthest) {
                        memProtein = moleculeList[index];
                        com = member.tmpComCoord;
                        iface = member.tmpICoords[i];
                        furthest = member.tmpICoords[i].length();
                    }
                }
            }
        }
        memProtein.comCoord = iface; // here targ is an ImplicitLipid, has only one interface
        const Vec3D ifaceToCom { iface - com };
        const double bond { ifaceToCom.length() };
        const double implicitLipidRadius { memProtein.comCoord.length() };
        memProtein.interfaceList[0].coord = (implicitLipidRadius - bond) / implicitLipidRadius * memProtein.comCoord;
    }
    if (memProtein.isLipid == false && memProtein.isImplicitLipid == false) {
        std::cout << "WRONG: failed to create memProtein, in the step to adjust complex's orientation on sphere. Exit..." << std::endl;
        exit(1);
    }
}

void find_lipid_sphere(const Complex& reactCom, Molecule& lipid,
    const std::vector<Molecule>& moleculeList, const Membrane& membraneObject)
{
    // for explicit lipid model; lipid is a member of reactCom.  Square roots
    // kept, as in set_memProtein_sphere() -- see \ref sphere_squared.
    double furthest { 0.0 };
    for (int mol : reactCom.memberList) {
        const Molecule& candidate = moleculeList[mol];
        if (!candidate.isLipid && !candidate.isImplicitLipid)
            continue;
        const double dist { candidate.tmpComCoord.length() };
        if (dist > furthest) {
            lipid = candidate;
            furthest = dist;
        }
    }
    // for implicit lipid model, on 2D->2D case, reactCom has no implicitlipid member.
    if (lipid.isLipid == false && lipid.isImplicitLipid == false && membraneObject.implicitLipid == true) {
        furthest = 0.0;
        Vec3D iface;
        Vec3D com;
        for (int mol : reactCom.memberList) {
            const Molecule& member = moleculeList[mol];
            for (size_t i = 0; i < member.interfaceList.size(); i++) {
                if (member.interfaceList[i].isBound == true) {
                    const int index = member.interfaceList[i].interaction.partnerIndex;
                    if (moleculeList[index].isImplicitLipid == true
                        && member.tmpICoords[i].length() > furthest) {
                        lipid = moleculeList[index];
                        com = member.tmpComCoord;
                        iface = member.tmpICoords[i];
                        furthest = member.tmpICoords[i].length();
                    }
                }
            }
        }
        lipid.set_tmp_association_coords();
        lipid.comCoord = iface; // here targ is an ImplicitLipid, has only one interface
        lipid.tmpComCoord = iface;
        lipid.interfaceList[0].coord = com;
        lipid.tmpICoords[0] = com;
    }

    if (lipid.isLipid == false && lipid.isImplicitLipid == false) {
        std::cout << "WRONG: failed to create memProtein, in the step to adjust complex's orientation on sphere. Exit..." << std::endl;
        exit(1);
    }
}
