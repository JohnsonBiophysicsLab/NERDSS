// Proves the closed-form sphere propagation agrees with the angle-based code it
// replaced, to rounding.  The new code makes a different rounding error, so a
// model run diverges from the old one and the suites cannot compare bitwise;
// this compares the two implementations move by move instead, on the same
// random draws.
//
// The old code is copied verbatim below (old_impl), from 205ee203, except that
// the step names old_impl::rotate_on_sphere(), which argument-dependent lookup
// would otherwise find twice:
//   * the step, create_complex_propagation_vectors_on_sphere()
//   * translate_on_sphere() with its Cramer's-rule coefficient solve
//   * rotate_on_sphere() with acos() and the |offsetI| == 1 guard
#include "classes/class_Molecule_Complex.hpp"
#include "math/rand_gsl.hpp"
#include "math/sincos.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "trajectory_functions/trajectory_functions.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_rng.h>
#include <random>

gsl_rng* r; // the generator GaussV() draws from
unsigned long totMatches = 0;

namespace old_impl {

const double sameCoordTolerance { 1E-8 };
const double sameCoordToleranceSquared { sameCoordTolerance * sameCoordTolerance };
const double unitProjectionTolerance { 1E-8 };

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


Vec3D step(const Parameters& params, Complex& targCom)
{
    Vec3D trajTrans;

    const double R_fixed = params.sphereR;
    const Vec3D COM_original = targCom.comCoord;
    const double COM_norm = COM_original.length();

    // project COM onto the surface
    Vec3D COM = Vec3D {
        R_fixed * COM_original.x / COM_norm,
        R_fixed * COM_original.y / COM_norm,
        R_fixed * COM_original.z / COM_norm,
    };

    double dx = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    double dy = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    double dl = sqrt(dx * dx + dy * dy); // propagation length
    double dangle = acos(dx / dl);
    if (dy < 0.0) {
        dangle = 2.0 * M_PI - dangle;
    } // propagation direction
    double rotangle = dangle - M_PI / 2.0;
    // Step along the polar angle by the arc length dl, then rotate that step
    // into the drawn direction below.  The named members are what keeps the
    // spherical triple from being read as a cartesian one.
    const SphericalCoord COMsphere = find_spherical_coords(COM);
    double dtheta = dl / COMsphere.r;
    Vec3D COMnewtmp = find_cartesian_coords(
        SphericalCoord { COMsphere.theta - dtheta, COMsphere.phi, COMsphere.r });
    // define the inner-coords-set
    Vec3D i = Vec3D { COM.x, COM.y, COM.z };
    i.normalize();
    Vec3D temp = Vec3D { 0.0, 0.0, COM.length() };
    if (std::abs(std::abs(COM.z) - COM.length()) < 1E-8) { // inner_coord_set()'s onAxisTolerance
        // COM is on the z axis, so the seed above is parallel to i and the cross
        // product below is the zero vector.  Vec3D::normalize() leaves a zero
        // vector alone rather than making NaNs, so j and k would come back as
        // (0, 0, 0), rotate_on_sphere() would build a zero tangential component,
        // and the complex would lose the whole tangential half of its step:
        // stuck on the axis for good, sinking dl^2/2R per step, with no warning.
        // The same guard, seed and tolerance as inner_coord_set()'s no-motion
        // branch, which this block is otherwise a copy of.  Calling that function
        // instead would not be a copy: it builds k from a j this block has
        // normalized once more, so the two frames differ in the last ulp at 44%
        // of surface positions.
        temp = Vec3D { -1.0, 0.0, 0.0 };
    }
    temp.normalize();
    Vec3D j = temp.unit_cross(i);
    j.normalize();
    Vec3D k = i.unit_cross(j);
    k.normalize();
    //std::vector<Vec3D> crdset;
    std::array<double, 9> crdset;
    crdset[0] = i.x;
    crdset[1] = i.y;
    crdset[2] = i.z;
    crdset[3] = j.x;
    crdset[4] = j.y;
    crdset[5] = j.z;
    crdset[6] = k.x;
    crdset[7] = k.y;
    crdset[8] = k.z;
    Vec3D COMnew = old_impl::rotate_on_sphere(COMnewtmp, COM, crdset, rotangle);

    trajTrans.x = COMnew.x - COM.x;
    trajTrans.y = COMnew.y - COM.y;
    trajTrans.z = COMnew.z - COM.z;
    return trajTrans;
}

} // namespace old_impl

namespace {

struct Worst {
    double value { 0.0 };
    long long trial { -1 };
    void see(double v, long long t)
    {
        if (!(v <= value)) { // NaN counts as worst
            value = v;
            trial = t;
        }
    }
};

Vec3D random_unit(std::mt19937_64& rng)
{
    std::normal_distribution<double> g(0.0, 1.0);
    Vec3D v { g(rng), g(rng), g(rng) };
    return v.normalized();
}

} // namespace

int main()
{
    r = gsl_rng_alloc(gsl_rng_mt19937);
    std::mt19937_64 rng(20261020);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    // ---- the step ----------------------------------------------------------
    // Relative to the step length dl, and to the sphere radius, over steps from
    // dl/R ~ 1e-4 (the sample models) up to dl/R ~ 1 (far past any physical
    // timestep, where the arc formula has to be right and not merely close).
    const long long nStep { 2000000 };
    Worst stepVsR, stepVsDl;
    // Invariants of an exact step, for each implementation: COM + step stays on
    // the sphere through COM, and the chord is 2 r sin(dl / 2r).
    Worst oldOffSphere, newOffSphere, oldChord, newChord;
    long long rngMismatch { 0 };
    for (long long t = 0; t < nStep; ++t) {
        Parameters params;
        params.sphereR = 20.0 + 480.0 * uni(rng);
        params.timeStep = std::pow(10.0, -2.0 + 2.5 * uni(rng));
        Complex com;
        // A complex's COM sits on or a little inside the surface.
        com.comCoord = params.sphereR * (1.0 - 0.05 * uni(rng)) * random_unit(rng);
        const double D { std::pow(10.0, -2.0 + 4.0 * uni(rng)) };
        com.D = Vec3D { D, D, 0.0 };

        const unsigned long seed { static_cast<unsigned long>(rng()) };
        gsl_rng_set(r, seed);
        const Vec3D a { old_impl::step(params, com) };
        const unsigned long nextOld { gsl_rng_get(r) };
        gsl_rng_set(r, seed);
        const Vec3D b { create_complex_propagation_vectors_on_sphere(params, com) };
        const unsigned long nextNew { gsl_rng_get(r) };
        if (nextOld != nextNew)
            ++rngMismatch;

        const double diff { (a - b).length() };
        stepVsR.see(diff / params.sphereR, t);
        stepVsDl.see(diff / std::max(b.length(), 1e-300), t);

        // Redraw the two Gaussians to recover dl.
        gsl_rng_set(r, seed);
        const double dx { std::sqrt(2.0 * params.timeStep * com.D.x) * GaussV() };
        const double dy { std::sqrt(2.0 * params.timeStep * com.D.y) * GaussV() };
        const double dl { std::sqrt(dx * dx + dy * dy) };
        const Vec3D proj { params.sphereR / com.comCoord.length() * com.comCoord };
        const double rr { proj.length() };
        const double chord { 2.0 * rr * std::sin(dl / (2.0 * rr)) };
        oldOffSphere.see(std::abs((proj + a).length() - rr) / rr, t);
        newOffSphere.see(std::abs((proj + b).length() - rr) / rr, t);
        oldChord.see(std::abs(a.length() - chord) / chord, t);
        newChord.see(std::abs(b.length() - chord) / chord, t);
    }
    std::printf("step: %lld trials, max |old-new|/R = %.3g (trial %lld), max |old-new|/|step| = %.3g (trial %lld), "
                "rng streams differing afterwards: %lld\n",
        nStep, stepVsR.value, stepVsR.trial, stepVsDl.value, stepVsDl.trial, rngMismatch);
    std::printf("  off the sphere, relative to r:    old %.3g, new %.3g\n", oldOffSphere.value, newOffSphere.value);
    std::printf("  chord vs 2r sin(dl/2r), relative: old %.3g, new %.3g\n", oldChord.value, newChord.value);

    // ---- translate + rotate of a member point --------------------------------
    // The frames come from inner_coord_set()/inner_coord_set_new(), unchanged.
    // Points are placed up to 20 nm from COM, radially, tangentially or both;
    // one step in ten does not move, which takes translate_on_sphere()'s
    // no-motion branch.
    const long long nPoint { 2000000 };
    Worst pointAbs;
    // A rigid move keeps the distance from the point to COM.
    Worst oldRigid, newRigid;
    long long radialOnly { 0 };
    for (long long t = 0; t < nPoint; ++t) {
        const double R { 20.0 + 480.0 * uni(rng) };
        const Vec3D com { R * (1.0 - 0.05 * uni(rng)) * random_unit(rng) };
        Vec3D comNew { com };
        if (t % 10 != 0) {
            const double stepLen { std::pow(10.0, -3.0 + 3.0 * uni(rng)) };
            const Vec3D n { com.normalized() };
            Vec3D tdir { random_unit(rng) };
            tdir = (tdir - n * tdir.dot(n)).normalized();
            comNew = com + stepLen * tdir - (stepLen * stepLen / (2.0 * R)) * n;
        }
        const std::array<double, 9> f { inner_coord_set(com, comNew) };
        const std::array<double, 9> fNew { inner_coord_set_new(com, comNew) };

        Vec3D targ { com };
        const int kind { static_cast<int>(t % 4) };
        const Vec3D n { com.normalized() };
        if (kind == 0) { // straight out along the radius
            targ = com + (40.0 * uni(rng) - 20.0) * n;
            ++radialOnly;
        } else {
            targ = com + 20.0 * uni(rng) * random_unit(rng);
        }
        const double angle { (2.0 * uni(rng) - 1.0) * M_PI };

        const Vec3D a { old_impl::rotate_on_sphere(
            old_impl::translate_on_sphere(targ, com, comNew, f, fNew), comNew, fNew, angle) };
        double s, c;
        sin_cos(angle, s, c);
        const Vec3D b { rotate_on_sphere(translate_on_sphere(targ, com, comNew, f, fNew), comNew, fNew, s, c) };
        pointAbs.see((a - b).length(), t);
        const Vec3D from { (t % 10 != 0) ? com : comNew };
        const double dist { (targ - from).length() };
        oldRigid.see(std::abs((a - comNew).length() - dist), t);
        newRigid.see(std::abs((b - comNew).length() - dist), t);
    }
    std::printf("translate+rotate: %lld points (%lld radial), max |old-new| = %.3g nm (trial %lld)\n", nPoint,
        radialOnly, pointAbs.value, pointAbs.trial);
    std::printf("  distance to COM not kept, nm:     old %.3g (trial %lld), new %.3g (trial %lld)\n", oldRigid.value,
        oldRigid.trial, newRigid.value, newRigid.trial);

    // ---- the rotation angle actually applied -----------------------------------
    // Measured with atan2, which is well conditioned at every angle, unlike the
    // acos() the old rotation recovered its starting angle with.
    const long long nRot { 2000000 };
    Worst oldAngle, newAngle;
    for (long long t = 0; t < nRot; ++t) {
        const double R { 20.0 + 480.0 * uni(rng) };
        const Vec3D com { R * random_unit(rng) };
        // An orthonormal frame by Gram-Schmidt.  inner_coord_set() builds k from
        // the cross product of two nearly parallel unit vectors, which leaves it
        // up to ~1e-10 off perpendicular for a short step, and that error would
        // show up here as an angle error of either implementation.
        const Vec3D i { com.normalized() };
        Vec3D seed { random_unit(rng) };
        const Vec3D j { (seed - i * seed.dot(i)).normalized() };
        const Vec3D k { i.cross(j) };
        const std::array<double, 9> f { { i.x, i.y, i.z, j.x, j.y, j.z, k.x, k.y, k.z } };
        const double phi0 { (2.0 * uni(rng) - 1.0) * M_PI };
        const double rho { 0.1 + 20.0 * uni(rng) };
        const Vec3D targ { com + (10.0 * uni(rng) - 5.0) * i + rho * std::cos(phi0) * j + rho * std::sin(phi0) * k };
        const double angle { (2.0 * uni(rng) - 1.0) * M_PI };
        double s, c;
        sin_cos(angle, s, c);
        const Vec3D a { old_impl::rotate_on_sphere(targ, com, f, angle) };
        const Vec3D b { rotate_on_sphere(targ, com, f, s, c) };
        auto applied = [&](const Vec3D& p) {
            const Vec3D o { p - com }, o0 { targ - com };
            const double before { std::atan2(o0.dot(k), o0.dot(j)) };
            const double after { std::atan2(o.dot(k), o.dot(j)) };
            return std::remainder(after - before - angle, 2.0 * M_PI);
        };
        oldAngle.see(std::abs(applied(a)), t);
        newAngle.see(std::abs(applied(b)), t);
    }
    std::printf("rotate: %lld points, rotation angle error, rad: old %.3g (trial %lld), new %.3g (trial %lld)\n", nRot,
        oldAngle.value, oldAngle.trial, newAngle.value, newAngle.trial);

    // What the numbers above show, and what is enforced:
    //   * both draw the same two Gaussians, so the random stream is untouched;
    //   * the new step is exact to rounding: on the sphere and with the right
    //     chord to ~1e-15, where the old one's arc was off by up to ~6e-10,
    //     because acos() near +-1 loses half the digits of its argument;
    //   * the new rotation turns by the requested angle to ~1e-13 rad, the old
    //     one by up to ~2e-9 rad off, for the same reason;
    //   * so old and new differ by the old code's error, which depends on the
    //     cases drawn: up to 7.9e-7 of a step and 1.8e-8 nm on a point with
    //     libstdc++'s distributions, 2.4e-7 and 1.9e-9 nm with libc++'s.  The
    //     bounds below leave a factor of ten over that and are still about six
    //     orders of magnitude under what a different move would give: O(1) of a
    //     step, O(0.1) nm on a point.
    //   * the new point transform keeps distances to 1e-9 nm.  What is left
    //     there is inner_coord_set()'s frame for a very short step, which is
    //     unchanged and shared by both.
    const bool ok { rngMismatch == 0 && newOffSphere.value < 1e-14 && newChord.value < 1e-13
        && stepVsDl.value < 1e-5 && pointAbs.value < 1e-7 && newRigid.value < 1e-9 && newAngle.value < 1e-11 };
    gsl_rng_free(r);
    if (!ok) {
        std::printf("sphere_propagation_check: FAILED\n");
        return 1;
    }
    std::printf("sphere_propagation_check: all passed\n");
    return 0;
}
