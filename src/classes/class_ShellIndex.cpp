/*! \file class_ShellIndex.cpp
 * \brief Latitude-longitude bands over a spherical membrane.
 *
 * The sizing rests on the haversine identity, which for two points at
 * colatitudes tA, tB separated by central angle g and longitude difference dp
 * reads
 *
 *     sin^2(g/2) = sin^2((tA - tB)/2) + sin(tA) sin(tB) sin^2(dp/2)
 *
 * Both terms are non-negative, so a pair within g of each other satisfies
 *
 *     |tA - tB| <= g            and     sin(|dp|/2) <= sin(g/2) / sqrt(sin tA sin tB)
 *
 * Bands a full g wide in colatitude therefore put such a pair at most one band
 * apart, and cells that subtend at least the second bound put it at most one
 * cell apart within a band.  That is what makes the plus-or-minus-one stencil
 * below complete.
 */

#include "classes/class_ShellIndex.hpp"

#include <algorithm>
#include <iostream>

// C++11 needs the out-of-line definition for a static constexpr member that is
// read at run time.
constexpr double ShellIndex::radiusFloorFraction;

namespace {

//! Smallest sin over a colatitude interval.  sin is concave on [0, pi], so its
//! minimum over any sub-interval sits at one of the two ends.
double min_sin_over(double thetaA, double thetaB)
{
    return std::max(0.0, std::min(std::sin(thetaA), std::sin(thetaB)));
}

/*! \brief Largest longitude difference a pair within gammaCut can have.
 *
 * From the haversine bound above with both colatitudes taken at the least
 * favourable point of the range.  Returns 2*pi when the bound does not bite,
 * which is the case near a pole.
 */
double max_delta_phi(double gammaCut, double sinMin)
{
    if (sinMin <= 0.0)
        return 2.0 * M_PI;
    const double arg { std::sin(0.5 * gammaCut) / sinMin };
    if (arg >= 1.0)
        return 2.0 * M_PI;
    return 2.0 * std::asin(arg);
}

int positive_mod(int value, int modulus)
{
    const int rem { value % modulus };
    return rem < 0 ? rem + modulus : rem;
}

} // namespace

void ShellIndex::build(double sphereR, double rMaxLimit, bool anyExplicitPairReaction)
{
    active = false;
    bandFirstCell.clear();
    bandLonCount.clear();
    neighborList.clear();
    memberMolList.clear();
    typeMask.clear();
    occupancyMask.clear();
    occupiedCells.clear();
    isShellBinned.clear();
    nBands = 0;
    totalCells = 0;

    if (!(sphereR > 0.0) || !(rMaxLimit > 0.0))
        return;

    // Nothing this index could offer can react, so offering it is pure cost.
    if (!anyExplicitPairReaction)
        return;

    radius = sphereR;
    cutoff = rMaxLimit;
    radiusFloor = radiusFloorFraction * sphereR;

    // Largest central angle a reacting pair can span, derived from the bound
    // the Cartesian grid already relies on rather than from an arc estimate.
    // rMaxLimit is a *chord*: it is built so that two molecules close enough to
    // react have |COM1 - COM2| <= rMaxLimit, arms included.  A chord c at
    // radius r subtends 2*asin(c / 2r), which grows as r shrinks, so the
    // largest angle any admitted pair can span is at the radius floor.
    //
    // The earlier form, cutoff / radiusFloor, is the small-angle limit of this
    // and is smaller, since asin(u) >= u.  It made the cells marginally too
    // narrow, and was safe only by the margin between where surface complexes
    // actually sit (0.971 and 0.995 of sphereR) and the floor at 0.9.  This
    // form needs no such margin.  It is also barely wider in practice: on the
    // R=70 sample it moves gammaCut from 0.30869 to 0.30994 rad, 0.4%, which
    // leaves the band and cell counts unchanged.
    const double halfChord{cutoff / (2.0 * radiusFloor)};
    gammaCut = (halfChord >= 1.0) ? M_PI : 2.0 * std::asin(halfChord);

    // One cutoff spans the sphere, so every surface molecule is a candidate
    // partner for every other one and there is nothing for an index to rule
    // out.  The Cartesian grid keeps the whole job.
    if (gammaCut >= M_PI)
        return;

    nBands = std::max(1, int(std::floor(M_PI / gammaCut)));
    bandWidth = M_PI / nBands;

    bandFirstCell.resize(nBands + 1);
    bandLonCount.resize(nBands);
    int cellCount { 0 };
    for (int bandItr { 0 }; bandItr < nBands; ++bandItr) {
        bandFirstCell[bandItr] = cellCount;
        const double sinMin { min_sin_over(bandItr * bandWidth, (bandItr + 1) * bandWidth) };
        const double dPhiMax { max_delta_phi(gammaCut, sinMin) };
        const int lonCount { (dPhiMax >= 2.0 * M_PI)
                ? 1
                : std::max(1, int(std::floor(2.0 * M_PI / dPhiMax))) };
        bandLonCount[bandItr] = lonCount;
        cellCount += lonCount;
    }
    bandFirstCell[nBands] = cellCount;
    totalCells = cellCount;

    neighborList.assign(totalCells, {});
    memberMolList.assign(totalCells, {});
    typeMask.assign(totalCells, 0);
    occupancyMask.assign((totalCells + 63) / 64, 0);

    for (int bandItr { 0 }; bandItr < nBands; ++bandItr) {
        const int lonCount { bandLonCount[bandItr] };
        const double cellPhi { 2.0 * M_PI / lonCount };

        for (int lonItr { 0 }; lonItr < lonCount; ++lonItr) {
            const int cellIndex { bandFirstCell[bandItr] + lonItr };
            std::vector<int>& neighbors = neighborList[cellIndex];

            // Forward along the band.  A ring of two cells would name the same
            // pair from both of them, so only the first one carries it; a ring
            // of one has no partner cell at all.
            if (lonCount >= 3)
                neighbors.push_back(bandFirstCell[bandItr] + positive_mod(lonItr + 1, lonCount));
            else if (lonCount == 2 && lonItr == 0)
                neighbors.push_back(bandFirstCell[bandItr] + 1);

            if (bandItr + 1 >= nBands)
                continue;

            // Forward into the next band.  Only the lower band carries the
            // pair, so each one is offered exactly once.  The longitude bound
            // is taken over both bands together, since the partner may sit
            // anywhere in the union of the two.
            const int nextCount { bandLonCount[bandItr + 1] };
            const double sinMinPair {
                min_sin_over(bandItr * bandWidth, (bandItr + 2) * bandWidth)
            };
            const double dPhiMax { max_delta_phi(gammaCut, sinMinPair) };
            const double loPhi { lonItr * cellPhi - dPhiMax };
            const double hiPhi { (lonItr + 1) * cellPhi + dPhiMax };

            if (hiPhi - loPhi >= 2.0 * M_PI) {
                for (int nextItr { 0 }; nextItr < nextCount; ++nextItr)
                    neighbors.push_back(bandFirstCell[bandItr + 1] + nextItr);
                continue;
            }

            const double nextPhi { 2.0 * M_PI / nextCount };
            const int firstNext { int(std::floor(loPhi / nextPhi)) };
            const int lastNext { int(std::floor(hiPhi / nextPhi)) };
            if (lastNext - firstNext + 1 >= nextCount) {
                for (int nextItr { 0 }; nextItr < nextCount; ++nextItr)
                    neighbors.push_back(bandFirstCell[bandItr + 1] + nextItr);
                continue;
            }
            // The span is shorter than the ring, so the wrapped indices are
            // all distinct and no de-duplication is needed.  Sorted only so
            // the walk over them runs forward through memory.
            const size_t spanStart { neighbors.size() };
            for (int nextItr { firstNext }; nextItr <= lastNext; ++nextItr)
                neighbors.push_back(
                    bandFirstCell[bandItr + 1] + positive_mod(nextItr, nextCount));
            std::sort(neighbors.begin() + spanStart, neighbors.end());
        }
    }

    active = true;
}

int ShellIndex::cell_of(const Vec3D& com) const
{
    const double len { std::sqrt(com.x * com.x + com.y * com.y + com.z * com.z) };
    if (!(len > 0.0))
        return 0;

    const double cosTheta { std::max(-1.0, std::min(1.0, com.z / len)) };
    const double theta { std::acos(cosTheta) };
    int bandItr { int(theta / bandWidth) };
    bandItr = std::max(0, std::min(nBands - 1, bandItr));

    double phi { std::atan2(com.y, com.x) };
    if (phi < 0.0)
        phi += 2.0 * M_PI;
    const int lonCount { bandLonCount[bandItr] };
    int lonItr { int(phi * lonCount / (2.0 * M_PI)) };
    lonItr = std::max(0, std::min(lonCount - 1, lonItr));

    return bandFirstCell[bandItr] + lonItr;
}

void ShellIndex::rebin(const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList)
{
    if (!active)
        return;

    for (int cellIndex : occupiedCells) {
        memberMolList[cellIndex].clear();
        typeMask[cellIndex] = 0;
    }
    std::fill(occupancyMask.begin(), occupancyMask.end(), 0);
    occupiedCells.clear();

    isShellBinned.assign(moleculeList.size(), 0);

    for (size_t molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
        const Molecule& mol = moleculeList[molItr];
        if (mol.isEmpty || mol.isImplicitLipid)
            continue;
        if (mol.myComIndex < 0 || size_t(mol.myComIndex) >= complexList.size())
            continue;
        if (!complexList[mol.myComIndex].OnSurface)
            continue;

        const double len { std::sqrt(mol.comCoord.x * mol.comCoord.x
            + mol.comCoord.y * mol.comCoord.y + mol.comCoord.z * mol.comCoord.z) };
        if (len < radiusFloor)
            continue; // the Cartesian grid keeps this one

        const int cellIndex { cell_of(mol.comCoord) };
        memberMolList[cellIndex].push_back(int(molItr));
        typeMask[cellIndex] |= (mol.molTypeIndex >= 0 && mol.molTypeIndex < 64)
            ? (uint64_t(1) << mol.molTypeIndex)
            : ~uint64_t(0);
        occupancyMask[cellIndex >> 6] |= uint64_t(1) << (cellIndex & 63);
        isShellBinned[molItr] = 1;
    }

    refresh_occupied();
}

void ShellIndex::refresh_occupied()
{
    occupiedCells.clear();
    for (size_t wordItr { 0 }; wordItr < occupancyMask.size(); ++wordItr) {
        uint64_t bits { occupancyMask[wordItr] };
        while (bits) {
#if defined(__GNUC__) || defined(__clang__)
            const int bit { __builtin_ctzll(bits) };
#else
            int bit { 0 };
            for (uint64_t probe { bits }; !(probe & 1); probe >>= 1)
                ++bit;
#endif
            occupiedCells.push_back(int(wordItr * 64) + bit);
            bits &= bits - 1;
        }
    }
}

void ShellIndex::display() const
{
    if (!active) {
        std::cout << "Spherical shell index: inactive.\n";
        return;
    }
    std::cout << "Spherical shell index:\n";
    std::cout << "\tShell radius: " << radius << ", cutoff arc: " << cutoff << '\n';
    std::cout << "\tBands: " << nBands << ", cells: " << totalCells << '\n';
    std::cout << "\tMax central angle for a reacting pair: " << gammaCut << " rad\n";
    std::cout << "\tMolecules below radius " << radiusFloor
              << " stay on the Cartesian grid\n";
}
