/*! \file class_ShellIndex.hpp
 * \brief A latitude-longitude neighbour index over a spherical membrane.
 *
 * SimulVolume lays one Cartesian grid over the (2R)^3 box that bounds a
 * spherical system, which is the right structure for everything moving through
 * the interior and the wrong one for anything pinned to the shell: only the
 * sub-volumes the shell passes through can ever hold a surface-bound molecule,
 * a fraction pi*h/(2R) of the grid, and the rest of the cell budget is spent on
 * volume that is empty by construction.  It is also the wrong metric.
 * get_distance() measures a surface pair by geodesic arc length while the grid
 * bins on Cartesian coordinates, so the grid over-includes by however much arc
 * exceeds chord.
 *
 * This index covers the surface pairs and nothing else.  Molecules stay in the
 * Cartesian grid as well, so surface-to-interior pairs are unaffected; the
 * pairwise search skips a pair there only when both of its molecules are in
 * this index, and picks it up here instead.
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"

struct ShellIndex {
    //! Cells are only useful while they are at least a cutoff wide, so a
    //! molecule has to be far enough out for its arc length to reach that.
    //! Surface complexes are pinned by the reflectors at sphereR - RS3D, which
    //! measures 0.971 and 0.995 of sphereR on the two spherical samples, so
    //! this admits every one of them with room to spare.
    //!
    //! Nothing about correctness rests on the value.  A molecule below the
    //! floor is simply left to the Cartesian grid, which handles it exactly as
    //! it did before this index existed; the constant trades one search for
    //! another and can never lose a pair.
    static constexpr double radiusFloorFraction { 0.9 };

    bool active { false };
    double radius { 0.0 };      //!< sphereR
    double cutoff { 0.0 };      //!< rMaxLimit, read as an arc length on the shell
    double radiusFloor { 0.0 }; //!< below this a molecule stays on the Cartesian path
    double gammaCut { 0.0 };    //!< the largest central angle a reacting pair can span
    double bandWidth { 0.0 };   //!< colatitude covered by one band
    int nBands { 0 };
    int totalCells { 0 };

    std::vector<int> bandFirstCell {};  //!< size nBands + 1
    std::vector<int> bandLonCount {};   //!< size nBands

    std::vector<std::vector<int>> neighborList {};  //!< forward-only, so each pair is offered once
    std::vector<std::vector<int>> memberMolList {};
    std::vector<uint64_t> typeMask {};      //!< molecule types present, as in SimulVolume
    std::vector<uint64_t> occupancyMask {}; //!< bit c set while cell c holds members
    std::vector<int> occupiedCells {};      //!< set bits of occupancyMask, ascending

    //! 1 while moleculeList[i] is in this index.  The Cartesian pass reads it
    //! to decide which pairs it no longer owns.
    std::vector<char> isShellBinned {};

    /*!
     * \brief Lays out the bands and their neighbour lists for one sphere.
     *
     * Leaves the index inactive, and the search unchanged, for a non-spherical
     * system, one whose cutoff spans the sphere, and one where no two explicit
     * molecules can react with each other at all.  That last case is not
     * hypothetical: the `sphere` sample binds A to an implicit lipid and
     * declares nothing else, so no pair this index could hold is ever worth
     * testing, and binning its 448 surface molecules every step -- an acos and
     * an atan2 apiece -- cost 9% of the run for nothing.
     */
    void build(double sphereR, double rMaxLimit, bool anyExplicitPairReaction);

    //! Cell holding a point in this direction.  Radius is not consulted.
    int cell_of(const Vec3D& com) const;

    //! Empties every cell, then bins every surface molecule that clears the floor.
    void rebin(const std::vector<Molecule>& moleculeList,
        const std::vector<Complex>& complexList);

    void display() const;

private:
    void refresh_occupied();
};
