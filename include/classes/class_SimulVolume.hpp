/*! \file class_simulbox.hpp

 * ### Created on 10/19/18 by Matthew Varga
 * ### Purpose Class for the simulation box cells
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

//#include "classes/class_coord.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Parameters.hpp"
#include "classes/class_Membrane.hpp"
#include "split.cpp"

/*! \ingroup SimulClasses
 * \brief Wrapper for arrays representing the simulation box
 *
 * TODO: Does this need to be variable as the number of proteins increases/decreases, due to maxPairs?
 * TODO: Create enum class for geometry
 */
struct SimulVolume {
    struct SubVolume {
        int absIndex{}; //!< absolute index of the SubBox in SimulBox::subBoxList
        int xIndex{}; //!< index of the SubBox in the x dimension
        int yIndex{}; //!< index of the SubBox in the y dimension
        int zIndex{}; //!< index of the SubBox in the z dimension

        std::vector<int> memberMolList; //!< list of Molecule indices in moleculeList which currently reside in the SubBox
        std::vector<int> neighborList; //!< list of SubBox absolute indices which are neighbors of this SubBox.
        /*! \brief Bit t is set while a Molecule of type t is a member.
         *
         * Lets the pairwise search drop a whole neighbouring SubBox when none
         * of the types in it can pair with the molecule being tested.  Derived
         * from memberMolList, so it is not serialized; the MPI ranks push into
         * memberMolList directly and their search does not read this.
         *
         * Molecule types past 63 all set every bit, which only costs the skip,
         * never correctness.
         */
        uint64_t typeMask{ 0 };

        void display();

        /*
        Function serialize serializes the SubVolume into array of bytes.
        */
        void serialize(unsigned char *arrayRank, int &nArrayRank) {
            PUSH(absIndex);
            PUSH(xIndex);
            PUSH(yIndex);
            PUSH(zIndex);
            serialize_primitive_vector<int>(memberMolList, arrayRank, nArrayRank);
            serialize_primitive_vector<int>(neighborList, arrayRank, nArrayRank);
        }
        /*
        Function deserialize deserializes the SubVolume from arrayRank of bytes.
        */
        void deserialize(unsigned char *arrayRank, int &nArrayRank) {
            POP(absIndex);
            POP(xIndex);
            POP(yIndex);
            POP(zIndex);
            deserialize_primitive_vector<int>(memberMolList, arrayRank, nArrayRank);
            deserialize_primitive_vector<int>(neighborList, arrayRank, nArrayRank);
        }
    };

    struct Dimensions {
        int x{ 0 }; //!< number of SubBoxes in the x dimension
        int y{ 0 }; //!< number of SubBoxes in the y dimension
        int z{ 0 }; //!< number of SubBoxes in the z dimension
        int tot{ 0 }; //!< total number of SubBoxes. For cubic, x*y*z = tot.

        /*! \func check_dimensions
         * \brief Checks the SubBoxes to make sure they are not too small
         */
        void check_dimensions(const Parameters& params, const Membrane &membraneObject);

        Dimensions() = default;
        explicit Dimensions(const Parameters& params, const Membrane &membraneObject);

        /*
        Function serialize serializes the SubVolume into arrayRank of bytes.
        */
        void serialize(unsigned char *arrayRank, int &nArrayRank) {
            PUSH(x);
            PUSH(y);
            PUSH(z);
            PUSH(tot);
        }
        /*
        Function deserialize deserializes the SubVolume from arrayRank of bytes.
        */
        void deserialize(unsigned char *arrayRank, int &nArrayRank) {
            POP(x);
            POP(y);
            POP(z);
            POP(tot);
        }
    };

    int maxNeighbors{ 13 }; //!< maximum number of neighbors a SubBox can have. Currently set to cubic
    Dimensions numSubCells{}; //!< number of SubBoxes in each dimension
    Vec3D subCellSize{}; //!< dimensions of each SubBox in nanometers
    std::vector<SubVolume> subCellList; //!< list of all the SubBoxes in the SimulBox. Size == numSubBoxes.tot
    /*! \brief Bit c is set while subCellList[c] holds members.
     *
     * The registry proper.  Setting a bit is idempotent, so add_member() does
     * not have to ask whether the SubBox was already registered, and the
     * registry cannot pick up a duplicate however the member lists are
     * manipulated between steps.
     *
     * INVARIANT: every SubBox with a non-empty memberMolList has its bit set.
     * clear_member_lists() relies on it -- a non-empty SubBox whose bit is
     * clear is never emptied, so its members survive into the next step and
     * are then added a second time by the re-binning pass.  Grow a
     * memberMolList through add_member(), never through subCellList directly.
     *
     * Derived from subCellList, so it is not part of the MPI wire format.  The
     * MPI-only sites in prepare.cpp and deserialize.cpp still push directly.
     * They are safe because no MPI path calls clear_member_lists() -- the
     * MpiContext overload of update_memberMolLists() sweeps all of subCellList
     * and clears the whole registry outright -- and they are left alone
     * because that path is untested here.
     */
    std::vector<uint64_t> occupancyMask {};

    /*! \brief The set bits of occupancyMask, ascending, as SubBox indices.
     *
     * Derived; refresh_occupied_cells() rebuilds it.  The pairwise search
     * walks this instead of all of subCellList, and ascending order makes it
     * exactly the non-empty subsequence of that walk, so the candidate pairs
     * come out in the same order as before.
     */
    std::vector<int> occupiedSubCells {};

    //! Index of the lowest set bit.  Only called on a non-zero word.
    static int lowest_set_bit(uint64_t word) {
#if defined(__GNUC__) || defined(__clang__)
        return __builtin_ctzll(word);
#else
        int bitItr { 0 };
        while (!(word & 1)) { word >>= 1; ++bitItr; }
        return bitItr;
#endif
    }

    /*!
     * \brief Puts one Molecule into a SubBox, registering the SubBox.
     */
    void add_member(int cellIndex, int molIndex, int molTypeIndex) {
        SubVolume& cell = subCellList[cellIndex];
        cell.memberMolList.push_back(molIndex);
        cell.typeMask |= (molTypeIndex >= 0 && molTypeIndex < 64)
            ? (uint64_t(1) << molTypeIndex)
            : ~uint64_t(0);
        occupancyMask[cellIndex >> 6] |= uint64_t(1) << (cellIndex & 63);
    }

    /*!
     * \brief Rebuilds occupiedSubCells from occupancyMask, in ascending order.
     *
     * Costs one pass over occupancyMask -- one word per 64 SubBoxes -- plus one
     * push_back per occupied SubBox.  This replaced sorting a list that
     * add_member() had filled in molecule order: on rev_3D, sorting the roughly
     * 1800 occupied SubBoxes cost more per step than the walk over all 27 000
     * that skipping them was supposed to save, and the case came out 2.3%
     * slower rather than faster.
     */
    void refresh_occupied_cells() {
        occupiedSubCells.clear();
        for (size_t wordItr{ 0 }; wordItr < occupancyMask.size(); ++wordItr) {
            uint64_t bits { occupancyMask[wordItr] };
            while (bits) {
                occupiedSubCells.push_back(
                    int(wordItr * 64) + lowest_set_bit(bits));
                bits &= bits - 1;
            }
        }
    }

    /*!
     * \brief Main function for the creation of the SubBoxes in the SimulBox.
     *
     * \param[in] params Parameters as given by the parameter file
     */
    void create_simulation_volume(const Parameters& params, const Membrane &membraneObject);

    /*!
     * \brief Set up the neighborLists for each SubBox.
     *
     * A SubBox only looks for neighbors forward and up. This prevents double counting in the pairwise interaction
     * search later in the main function
     */
    void create_cell_neighbor_list_cubic();

    /*!
     * \brief Update the lists of Molecule members in each SubVolume.
     * \param[in] params Parameters as provided by user.
     * \param[in] moleculeList List of all Molecules in the system.
     * \param[in] complexList List of all Complexes in the system.
     * \param[in] molTemplateList List of all provided MolTemplates.
     *
     * Replaces get_bin2.cpp. Also checks to make sure they're still in the confines of the SimulVolume.
     * TODO: I think this can be made more efficient -- it restarts the search for member molecules every time a
     * Molecule doesn't fit.
     */
    void update_memberMolLists(const Parameters& params, std::vector<Molecule>& moleculeList,
			        std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList, const Membrane &membraneObject, int simItr);

    void update_memberMolLists(const Parameters &params, std::vector<Molecule> &moleculeList,
                    std::vector<Complex> &complexList, std::vector<MolTemplate> &molTemplateList, const Membrane &membraneObject, int simItr,
                    MpiContext &mpiContext);  // For parallel programming

    /*!
     * \brief Empties every non-empty memberMolList, using occupiedSubCells.
     *
     * Equivalent to clearing every SubBox in subCellList, but touches only the
     * SubBoxes that actually hold members.  Leaves occupiedSubCells empty, so
     * the "occupiedSubCells lists exactly the non-empty SubBoxes" invariant
     * still holds afterwards.
     */
    void clear_member_lists();

    void display();

    /*
    Function serialize serializes the Vec3D into arrayRank of bytes.
    */
    void serialize(unsigned char *arrayRank, int &nArrayRank) {
        PUSH(maxNeighbors);
        numSubCells.serialize(arrayRank, nArrayRank);
        subCellSize.serialize(arrayRank, nArrayRank);
        serialize_abstract_vector<SubVolume>(subCellList, arrayRank, nArrayRank);
    }
    /*
    Function deserialize deserializes the Vec3D from arrayRank of bytes.
    */
    void deserialize(unsigned char *arrayRank, int &nArrayRank) {
        POP(maxNeighbors);
        numSubCells.deserialize(arrayRank, nArrayRank);
        subCellSize.deserialize(arrayRank, nArrayRank);
        deserialize_abstract_vector<SubVolume>(subCellList, arrayRank, nArrayRank);
    }
};
