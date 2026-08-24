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
    /*! \brief Indices of the SubBoxes whose memberMolList is currently non-empty.
     *
     * Only ever a few hundred entries even when subCellList holds thousands of
     * SubBoxes, which lets clear_member_lists() empty the member lists in time
     * proportional to the number of occupied SubBoxes rather than to the total.
     * Derived from subCellList, so it is not part of the MPI wire format.
     *
     * INVARIANT: every SubBox with a non-empty memberMolList appears here
     * exactly once.  clear_member_lists() relies on it -- a non-empty SubBox
     * that is missing from this list is never emptied, so its members survive
     * into the next step and are then added a second time by the re-binning
     * pass.  Grow a memberMolList through add_member(), never through
     * subCellList directly: pushing directly is what broke the invariant
     * before, because molecules created by zeroth-order and unimolecular
     * creation reactions are binned before update_memberMolLists() runs.
     *
     * The MPI-only sites in prepare.cpp and deserialize.cpp still push
     * directly.  They are safe because no MPI path calls clear_member_lists()
     * -- the MpiContext overload of update_memberMolLists() sweeps all of
     * subCellList and empties this list outright -- and they are left alone
     * because that path is untested here.
     */
    std::vector<int> occupiedSubCells {};

    /*!
     * \brief Puts one Molecule into a SubBox, keeping occupiedSubCells right.
     *
     * The registration test is "was this SubBox empty", which both keeps the
     * list complete and keeps it free of duplicates: the second and later
     * members of a SubBox find it already non-empty.
     */
    void add_member(int cellIndex, int molIndex) {
        if (subCellList[cellIndex].memberMolList.empty())
            occupiedSubCells.push_back(cellIndex);
        subCellList[cellIndex].memberMolList.push_back(molIndex);
    }

    /*!
     * \brief Puts occupiedSubCells into ascending SubBox order, without repeats.
     *
     * The pairwise search used to walk all of subCellList, which visits
     * SubBoxes in ascending absIndex order.  Sorted, occupiedSubCells is
     * exactly the non-empty subsequence of that walk, so the candidate pairs
     * come out in the same order as before and the random stream -- and with it
     * the trajectory -- is unchanged.
     *
     * add_member() appends in the order molecules are binned rather than in
     * SubBox order, hence the sort.  The uniquing covers one sequence inside a
     * single step: a SubBox emptied by a dissociation stays on the list, and a
     * creation reaction that then bins a molecule into it finds it empty and
     * registers it a second time.  clear_member_lists() tolerates that repeat;
     * visiting the SubBox twice in the pairwise search would not.
     */
    void sort_occupied_cells() {
        std::sort(occupiedSubCells.begin(), occupiedSubCells.end());
        occupiedSubCells.erase(
            std::unique(occupiedSubCells.begin(), occupiedSubCells.end()),
            occupiedSubCells.end());
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
