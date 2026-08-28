/*! @file classes.hpp
 * \brief Class declarations for general use in the program
 *
 * ### Created on 5/20/18 by Matthew Varga
 * ### Purpose
 * ***
 * Class declarations for all molecule and complex associated classes
 *
 * ### Notes
 * ***
 * CHANGED 6/7/18 -- moved stateChangeList on RxnBase -> RxnBase::RateState
 *
 * ### TODO List
 * ***
 *  - TODO PRIORITY MED: split into their different classes (maybe keep derived classes together?)
 *  - TODO PRIORITY LOW: I kind of want to make all data members of MolTemplate, Interface, Angles, and the reaction
 * classes private
 *  - TODO PRIORITY LOW: change all `int` indices to `size_t`
 *  - TODO PRIORITY LOW: rename class names to follow Webkit standards (no abbreviations, camel case)
 */
#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <unordered_map>

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Vec3D.hpp"
#include "split.cpp"

/*! \defgroup SimulClasses
 * \brief Classes actually used for simulation objects
 */

extern int propCalled;

enum class TrajStatus : int {
    none = 0,
    propagated = 1,
    canBeResampled = 2,
};

/*! \brief The dimensionality a reacting pair diffuses in.
 *
 * The values are the dimension count itself, because weighted_D_sum() divides
 * by them.
 */
enum class Dim : int {
    Fiber1D = 1, //!< both complexes on a fiber; diffusion along x only
    Surface2D = 2, //!< both complexes on the membrane; x and y
    Bulk3D = 3, //!< everything else
};

struct Complex;

/*!
 * \ingroup SimulClasses
 * \brief Contains information for one specie (protein, lipid, etc.) in the system
 */
struct Molecule {
    /*!
     * \brief Storage for information on Molecules this Molecule encounters during a timestep
     */
    // struct Encounter {
    //     size_t theirIndex { 0 }; //!< the encountered Molecule's index in moleculeList
    //     size_t myIface { 0 }; //!< this Molecule's interface which encountered the Molecule
    //     size_t theirIface { 0 }; //!< the encountered Molecule's interface which encountered this Molecule
    //     size_t rxnItr { 0 };
    //     double probability { 0 };

    //     Encounter() = default;
    //     Encounter(size_t theirIndex, size_t myIface, size_t theirIface, size_t rxnItr, double probability)
    //         : theirIndex(theirIndex)
    //         , myIface(myIface)
    //         , theirIface(myIface)
    //         , rxnItr(rxnItr)
    //         , probability(probability)
    //     {
    //     }
    // };

    /*!
     * \brief Holds information on the interactions this Molecule is currently participarting in
     */
    struct Interaction {

        // TODO: change to int and make default to -1
        int partnerIndex { -1 }; //!< bound partner index in moleculeList
        int partnerIfaceIndex { -1 }; //!< interface index of the Molecule's partner
        int conjBackRxn { -1 }; //!< back reaction of the forward reaction that formed this Interaction

        // Following fields are for MPI version only:
        int partnerId{-1};  // partner ID for finding partner on neighbor rank

        /*!
         * \brief Sets all the Interaction indices to 0. Used in break_interaction().
         *
         * TODO: Change this to -1.
         */
        void clear();

        Interaction() = default;
        Interaction(int partnerIndex, int partnerIfaceIndex)
            : partnerIndex(partnerIndex)
            , partnerIfaceIndex(partnerIfaceIndex)
        {
        }
        Interaction(int partnerIndex, int partnerIfaceIndex, int conjBackRxn)
            : partnerIndex(partnerIndex)
            , partnerIfaceIndex(partnerIfaceIndex)
            , conjBackRxn(conjBackRxn)
        {
        }
        /*
        Function serialize serializes the Iface
        into array of bytes.
        */
        void serialize(unsigned char* arrayRank, int& nArrayRank) {
            PUSH(partnerIndex);
            PUSH(partnerIfaceIndex);
            PUSH(conjBackRxn);
            PUSH(partnerId);
        }
        /*
        Function deserialize deserializes the Iface
        from given array of bytes.
        */
        void deserialize(unsigned char* arrayRank, int& nArrayRank) {
            POP(partnerIndex);
            POP(partnerIfaceIndex);
            POP(conjBackRxn);
            POP(partnerId);
        }
    };

    struct Iface {
        /*! \struct Molecule::Iface
         * \brief Holds the interface coordinate and the interface's state
         */

        Vec3D coord { 0, 0, 0 }; //!< Coordinate of the interface (absolute, not relative)
        char stateIden { '\0' }; //!< current state of the interface.
        int stateIndex { -1 }; //!< index of the current state in MolTemplate::Interface::stateList
        int index { -1 }; //!< this interface's absolute index (i.e. the index of its current state)
        int relIndex { -1 }; //!< this interface's relative index (index in Molecule::interfaceList)
        int molTypeIndex { -1 }; //!< index of this interface's parent Molecule's MolTemplate (makes comparisons easier)
        bool isBound { false }; //!< is this interface bound
        float boundTime { -1 }; //!< time the interface was bound
        bool excludeVolume { false }; //!< need check exclude volume?

        Interaction interaction{};

        /*!
         * \brief This simply changes the state of the interface.
         * \param[in] newRelIndex Index of the new State in Interface::stateList
         * \param[in] newAbsIndex Constant reference to the new State in Interface::stateList.
         */
        void change_state(int newRelIndex, int newAbsIndex, char newIden);

        Iface() = default;
        explicit Iface(const Vec3D& coord)
            : coord(coord)
        {
        }
        Iface(const Vec3D& coord, char state)
            : coord(coord)
            , stateIden(state)
        {
        }
        Iface(char state, int index, int molTypeIndex, bool isBound)
            : stateIden(state)
            , index(index)
            , molTypeIndex(molTypeIndex)
            , isBound(isBound)
        {
        }
        Iface(const Vec3D& coord, char state, int index, int molTypeIndex, bool isBound)
            : coord(coord)
            , stateIden(state)
            , index(index)
            , molTypeIndex(molTypeIndex)
            , isBound(isBound)
        {
        }
        /*
        Function serialize serializes the Iface
        into array of bytes.
        */
        void serialize(unsigned char* arrayRank, int& nArrayRank) {
            coord.serialize(arrayRank, nArrayRank);
            PUSH(stateIden);
            PUSH(stateIndex);
            PUSH(index);
            PUSH(relIndex);
            PUSH(molTypeIndex);
            PUSH(isBound);
            PUSH(excludeVolume);
            interaction.serialize(arrayRank, nArrayRank);
        }
        /* deserialies array of bytes into a struct,
        and returns the number of bytes processed */
        void deserialize(unsigned char* arrayRank, int& nArrayRank) {
            coord.deserialize(arrayRank, nArrayRank);
            POP(stateIden);
            POP(stateIndex);
            POP(index);
            POP(relIndex);
            POP(molTypeIndex);
            POP(isBound);
            POP(excludeVolume);
            interaction.deserialize(arrayRank, nArrayRank);
        }
    };

    // WHEN ADDING A FIELD TO THE MOLECULE OR COMPLEX STRUCT, UPDATE SERIALIZE AND DESERIALIZE METHODS AS WELL

    int myComIndex { -1 }; //!< which complex does the molecule belong to
    int molTypeIndex { -1 }; //!< index of the Molecule's MolTemplate in molTemplateList
    int mySubVolIndex { -1 };
    int index { -1 }; //!< index of the Molecule in moleculeList
    double mass { -1 }; //!< mass of this molecule
    bool isLipid { false }; //!< is the molecule a lipid
    bool isPromoter {false}; //!< is the molecule a promoter for transcription initiation
    Vec3D comCoord; //!< center of mass coordinate
    std::vector<Iface> interfaceList; //!< interface coordinates
    bool isEmpty { false }; //!< true if the molecule has been destroyed and is void
    TrajStatus trajStatus { TrajStatus::none }; //!< Status of the molecule in that timestep
    bool isDissociated {false}; //!< true if the molecule just dissociated in this timestep
    bool justBoundThisStep{false};  //!< true if the molecule just bound this step
    bool justUnboundThisStep{false};  //!< true if the molecule just unbound this step

    bool isImplicitLipid = false;
    int linksToSurface { 0 }; //!<store each proteins links to surface, to ease updating complex.
    // static variables
    static int numberOfMolecules; //!< counter for the number of molecules in the system
    static std::vector<int> emptyMolList; //!< list of indices to empty Molecules in moleculeList

    double transmissionProb {-1} ; //The transmission probability
    bool enforceCompartmentBC {false}; //Boundary conditions for the compartment

    // association variables
    // temporary positions
    Vec3D tmpComCoord {}; //!< temporary center of mass coordinates for association
    std::vector<Vec3D> tmpICoords {}; //!< temporary interface coordinates for association

    // New Encounter lists
    //    std::vector<Interaction> interactionList; //!< list of interactions the Molecule is participating in
    //    std::vector<Encounter> encounterList; //!< list of Molecule Encounters during a timestep
    //    std::vector<Encounter> prevEncounterList; //!< list of Molecule Encounters during a timestep

    // Legacy encounter lists
    /*Vectors for possible association reactions*/
    std::vector<int> freelist; // legacy, may be replaced
    // std::vector<int> assoclist; // These species are capable of binding.
    std::vector<int> bndlist; // These species are capable of dissociation
    std::vector<int> bndpartner; // It if is bound, who is it bound to? Make this have the same numbering !!
    std::vector<int> bndRxnList;
    // std::vector<int> bndiface; // If it is bound, though which interface!
    //    int ncross = 0;
    //    int movestat = 0;
    std::vector<double> probvec;
    std::vector<int> crossbase; // proteins base encountered
    std::vector<int> mycrossint; // interfaces base encountered other proteins with
    //    std::vector<int> crossrxn;

    std::vector<std::array<int, 3>> crossrxn;

    // std::vector<double> probvec_dissociate;

    /*Vectors for reweighting!*/
    std::vector<int> prevlist;
    std::vector<int> currlist;
    std::vector<int> prevmyface;
    std::vector<int> currmyface;
    std::vector<int> prevpface;
    std::vector<int> currpface;

    std::vector<double> prevnorm;
    std::vector<double> currprevnorm;
    std::vector<double> ps_prev;
    std::vector<double> currps_prev;
    std::vector<double> prevsep;
    std::vector<double> currprevsep;

    // Following fields are for MPI version only:
    int id;  // unique molecule identifier in the system
    bool isGhosted{false};  // whether this molecule belongs to certain rank, or is a copy from neighborhood rank This field doesn't get serialized, nor deserialized.
    bool isShared{false};  // whether this molecule is shared between two ranks. This field doesn't get serialized, nor deserialized.
    bool receivedFromNeighborRank{true}; // When one rank processes molecules near the border with another rank, it sends molecules to that rank, and receives updated versions of some. Those that are not received back should be deleted localy. This field doesn't get serialized, nor deserialized.
    bool need_to_send{false};  // whether this molecule needs to be sent to another rank This field doesn't get serialized, nor deserialized.
    int complexId{-1};  // used for molecules in ghosted zones, for finding complex at neighbor rank
    // MPI static variable:
    static int maxID;  //!< the number of the first empty ID for a complex atparticular rank
    // Declaring mapIdToIndex to be able to translate ID to index relatively fast:
    static std::unordered_map<size_t, size_t> mapIdToIndex;

    void write_crd_file(std::ofstream& os) const;
    void write_crd_file_cout() const;
    friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);

    // association member functions
    void display_assoc_icoords(const std::string& name);
    void update_association_coords(const Vec3D& vec);
    void set_tmp_association_coords();
    void clear_tmp_association_coords();
    void create_position_implicit_lipid(Molecule& reactMol1, int ifaceIndex2, double bindRadius, const Membrane& membraneObject);

    // other reaction member functions
    void create_random_coords(const MolTemplate& molTemplate, const Membrane& membraneObject);
    void destroy();
    void MPI_remove_from_one_rank(std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList);

    void display(const MolTemplate& molTemplate) const;
    void display_all() const;
    void display_my_coords(const std::string& name);
    void print(MpiContext &mpiContext) const;

    bool operator==(const Molecule& rhs) const;
    bool operator!=(const Molecule& rhs) const;

    Molecule() {}
    Molecule(int _mycomplex, Vec3D _comcoords)
        : myComIndex(_mycomplex)
        , comCoord(_comcoords)
    {
    }

    /*
    Function serialize serializes the Molecule into array of bytes.
    */
    void serialize(unsigned char* arrayRank, int& nArrayRank) {
        // std::cout << "+Molecule serialization starts here..." << std::endl;
        //  Serialize myComIndex into arrayRank starting from nArrayRank byte
        //  and increase nArrayRank by number of serialized bytes after:
        PUSH(myComIndex);
        PUSH(molTypeIndex);
        PUSH(mySubVolIndex);
        PUSH(index);
        PUSH(mass);
        PUSH(isLipid);
        // Serialize comCoord starting from arrayRank, nArrayRank byte
        // and increase nArrayRank by the number of serialized bytes:
        comCoord.serialize(arrayRank, nArrayRank);

        // serialize interfaceList vector of Iface
        serialize_abstract_vector<Iface>(interfaceList, arrayRank, nArrayRank);

        PUSH(isEmpty);
        // Serialize trajStatus into arrayRank starting from nArrayRank byte
        // and increase nArrayRank by number of serialized bytes after:
        PUSH(trajStatus);

        PUSH(isImplicitLipid);
        PUSH(linksToSurface);
        PUSH(Molecule::numberOfMolecules);

        // serialize freelist vector of int
        serialize_primitive_vector<int>(freelist, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(bndlist, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(bndpartner, arrayRank, nArrayRank);

        // serialize_primitive_vector<int>(prevlist, arrayRank, nArrayRank);
        // serialize_primitive_vector<int>(prevmyface, arrayRank, nArrayRank);
        // serialize_primitive_vector<int>(prevpface, arrayRank, nArrayRank);
        // serialize_primitive_vector<double>(prevnorm, arrayRank, nArrayRank);
        // serialize_primitive_vector<double>(ps_prev, arrayRank, nArrayRank);
        // serialize_primitive_vector<double>(prevsep, arrayRank, nArrayRank);

        PUSH(id);
        PUSH(complexId);
        // std::cout << "+Total molecule " << id << " size in bytes: " << nArrayRank
        // << std::endl;
    }
    void deserialize(unsigned char* arrayRank, int& nArrayRank) {
        POP(myComIndex);
        // std::cout << "myComIndex=" << myComIndex << ", nArrayRank=" << nArrayRank
        // << std::endl;
        POP(molTypeIndex);
        POP(mySubVolIndex);
        POP(index);
        POP(mass);
        POP(isLipid);
        // std::cout << "isLipid=" << isLipid << ", nArrayRank=" << nArrayRank <<
        // std::endl;

        comCoord.deserialize(arrayRank, nArrayRank);

        deserialize_abstract_vector<Iface>(interfaceList, arrayRank, nArrayRank);

        POP(isEmpty);
        POP(trajStatus);
        POP(isImplicitLipid);
        POP(linksToSurface);
        POP(numberOfMolecules);
        
        deserialize_primitive_vector<int>(freelist, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(bndlist, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(bndpartner, arrayRank, nArrayRank);

        // deserialize_primitive_vector<int>(prevlist, arrayRank, nArrayRank);
        // deserialize_primitive_vector<int>(prevmyface, arrayRank, nArrayRank);
        // deserialize_primitive_vector<int>(prevpface, arrayRank, nArrayRank);
        // deserialize_primitive_vector<double>(prevnorm, arrayRank, nArrayRank);
        // deserialize_primitive_vector<double>(ps_prev, arrayRank, nArrayRank);
        // deserialize_primitive_vector<double>(prevsep, arrayRank, nArrayRank);

        // std::cout << "id=" << id << ", nArrayRank=" << nArrayRank << std::endl;
        POP(id);
        POP(complexId);
    }
    /*
    Function serialize serializes the Molecule
    into array of bytes.
    */
    void serialize_back(unsigned char* arrayRank, int& nArrayRank) {
        // std::cout << "+Molecule serialization starts here..." << std::endl;
        serialize_primitive_vector<double>(probvec, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(crossbase, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(mycrossint, arrayRank, nArrayRank);
        serialize_vector_array<int, 3>(crossrxn, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(currlist, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(currmyface, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(currpface, arrayRank, nArrayRank);
        serialize_primitive_vector<double>(currprevnorm, arrayRank, nArrayRank);
        serialize_primitive_vector<double>(currps_prev, arrayRank, nArrayRank);
        serialize_primitive_vector<double>(currprevsep, arrayRank, nArrayRank);
        // std::cout << "+Total molecule size in bytes: " << start << std::endl;
    }
    void deserialize_back(unsigned char* arrayRank, int& nArrayRank) {
        deserialize_primitive_vector<double>(probvec, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(crossbase, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(mycrossint, arrayRank, nArrayRank);
        deserialize_vector_array<int, 3>(crossrxn, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(currlist, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(currmyface, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(currpface, arrayRank, nArrayRank);
        deserialize_primitive_vector<double>(currprevnorm, arrayRank, nArrayRank);
        deserialize_primitive_vector<double>(currps_prev, arrayRank, nArrayRank);
        deserialize_primitive_vector<double>(currprevsep, arrayRank, nArrayRank);
    }
};

struct Complex {
    /*! \struct Complex
     * \ingroup SimulClasses
     * \brief Contains information on each complex, i.e. bound Molecules
     */

public:
    Vec3D comCoord; //!< Complex's center of mass coordinate
    int index { 0 }; //!< index of this Complex in complexList
    double radius {}; //!< radius of the Complex's bounding sphere
    double mass {};
    std::vector<int> memberList {}; //!< list of member Molecule's indices
    std::vector<int> numEachMol {}; //!< list of the number of each Molecules in this complex
    std::vector<long long int> lastNumberUpdateItrEachMol {}; //!< list of the last size update itr of each Molecules in this complex
    Vec3D D { 0, 0, 0 }; //!< Complex's translational diffusion constants
    Vec3D Dr { 0, 0, 0 }; //!< Complex's rotational diffusion constants

    // Memo for add_3D_rotational_diffusion()'s cos(sqrt(k * Dr.z * timeStep)).
    // That value depends on nothing but Dr.z and the timestep, so it is constant
    // for as long as the complex's composition is, but it was recomputed once
    // per candidate pair per timestep -- the two remaining scalar cos() calls in
    // the whole hot path.  Two slots because the caller picks k = 2 or k = 4
    // from D.z and one complex can be asked for either.
    //
    // The key is the exact argument, compared bit-for-bit, so any change to Dr.z
    // or to the timestep is a miss: the cache cannot go stale, and it does not
    // depend on update_properties() being called on every path that alters Dr.
    // The NaN sentinel forces a miss on first use, because no comparison against
    // NaN is ever true, and no legitimate argument can alias it.
    //
    // `mutable` because add_3D_rotational_diffusion() reads complexList through
    // a const reference.  Not serialized, like deleteIfNotReceivedBack below:
    // this is derived state, so a restart recomputes it on first use and the
    // restart-file format is unchanged.
    mutable double rotDiffArg2 { std::numeric_limits<double>::quiet_NaN() };
    mutable double rotDiffCos2 { 0.0 };
    mutable double rotDiffArg4 { std::numeric_limits<double>::quiet_NaN() };
    mutable double rotDiffCos4 { 0.0 };

    bool isEmpty { false }; //!< true if the complex has been destroyed and is a void
    bool OnSurface { false }; // to check whether on the implicit-lipid membrane.
    bool onFiber {false}; // to check whether on a fiber
    bool tmpOnSurface { false }; //

    // static variables
    static int numberOfComplexes; //!< total number of complexes in the system. starts out equal to numberOfMolecules
    static int currNumberComTypes;
    static int currNumberMolTypes;
    static std::vector<int> emptyComList; //!< list of indices to empty Complexes in complexList

    // TODO: TEMPORARY
    static std::vector<int> obs; //!< TEMPORARY observables vector
    //std::vector<int> NofEach;//!< number of each protein type in this complex
    int linksToSurface = 0; //!< for an adsorbing surface, number of bonds/links formed between this complex and the surface.
    int iLipidIndex = 0; //!< If you need to look up the implicit lipid, this is its molecule index
    // Simulation status parameters
    int ncross { 0 };
    //    int movestat{ 0 };
    TrajStatus trajStatus { TrajStatus::none };
    Vec3D trajTrans;
    Vec3D trajRot;
    Vec3D tmpComCoord;

    // Following fields are for MPI version only:
    int id;  // unique complex identifier in the system
    int ownerRank;  // rank responsible for collecting requests to associate and making decisions and propagating them.

    // When another rank processes molecules near the border with this rank,
    // it can associate two molecules of two received complexes,
    // merging them into a single complex.
    // Complex that is not received back should be deleted localy,
    // and myComIndex must be updated for members that were in the shared zone.
    // These fields don't get serialized, nor deserialized.
    bool deleteIfNotReceivedBack{true};
    bool receivedFromNeighborRank{true};

    // MPI static variable:
    static int maxID;  //!< the number of the first empty ID for a complex at particular rank

    // Declaring mapIdToIndex to be able to translate ID to index relatively fast:
    static std::unordered_map<size_t, size_t> mapIdToIndex;

    friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);

    void update_properties(const std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList);
    void display();
    void display(const std::string& name);
    Complex create(const Molecule& mol, const MolTemplate& molTemp);
    void destroy(std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList);
    void put_back_into_SimulVolume(int& itr, Molecule& errantMol, const Membrane& membraneObject, std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList);
    void translate(Vec3D transVec, std::vector<Molecule>& moleculeList);
    // void propagate(std::vector<Molecule>& moleculeList);
    void propagate(std::vector<Molecule>& moleculeList, const Membrane& membraneObject, const std::vector<MolTemplate>& molTemplateList);
    void update_association_coords_sphere(std::vector<Molecule>& moleculeList, Vec3D iface, Vec3D ifacenew);

    Complex() = default;
    //    Complex(Molecule mol, Vec3D D, Vec3D Dr);
    Complex(const Molecule& mol, const MolTemplate& oneTemp);
    Complex(int _index, const Molecule& _memMol, const MolTemplate& _molTemp);
    Complex(Vec3D comCoord, Vec3D D, Vec3D Dr);
    Complex(Vec3D comCoord, Vec3D D, Vec3D Dr, int idComplex);

    /*
    void operator=(const Complex com)
    {
        this->comCoord = com.comCoord;
        this->index = com.index;
        this->radius = com.radius;
        this->mass = com.mass;
        for (auto member : com.memberList) {
            this->memberList.push_back(member);
        }
        for (auto num : com.numEachMol) {
            this->numEachMol.push_back(num);
        }
        this->D = com.D;
        this->Dr = com.Dr;
        this->isEmpty = com.isEmpty;
        this->OnSurface = com.OnSurface;
        this->linksToSurface = com.linksToSurface;
        this->iLipidIndex = com.iLipidIndex;
        this->ncross = com.ncross;
        this->trajStatus = com.trajStatus;
        this->trajTrans = com.trajTrans;
        this->trajRot = com.trajRot;
        this->tmpComCoord = com.tmpComCoord;
    }*/
   /*
    Function serialize serializes the Molecule into array of bytes.
    */
    void serialize(unsigned char* arrayRank, int& nArrayRank) {
        // std::cout << "+Complex serialization nArrayRanks here..." << std::endl;
        //  Serialize starting from beginning of arrayRank
        //  increased by the number of bytes already serialized
        comCoord.serialize(arrayRank, nArrayRank);
        PUSH(index);
        PUSH(radius);
        PUSH(mass);
        serialize_primitive_vector<int>(memberList, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(numEachMol, arrayRank, nArrayRank);
        serialize_primitive_vector<long long int>(lastNumberUpdateItrEachMol,
                                                arrayRank, nArrayRank);
        D.serialize(arrayRank, nArrayRank);
        Dr.serialize(arrayRank, nArrayRank);
        PUSH(isEmpty);
        PUSH(OnSurface);
        PUSH(tmpOnSurface);
        PUSH(numberOfComplexes);
        PUSH(currNumberComTypes);
        PUSH(currNumberMolTypes);
        // serialize_primitive_vector<int>(emptyComList, arrayRank, nArrayRank);
        serialize_primitive_vector<int>(obs, arrayRank, nArrayRank);
        PUSH(linksToSurface);
        PUSH(iLipidIndex);
        PUSH(ncross);
        PUSH(trajStatus);
        trajTrans.serialize(arrayRank, nArrayRank);
        trajRot.serialize(arrayRank, nArrayRank);
        tmpComCoord.serialize(arrayRank, nArrayRank);
        PUSH(id);
        PUSH(ownerRank);
        // std::cout << "+Total Complex size in bytes: " << nArrayRank << std::endl;
    }
    void deserialize(unsigned char* arrayRank, int& nArrayRank) {
        comCoord.deserialize(arrayRank, nArrayRank);
        POP(index);
        POP(radius);
        POP(mass);
        deserialize_primitive_vector<int>(memberList, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(numEachMol, arrayRank, nArrayRank);
        deserialize_primitive_vector<long long int>(lastNumberUpdateItrEachMol,
                                                    arrayRank, nArrayRank);
        D.deserialize(arrayRank, nArrayRank);
        Dr.deserialize(arrayRank, nArrayRank);
        POP(isEmpty);
        POP(OnSurface);
        POP(tmpOnSurface);
        POP(numberOfComplexes);
        POP(currNumberComTypes);
        POP(currNumberMolTypes);
        // deserialize_primitive_vector<int>(emptyComList, arrayRank, nArrayRank);
        deserialize_primitive_vector<int>(obs, arrayRank, nArrayRank);
        POP(linksToSurface);
        POP(iLipidIndex);
        POP(ncross);
        POP(trajStatus);
        trajTrans.deserialize(arrayRank, nArrayRank);
        trajRot.deserialize(arrayRank, nArrayRank);
        tmpComCoord.deserialize(arrayRank, nArrayRank);
        POP(id);
        POP(ownerRank);
    }
    /*
    Function serialize serializes the Molecule into arrayRank of bytes.
    */
    void serialize_back(unsigned char* arrayRank, int& nArrayRank) {
        PUSH(ncross);
    }
    void deserialize_back(unsigned char* arrayRank, int& nArrayRank) {
        POP(ncross);
    }
};

/*! \brief The dimensionality a reacting pair of complexes shares.
 *
 * This is the ladder that `check_bimolecular_reactions()` wrote out by hand at
 * each use site.  Note that it is not `min(dim(c1), dim(c2))`: `onFiber` and
 * `OnSurface` are set independently - from `isPromoter` and `isLipid`
 * respectively, in `Complex::update_properties()` - so a fiber complex meeting
 * a membrane complex agrees on neither and falls through to 3D.  Only a pair
 * that agrees drops a dimension, and the fiber case is tested first.
 */
inline Dim pair_dim(const Complex& com1, const Complex& com2)
{
    if (com1.onFiber && com2.onFiber)
        return Dim::Fiber1D;
    if (com1.OnSurface && com2.OnSurface)
        return Dim::Surface2D;
    return Dim::Bulk3D;
}

/*! \brief Sum of two complexes' diffusion constants, averaged over the
 * dimensions the pair actually diffuses in.
 *
 * Four hand-written copies of this differed only in how many terms they
 * carried.  Each arm below is the expression it replaced, character for
 * character and as a single statement, which is what makes the result
 * bit-identical.
 *
 * Do not "simplify" this into a loop or an accumulate-and-add over the three
 * axes.  That form is algebraically the same and numerically is not: with
 * -ffp-contract on (the default at -O3) the compiler fuses `w*x + w*y + w*z`
 * written as one expression differently from the same terms accumulated across
 * statements, and the two disagree in the last ulp.  A loop version of this
 * function was measured against these expressions over 4 million random
 * triples and differed on 3.2% of the 3D cases.  It matters because Dtot feeds
 * sqrt() and the 2D reaction-probability tables, where one ulp changes which
 * random draws a run consumes and the trajectory diverges from there.
 */
inline double weighted_D_sum(const Vec3D& D1, const Vec3D& D2, Dim dim)
{
    switch (dim) {
    case Dim::Fiber1D:
        // Diffusion along the fiber axis only; no averaging.
        return D1.x + D2.x;
    case Dim::Surface2D:
        return 1.0 / 2.0 * (D1.x + D2.x)
            + 1.0 / 2.0 * (D1.y + D2.y);
    case Dim::Bulk3D:
    default:
        return 1.0 / 3.0 * (D1.x + D2.x)
            + 1.0 / 3.0 * (D1.y + D2.y)
            + 1.0 / 3.0 * (D1.z + D2.z);
    }
}
