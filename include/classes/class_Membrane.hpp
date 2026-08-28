/*! \file class_Membrane.hpp
 *
 */

#pragma once

#include <string>
#include <vector>

#include "classes/mpi_functions.hpp"

/*! \enum BoundaryKeywords
 * \ingroup Parser
 * \brief Boundary parameters read in from the command line
 */
enum class BoundaryKeyword : int {
    implicitLipid = 0, //!< use implicit lipid model
    waterBox = 1, //!< use a rectangular box, specify x y and z lengths
    xBCtype = 2, //!< reflecting or periodic?
    yBCtype = 3, //!< reflecting or periodic?
    zBCtype = 4, //!< reflecting or periodic?
    isSphere = 5, //!< use a sphere boundary, if provide sphereR, this will be set to true.
	sphereR = 6, //!< sphere radius
	hasCompartment = 7, //!< Include a compartment centered at the origin
    compartmentR = 8, //!< Radius of the compartment
    compartmentSiteD = 9, //!< The diffusion constant of the surface binding sites
    compartmentSiteRho = 10, //!< The density of the surface binding sites
};

/*! \brief Which boundary the system is enclosed by.
 *
 * This replaces the `isBox` / `isSphere` boolean pair, which could represent
 * two shapes at once and had no way to say "neither".  `Unspecified` is that
 * fourth state made explicit: it is what a run has before the input names a
 * boundary, and it is what the old pair meant when both were false.
 *
 * This is the geometry alone.  Whether the input supplied a `waterBox` is a
 * separate fact - see `waterBoxGiven` - because a sphere run gets a water box
 * too, fabricated by create_water_box() as a bounding box.  The old `isBox`
 * flag recorded that provenance, not a shape: only the `waterBox` keyword ever
 * set it and nothing ever cleared it.
 *
 * The underlying values are not depended on anywhere; serialize() writes the
 * shape as a byte and write_restart() writes the two flags through the
 * accessors.  They are given explicitly only so the wire byte is readable.
 */
enum class BoundaryShape : int {
    Unspecified = 0, //!< no boundary keyword was read
    Box = 1, //!< rectangular water box; dimensions in `waterBox`
    Sphere = 2, //!< spherical boundary of radius `sphereR`
};

struct Membrane {
  //public:

  struct Droplet {
    double D {0};
    double rho {0};

    Droplet() = default;
    /*
    Function serialize serializes the Droplet
    into array of bytes.
    */
    void serialize(unsigned char *arrayRank, int &nArrayRank) {
      PUSH(D);
      PUSH(rho);
    }
    /*
    Function deserialize deserializes the Droplet
    from array of bytes.
    */
    void deserialize(unsigned char *arrayRank, int &nArrayRank) {
      POP(D);
      POP(rho);
    }
  };

  struct WaterBox {
	/*!
	 * \brief Just a container for the water box dimensions
	 * Only cubic at the moment. Not a Vec3D because then it'll be a circular include (since Vec3D needs Parameters)
	 */

	      double x { 0 };
	      double y { 0 };
        double z { 0 };
        double xLeft{0.0};   // left bound of the rank
        double xRight{0.0};  // right bound of the rank
        double volume { 0 };
        WaterBox() = default;
        explicit WaterBox(std::vector<double> vals)
            : x(vals[0])
            , y(vals[1])
            , z(vals[2])
        {
            volume = x * y * z;
            xLeft = -x / 2.0;
            xRight = x / 2.0;
        }
        /*
        Function serialize serializes the WaterBox
        into array of bytes.
        */
        void serialize(unsigned char *arrayRank, int &nArrayRank) {
          PUSH(x);
          PUSH(y);
          PUSH(z);
          PUSH(xLeft);
          PUSH(xRight);
          PUSH(volume);
        }
        /*
        Function deserialize deserializes the WaterBox
        from array of bytes.
        */
        void deserialize(unsigned char *arrayRank, int &nArrayRank) {
          POP(x);
          POP(y);
          POP(z);
          POP(xLeft);
          POP(xRight);
          POP(volume);
        }
    };
    Droplet droplet;
    WaterBox waterBox; //!< water box x, y, z. used to be xboxl, yboxl, zboxl
    double sphereR = 0; //!< for sphere, value of radius in nm.
    double sphereVol = 0;
    // The initializers below matter beyond hygiene: write_restart() prints
    // nSites, No_free_lipids, No_protein and totalSA unconditionally, so
    // without them a run that never sets these members writes indeterminate
    // values into DATA/restart.dat.  That made restart files differ between any
    // two builds of the same source, which in turn made bitwise output
    // comparison useless for every model that does not use implicit lipids.
    int nSites { 0 };
    int nStates { 0 }; // number of the states of implict lipid
    int No_free_lipids { 0 };
    std::vector<int> numberOfFreeLipidsEachState {}; // record the free lipids of each state for IL, updated each step in main function
    int No_protein { 0 }; // use for implicit-lipid model;
    std::vector<int> numberOfProteinEachState {}; // record the number of proteins that can bound to each state for IL
    int implicitlipidIndex { -1 };
    std::vector<double> RS3Dvect; //this is the look-up table for RS3D, which is the reflecting-surface for 3D-->2D reaction of implicit-lipid case

    //    double RD2D = 0; // block-distance for 2D->2D reaction of implicit-lipid case
    double totalSA { 0 };
    double Dx { 0 };
    double Dy { 0 };
    double Dz { 0 };
    double Drx { 0 };
    double Dry { 0 };
    double Drz { 0 };
    double offset { 0 };
    double lipidLength { 0.0 };
    bool implicitLipid = false;
    bool TwoD = false;
    //! The enclosing boundary.  Read it through isSphere(); that is the only
    //! question the other 65 sites ask of it.
    BoundaryShape shape { BoundaryShape::Unspecified };
    //! True once a `waterBox` keyword has been read.  This is the old `isBox`
    //! flag under a name that says what it meant: it is parse provenance, not
    //! geometry, and a sphere run can carry both it and a Sphere `shape`.
    bool waterBoxGiven { false };
    std::string xBCtype; //allow reflect, or pbc
    std::string yBCtype;
    std::string zBCtype;

    bool hasCompartment = false;
    double compartmentR = 0.0;

  /*set_value_BC is defined in src/parser/parse_input.cpp
      And the map to BoundaryKeyword keywords is also defined in that file.
      BoundaryKeyword keywords are defined above.
      ParameterKeywords are in include/classes/class_Parameters.hpp
     */

    //! \brief True when the input supplied a `waterBox`.  NOT the complement
    //! of isSphere(): a sphere run that also names a waterBox answers true to
    //! both, exactly as the old flag pair did.  Only display() and
    //! write_restart() ask this.
    bool hasWaterBox() const { return waterBoxGiven; }
    //! \brief True when the system is enclosed by a sphere of radius `sphereR`.
    bool isSphere() const { return shape == BoundaryShape::Sphere; }

    //! \brief Abort unless the input named a boundary.  See parse_input.cpp.
    void require_boundary() const;

    void set_value_BC(std::string value, BoundaryKeyword keywords);
    /*In here, we could also store coordinate vector
      for a single representative lipid
    */

    void display(); // display the information for the boundary, define in the src/parse/parse_input.cpp

    void create_water_box(); // create box for sphere boundary, define in the src/parse/parse_input.cpp

    /*
  Function serialize serializes the Molecule
  into array of bytes.
  */
  void serialize(unsigned char *arrayRank, int &nArrayRank) {
    // std::cout << "+Membrane serialization starts here..." << std::endl;
    droplet.serialize(arrayRank, nArrayRank);
    waterBox.serialize(arrayRank, nArrayRank);  // serialize starting
    // from beginning of arrayRank
    // increased by the number of bytes already serialized
    PUSH(sphereR);
    PUSH(sphereVol);
    PUSH(nSites);
    PUSH(nStates);
    PUSH(No_free_lipids);
    serialize_primitive_vector<int>(numberOfFreeLipidsEachState, arrayRank,
                                    nArrayRank);
    PUSH(No_protein);
    serialize_primitive_vector<int>(numberOfProteinEachState, arrayRank,
                                    nArrayRank);
    PUSH(implicitlipidIndex);
    serialize_primitive_vector<double>(RS3Dvect, arrayRank, nArrayRank);
    //    double RD2D = 0; // block-distance for 2D->2D reaction of
    //    implicit-lipid case
    PUSH(totalSA);
    PUSH(Dx);
    PUSH(Dy);
    PUSH(Dz);
    PUSH(Drx);
    PUSH(Dry);
    PUSH(Drz);
    PUSH(offset);
    PUSH(lipidLength);
    PUSH(implicitLipid);
    PUSH(TwoD);
    // Serialized as a byte, not as the 4-byte enum: everything before this
    // point leaves the cursor at an offset of 2 (mod 4), so a 4-byte store
    // through PUSH's cast pointer would be misaligned, which is undefined
    // behavior and a SIGBUS on a strict-alignment target.  A byte also keeps
    // the record the same size as the two bools this replaced.
    unsigned char shapeByte { static_cast<unsigned char>(shape) };
    PUSH(shapeByte);
    PUSH(waterBoxGiven);
    serialize_string(xBCtype, arrayRank, nArrayRank);
    serialize_string(yBCtype, arrayRank, nArrayRank);
    serialize_string(zBCtype, arrayRank, nArrayRank);
    // std::cout << "+Total Membrane size in bytes: " << nArrayRank <<
    // std::endl;
  }
  void deserialize(unsigned char *arrayRank, int &nArrayRank) {
    droplet.deserialize(arrayRank, nArrayRank);
    waterBox.deserialize(arrayRank, nArrayRank);
    POP(sphereR);
    POP(sphereVol);
    POP(nSites);
    POP(nStates);
    POP(No_free_lipids);
    deserialize_primitive_vector<int>(numberOfFreeLipidsEachState, arrayRank,
                                      nArrayRank);
    POP(No_protein);
    deserialize_primitive_vector<int>(numberOfProteinEachState, arrayRank,
                                      nArrayRank);
    POP(implicitlipidIndex);
    deserialize_primitive_vector<double>(RS3Dvect, arrayRank, nArrayRank);
    //    double RD2D = 0; // block-distance for 2D->2D reaction of
    //    implicit-lipid case
    POP(totalSA);
    POP(Dx);
    POP(Dy);
    POP(Dz);
    POP(Drx);
    POP(Dry);
    POP(Drz);
    POP(offset);
    POP(lipidLength);
    POP(implicitLipid);
    POP(TwoD);
    unsigned char shapeByte { 0 };
    POP(shapeByte);
    shape = static_cast<BoundaryShape>(shapeByte);
    POP(waterBoxGiven);
    deserialize_string(xBCtype, arrayRank, nArrayRank);
    deserialize_string(yBCtype, arrayRank, nArrayRank);
    deserialize_string(zBCtype, arrayRank, nArrayRank);
  }
};
