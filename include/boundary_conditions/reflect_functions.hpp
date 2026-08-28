/*! \file old_funcs.hpp
 * \brief Functions for enforcing reflecting boundary conditions
 * ### Created on 11/1/18 by Matthew Varga
 * ### TODO List
 * ***
 */
#pragma once

#include "classes/class_Membrane.hpp"
#include "classes/class_Molecule_Complex.hpp"

/*! \defgroup BoundaryConditions
 * \brief Functions associated with enforcing various boundary conditions
 */

/*! \ingroup BoundaryConditions
 * \brief Which side of a spherical boundary a complex is required to stay on.
 *
 * A sphere boundary contains: the complex must stay at radius below it.  A
 * compartment excludes: the complex must stay at radius above it.  The two are
 * otherwise the same reflection, which is why the three `*_compartment` files
 * were character-for-character their `*_sphere` twins with a handful of signs
 * flipped.  The values are +1 / -1 so they can multiply a radial distance
 * directly.
 */
enum class RadialSide : int {
    Inside = 1, //!< stay within the radius; a sphere or membrane boundary
    Outside = -1, //!< stay beyond the radius; a compartment
};

//! \brief `side` as the factor +1 / -1 it stands for.
inline double radial_sign(RadialSide side) { return static_cast<double>(static_cast<int>(side)); }

/* DISSOCIATION REFLECTION */

/*! \ingroup BoundaryConditions
 * \brief Enforces reflecting boundary conditions during dissociation reactions.
 */
// void reflect_complex_rad_rot(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput);
void reflect_complex_rad_rot_inx(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList);
void reflect_complex_rad_rot(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput, bool isInsideCompartment);
void reflect_complex_rad_rot_box(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput);
void reflect_complex_rad_rot_sphere(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double radius, double RS3Dinput);
void reflect_complex_compartment(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput);

/* TRAJ PROPAGATION FUNCTIONS */

/*! \ingroup BoundaryConditions
 * \brief Enforces reflecting boundary conditions during complex propagation.
 */
void reflect_traj_complex_rad_rot(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput, bool isInsideCompartment);
void reflect_traj_complex_rad_rot_box(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput);
//! \brief Reflects the sampled translation back across a spherical boundary; see reflect_traj_complex_radial.cpp.
void reflect_traj_complex_radial(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom,
    const Membrane& membraneObject, double radius, RadialSide side, double RS3Dinput, bool skipOnSurface,
    bool recheckSpan);
// void reflect_traj_complex_rad_rot_new(
//     const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 9>& M, const Membrane& membraneObject, double RS3Dinput);

/*! \ingroup BoundaryConditions
 * \brief Checks to make sure the complex doesn't span the box (go out of the box in both sides).
 *
 * Child function of reflect_traj_complex_rad_rot
 */
// void reflect_traj_check_span(double xtot, double ytot, double ztot, const Parameters& params, Complex& targCom,
//     std::vector<Molecule>& moleculeList, std::array<double, 9>& M, const Membrane& membraneObject, double RS3Dinput);
void reflect_traj_check_span(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput);
void reflect_traj_check_span_box(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput);
void reflect_traj_check_span_sphere(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double radius, double RS3Dinput);

/*!
 * \brief Enforces reflecting boundary conditions by placing Complex and component Molecules back into the box without
 * checking bounds after doing so.
 *
 * Child function of reflect_traj_check_span
 */
// void reflect_traj_complex_rad_rot_nocheck(
//     const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, std::array<double, 9>& M, const Membrane& membraneObject, double RS3Dinput);
void reflect_traj_complex_rad_rot_nocheck(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput);
void reflect_traj_complex_rad_rot_nocheck_box(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput);
void reflect_traj_complex_rad_rot_nocheck_sphere(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput);

/*!
 * \brief Checks if the newly associated complex spans the simulation volume.
 *
 * If it does, cancels association. If complex is out of the box in either the positive or negative direction (of any x,
 * y, z dimensions), it puts the complex back into the box at the edge, not bounced off, as in
 * reflect_traj_complex_rad_rot, and its children.
 */

// void check_if_spans_box(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
//     std::vector<Molecule>& moleculeList, const Membrane& membraneObject);
void check_if_spans(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList, const Membrane& membraneObject);
void check_if_spans_box(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList, const Membrane& membraneObject);
void check_if_spans_sphere(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double radius);

/*! \ingroup BoundaryConditions
 * \brief evaluates size of reflection off of walls during association, stores in temporary vector traj.
 *based on tmpCoords, does not update complex.traj vectors..
 */
// void reflect_traj_tmp_crds(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput);
void reflect_traj_tmp_crds(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput, bool isInsideCompartment);
void reflect_traj_tmp_crds_box(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput);
//! \brief Reflects a trial translation back across a spherical boundary; see reflect_traj_tmp_crds_radial.cpp.
void reflect_traj_tmp_crds_radial(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom,
    std::array<double, 3>& traj, double radius, RadialSide side, double RS3Dinput);

// function to calculate the position of one interface after translation and rotation on sphere surface
Vec3D calculate_update_position_interface(const Complex& targCom, const Vec3D ifacecrds); // iface is cardesian coords
