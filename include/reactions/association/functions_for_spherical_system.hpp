/*
 * ### Created on 2/05/2020 by Yiben Fu
 * ### Purpose
 * ***
 * all functions that are used in spherical systems
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vec3D.hpp"

#include <array>
#include <cmath>

double radius(Vec3D mol);

Vec3D find_spherical_coords(Vec3D mol); // mol: cardesian coords, output spherical coords

Vec3D find_cardesian_coords(Vec3D mol); // mol: spherical coords, output cardesian coords

Vec3D find_position_after_association(double alpha1, Vec3D Iface1, Vec3D Iface2, double alpha_total, double bindRadius); // on sphere, when association, the new position of Iface1. alpha1 is the geodesic angle that Iface1 moves.

std::array<double, 9> inner_coord_set(Vec3D com, Vec3D comnew);
std::array<double, 9> inner_coord_set_new(Vec3D com, Vec3D comnew);
std::array<double, 3> calculate_inner_coord_coefficients(Vec3D TARG, Vec3D COM, std::array<double, 9> crdset);
Vec3D translate_on_sphere(Vec3D targ, Vec3D COM, Vec3D COMnew, std::array<double, 9> crdset, std::array<double, 9> crdsetnew);
Vec3D rotate_on_sphere(Vec3D Targ, Vec3D COM, std::array<double, 9> crdset, double dangle);

double calc_bindRadius2D(double bindRadius, Vec3D iFace);

void set_memProtein_sphere(Complex reactCom, Molecule& memProtein, std::vector<Molecule> moleculeList, const Membrane membraneObject);
void find_Lipid_sphere(Complex reactCom, Molecule& Lipid, std::vector<Molecule> moleculeList, const Membrane membraneObject);