/*! \file calculate_update_position_interface.cpp
 * ### Created on 02/28/2020 by Yiben Fu
 * ### Purpose: to calculate the new position after the translation and rotation on the sphere surface
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "reactions/association/functions_for_spherical_system.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

Vec3D calculate_update_position_interface(const Complex& targCom, const Vec3D ifacecrds) // iface is cardesian coords
{
    Vec3D finalcrds; // for output
    Vec3D trajTrans = targCom.trajTrans;
    double dangle = targCom.trajRot.x;
    if (trajTrans.length() < 1E-14) {
        finalcrds = ifacecrds;
        return finalcrds;
    }
    Vec3D COM = targCom.comCoord;
    Vec3D COMnew = targCom.comCoord + targCom.trajTrans;
    std::array<double, 9> Crdset = inner_coord_set(COM, COMnew);
    std::array<double, 9> Crdsetnew = inner_coord_set_new(COM, COMnew);
    finalcrds = translate_on_sphere(ifacecrds, COM, COMnew, Crdset, Crdsetnew);
    finalcrds = rotate_on_sphere(finalcrds, COMnew, Crdsetnew, dangle);
    return finalcrds;
}