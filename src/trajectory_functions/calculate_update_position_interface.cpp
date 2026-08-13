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
    /*
    // define reference points
        Vec3D COM = targCom.comCoord;
        Vec3D COMnew = targCom.comCoord + targCom.trajTrans;
	     Vec3D REF = COMnew;
	     Vec3D dist = COMnew - COM;
	     double l = dist.length();
	     double R = COM.length();
	     double lnew = l + l*R*R/(R*R - l*l);
	     Vec3D dREFnew = (lnew/l) * dist;
	     Vec3D REFnew = COM + dREFnew;
	     REFnew = (R/REFnew.length()) * REFnew;
	
	     // trans into spherical coords
	     COM = find_spherical_coords(COM);
	     COMnew = find_spherical_coords(COMnew);
	     REF = find_spherical_coords(REF);
	     REFnew = find_spherical_coords(REFnew);
    // get the rotation angle: dangle
    double dangle = targCom.trajRot.x; // we select the rotation angle x as the angle that the complex rotate on the sphere
  
    Vec3D temp = find_spherical_coords( ifacecrds );
    Vec3D targ = translate_on_sphere(temp, COM, REF, COMnew, REFnew);
    Vec3D finalcrds1 = rotate_on_sphere(targ, COMnew, dangle);
    finalcrds = find_cardesian_coords(finalcrds1);
    */
    Vec3D COM = targCom.comCoord;
    Vec3D COMnew = targCom.comCoord + targCom.trajTrans;
    std::array<double, 9> Crdset = inner_coord_set(COM, COMnew);
    std::array<double, 9> Crdsetnew = inner_coord_set_new(COM, COMnew);
    finalcrds = translate_on_sphere(ifacecrds, COM, COMnew, Crdset, Crdsetnew);
    finalcrds = rotate_on_sphere(finalcrds, COMnew, Crdsetnew, dangle);
    return finalcrds;
}