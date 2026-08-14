/*! \file bimolecular_reactions.hpp
 *
 * \brief
 *
 * ### Created on 2019-02-12 by Matthew Varga
 */
#pragma once

#include "classes/class_Rxns.hpp"
#include <classes/class_copyCounters.hpp>
#include <cmath>
#include <gsl/gsl_matrix.h>

struct BiMolData {
    int pro1Index { 0 };
    int pro2Index { 0 };
    int com1Index { 0 };
    int com2Index { 0 };
    int relIface1 { 0 };
    int relIface2 { 0 };
    int absIface1 { 0 };
    int absIface2 { 0 };
    double Dtot { 0 };
    double magMol1 { 0 };
    double magMol2 { 0 };

    BiMolData() = default;
    BiMolData(int pro1Index, int pro2Index, int com1Index, int com2Index, int relIface1, int relIface2, int absIface1,
        int absIface2, double Dtot, double magMol1, double magMol2)
        : pro1Index(pro1Index)
        , pro2Index(pro2Index)
        , com1Index(com1Index)
        , com2Index(com2Index)
        , relIface1(relIface1)
        , relIface2(relIface2)
        , absIface1(absIface1)
        , absIface2(absIface2)
        , Dtot(Dtot)
        , magMol1(magMol1)
        , magMol2(magMol2)
    {
    }
};

/*! \defgroup RotationalDiffusion Rotational contribution to the pair diffusion constant
 *
 * Before a candidate pair can be tested for reaction, the rotation each partner
 * undergoes during one timestep has to be folded into `BiMolData::Dtot`: an
 * interface swings about its complex COM, and that motion closes distance just
 * as translation does.  The two functions below are that step.
 *
 * Each was written out by hand at four call sites - twice inside
 * `check_bimolecular_reactions()`, which has its own copy of the 2D and the 3D
 * path, and once in each `determine_*_reaction_probability()`.  The arithmetic
 * here is character for character what those sites computed, in the same order,
 * so `Dtot` rounds and contracts exactly as it did.
 * @{
 */

/*! \brief Adds the rotational contribution to `Dtot` for a pair on the membrane.
 *
 * Both partners rotate about the membrane normal only, so both get the same 2D
 * treatment and the two contributions are summed before the single divide.
 */
inline void add_2D_rotational_diffusion(
    BiMolData& biMolData, const std::vector<Complex>& complexList, const Parameters& params)
{
    double Dr1 {};
    {
        double cf { cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep)) };
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
    }

    double Dr2 {};
    {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
    }

    biMolData.Dtot += (Dr1 + Dr2) / (4.0 * params.timeStep); // add in contributions from rotation
}

/*! \brief Snaps `Dtot` onto the coarse grid the 2D lookup tables are keyed on.
 *
 * Only allow 2D diffusion at certain intervals, to avoid generating too many 2D
 * tables.  Keep only one significant figure for <0.1, 2 for 0.1<d<10, 3 for
 * 10<d<100, and so on.  Always paired with \ref add_2D_rotational_diffusion:
 * the value being rounded is the one that function just finished.
 */
inline void discretize_2D_Dtot(BiMolData& biMolData)
{
    double dtmp;
    if (biMolData.Dtot < 0.0001)
        dtmp = biMolData.Dtot * 100000;
    else if (biMolData.Dtot < 0.001)
        dtmp = biMolData.Dtot * 10000;
    else if (biMolData.Dtot < 0.01)
        dtmp = biMolData.Dtot * 1000;
    else if (biMolData.Dtot < 0.1)
        dtmp = biMolData.Dtot * 100;
    else
        dtmp = biMolData.Dtot * 100;

    int d_ones = int(round(dtmp));

    if (biMolData.Dtot < 0.0001)
        biMolData.Dtot = d_ones * 0.00001;
    else if (biMolData.Dtot < 0.001)
        biMolData.Dtot = d_ones * 0.0001;
    else if (biMolData.Dtot < 0.01)
        biMolData.Dtot = d_ones * 0.001;
    else if (biMolData.Dtot < 0.1)
        biMolData.Dtot = d_ones * 0.01;
    else
        biMolData.Dtot = d_ones * 0.01;

    if (biMolData.Dtot < 1E-50)
        biMolData.Dtot = 0;
}

/*! \brief Adds the rotational contribution to `Dtot` for a pair in solution.
 *
 * Each partner is treated separately, because one of the two may be stuck to
 * the membrane while the other is free: a partner whose complex cannot move in
 * z rotates about one axis (2 Dr t, averaged over 4 t), a free one about three
 * (4 Dr t, averaged over 6 t).
 *
 * \param[in] flatCutoff how small `D.z` has to be for a complex to count as
 * membrane-bound.  A parameter rather than a constant because the callers do
 * not agree: the explicit-lipid paths use 1E-10 and
 * `determine_3D_implicitlipid_reaction_probability` uses 1E-15.  Whether that
 * disagreement is deliberate is not recorded, so it is preserved.
 */
inline void add_3D_rotational_diffusion(
    BiMolData& biMolData, const std::vector<Complex>& complexList, const Parameters& params, double flatCutoff)
{
    double Dr1 {};
    if (std::abs(complexList[biMolData.com1Index].D.z - 0) < flatCutoff) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (6.0 * params.timeStep);
    }

    double Dr2;
    if (std::abs(complexList[biMolData.com2Index].D.z - 0) < flatCutoff) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (6.0 * params.timeStep);
    }
}

/*! @} */

/*!
 * \brief Gets the distance between two Molecule's Interfaces and determines if they are within Rmax, and can therefore
 * react.
 */
bool get_distance(int pro1, int pro2, int iface1, int iface2, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    double& sep, double& R1, double Rmax, std::vector<Complex>& complexList, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, bool isSphere);

double passocF(double r0, double tCurr, double Dtot, double bindRadius, double alpha, double cof);

double passocF_1D(double r0, double tCurr, double Dtot, double bindRadius, double ka);

double pirr_pfree_ratio_psF_1D(double rCurr, double r0, double tCurr,
                               double Dtot, double bindrad, double ka,
                               double ps_prev);

void determine_1D_bimolecular_reaction_probability(
    int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    BiMolData &biMolData, const Parameters &params,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    const std::vector<ForwardRxn> &forwardRxns,
    const std::vector<BackRxn> &backRxns);

void determine_2D_bimolecular_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    unsigned& DDTableIndex, double* tableIDs, BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices, std::vector<gsl_matrix*>& pirMatrices);

void determine_3D_bimolecular_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns);

void perform_bimolecular_state_change(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_bimolecular_state_change_box(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_bimolecular_state_change_sphere(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);

void perform_implicitlipid_state_change(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_implicitlipid_state_change_box(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_implicitlipid_state_change_sphere(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);