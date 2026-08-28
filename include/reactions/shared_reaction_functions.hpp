/*! \file shared_reaction_functions.hpp

 * ### Created on 11/6/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#pragma once

#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "gsl/gsl_matrix.h"
#include "math/constants.hpp"

#include <algorithm>
#include <cmath>

/*!
 * \brief Determines if the Interface::State of the reactant is equivalent to the reactant as contained in the
 * reaction
 *
 * \param reactant Interface::State of the reacting interface
 * \param tempReactant a reactant of the ForwardRxn/BackRxn/CreateDestructRxn being evaluated
 */
bool isReactant(const Molecule::Iface& reactIface, const Molecule& reactMol, const RxnIface& tempReactant);

/*!
 * \brief Finds a bond in `bndlist` by the interface that carries it.
 *
 * `Molecule::bndlist` holds the relative indices of this molecule's bound
 * interfaces and `Molecule::bndpartner` the molecule each one is bound to, so
 * entry *k* of one describes the same bond as entry *k* of the other.  An
 * interface can be bound to at most one partner, which makes `bndlist` the only
 * key that identifies a bond uniquely -- `bndpartner` can hold the same
 * molecule more than once, when two interfaces bind the same neighbour.
 *
 * Every erase must therefore locate the bond through `bndlist` and then act on
 * both vectors at that position.  Three sites used to do it three different
 * ways, and one of them left the two permuted; see
 * `docs/bond_bookkeeping_defects.md`.
 *
 * \return the position of the bond, or `bndlist.size()` if the interface is not
 * bound.  Callers must check: erasing at a not-found position is what the
 * previous implicit-lipid code did, by calling `erase()` on `end()`.
 */
size_t find_bond_slot(const Molecule& mol, int relIface);

/*!
 * \brief Erases one bond from `bndlist` and `bndpartner` together.
 *
 * Does nothing if `relIface` is not bound, which is the guard the previous
 * `erase(std::find_if(...))` call sites lacked.
 *
 * \return true if a bond was erased.
 */
bool erase_bond(Molecule& mol, int relIface);

/*!
 * \brief Records that two molecules are close enough to interact this timestep.
 *
 * Each of them appends one `Molecule::CrossEntry` naming the other, the
 * interface it presents and the reaction they might do; each parent complex
 * bumps `ncross`.  Six sites wrote this out: `get_distance()`, and five in the
 * volume-exclusion half of `check_bimolecular_reactions()`.
 *
 * The entry's probability starts at zero. Callers that evaluate one overwrite it
 * through `crossings.back().prob`; the exclusion sites leave it at zero, which
 * is what they used to push explicitly.
 *
 * \param[in] crossRxn the `{rxnIndex, rateIndex, isStateChangeBackRxn}` triple.
 * Every copy pushed the same value onto both molecules, so it is one parameter.
 * The exclusion sites pass `{rxnIndex, 0, false}`, because no probability is
 * evaluated for them.
 */
void record_crossing_pair(int pro1, int pro2, int relIface1, int relIface2,
    const std::array<int, 3>& crossRxn, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList);

/*!
 * \brief Zeroes the reaction probability every partner recorded for `mol`.
 *
 * Called once a reaction involving `mol` has fired, so that `mol` cannot react
 * again this timestep while its partners still avoid overlapping it.
 *
 * \param[in] skipImplicitLipidPartners leave the implicit lipid's own entry
 * alone.  The explicit-lipid association paths zero every partner; the
 * implicit-lipid ones skip it, and that difference is theirs to state.
 */
void zero_partner_probvec(
    const Molecule& mol, std::vector<Molecule>& moleculeList, bool skipImplicitLipidPartners = false);

/*!
 * \brief Moves the observable counter for a bimolecular state change.
 *
 * Up for a forward reaction, down for the back reaction; nothing at all if the
 * reaction is not observed or the label was never declared.
 */
void update_state_change_observable(bool isStateChangeBackRxn, int rxnIndex, int backRxnIndex,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    std::map<std::string, int>& observablesList);

/*!
 * \brief Counts an observed unimolecular reaction, in whichever direction fired.
 *
 * Distinct from \ref update_state_change_observable: this indexes both reaction
 * lists with the same `rxnIndex`, and counts up for the back reaction as well
 * as the forward one.  That is what the four unimolecular sites have always
 * done, so it is kept as its own function rather than merged with the other.
 */
void count_unimolecular_observable(bool isStateChangeBackRxn, int rxnIndex,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    std::map<std::string, int>& observablesList);

/*!
 * \brief Determines if the Molecule is a reactant for a Creation/Destruction reaction.
 */
bool isReactant(const Molecule& currMol, const Complex& currCom, const CreateDestructRxn& currRxn,
    const std::vector<Molecule>& moleculeList);

/*!
 * \brief Determines if the associating Molecules have the ancillary interfaces required by this reaction.
 */
bool hasIntangibles(int reactIndex1, int reactIndex2, const Molecule& reactMol1, const Molecule& reactMol2,
    const RxnBase::RateState& currRxnState);

/*!
 * \brief Determines if the Molecule in a unimolecular reaction has he ancillary interaces required
 */
bool hasIntangibles(int reactantIndex, const Molecule& reactMol, const RxnBase::RateState& currRxnState);

/*!
 * \brief The single best matching RxnBase::RateState found in a rate list.
 *
 * \var matches how many rate states accepted the reactant pair
 * \var bestRateIndex index into the scanned rate list, -1 when nothing matched
 * \var mostAncillaryIfaces ancillary-interface count of the selected match
 */
struct RateMatch {
    int matches { 0 };
    int bestRateIndex { -1 };
    std::size_t mostAncillaryIfaces { 0 };
};

/*!
 * \brief Scans a rate list once and keeps only the best matching rate state.
 *
 * When several rate states accept the reactant pair, the one requiring the most
 * ancillary interfaces wins, because it is the most specific description of the
 * pair.  The comparison is strict, so the first of several equally specific
 * matches is kept.
 *
 * This replaces the temporary match vector that find_which_reaction() and
 * find_reaction_rate_state() used to build and then rescan.  Both are called
 * once per candidate interface pair per timestep, so the vector cost a heap
 * allocation on a hot path to hold information that a running best already
 * carries.
 */
inline RateMatch best_matching_rate(const std::vector<RxnBase::RateState>& rateList, int reactIndex1, int reactIndex2,
    const Molecule& reactMol1, const Molecule& reactMol2)
{
    RateMatch result {};
    for (std::size_t rateItr { 0 }; rateItr < rateList.size(); ++rateItr) {
        const RxnBase::RateState& oneRate = rateList[rateItr];
        if (!hasIntangibles(reactIndex1, reactIndex2, reactMol1, reactMol2, oneRate))
            continue;

        ++result.matches;
        const std::size_t ancillaryCount { oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size() };
        if (result.bestRateIndex == -1 || ancillaryCount > result.mostAncillaryIfaces) {
            result.bestRateIndex = static_cast<int>(rateItr);
            result.mostAncillaryIfaces = ancillaryCount;
        }
    }
    return result;
}

/*!
 * \brief Probability that a first-order event of this rate occurs in one timestep.
 *
 * \param rate as stored in `RxnBase::RateState::rate`, in s^-1.
 * \param timeStep `Parameters::timeStep`, in us.
 *
 * An infinite rate returns certainty directly rather than arriving at it through
 * `exp()` underflow.  The value is the same either way -- `1 - exp(-INF)` is also
 * 1 -- but this way the product `rate * timeStep` is never formed, and that
 * product is where an infinite rate turns dangerous: `INF * 0` is NaN, every
 * caller decides either with `prob > rand_gsl()` or by passing `prob` to
 * `gsl_ran_binomial()`, and both read NaN as "no event".  A rate meaning
 * "always" would silently become one meaning "never".  A NaN rate lands in the
 * same trap, and `-INF` would reach `gsl_ran_binomial()` as a probability of
 * `-INF`; both return 0 here instead.
 *
 * A *finite* negative rate is left to the arithmetic, because it is a live
 * sentinel: `check_for_unimolecular_reactions_population.cpp` uses `rate == -1`
 * to select strong destruction and overwrites the event count afterwards.
 *
 * For every finite rate this is the arithmetic the eleven call sites performed
 * inline, unchanged, so existing models are unaffected.
 */
inline double first_order_probability(double rate, double timeStep)
{
    if (!std::isfinite(rate))
        return rate > 0.0 ? 1.0 : 0.0;

    return 1.0 - std::exp(-rate * timeStep * Constants::usToSeconds);
}

/* FUNCTIONS TO DETERMINE WHICH REACTION TO PERFORM */

/*!
 * \brief Determines which RxnState to use for the association/dissociation reaction.
 */
int find_reaction_rate_state(int simItr, int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1,
    const Molecule& reactMol2, const BackRxn& backRxn,
    const std::vector<MolTemplate>& molTemplateList);

size_t find_reaction_rate_state(const Molecule& reactMol1, const Molecule& reactMol2, const ForwardRxn& forwardRxn,
    const std::vector<MolTemplate>& molTemplateList);

/*!
 * \brief This function determines which reaction to use based on the identities of the reacting Molecules and
 * interfaces, and returns its reaction index and rate index
 *
 * \param[out] std::array<size_t, 2> array of two indices, [index of ForwardRxn, index of RxnBase::RateState]
 */
void find_which_reaction(int ifaceIndex1, int ifaceIndex2, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Interface::State& currState, const Molecule& reactMol1, const Molecule& reactMol2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<MolTemplate>& molTemplateList);

/*!
 * \brief Determines which state change reaction to use based on the identity of the current Interface::State of the
 * target Interface.
 *
 * \params[out] std::array<int, 3> [rxnIndex, rateIndex, rxnType], where rxnType = 0 if forward, 1 if back
 */
void find_which_state_change_reaction(int ifaceIndex, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Molecule& reactMol, const Interface::State& currState,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

std::array<int, 2> find_which_reaction(int ifaceIndex, const Molecule& reactMol, const Interface::State& currState,
    const std::vector<BackRxn>& forwardRxns);

/*!
 * \brief
 */
// double passocF(double r0, double tCurr, double Dtot, double bindRadius, double alpha, double cof);

extern int evalBindNum;
/*!
 * \brief Main function for evaluating the potential interactions between two Molecules.
 *
 * \param[in] pro1Index Index of first protein in moleculeList.
 * \param[in] pro2Index Index of second protein in moleculeList.
 * \param[in] simItr Current simulation iteration.
 * \param[in] tableIDs C-style array containing the rates and their total diffusion constant.
 * \param[in] DDTableIndex Current index in tableIDs
 * \param[in] params User-provided simulation Parameters.
 * \param[in] normMatrices List of all previously calculated normMatrix
 * \param[in] survMatrices List of all previously calculated survMatrix
 * \param[in] pirMatrices List of all previously calculated pirMatrix
 * \param[in] moleculeList List of all Molecules in the system.
 * \param[in] complexList List of all Complexes in the system.
 * \param[in] molTemplateList List of all user-provided MolTemplates.
 * \param[in] forwardRxns List of all user-provided ForwardRxns.
 *
 * If 2D, also calculates/looks up values for 2D reaction tables (normMatrix, survMatrix, and pirMatrix)
 */
void check_bimolecular_reactions(int pro1Index, int pro2Index, int simItr, double* tableIDs, unsigned& DDTableIndex,
    const Parameters& params, std::vector<gsl_matrix*>& normMatrices, std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject);

/*!
 * \brief Determines if binding of two molecules within the same complex can occur.
 *
 * If the two molecules are within Rmax, probability of reaction is set to 1.0. If not, it is set to 0.
 */

void evaluate_binding_within_complex(int pro1Index, int pro2Index, int iface1Index, int iface2Index, int rxnIndex,
    int rateIndex, bool isBiMolStateChange, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const ForwardRxn& oneRxn,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, copyCounters& counterArrays);

bool determine_if_reaction_occurs(int& crossIndex1, int& crossIndex2, const double maxRandInt, Molecule& mol,
    std::vector<Molecule>& moleculeList, const std::vector<ForwardRxn>& forwardRxns);

void update_Nboundpairs(int ptype1, int ptype2, int chg, const Parameters& params, copyCounters& counterArrays);

void check_implicit_reactions(int pro1Index, int pro2Index, int simItr,
    const Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs);

void check_perform_zeroth_first_order_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<double>& IL2DbindingVec,
    std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs,
    MpiContext& mpiContext);

void measure_separations_to_identify_possible_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<double>& IL2DbindingVec,
    std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs,
    std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, int implicitlipidIndex,
    double* tableIDs, unsigned& DDTableIndex);

void perform_bimolecular_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<int>& region);

void check_overlap(std::vector<int>& region, unsigned simItr,
                   Parameters& params, std::vector<Molecule>& moleculeList,
                   std::vector<Complex>& complexList, SimulVolume& simulVolume,
                   std::vector<ForwardRxn>& forwardRxns,
                   std::vector<BackRxn>& backRxns,
                   std::vector<CreateDestructRxn>& createDestructRxns,
                   std::vector<MolTemplate>& molTemplateList,
                   std::map<std::string, int>& observablesList,
                   copyCounters& counterArrays, Membrane& membraneObject,
                   MpiContext& mpiContext);

void remove_empty_slots(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, MpiContext& mpiContext);
    