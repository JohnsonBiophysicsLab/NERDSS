#include "reactions/shared_reaction_functions.hpp"

namespace {

/*!
 * \brief Tests one ancillary-interface requirement against a Molecule.
 *
 * RxnIface::relIfaceIndex is documented as an index into MolTemplate's (and
 * therefore Molecule's) interfaceList, and every site that populates a
 * Molecule sets interfaceList[i].relIndex = i.  So the requirement can be
 * looked up directly instead of scanning the whole interface list for it,
 * which turns the cost of hasIntangibles() from
 * O(requirements x interfaces) into O(requirements).
 *
 * The direct lookup is guarded two ways so a malformed model behaves exactly
 * as it did before: an out-of-range relIfaceIndex simply fails to match, and
 * if a Molecule ever violates the relIndex == position invariant the original
 * linear scan is used instead.
 */
bool hasAncillaryIface(const Molecule& reactMol, int molTypeIndex, int relIfaceIndex, bool requiresInteraction,
    char requiresState)
{
    if (relIfaceIndex >= 0 && static_cast<std::size_t>(relIfaceIndex) < reactMol.interfaceList.size()) {
        const Molecule::Iface& oneIface = reactMol.interfaceList[relIfaceIndex];
        if (oneIface.relIndex == relIfaceIndex) {
            // relIndex is unique within a Molecule, so this is the only
            // interface the scan below could ever have matched.
            return molTypeIndex == oneIface.molTypeIndex && requiresInteraction == oneIface.isBound
                && requiresState == oneIface.stateIden;
        }
    }

    for (const auto& oneIface : reactMol.interfaceList) {
        if (molTypeIndex == oneIface.molTypeIndex && relIfaceIndex == oneIface.relIndex
            && requiresInteraction == oneIface.isBound && requiresState == oneIface.stateIden) {
            return true;
        }
    }
    return false;
}

bool hasAncillaryIface(const Molecule& reactMol, const RxnIface& anccIface)
{
    return hasAncillaryIface(reactMol, anccIface.molTypeIndex, anccIface.relIfaceIndex,
        anccIface.requiresInteraction, anccIface.requiresState);
}

/*!
 * \brief The fields an ancillary requirement compares, packed so that one
 * integer compare tests all three.  Each field has bits of its own, so two keys
 * are equal exactly when all three fields are.
 */
std::uint64_t requirement_key(int molTypeIndex, bool isBound, char state)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(molTypeIndex)) << 16)
        | (static_cast<std::uint64_t>(isBound) << 8) | static_cast<std::uint64_t>(static_cast<unsigned char>(state));
}

} // namespace

bool hasIntangibles(int reactantIndex, const Molecule& reactMol, const RxnBase::RateState& currRxnState)
{
    for (const auto& anccIface : currRxnState.otherIfaceLists[reactantIndex]) {
        // if the ancillary interface was found in either of the two reactants, iterate the number of matches
        // TODO: Does it matter if there are matches in both reactants?
        if (!hasAncillaryIface(reactMol, anccIface))
            return false;
    }
    return true;
}

bool hasIntangibles(int reactIndex1, int reactIndex2, const Molecule& reactMol1, const Molecule& reactMol2,
    const RxnBase::RateState& currRxnState)
{
    if (currRxnState.otherIfaceLists[0].empty() && currRxnState.otherIfaceLists[1].empty()) {
        return true;
    }

    for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex1]) {
        if (!hasAncillaryIface(reactMol1, anccIface))
            return false;
    }
    for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex2]) {
        if (!hasAncillaryIface(reactMol2, anccIface))
            return false;
    }

    return true;
}

void compile_rate_list(const RxnBase& rxn)
{
    RxnBase::CompiledRates& compiled = rxn.compiledRates;
    compiled.reset();
    compiled.built = true;
    compiled.stateCount = rxn.rateList.size();

    // The scan indexes otherIfaceLists[0] and [1] unconditionally; a list
    // shaped any other way is left to it.
    for (const auto& oneRate : rxn.rateList) {
        if (oneRate.otherIfaceLists.size() != 2)
            return;
    }

    // Distinct requirements, grouped by the interface they read, so that a
    // selection loads each interface once.  Probes and tests keep the order of
    // first appearance; any order would do, since a mask is order-free.
    using Test = RxnBase::CompiledRates::Test;
    using Probe = RxnBase::CompiledRates::Probe;
    std::vector<std::vector<Test>> testsByProbe;
    for (const auto& oneRate : rxn.rateList) {
        for (int side { 0 }; side < 2; ++side) {
            for (const RxnIface& req : oneRate.otherIfaceLists[side]) {
                std::size_t p { 0 };
                while (p < compiled.probes.size()
                    && !(compiled.probes[p].side == side && compiled.probes[p].relIfaceIndex == req.relIfaceIndex))
                    ++p;
                if (p == compiled.probes.size()) {
                    Probe probe;
                    probe.side = side;
                    probe.relIfaceIndex = req.relIfaceIndex;
                    compiled.probes.push_back(probe);
                    testsByProbe.emplace_back();
                }
                const std::uint64_t key { requirement_key(
                    req.molTypeIndex, req.requiresInteraction, req.requiresState) };
                bool seen { false };
                for (const Test& t : testsByProbe[p])
                    seen = seen || t.key == key;
                if (!seen) {
                    Test test;
                    test.key = key;
                    test.molTypeIndex = req.molTypeIndex;
                    test.relIfaceIndex = req.relIfaceIndex;
                    test.requiresState = req.requiresState;
                    test.requiresInteraction = req.requiresInteraction;
                    testsByProbe[p].push_back(test);
                }
            }
        }
    }
    for (std::size_t p { 0 }; p < compiled.probes.size(); ++p) {
        compiled.probes[p].firstTest = static_cast<int>(compiled.tests.size());
        compiled.tests.insert(compiled.tests.end(), testsByProbe[p].begin(), testsByProbe[p].end());
        compiled.probes[p].endTest = static_cast<int>(compiled.tests.size());
    }
    if (compiled.tests.size() > 64) {
        compiled.probes.clear();
        compiled.tests.clear();
        return;
    }

    for (const auto& oneRate : rxn.rateList) {
        std::uint64_t mask { 0 };
        for (int side { 0 }; side < 2; ++side) {
            for (const RxnIface& req : oneRate.otherIfaceLists[side]) {
                const std::uint64_t key { requirement_key(
                    req.molTypeIndex, req.requiresInteraction, req.requiresState) };
                for (const Probe& probe : compiled.probes) {
                    if (probe.side != side || probe.relIfaceIndex != req.relIfaceIndex)
                        continue;
                    for (int t { probe.firstTest }; t < probe.endTest; ++t)
                        if (compiled.tests[t].key == key)
                            mask |= std::uint64_t { 1 } << t;
                }
            }
        }
        compiled.stateMasks.push_back(mask);
        compiled.ancillaryCounts.push_back(oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size());
    }
    compiled.usable = true;
}

RateMatch best_matching_rate_compiled(const RxnBase& rxn, int reactIndex1, int reactIndex2,
    const Molecule& reactMol1, const Molecule& reactMol2)
{
    const RxnBase::CompiledRates& compiled = rxn.compiledRates;
    if (!compiled.built || compiled.stateCount != rxn.rateList.size())
        compile_rate_list(rxn);

    // The callers pass the two reactant positions of a bimolecular reaction,
    // one each.  Anything else -- both molecules tested against one side --
    // is the scan's to answer.
    const bool inOrder { reactIndex1 == 0 && reactIndex2 == 1 };
    if (!compiled.usable || !(inOrder || (reactIndex1 == 1 && reactIndex2 == 0)))
        return best_matching_rate_scan(rxn.rateList, reactIndex1, reactIndex2, reactMol1, reactMol2);

    // otherIfaceLists[side] is tested against this molecule.
    const Molecule* sideMol[2] { inOrder ? &reactMol1 : &reactMol2, inOrder ? &reactMol2 : &reactMol1 };

    std::uint64_t satisfied { 0 };
    for (const auto& probe : compiled.probes) {
        const Molecule& mol = *sideMol[probe.side];
        // The same direct lookup and fallback as hasAncillaryIface(), with the
        // interface loaded once for all the requirements on it.
        if (probe.relIfaceIndex >= 0 && static_cast<std::size_t>(probe.relIfaceIndex) < mol.interfaceList.size()
            && mol.interfaceList[probe.relIfaceIndex].relIndex == probe.relIfaceIndex) {
            const Molecule::Iface& oneIface = mol.interfaceList[probe.relIfaceIndex];
            const std::uint64_t key { requirement_key(oneIface.molTypeIndex, oneIface.isBound, oneIface.stateIden) };
            for (int t { probe.firstTest }; t < probe.endTest; ++t)
                satisfied |= static_cast<std::uint64_t>(compiled.tests[t].key == key) << t;
        } else {
            for (int t { probe.firstTest }; t < probe.endTest; ++t) {
                const auto& test = compiled.tests[t];
                if (hasAncillaryIface(
                        mol, test.molTypeIndex, test.relIfaceIndex, test.requiresInteraction, test.requiresState))
                    satisfied |= std::uint64_t { 1 } << t;
            }
        }
    }

    // The scan's selection, unchanged: count every match, keep the first of
    // those with the most ancillary interfaces.
    RateMatch result {};
    for (std::size_t rateItr { 0 }; rateItr < compiled.stateCount; ++rateItr) {
        if ((compiled.stateMasks[rateItr] & ~satisfied) != 0)
            continue;

        ++result.matches;
        const std::size_t ancillaryCount { compiled.ancillaryCounts[rateItr] };
        if (result.bestRateIndex == -1 || ancillaryCount > result.mostAncillaryIfaces) {
            result.bestRateIndex = static_cast<int>(rateItr);
            result.mostAncillaryIfaces = ancillaryCount;
        }
    }
    return result;
}
