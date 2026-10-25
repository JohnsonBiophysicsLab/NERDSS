#include "reactions/shared_reaction_functions.hpp"
#include "classes/iface_state_bits.hpp"

#ifdef NERDSS_CHECK_STATE_BITS
#include <cstdio>
#include <cstdlib>
#endif

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

/*!
 * \brief Compiles each rate state's requirements into masks on Molecule::ifaceStateBits.
 *
 * Every requirement on a side has to name one molecule type, the one whose
 * words that side is tested against; selection falls back to the interfaces
 * for a molecule of any other type.  A requirement naming another type, or an
 * interface the type does not have, can never hold for such a molecule, and
 * neither can two requirements on one interface that differ, so those states
 * are compiled to never match.  A state name the template does not list is the
 * one thing a word cannot decide, and leaves the list without bit states.
 */
void compile_bit_states(const RxnBase& rxn, RxnBase::CompiledRates& compiled)
{
    if (compiled.bitsGeneration == 0)
        return; // no layouts yet
    for (int side { 0 }; side < 2; ++side) {
        for (const auto& oneRate : rxn.rateList) {
            if (!oneRate.otherIfaceLists[side].empty()) {
                compiled.sideMolType[side] = oneRate.otherIfaceLists[side].front().molTypeIndex;
                break;
            }
        }
    }

    using BitState = RxnBase::CompiledRates::BitState;
    using iface_state_bits::Requirement;
    for (const auto& oneRate : rxn.rateList) {
        BitState bits;
        bool never { false };
        for (int side { 0 }; side < 2; ++side) {
            for (const RxnIface& req : oneRate.otherIfaceLists[side]) {
                if (req.molTypeIndex != compiled.sideMolType[side]) {
                    never = true;
                    continue;
                }
                std::uint64_t mask { 0 };
                std::uint64_t value { 0 };
                const Requirement compiledReq { iface_state_bits::compile_requirement(
                    req.molTypeIndex, req.relIfaceIndex, req.requiresInteraction, req.requiresState, mask, value) };
                if (compiledReq == Requirement::unlisted) {
                    compiled.bitStates.clear();
                    compiled.sideMolType[0] = -1;
                    compiled.sideMolType[1] = -1;
                    return;
                }
                if (compiledReq == Requirement::unsatisfiable) {
                    never = true;
                    continue;
                }
                if ((bits.mask[side] & mask) != 0 && (bits.value[side] & mask) != value)
                    never = true;
                bits.mask[side] |= mask;
                bits.value[side] |= value;
            }
        }
        if (never) {
            // A word never has the invalid bit set when it is tested, so this
            // state fails whatever the word holds.
            bits.mask[0] = iface_state_bits::invalid;
            bits.value[0] = iface_state_bits::invalid;
        }
        compiled.bitStates.push_back(bits);
    }
    compiled.bitsUsable = true;
}

/*!
 * \brief The selection from the two molecules' words; false when a word cannot
 * answer (an invalid word, or a molecule of another type than the masks).
 */
bool best_matching_rate_from_words(
    const RxnBase::CompiledRates& compiled, const Molecule* const sideMol[2], RateMatch& result)
{
    std::uint64_t word[2] { 0, 0 };
    for (int side { 0 }; side < 2; ++side) {
        if (compiled.sideMolType[side] < 0)
            continue; // no rate state asks anything of this molecule
        const Molecule& mol = *sideMol[side];
        if (mol.molTypeIndex != compiled.sideMolType[side] || (mol.ifaceStateBits & iface_state_bits::invalid) != 0)
            return false;
        word[side] = mol.ifaceStateBits;
    }

    // The scan's selection, unchanged: count every match, keep the first of
    // those with the most ancillary interfaces.
    for (std::size_t rateItr { 0 }; rateItr < compiled.stateCount; ++rateItr) {
        const RxnBase::CompiledRates::BitState& bits = compiled.bitStates[rateItr];
        if ((((word[0] ^ bits.value[0]) & bits.mask[0]) | ((word[1] ^ bits.value[1]) & bits.mask[1])) != 0)
            continue;

        ++result.matches;
        const std::size_t ancillaryCount { compiled.ancillaryCounts[rateItr] };
        if (result.bestRateIndex == -1 || ancillaryCount > result.mostAncillaryIfaces) {
            result.bestRateIndex = static_cast<int>(rateItr);
            result.mostAncillaryIfaces = ancillaryCount;
        }
    }
    return true;
}

/*! \brief The selection from the molecules' interfaces, read once per probe. */
RateMatch best_matching_rate_from_probes(const RxnBase::CompiledRates& compiled, const Molecule* const sideMol[2])
{
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

#ifdef NERDSS_CHECK_STATE_BITS
// Tallies for the checking build, reported when the program exits.
struct StateBitsTally {
    unsigned long long fromWords { 0 }, fromProbes { 0 };
    ~StateBitsTally()
    {
        std::fprintf(stderr, "STATE BITS CHECK: %llu selections from words, %llu from interfaces, 0 mismatches\n",
            fromWords, fromProbes);
    }
} stateBitsTally;
#endif

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
    compiled.bitsGeneration = iface_state_bits::generation();

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
    compile_bit_states(rxn, compiled);
}

RateMatch best_matching_rate_compiled(const RxnBase& rxn, int reactIndex1, int reactIndex2,
    const Molecule& reactMol1, const Molecule& reactMol2)
{
    const RxnBase::CompiledRates& compiled = rxn.compiledRates;
    if (!compiled.built || compiled.stateCount != rxn.rateList.size()
        || compiled.bitsGeneration != iface_state_bits::generation())
        compile_rate_list(rxn);

    // The callers pass the two reactant positions of a bimolecular reaction,
    // one each.  Anything else -- both molecules tested against one side --
    // is the scan's to answer.
    const bool inOrder { reactIndex1 == 0 && reactIndex2 == 1 };
    if (!compiled.usable || !(inOrder || (reactIndex1 == 1 && reactIndex2 == 0)))
        return best_matching_rate_scan(rxn.rateList, reactIndex1, reactIndex2, reactMol1, reactMol2);

    // otherIfaceLists[side] is tested against this molecule.
    const Molecule* const sideMol[2] { inOrder ? &reactMol1 : &reactMol2, inOrder ? &reactMol2 : &reactMol1 };

#ifdef NERDSS_CHECK_STATE_BITS
    iface_state_bits::verify(reactMol1, "best_matching_rate_compiled");
    iface_state_bits::verify(reactMol2, "best_matching_rate_compiled");
    const RateMatch fromProbes { best_matching_rate_from_probes(compiled, sideMol) };
    RateMatch fromWords {};
    if (compiled.bitsUsable && best_matching_rate_from_words(compiled, sideMol, fromWords)) {
        ++stateBitsTally.fromWords;
        if (fromWords.matches != fromProbes.matches || fromWords.bestRateIndex != fromProbes.bestRateIndex
            || fromWords.mostAncillaryIfaces != fromProbes.mostAncillaryIfaces) {
            std::fprintf(stderr,
                "STATE BITS MISMATCH in rate selection: rxn %d, positions (%d, %d), molecules %d and %d: words (%d, %d, "
                "%zu), interfaces (%d, %d, %zu)\n",
                rxn.absRxnIndex, reactIndex1, reactIndex2, reactMol1.index, reactMol2.index, fromWords.matches,
                fromWords.bestRateIndex, fromWords.mostAncillaryIfaces, fromProbes.matches, fromProbes.bestRateIndex,
                fromProbes.mostAncillaryIfaces);
            std::abort();
        }
    } else {
        ++stateBitsTally.fromProbes;
    }
    return fromProbes;
#else
    RateMatch result {};
    if (compiled.bitsUsable && best_matching_rate_from_words(compiled, sideMol, result))
        return result;
    return best_matching_rate_from_probes(compiled, sideMol);
#endif
}

bool best_matching_rate_from_state_bits(const RxnBase& rxn, int reactIndex1, int reactIndex2,
    const Molecule& reactMol1, const Molecule& reactMol2, RateMatch& result)
{
    const RxnBase::CompiledRates& compiled = rxn.compiledRates;
    if (!compiled.built || compiled.stateCount != rxn.rateList.size()
        || compiled.bitsGeneration != iface_state_bits::generation())
        compile_rate_list(rxn);
    const bool inOrder { reactIndex1 == 0 && reactIndex2 == 1 };
    if (!compiled.usable || !compiled.bitsUsable || !(inOrder || (reactIndex1 == 1 && reactIndex2 == 0)))
        return false;
    const Molecule* const sideMol[2] { inOrder ? &reactMol1 : &reactMol2, inOrder ? &reactMol2 : &reactMol1 };
    result = RateMatch {};
    return best_matching_rate_from_words(compiled, sideMol, result);
}
