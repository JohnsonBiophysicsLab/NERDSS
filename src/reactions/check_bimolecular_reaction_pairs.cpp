#include "reactions/shared_reaction_functions.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>

#ifdef NERDSS_USE_OPENMP
#include <omp.h>
#endif

namespace {

bool pair_is_parallel_3d_safe(
    const BimolecularReactionPair& pair,
    const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList) {
  if (pair.pro1Index < 0 || pair.pro2Index < 0 ||
      static_cast<std::size_t>(pair.pro1Index) >= moleculeList.size() ||
      static_cast<std::size_t>(pair.pro2Index) >= moleculeList.size()) {
    return false;
  }

  const Molecule& molecule1 = moleculeList[pair.pro1Index];
  const Molecule& molecule2 = moleculeList[pair.pro2Index];
  if (molecule1.isEmpty || molecule2.isEmpty || molecule1.isImplicitLipid ||
      molecule2.isImplicitLipid || molecule1.isGhosted ||
      molecule2.isGhosted) {
    return false;
  }

  if (molecule1.myComIndex < 0 || molecule2.myComIndex < 0 ||
      static_cast<std::size_t>(molecule1.myComIndex) >= complexList.size() ||
      static_cast<std::size_t>(molecule2.myComIndex) >= complexList.size() ||
      molecule1.myComIndex == molecule2.myComIndex) {
    return false;
  }

  const Complex& complex1 = complexList[molecule1.myComIndex];
  const Complex& complex2 = complexList[molecule2.myComIndex];
  if ((complex1.onFiber && complex2.onFiber) ||
      (complex1.OnSurface && complex2.OnSurface)) {
    return false;
  }

  const MolTemplate& template1 = molTemplateList[molecule1.molTypeIndex];
  const MolTemplate& template2 = molTemplateList[molecule2.molTypeIndex];
  const bool hasExcludeVolumeMutation =
      (!molecule1.bndlist.empty() && template1.excludeVolumeBound) ||
      (!molecule2.bndlist.empty() && template2.excludeVolumeBound);
  return !hasExcludeVolumeMutation;
}

std::size_t minimum_parallel_pair_count() {
  static const std::size_t value = [] {
    const char* text = std::getenv("NERDSS_OMP_MIN_PAIRS");
    if (text == nullptr || *text == '\0') return std::size_t{64};
    const unsigned long parsed = std::strtoul(text, nullptr, 10);
    return std::max<std::size_t>(1, parsed);
  }();
  return value;
}

}  // namespace

void check_bimolecular_reaction_pairs(
    const std::vector<BimolecularReactionPair>& reactionPairs, int simItr,
    double* tableIDs, unsigned& DDTableIndex, const Parameters& params,
    std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays,
    Membrane& membraneObject) {
  auto evaluatePair = [&](const BimolecularReactionPair& pair) {
    check_bimolecular_reactions(
        pair.pro1Index, pair.pro2Index, simItr, tableIDs, DDTableIndex, params,
        normMatrices, survMatrices, pirMatrices, moleculeList, complexList,
        molTemplateList, forwardRxns, backRxns, counterArrays,
        membraneObject);
  };

#ifndef NERDSS_USE_OPENMP
  for (const auto& pair : reactionPairs) evaluatePair(pair);
#else
  if (reactionPairs.size() < minimum_parallel_pair_count() ||
      omp_get_max_threads() <= 1) {
    for (const auto& pair : reactionPairs) evaluatePair(pair);
    return;
  }

  std::vector<int> lastMoleculeWave(moleculeList.size(), -1);
  std::vector<int> lastComplexWave(complexList.size(), -1);
  std::vector<std::vector<BimolecularReactionPair>> waves;
  int lastGloballyOrderedWave = -1;

  for (const auto& pair : reactionPairs) {
    const int complex1 = moleculeList[pair.pro1Index].myComIndex;
    const int complex2 = moleculeList[pair.pro2Index].myComIndex;
    int wave = 0;
    wave = std::max(wave, lastMoleculeWave[pair.pro1Index] + 1);
    wave = std::max(wave, lastMoleculeWave[pair.pro2Index] + 1);
    if (complex1 >= 0 && static_cast<std::size_t>(complex1) < complexList.size())
      wave = std::max(wave, lastComplexWave[complex1] + 1);
    if (complex2 >= 0 && static_cast<std::size_t>(complex2) < complexList.size())
      wave = std::max(wave, lastComplexWave[complex2] + 1);

    const bool globallyOrdered = !pair_is_parallel_3d_safe(
        pair, moleculeList, complexList, molTemplateList);
    if (globallyOrdered)
      wave = std::max(wave, lastGloballyOrderedWave + 1);

    if (static_cast<std::size_t>(wave) >= waves.size())
      waves.resize(static_cast<std::size_t>(wave) + 1);
    waves[wave].push_back(pair);

    lastMoleculeWave[pair.pro1Index] = wave;
    lastMoleculeWave[pair.pro2Index] = wave;
    if (complex1 >= 0 && static_cast<std::size_t>(complex1) < complexList.size())
      lastComplexWave[complex1] = wave;
    if (complex2 >= 0 && static_cast<std::size_t>(complex2) < complexList.size())
      lastComplexWave[complex2] = wave;
    if (globallyOrdered) lastGloballyOrderedWave = wave;
  }

#pragma omp parallel
  {
    for (std::size_t wave = 0; wave < waves.size(); ++wave) {
#pragma omp for schedule(static)
      for (std::int64_t index = 0;
           index < static_cast<std::int64_t>(waves[wave].size()); ++index) {
        evaluatePair(waves[wave][static_cast<std::size_t>(index)]);
      }
    }
  }
#endif
}
