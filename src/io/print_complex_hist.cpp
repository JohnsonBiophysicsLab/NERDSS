#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "debug/debug.hpp"
#include "io/io.hpp"
#include "mpi/mpi_function.hpp"
#include "tracing.hpp"

std::string complex_composition_string(const Complex& com,
                                       const std::vector<MolTemplate>& molTemplateList,
                                       const std::vector<Molecule>& moleculeList,
                                       bool byState) {
  const int nTypes = molTemplateList.size();
  const bool hasImplicitLipid = nTypes > 0 && molTemplateList[0].isImplicitLipid;
  std::string line;

  // One "name: n. " entry.  states holds the state identity of each interface
  // that has more than one state, in interface order, and is empty when the
  // entry is by type alone.  Named like the species columns of
  // copy_numbers_time.dat: A(ser~P), or A(ser~P,thr~U) with two such interfaces.
  auto emit = [&](int typeIndex, const std::vector<char>& states, int n) {
    const MolTemplate& molTemplate = molTemplateList[typeIndex];
    line += molTemplate.molName;
    if (!states.empty()) {
      line += '(';
      size_t stateItr = 0;
      for (const auto& iface : molTemplate.interfaceList) {
        if (iface.stateList.size() <= 1) continue;
        if (stateItr > 0) line += ',';
        line += iface.name + '~' + states[stateItr++];
      }
      line += ')';
    }
    line += ": " + std::to_string(n) + ". ";
  };

  // An implicit lipid is one Molecule standing for every lipid; the complexes
  // bound to it count it through linksToSurface, not through memberList, and
  // its own complex (the one member it is) reports nothing.  The links carry
  // no state, so its entry is the plain type count either way.
  if (!byState) {
    // The tally the histogram always wrote: one count per type, in type order.
    std::vector<int> numEachMol(nTypes, 0);
    for (int memberIndex : com.memberList) numEachMol[moleculeList[memberIndex].molTypeIndex] += 1;
    if (hasImplicitLipid) numEachMol[0] = com.linksToSurface;
    for (int t = 0; t < nTypes; ++t) {
      if (numEachMol[t] != 0) emit(t, {}, numEachMol[t]);
    }
    return line;
  }

  // One entry per distinct member kind: the type and the states of its
  // multi-state interfaces.  std::map orders the entries by type index first
  // and then by those states, so the line is the same for the same composition
  // whatever order the members were added in.
  using Kind = std::pair<int, std::vector<char>>;
  std::map<Kind, int> counts;
  for (int memberIndex : com.memberList) {
    const Molecule& mol = moleculeList[memberIndex];
    if (hasImplicitLipid && mol.molTypeIndex == 0) continue;
    Kind kind{mol.molTypeIndex, {}};
    const auto& ifaceTemplates = molTemplateList[mol.molTypeIndex].interfaceList;
    for (size_t k = 0; k < ifaceTemplates.size(); ++k) {
      if (ifaceTemplates[k].stateList.size() > 1)
        kind.second.push_back(mol.interfaceList[k].stateIden);
    }
    counts[kind] += 1;
  }
  if (hasImplicitLipid && com.linksToSurface != 0) emit(0, {}, com.linksToSurface);
  for (const auto& entry : counts) emit(entry.first.first, entry.first.second, entry.second);
  return line;
}

double print_complex_hist(std::vector<Complex>& complexList,
                          std::ofstream& outfile, int it, Parameters params,
                          std::vector<MolTemplate>& molTemplateList,
                          std::vector<Molecule>& moleculeList,
                          int nImplicitLipids) {
  int n_complexes = complexList.size();
  std::unordered_map<std::string, int> complex_count_dict;

  for (auto i = 0; i < n_complexes; i++) {
    auto& com = complexList[i];
    if (com.isEmpty) {
      continue;
    }

    std::string components_string = complex_composition_string(
        com, molTemplateList, moleculeList, params.histogramMolStates);

    // update the count of this complex
    if(components_string!=""){
      if (complex_count_dict.find(components_string) == complex_count_dict.end()) {
        complex_count_dict[components_string] = 1;
      } else {
        complex_count_dict[components_string] += 1;
      }
    }
  }

  // write the complex counts to the output file
  outfile << "Time (s): "
          << (it - params.itrRestartFrom) * params.timeStep * 1E-6 +
                 params.timeRestartFrom
          << "\n";

  if(molTemplateList[0].isImplicitLipid){
    outfile << nImplicitLipids << '\t' << molTemplateList[0].molName << ": 1." << "\n";
  }

  for (auto& pair : complex_count_dict) {
    auto& components_string = pair.first;
    auto& count = pair.second;
    outfile << count << '\t' << components_string << "\n";
  }

  outfile << std::flush;

  return 0.0;
}
