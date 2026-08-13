#include "system_setup/system_setup.hpp"

bool areInVicinity(const Molecule& mol1, const Molecule& mol2,
                   const std::vector<MolTemplate>& molTemplateList) {
  Vec3D tmpVec{mol2.comCoord - mol1.comCoord};
  /*this does not work for points that still need to avoid each other.
    To be effective, it should replace molecule radius with the largest binding
    radius for that interface.
   */
  return tmpVec.length() < (molTemplateList[mol1.molTypeIndex].radius +
                             molTemplateList[mol2.molTypeIndex].radius + 2.0);
}
