#include "system_setup/system_setup.hpp"

void set_rMaxLimit(Parameters &params,
                   const std::vector<MolTemplate> &molTemplateList,
                   const std::vector<ForwardRxn> &forwardRxns) {
  /*For each reaction, need distance from the interface to the COM for both
   * partners, plus the bindrad+sqrt(6*Dtot*deltat)
   */
  double rMaxTot{0};
  params.rMaxLimit = 0.0;
  for (auto &oneRxn : forwardRxns) {
    if (oneRxn.rxnType == ReactionType::bimolecular) {
      const RxnIface &rxnIface1 = oneRxn.reactantListNew.at(0);
      const RxnIface &rxnIface2 = oneRxn.reactantListNew.at(1);

      // MolTemplates of the interfaces involved in the reaction
      const MolTemplate &pro1Temp = molTemplateList.at(rxnIface1.molTypeIndex);
      const MolTemplate &pro2Temp = molTemplateList.at(rxnIface2.molTypeIndex);

      // A reactant's absIfaceIndex is the index of one of its interface's
      // states.  An add file numbers its templates' states after the existing
      // ones and renumbers only products, so absToRelIface maps a reactant
      // straight to its interface, whichever file its template came from.
      const Interface &iface1 = pro1Temp.interfaceList.at(
          MolTemplate::absToRelIface.at(rxnIface1.absIfaceIndex));
      const Interface &iface2 = pro2Temp.interfaceList.at(
          MolTemplate::absToRelIface.at(rxnIface2.absIfaceIndex));

      double scal{1.0 / 3.0};
      double Dtot{(scal * (pro1Temp.D.x + pro2Temp.D.x)) +
                  (scal * (pro1Temp.D.y + pro2Temp.D.y)) +
                  (scal * (pro1Temp.D.z + pro2Temp.D.z))};
      // TODO: add rotational diffusion contribution

      /*Now calculate distance from the interface to the protein COM.*/
      double pro1R1{0};
      double pro2R1{0};
      // pro1
      if (pro1Temp.isPromoter) {
        double R2 = (iface1.iCoord.x * iface1.iCoord.x);
        pro1R1 = sqrt(R2);
      }
      // if (std::abs(pro1Temp.D.z) < 1E-10) {
      if (pro1Temp.isImplicitLipid || pro1Temp.isLipid) {
        double R2 = (iface1.iCoord.x * iface1.iCoord.x) +
                    (iface1.iCoord.y * iface1.iCoord.y);
        pro1R1 = sqrt(R2);
        double einsStks = cos(sqrt(2.0 * pro1Temp.Dr.z * params.timeStep));
        double Dr = 2.0 * R2 * (1.0 - einsStks);
        Dtot += Dr / (4.0 * params.timeStep);
      } else {
        double R2 = (iface1.iCoord.x * iface1.iCoord.x) +
                    (iface1.iCoord.y * iface1.iCoord.y) +
                    (iface1.iCoord.z * iface1.iCoord.z);
        pro1R1 = sqrt(R2);
        double einsStks = cos(sqrt(4.0 * pro1Temp.Dr.z * params.timeStep));
        double Dr = 2.0 * R2 * (1.0 - einsStks);
        Dtot += Dr / (6.0 * params.timeStep);
      }

      // pro2
      if (pro2Temp.isPromoter) {
        double R2 = (iface2.iCoord.x * iface2.iCoord.x);
        pro2R1 = sqrt(R2);
      }
      // if (std::abs(pro2Temp.D.z) < 1E-10) {
      else if (pro2Temp.isImplicitLipid || pro2Temp.isLipid) {
        double R2 = (iface2.iCoord.x * iface2.iCoord.x) +
                    (iface2.iCoord.y * iface2.iCoord.y);
        pro2R1 = sqrt(R2);
        double einsStks = cos(sqrt(2.0 * pro2Temp.Dr.z * params.timeStep));
        double Dr = 2.0 * R2 * (1.0 - einsStks);
        Dtot += Dr / (4.0 * params.timeStep);
      } else {
        double R2 = (iface2.iCoord.x * iface2.iCoord.x) +
                    (iface2.iCoord.y * iface2.iCoord.y) +
                    (iface2.iCoord.z * iface2.iCoord.z);
        pro2R1 = sqrt(R2);
        double einsStks = cos(sqrt(4.0 * pro2Temp.Dr.z * params.timeStep));
        double Dr = 2.0 * R2 * (1.0 - einsStks);
        Dtot += Dr / (6.0 * params.timeStep);
      }

      double RmaxDiff =
          3.0 * sqrt(6.0 * Dtot * params.timeStep) + oneRxn.bindRadius;
      // A pair of fiber molecules (both templates promoters; pair_dim() gives
      // Fiber1D for their complexes) is evaluated by the 1D rules, which reach
      // 4 sqrt(2 (D1.x + D2.x) dt) + sigma along the fiber: with D.y = D.z = 0
      // that is 4/3 of the three-axis estimate above.  The cell grid and the
      // pre-cut in check_bimolecular_reactions() both rest on rMaxLimit being
      // the longest reach, so take the larger of the two.  Reactions without
      // promoters are untouched.
      if (pro1Temp.isPromoter && pro2Temp.isPromoter) {
        const double RmaxDiff1D =
            4.0 * sqrt(2.0 * (pro1Temp.D.x + pro2Temp.D.x) * params.timeStep) +
            oneRxn.bindRadius;
        if (RmaxDiff1D > RmaxDiff)
          RmaxDiff = RmaxDiff1D;
      }
      rMaxTot = RmaxDiff + pro1R1 + pro2R1;
      if (rMaxTot > params.rMaxLimit) {
        params.rMaxLimit = rMaxTot;
        params.rMaxRadius = pro1R1 + pro2R1;
      }
    }
  }
  if (forwardRxns.size() == 0) {
    params.rMaxLimit = 40.0;
  }
  std::cout << "Rmaxlimit: " << params.rMaxLimit << std::endl;
}
