#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

#include "json.hpp"

#include <cmath>

using json = nlohmann::json;

/* The restart file is JSON: one object with the top-level keys
 *
 *   parameters, molTemplates, reactions, molecules, complexes, observables,
 *   counterArrays, and implicitLipid (implicit-lipid models only).
 *
 * The layout follows the json-restarts branch of nerdss_development, so a file
 * written by either build reads on the other as far as the state they share.
 * The state this branch's restart carries beyond it is written under keys of
 * its own, which that build's reader never looks at:
 *
 *   parameters.membrane.compartmentSiteD, .compartmentSiteRho
 *   parameters.numerics
 *   molTemplates.numEachMolType, .totalNumOfStates
 *   reactions.back[].stateChangeIface
 *   reactions.transmission[].bindRadius, .reactantList, .productList
 *   molecules.list[].enforceCompartmentBC
 *   implicitLipid
 *
 * Nothing in the file is positional, so a field can be added, or written for
 * the models that need it alone, without changing what any other model's file
 * says; the .dat format of LEGACY_write_restart() below could not do that (see
 * expect_restart_key() in read_restart.cpp).  Doubles are written with enough
 * digits to read back exactly.  NaN is not a JSON number and is written as
 * null; read_restart() turns null back into NaN.
 */
namespace {

json maybe_null(double val)
{
    return std::isnan(val) ? json(nullptr) : json(val);
}

json serialize_vec3(const Vec3D& v)
{
    return { v.x, v.y, v.z };
}

json serialize_iface_ref(const RxnIface& iface)
{
    // requiresState is a char.  The .dat format writes it raw, a NUL byte for
    // "no state required"; here it is the character's integer value.
    return {
        { "molTypeIndex", iface.molTypeIndex },
        { "ifaceName", iface.ifaceName },
        { "absIfaceIndex", iface.absIfaceIndex },
        { "relIfaceIndex", iface.relIfaceIndex },
        { "requiresState", iface.requiresState },
        { "requiresInteraction", iface.requiresInteraction }
    };
}

json serialize_iface_list(const std::vector<RxnIface>& ifaceList)
{
    json j = json::array();
    for (const auto& iface : ifaceList)
        j.push_back(serialize_iface_ref(iface));
    return j;
}

json serialize_state_change_iface(const std::pair<RxnIface, RxnIface>& stateChangeIface)
{
    json j = json::array();
    j.push_back(serialize_iface_ref(stateChangeIface.first));
    j.push_back(serialize_iface_ref(stateChangeIface.second));
    return j;
}

json serialize_rate_list(const std::vector<RxnBase::RateState>& rateList)
{
    json jrates = json::array();
    for (const auto& rate : rateList) {
        json jr;
        jr["rate"] = maybe_null(rate.rate);
        jr["otherIfaceLists"] = json::array();
        for (const auto& list : rate.otherIfaceLists)
            jr["otherIfaceLists"].push_back(serialize_iface_list(list));
        jrates.push_back(jr);
    }
    return jrates;
}

json serialize_coupled_rxn(const RxnBase::CoupledRxn& coupled)
{
    return {
        { "absRxnIndex", coupled.absRxnIndex },
        { "relRxnIndex", coupled.relRxnIndex },
        { "rxnType", static_cast<std::underlying_type<ReactionType>::type>(coupled.rxnType) },
        { "label", coupled.label },
        { "probCoupled", coupled.probCoupled }
    };
}

// CreateDestructRxn::CreateDestructMol and TransmissionRxn::TransmissionMol.
template <typename MolWithIfaces>
json serialize_mol_with_ifaces(const MolWithIfaces& mol)
{
    json jm;
    jm["molTypeIndex"] = mol.molTypeIndex;
    jm["molName"] = mol.molName;
    jm["interfaceList"] = serialize_iface_list(mol.interfaceList);
    return jm;
}

json serialize_molecule_interface(const Molecule::Iface& iface)
{
    json ji;
    ji["index"] = iface.index;
    ji["relIndex"] = iface.relIndex;
    ji["molTypeIndex"] = iface.molTypeIndex;
    ji["stateIndex"] = iface.stateIndex;
    ji["stateIden"] = iface.stateIden;
    ji["isBound"] = iface.isBound;
    ji["coord"] = serialize_vec3(iface.coord);
    if (iface.isBound) {
        ji["interaction"] = {
            { "partnerIndex", iface.interaction.partnerIndex },
            { "partnerIfaceIndex", iface.interaction.partnerIfaceIndex },
            { "conjBackRxn", iface.interaction.conjBackRxn }
        };
    }
    return ji;
}

// The same names as the #NumericalSettings block of the .dat format.
json serialize_numerics(const NumericalSettings& numerics)
{
    return {
        { "integrationAbsError", numerics.integration.tableAbsoluteError },
        { "integrationRelError", numerics.integration.tableRelativeError },
        { "integrationFallbackError", numerics.integration.fallbackError },
        { "integrationTailCutoff", numerics.integration.tailCutoff },
        { "normalizationAbsError", numerics.integration.normalizationAbsoluteError },
        { "normalizationRelError", numerics.integration.normalizationRelativeError },
        { "tableRateAbsTolerance", numerics.tableLookup.reactionRate.absolute },
        { "tableRateRelTolerance", numerics.tableLookup.reactionRate.relative },
        { "tableDiffusionAbsTolerance", numerics.tableLookup.diffusionCoefficient.absolute },
        { "tableDiffusionRelTolerance", numerics.tableLookup.diffusionCoefficient.relative },
        { "explicitLipidFlatDiffusion", numerics.classification.explicitLipidFlatDiffusion },
        { "implicitLipidFlatDiffusion", numerics.classification.implicitLipidFlatDiffusion },
        { "associationSameAngleAbsTolerance", numerics.associationAngles.sameAngle.absolute },
        { "associationSameAngleRelTolerance", numerics.associationAngles.sameAngle.relative },
        { "associationRotationTolerance", numerics.associationAngles.rotationConvergenceTolerance },
        { "associationEndpointSignTolerance", numerics.associationAngles.endpointSignTolerance },
        { "vec3DCoordinatePrecision", numerics.vec3D.coordinateEqualityPrecision }
    };
}

} // namespace

void write_json_restart(long long int simItr, std::ofstream& restartFile, const Parameters& params, const SimulVolume& simulVolume,
    const std::vector<Molecule>& moleculeList, const std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, const std::vector<CreateDestructRxn>& createDestructRxns,
    const std::vector<TransmissionRxn>& transmissionRxns,
    const std::map<std::string, int>& observablesList, const Membrane& membraneObject, const copyCounters& counterArrays)
{
    // TRACE();
    json j;

    // parameters
    {
        json& jp = j["parameters"];
        jp["nItr"] = params.nItr;
        jp["simItr"] = simItr;
        jp["currSimTime"] = (simItr - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom;
        jp["numMolTypes"] = params.numMolTypes;
        jp["numTotalSpecies"] = params.numTotalSpecies;
        jp["numTotalComplex"] = params.numTotalComplex;
        jp["numTotalUnits"] = params.numTotalUnits;
        jp["numLipids"] = params.numLipids;
        jp["timeStep"] = params.timeStep;
        jp["max2DRxns"] = params.max2DRxns;
        jp["waterBox"] = { membraneObject.waterBox.x, membraneObject.waterBox.y, membraneObject.waterBox.z };
        jp["membrane"] = {
            { "implicitlipidIndex", membraneObject.implicitlipidIndex },
            { "nSites", membraneObject.nSites },
            { "nStates", membraneObject.nStates },
            { "No_free_lipids", membraneObject.No_free_lipids },
            { "No_protein", membraneObject.No_protein },
            { "totalSA", membraneObject.totalSA },
            // The diffusion constant and density of the compartment's surface
            // sites enter every transmission probability, and nothing else in
            // the file records them.  Zero for a model without a compartment.
            { "compartmentSiteD", membraneObject.droplet.D },
            { "compartmentSiteRho", membraneObject.droplet.rho }
        };
        jp["numberOfFreeLipidsEachState"] = membraneObject.numberOfFreeLipidsEachState;
        // isBox and isSphere keep the names of the two flags the .dat format
        // carries: the box flag is the waterBox provenance and the sphere flag
        // the geometry; see LEGACY_read_restart().
        jp["implicitLipidParams"] = {
            { "implicitLipid", membraneObject.implicitLipid },
            { "TwoD", membraneObject.TwoD },
            { "isBox", membraneObject.hasWaterBox() },
            { "isSphere", membraneObject.isSphere() },
            { "sphereR", membraneObject.sphereR },
            { "hasCompartment", membraneObject.hasCompartment },
            { "compartmentR", membraneObject.compartmentR }
        };
        jp["overlapSepLimit"] = params.overlapSepLimit;
        jp["rMaxLimit"] = params.rMaxLimit;
        jp["timeWrite"] = params.timeWrite;
        jp["trajWrite"] = params.trajWrite;
        jp["restartWrite"] = params.restartWrite;
        jp["pdbWrite"] = params.pdbWrite;
        jp["assocDissocWrite"] = params.assocDissocWrite;
        jp["checkPoint"] = params.checkPoint;
        jp["scaleMaxDisplace"] = params.scaleMaxDisplace;
        jp["transitionWrite"] = params.transitionWrite;
        jp["clusterOverlapCheck"] = params.clusterOverlapCheck;
        jp["rngwrite"] = params.rngwrite;
        jp["bondedComplexWrite"] = params.bondedComplexWrite;
        jp["lastUpdateTransition"] = Parameters::lastUpdateTransition;
        jp["numerics"] = serialize_numerics(params.numerics);
    }

    // molecule templates
    {
        json templates = json::array();
        for (const auto& oneTemp : molTemplateList) {
            json jt;
            jt["molTypeIndex"] = oneTemp.molTypeIndex;
            jt["molName"] = oneTemp.molName;
            jt["copies"] = oneTemp.copies;
            jt["mass"] = oneTemp.mass;
            jt["radius"] = oneTemp.radius;
            jt["isLipid"] = oneTemp.isLipid;
            jt["isImplicitLipid"] = oneTemp.isImplicitLipid;
            jt["isRod"] = oneTemp.isRod;
            jt["isPoint"] = oneTemp.isPoint;
            jt["checkOverlap"] = oneTemp.checkOverlap;
            jt["countTransition"] = oneTemp.countTransition;
            jt["transitionMatrixSize"] = oneTemp.transitionMatrixSize;
            jt["outsideCompartment"] = oneTemp.outsideCompartment;
            jt["insideCompartment"] = oneTemp.insideCompartment;
            jt["crossesCompartment"] = oneTemp.crossesCompartment;
            jt["transmissionRxnIndex"] = oneTemp.transmissionRxnIndex;
            jt["comCoord"] = serialize_vec3(oneTemp.comCoord);
            jt["D"] = serialize_vec3(oneTemp.D);
            jt["Dr"] = serialize_vec3(oneTemp.Dr);
            jt["rxnPartners"] = oneTemp.rxnPartners;

            jt["bondList"] = json::array();
            for (const auto& bond : oneTemp.bondList)
                jt["bondList"].push_back({ bond[0], bond[1] });

            jt["interfaces"] = json::array();
            for (const auto& oneIface : oneTemp.interfaceList) {
                json ji;
                ji["index"] = oneIface.index;
                ji["name"] = oneIface.name;
                ji["coord"] = serialize_vec3(oneIface.iCoord);
                ji["states"] = json::array();
                for (const auto& oneState : oneIface.stateList) {
                    json js;
                    js["index"] = oneState.index;
                    js["iden"] = oneState.iden;
                    js["rxnPartners"] = oneState.rxnPartners;
                    js["myForwardRxns"] = oneState.myForwardRxns;
                    js["myCreateDestructRxns"] = oneState.myCreateDestructRxns;
                    js["stateChangeRxns"] = json::array();
                    for (const auto& elem : oneState.stateChangeRxns)
                        js["stateChangeRxns"].push_back({ elem.first, elem.second });
                    ji["states"].push_back(js);
                }
                jt["interfaces"].push_back(ji);
            }

            jt["ifacesWithStates"] = oneTemp.ifacesWithStates;
            jt["monomerList"] = oneTemp.monomerList;

            if (oneTemp.countTransition) {
                jt["lifeTime"] = json::array();
                for (int i = 0; i < oneTemp.transitionMatrixSize; ++i)
                    jt["lifeTime"].push_back(oneTemp.lifeTime[i]);
                jt["transitionMatrix"] = json::array();
                for (int i = 0; i < oneTemp.transitionMatrixSize; ++i)
                    jt["transitionMatrix"].push_back(oneTemp.transitionMatrix[i]);
            }

            templates.push_back(jt);
        }
        json& jtemplates = j["molTemplates"];
        jtemplates["templates"] = templates;
        jtemplates["absToRelIface"] = MolTemplate::absToRelIface;
        // The count of each type is kept by the reactions that create and
        // destroy molecules, so it cannot be recounted from the molecule list
        // once any molecule has been destroyed: a destroyed molecule keeps its
        // slot, and its type, in the list.
        jtemplates["numEachMolType"] = MolTemplate::numEachMolType;
        // The states an add file's templates bring are numbered from here
        // (parse_molFile.cpp), so a restart with an add file needs it.
        jtemplates["totalNumOfStates"] = Interface::State::totalNumOfStates;
    }

    // reactions
    {
        json& jr = j["reactions"];
        jr["numberOfRxns"] = RxnBase::numberOfRxns;
        jr["numForward"] = forwardRxns.size();
        jr["numBack"] = backRxns.size();
        jr["numCreateDestruct"] = createDestructRxns.size();
        jr["numTransmission"] = transmissionRxns.size();
        jr["totRxnSpecies"] = RxnBase::totRxnSpecies;

        jr["forward"] = json::array();
        for (const auto& oneRxn : forwardRxns) {
            json jx;
            jx["absRxnIndex"] = oneRxn.absRxnIndex;
            jx["relRxnIndex"] = oneRxn.relRxnIndex;
            jx["rxnLabel"] = oneRxn.rxnLabel;
            jx["rxnType"] = static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType);
            jx["isSymmetric"] = oneRxn.isSymmetric;
            jx["isOnMem"] = oneRxn.isOnMem;
            jx["hasStateChange"] = oneRxn.hasStateChange;
            jx["isObserved"] = oneRxn.isObserved;
            jx["observeLabel"] = oneRxn.observeLabel;
            jx["productName"] = oneRxn.productName;
            jx["isReversible"] = oneRxn.isReversible;
            jx["conjBackRxnIndex"] = oneRxn.conjBackRxnIndex;
            jx["irrevRingClosure"] = oneRxn.irrevRingClosure;
            jx["bindRadSameCom"] = oneRxn.bindRadSameCom;
            jx["loopCoopFactor"] = oneRxn.loopCoopFactor;
            jx["length3Dto2D"] = oneRxn.length3Dto2D;
            jx["area3Dto1D"] = oneRxn.area3Dto1D;
            jx["bindRadius"] = oneRxn.bindRadius;
            jx["assocAngles"] = {
                { "theta1", maybe_null(oneRxn.assocAngles.theta1) },
                { "theta2", maybe_null(oneRxn.assocAngles.theta2) },
                { "phi1", maybe_null(oneRxn.assocAngles.phi1) },
                { "phi2", maybe_null(oneRxn.assocAngles.phi2) },
                { "omega", maybe_null(oneRxn.assocAngles.omega) }
            };
            jx["norm1"] = serialize_vec3(oneRxn.norm1);
            jx["norm2"] = serialize_vec3(oneRxn.norm2);
            jx["excludeVolumeBound"] = oneRxn.excludeVolumeBound;
            jx["isCoupled"] = oneRxn.isCoupled;
            if (oneRxn.isCoupled)
                jx["coupledRxn"] = serialize_coupled_rxn(oneRxn.coupledRxn);
            jx["intReactantList"] = oneRxn.intReactantList;
            jx["intProductList"] = oneRxn.intProductList;
            jx["reactantList"] = serialize_iface_list(oneRxn.reactantListNew);
            jx["productList"] = serialize_iface_list(oneRxn.productListNew);
            jx["rateList"] = serialize_rate_list(oneRxn.rateList);
            // Not in the .dat format.
            jx["stateChangeIface"] = serialize_state_change_iface(oneRxn.stateChangeIface);
            jr["forward"].push_back(jx);
        }

        jr["back"] = json::array();
        for (const auto& oneRxn : backRxns) {
            json jx;
            jx["absRxnIndex"] = oneRxn.absRxnIndex;
            jx["relRxnIndex"] = oneRxn.relRxnIndex;
            jx["rxnType"] = static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType);
            jx["isSymmetric"] = oneRxn.isSymmetric;
            jx["isOnMem"] = oneRxn.isOnMem;
            jx["hasStateChange"] = oneRxn.hasStateChange;
            jx["isObserved"] = oneRxn.isObserved;
            jx["observeLabel"] = oneRxn.observeLabel;
            jx["conjForwardRxnIndex"] = oneRxn.conjForwardRxnIndex;
            jx["isCoupled"] = oneRxn.isCoupled;
            if (oneRxn.isCoupled)
                jx["coupledRxn"] = serialize_coupled_rxn(oneRxn.coupledRxn);
            jx["intReactantList"] = oneRxn.intReactantList;
            jx["intProductList"] = oneRxn.intProductList;
            jx["reactantList"] = serialize_iface_list(oneRxn.reactantListNew);
            jx["productList"] = serialize_iface_list(oneRxn.productListNew);
            jx["rateList"] = serialize_rate_list(oneRxn.rateList);
            // A back reaction carries its forward reaction's state change with
            // the two interfaces swapped (BackRxn::BackRxn()).  nerdss_development
            // writes this for forward reactions only and then requires it of
            // every back reaction that has a state change, so its own files
            // with a reversible state change do not read back; written for
            // both here.
            jx["stateChangeIface"] = serialize_state_change_iface(oneRxn.stateChangeIface);
            jr["back"].push_back(jx);
        }

        jr["createDestruct"] = json::array();
        for (const auto& oneRxn : createDestructRxns) {
            json jx;
            jx["absRxnIndex"] = oneRxn.absRxnIndex;
            jx["relRxnIndex"] = oneRxn.relRxnIndex;
            jx["rxnType"] = static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType);
            jx["isOnMem"] = oneRxn.isOnMem;
            jx["isObserved"] = oneRxn.isObserved;
            jx["observeLabel"] = oneRxn.observeLabel;
            jx["creationRadius"] = oneRxn.creationRadius;
            jx["intReactantList"] = oneRxn.intReactantList;
            jx["intProductList"] = oneRxn.intProductList;
            jx["reactantMolList"] = json::array();
            for (const auto& mol : oneRxn.reactantMolList)
                jx["reactantMolList"].push_back(serialize_mol_with_ifaces(mol));
            jx["productMolList"] = json::array();
            for (const auto& mol : oneRxn.productMolList)
                jx["productMolList"].push_back(serialize_mol_with_ifaces(mol));
            jx["rateList"] = serialize_rate_list(oneRxn.rateList);
            jr["createDestruct"].push_back(jx);
        }

        jr["transmission"] = json::array();
        for (const auto& oneRxn : transmissionRxns) {
            json jx;
            jx["absRxnIndex"] = oneRxn.absRxnIndex;
            jx["relRxnIndex"] = oneRxn.relRxnIndex;
            jx["rxnType"] = static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType);
            jx["isOnMem"] = oneRxn.isOnMem;
            jx["isObserved"] = oneRxn.isObserved;
            jx["observeLabel"] = oneRxn.observeLabel;
            jx["intReactantList"] = oneRxn.intReactantList;
            jx["intProductList"] = oneRxn.intProductList;
            jx["reactantMolList"] = json::array();
            for (const auto& mol : oneRxn.reactantMolList)
                jx["reactantMolList"].push_back(serialize_mol_with_ifaces(mol));
            jx["productMolList"] = json::array();
            for (const auto& mol : oneRxn.productMolList)
                jx["productMolList"].push_back(serialize_mol_with_ifaces(mol));
            jx["rateList"] = serialize_rate_list(oneRxn.rateList);
            // The compartment code takes the crossing interface from
            // reactantListNew and the capture radius from bindRadius, and
            // neither can be rebuilt from the lists above; see
            // LEGACY_write_restart().
            jx["bindRadius"] = oneRxn.bindRadius;
            jx["reactantList"] = serialize_iface_list(oneRxn.reactantListNew);
            jx["productList"] = serialize_iface_list(oneRxn.productListNew);
            jr["transmission"].push_back(jx);
        }
    }

    // molecules
    {
        json& jmols = j["molecules"];
        jmols["count"] = moleculeList.size();
        jmols["numberOfMolecules"] = Molecule::numberOfMolecules;
        jmols["list"] = json::array();
        for (const auto& oneMol : moleculeList) {
            json jm;
            jm["index"] = oneMol.index;
            jm["isEmpty"] = oneMol.isEmpty;
            jm["myComIndex"] = oneMol.myComIndex;
            jm["molTypeIndex"] = oneMol.molTypeIndex;
            jm["mySubVolIndex"] = oneMol.mySubVolIndex;
            jm["mass"] = oneMol.mass;
            jm["isLipid"] = oneMol.isLipid;
            jm["isImplicitLipid"] = oneMol.isImplicitLipid;
            jm["linksToSurface"] = oneMol.linksToSurface;
            jm["isPromoter"] = oneMol.isPromoter;
            // Set the first time a molecule is kept from crossing the
            // compartment and never cleared, so it is state, not a per-step
            // flag.
            jm["enforceCompartmentBC"] = oneMol.enforceCompartmentBC;
            jm["comCoord"] = serialize_vec3(oneMol.comCoord);
            jm["freelist"] = oneMol.freelist;
            jm["bndlist"] = oneMol.bndlist;
            jm["bndpartner"] = oneMol.bndpartner;
            jm["interfaceList"] = json::array();
            for (const auto& iface : oneMol.interfaceList)
                jm["interfaceList"].push_back(serialize_molecule_interface(iface));

            // Reweighting lists.  One std::vector<ReweightEntry> in memory,
            // the six parallel arrays the .dat format carries in the file,
            // under the names the entries' fields used to have.
            json prevlist = json::array();
            json prevmyface = json::array();
            json prevpface = json::array();
            json prevnorm = json::array();
            json ps_prev = json::array();
            json prevsep = json::array();
            for (const auto& oneEntry : oneMol.prevReweight) {
                prevlist.push_back(oneEntry.partner);
                prevmyface.push_back(oneEntry.myFace);
                prevpface.push_back(oneEntry.partnerFace);
                prevnorm.push_back(maybe_null(oneEntry.norm));
                ps_prev.push_back(maybe_null(oneEntry.survProb));
                prevsep.push_back(maybe_null(oneEntry.sep));
            }
            jm["prevlist"] = prevlist;
            jm["prevmyface"] = prevmyface;
            jm["prevpface"] = prevpface;
            jm["prevnorm"] = prevnorm;
            jm["ps_prev"] = ps_prev;
            jm["prevsep"] = prevsep;

            jmols["list"].push_back(jm);
        }
        jmols["emptyMolList"] = Molecule::emptyMolList;
    }

    // complexes
    {
        json& jcoms = j["complexes"];
        jcoms["count"] = complexList.size();
        jcoms["numberOfComplexes"] = Complex::numberOfComplexes;
        jcoms["list"] = json::array();
        for (const auto& oneCom : complexList) {
            json jc;
            jc["index"] = oneCom.index;
            jc["isEmpty"] = oneCom.isEmpty;
            jc["radius"] = oneCom.radius;
            jc["mass"] = oneCom.mass;
            jc["linksToSurface"] = oneCom.linksToSurface;
            jc["iLipidIndex"] = oneCom.iLipidIndex;
            jc["OnSurface"] = oneCom.OnSurface;
            jc["onFiber"] = oneCom.onFiber;
            jc["comCoord"] = serialize_vec3(oneCom.comCoord);
            jc["D"] = serialize_vec3(oneCom.D);
            jc["Dr"] = serialize_vec3(oneCom.Dr);
            jc["memberList"] = oneCom.memberList;
            jc["numEachMol"] = oneCom.numEachMol;
            jc["lastNumberUpdateItrEachMol"] = oneCom.lastNumberUpdateItrEachMol;
            jcoms["list"].push_back(jc);
        }
        jcoms["emptyComList"] = Complex::emptyComList;
    }

    // observables
    j["observables"] = json::array();
    for (const auto& observable : observablesList)
        j["observables"].push_back({ { "name", observable.first }, { "value", observable.second } });

    // counter arrays
    {
        json& jc = j["counterArrays"];
        jc["nLoops"] = counterArrays.nLoops;
        jc["nCancelOverlapPartner"] = counterArrays.nCancelOverlapPartner;
        jc["nCancelOverlapSystem"] = counterArrays.nCancelOverlapSystem;
        jc["nCancelDisplace2D"] = counterArrays.nCancelDisplace2D;
        jc["nCancelDisplace3D"] = counterArrays.nCancelDisplace3D;
        jc["nCancelDisplace3Dto2D"] = counterArrays.nCancelDisplace3Dto2D;
        jc["nCancelSpanBox"] = counterArrays.nCancelSpanBox;
        jc["nAssocSuccess"] = counterArrays.nAssocSuccess;
        jc["eventArraySize"] = counterArrays.eventArraySize;
        jc["events3D"] = counterArrays.events3D;
        jc["events3Dto2D"] = counterArrays.events3Dto2D;
        jc["events2D"] = counterArrays.events2D;
        jc["bindPairList"] = counterArrays.bindPairList;
    }

    // The implicit lipid's 2D binding table, and the protein counts its entries
    // are built from.  Each entry is computed the first time it is needed,
    // from the free-lipid count at that step, and kept for the rest of the
    // run; a restart that rebuilt it would use the count at the restart
    // instead and bind with different probabilities from then on.  The protein
    // counts are those at step 0, which a restart cannot recount from its own
    // molecules once any have been created, destroyed or changed state.
    if (membraneObject.implicitLipid) {
        json& jil = j["implicitLipid"];
        jil["numberOfProteinEachState"] = membraneObject.numberOfProteinEachState;
        jil["binding2DTable"] = json::array();
        for (std::size_t entry { 0 }; entry < membraneObject.IL2DbindingVec.size(); ++entry) {
            json row = {
                { "ka", membraneObject.ILTableIDs[3 * entry] },
                { "Dtot", membraneObject.ILTableIDs[3 * entry + 1] },
                { "kb", membraneObject.ILTableIDs[3 * entry + 2] },
                { "probability", membraneObject.IL2DbindingVec[entry] }
            };
            jil["binding2DTable"].push_back(row);
        }
    }

    restartFile << j.dump() << std::endl;
}

void write_restart(long long int simItr, std::ofstream& restartFile, const Parameters& params, const SimulVolume& simulVolume,
    const std::vector<Molecule>& moleculeList, const std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, const std::vector<CreateDestructRxn>& createDestructRxns,
    const std::vector<TransmissionRxn>& transmissionRxns,
    const std::map<std::string, int>& observablesList, const Membrane& membraneObject, const copyCounters& counterArrays)
{
    // JSON unless the run asked for the .dat format of earlier builds.
    if (params.legacyRestartFormat) {
        LEGACY_write_restart(simItr, restartFile, params, simulVolume, moleculeList, complexList, molTemplateList,
            forwardRxns, backRxns, createDestructRxns, transmissionRxns, observablesList, membraneObject,
            counterArrays);
    } else {
        write_json_restart(simItr, restartFile, params, simulVolume, moleculeList, complexList, molTemplateList,
            forwardRxns, backRxns, createDestructRxns, transmissionRxns, observablesList, membraneObject,
            counterArrays);
    }
}

/* The .dat format of builds before the JSON format.  A run writes it when
 * legacyRestartFormat is set (a restart from a .dat file keeps it), the
 * regression harness converts between the two formats with it, and it is
 * the reference for what LEGACY_read_restart() expects.  Its layout must not
 * change: files it writes are read by the builds it comes from.
 */

void LEGACY_write_restart(long long int simItr, std::ofstream& restartFile, const Parameters& params, const SimulVolume& simulVolume,
    const std::vector<Molecule>& moleculeList, const std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, const std::vector<CreateDestructRxn>& createDestructRxns,
    const std::vector<TransmissionRxn>& transmissionRxns,
    const std::map<std::string, int>& observablesList, const Membrane& membraneObject, const copyCounters& counterArrays)
{
    // TRACE();
    // Write parameters
    restartFile.precision(20);
    {
        restartFile << "#Parameters--update these for restart \n";
        restartFile << "numItr = " << params.nItr << '\n';
        restartFile << "currItr = " << simItr << '\n';
        restartFile << "currSimTime (s) = " << std::scientific << (simItr - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom << '\n';
        restartFile << "numMolTypes = " << params.numMolTypes << '\n';
        restartFile << "numTotalSpecies = " << params.numTotalSpecies << '\n';
        restartFile << "numComplexs = " << params.numTotalComplex << '\n';
        restartFile << "numTotalUnits = " << params.numTotalUnits << '\n';
        restartFile << "numLipids = " << params.numLipids << '\n';
        restartFile << "timestep = " << std::fixed << params.timeStep << '\n';
        restartFile << "max2DRxns = " << params.max2DRxns << '\n';
        restartFile << "simulDimensions = " << membraneObject.waterBox.x << ' ' << membraneObject.waterBox.y << ' ' << membraneObject.waterBox.z
                    << '\n';
        restartFile << "membrane = " << membraneObject.implicitlipidIndex << ' ' << membraneObject.nSites << ' ' << membraneObject.nStates << ' ' << membraneObject.No_free_lipids << ' ' << membraneObject.No_protein << ' ' << membraneObject.totalSA << '\n';
        restartFile << "implicitLipidStates = ";
        for (int i = 0; i < membraneObject.nStates; i++) {
            restartFile << membraneObject.numberOfFreeLipidsEachState[i];
            if (i != membraneObject.nStates - 1) {
                restartFile << ' ';
            }
        }
        restartFile << '\n';
        restartFile << "implicitLipidsParams = " << membraneObject.implicitLipid << ' ' << membraneObject.TwoD << ' ' << membraneObject.hasWaterBox() << ' ' << membraneObject.isSphere() << ' ' << membraneObject.sphereR <<' ' << membraneObject.hasCompartment << ' ' << membraneObject.compartmentR << '\n';
        // The diffusion constant and density of the compartment's surface sites
        // enter every transmission probability, and nothing else in the file
        // records them, so a restart used to run with both at zero.  Written only
        // when there is a compartment, so no other model's file changes.
        // Scientific, unlike the fixed-point lines around it: a site density in
        // nm^-2 can be small enough that twenty decimal places drop digits.
        if (membraneObject.hasCompartment) {
            const std::ios_base::fmtflags savedFlags { restartFile.flags() };
            restartFile << std::scientific;
            restartFile << "compartmentSiteD = " << membraneObject.droplet.D << '\n';
            restartFile << "compartmentSiteRho = " << membraneObject.droplet.rho << '\n';
            restartFile.flags(savedFlags);
        }
        restartFile << "ifaceOverlapSepLimit = " << params.overlapSepLimit << '\n';
        restartFile << "rMaxLimit = " << params.rMaxLimit << '\n';
        restartFile << "timeWrite = " << params.timeWrite << '\n';
        restartFile << "trajWrite = " << params.trajWrite << '\n';
        restartFile << "restartWrite = " << params.restartWrite << '\n';
        restartFile << "pdbWrite = " << params.pdbWrite << '\n';
        restartFile << "accocDissocWrite = " << params.assocDissocWrite << '\n';
        restartFile << "checkPoint = " << params.checkPoint << '\n';
        restartFile << "scaleMaxDisplace = " << params.scaleMaxDisplace << '\n';
        restartFile << "transitionWrite = " << params.transitionWrite << '\n';
        restartFile << "clusterOverlapCheck = " << params.clusterOverlapCheck << '\n';
        restartFile << "RNGwrite = " << params.rngwrite << '\n';
        restartFile << "bondedComplexWrite = " << params.bondedComplexWrite << '\n';

        restartFile << Parameters::lastUpdateTransition.size();
        for (auto& index : Parameters::lastUpdateTransition)
            restartFile << ' ' << index;
        restartFile << '\n';

        restartFile << "#NumericalSettings version = 3\n" << std::scientific;
        restartFile << "integrationAbsError = " << params.numerics.integration.tableAbsoluteError << '\n';
        restartFile << "integrationRelError = " << params.numerics.integration.tableRelativeError << '\n';
        restartFile << "integrationFallbackError = " << params.numerics.integration.fallbackError << '\n';
        restartFile << "integrationTailCutoff = " << params.numerics.integration.tailCutoff << '\n';
        restartFile << "normalizationAbsError = "
                    << params.numerics.integration.normalizationAbsoluteError << '\n';
        restartFile << "normalizationRelError = "
                    << params.numerics.integration.normalizationRelativeError << '\n';
        restartFile << "tableRateAbsTolerance = " << params.numerics.tableLookup.reactionRate.absolute << '\n';
        restartFile << "tableRateRelTolerance = " << params.numerics.tableLookup.reactionRate.relative << '\n';
        restartFile << "tableDiffusionAbsTolerance = "
                    << params.numerics.tableLookup.diffusionCoefficient.absolute << '\n';
        restartFile << "tableDiffusionRelTolerance = "
                    << params.numerics.tableLookup.diffusionCoefficient.relative << '\n';
        restartFile << "explicitLipidFlatDiffusion = "
                    << params.numerics.classification.explicitLipidFlatDiffusion << '\n';
        restartFile << "implicitLipidFlatDiffusion = "
                    << params.numerics.classification.implicitLipidFlatDiffusion << '\n';
        restartFile << "associationSameAngleAbsTolerance = "
                    << params.numerics.associationAngles.sameAngle.absolute << '\n';
        restartFile << "associationSameAngleRelTolerance = "
                    << params.numerics.associationAngles.sameAngle.relative << '\n';
        restartFile << "associationRotationTolerance = "
                    << params.numerics.associationAngles.rotationConvergenceTolerance << '\n';
        restartFile << "associationEndpointSignTolerance = "
                    << params.numerics.associationAngles.endpointSignTolerance << '\n';
        restartFile << "vec3DCoordinatePrecision = "
                    << params.numerics.vec3D.coordinateEqualityPrecision << '\n';
        restartFile << "#EndNumericalSettings\n";
    }
    /*
    // write Simulation Volume
    {
        restartFile << simulVolume.numSubCells.x << ' ' << simulVolume.numSubCells.y << ' ' << simulVolume.numSubCells.z
                    << ' ' << simulVolume.numSubCells.tot << '\n';
        restartFile << simulVolume.subCellSize.x << ' ' << simulVolume.subCellSize.y << ' ' << simulVolume.subCellSize.z
                    << '\n';
        for (auto& subCell : simulVolume.subCellList) {
            restartFile << subCell.absIndex << ' ' << subCell.xIndex << ' ' << subCell.yIndex << ' ' << subCell.zIndex
                        << '\n';
            restartFile << subCell.neighborList.size();
            for (const auto& neighbor : subCell.neighborList)
                restartFile << ' ' << neighbor;
            restartFile << '\n';

            restartFile << subCell.memberMolList.size();
            for (const auto& memMol : subCell.memberMolList)
                restartFile << ' ' << memMol;
            restartFile << '\n';
        }
    }
    */
    // write MolTemplates
    {
        restartFile << "#MolTemplates \n";
        restartFile << MolTemplate::numMolTypes;
        for (auto& num : MolTemplate::numEachMolType)
            restartFile << ' ' << num;
        restartFile << '\n';
        restartFile << MolTemplate::absToRelIface.size();
        for (auto& iface : MolTemplate::absToRelIface)
            restartFile << ' ' << iface;
        restartFile << '\n';
        restartFile << Interface::State::totalNumOfStates << '\n';

        for (const auto& oneTemp : molTemplateList) {
            restartFile << oneTemp.molTypeIndex << ' ' << oneTemp.molName << '\n';
            restartFile << oneTemp.copies << ' ' << oneTemp.mass << ' ' << oneTemp.radius << '\n';
            restartFile << oneTemp.isLipid << ' ' << oneTemp.isImplicitLipid << ' ' << oneTemp.isRod << ' ' << oneTemp.isPoint << ' '
                        << oneTemp.checkOverlap << ' ' << oneTemp.countTransition << ' ' << oneTemp.transitionMatrixSize << ' ' 
                        << oneTemp.outsideCompartment << ' ' << oneTemp.insideCompartment << ' ' << oneTemp.crossesCompartment << ' ' << oneTemp.transmissionRxnIndex << '\n';
            restartFile << oneTemp.comCoord.x << ' ' << oneTemp.comCoord.y << ' ' << oneTemp.comCoord.z
                        << '\n';
            restartFile << oneTemp.D.x << ' ' << oneTemp.D.y << ' ' << oneTemp.D.z << '\n';
            restartFile << oneTemp.Dr.x << ' ' << oneTemp.Dr.y << ' ' << oneTemp.Dr.z << '\n';

            // reaction partners
            restartFile << oneTemp.rxnPartners.size();
            for (const auto& partner : oneTemp.rxnPartners)
                restartFile << ' ' << partner;
            restartFile << '\n';

            // optional bonds
            restartFile << oneTemp.bondList.size() << '\n';
            for (const auto& bond : oneTemp.bondList)
                restartFile << bond[0] << ' ' << bond[1] << '\n';

            // write interfaces
            restartFile << oneTemp.interfaceList.size() << '\n';
            for (const auto& oneIface : oneTemp.interfaceList) {
                restartFile << oneIface.index << ' ' << oneIface.name << '\n';
                restartFile << std::fixed << oneIface.iCoord.x << ' ' << oneIface.iCoord.y << ' ' << oneIface.iCoord.z
                            << '\n';
                restartFile << oneIface.stateList.size() << '\n';
                for (auto& oneState : oneIface.stateList) {
                    restartFile << oneState.index << ' ' << oneState.iden << '\n';

                    // partner list
                    restartFile << oneState.rxnPartners.size();
                    for (auto elem : oneState.rxnPartners)
                        restartFile << ' ' << elem;
                    restartFile << '\n';

                    // reaction lists
                    restartFile << oneState.myForwardRxns.size();
                    for (auto elem : oneState.myForwardRxns)
                        restartFile << ' ' << elem;
                    restartFile << '\n';
                    restartFile << oneState.myCreateDestructRxns.size();
                    for (auto elem : oneState.myCreateDestructRxns)
                        restartFile << ' ' << elem;
                    restartFile << '\n';
                    restartFile << oneState.stateChangeRxns.size();
                    for (auto elem : oneState.stateChangeRxns)
                        restartFile << ' ' << elem.first << ' ' << elem.second;
                    restartFile << '\n';
                }
            }

            //write ifacesWithStates
            restartFile << oneTemp.ifacesWithStates.size();
            for (auto elem : oneTemp.ifacesWithStates) {
                restartFile << ' ' << elem;
            }
            restartFile << '\n';

            //write monomerList
            restartFile << oneTemp.monomerList.size();
            for (auto elem : oneTemp.monomerList) {
                restartFile << ' ' << elem;
            }
            restartFile << '\n';

            //write lifetime
            if(oneTemp.countTransition == true){
                for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne){
                    restartFile << oneTemp.lifeTime[indexOne].size();
                    for (auto elem : oneTemp.lifeTime[indexOne]) {
                        restartFile << ' ' << elem;
                    }
                    restartFile << '\n';
                }
                restartFile << '\n';
            }

            //write transition matrix
            if(oneTemp.countTransition == true){
                for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne){
                    for (int indexTwo = 0; indexTwo < oneTemp.transitionMatrixSize; ++indexTwo){
                        restartFile << ' ' << oneTemp.transitionMatrix[indexOne][indexTwo];
                    }
                    restartFile << '\n';
                }
                restartFile << '\n';
            }
        }
    }

    // write Reactions
    {
        restartFile << "#Reactions \n";
        restartFile << RxnBase::numberOfRxns << ' ' << forwardRxns.size() << ' ' << backRxns.size() << ' '
                    << createDestructRxns.size() << ' ' << transmissionRxns.size() << ' ' << RxnBase::totRxnSpecies << '\n';

        // forward reactions
        for (const auto& oneRxn : forwardRxns) {
            restartFile << oneRxn.absRxnIndex << ' ' << oneRxn.relRxnIndex << ' ' << oneRxn.rxnLabel << '\n';
            restartFile << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType) << ' '
                        << oneRxn.isSymmetric << ' ' << oneRxn.isOnMem << ' ' << oneRxn.hasStateChange << '\n';
            restartFile << oneRxn.isObserved << ' ' << oneRxn.observeLabel << ' ' << oneRxn.productName << '\n';
            restartFile << oneRxn.isReversible << ' ' << oneRxn.conjBackRxnIndex << ' ' << oneRxn.irrevRingClosure << ' ' << oneRxn.bindRadSameCom << ' '
                        << oneRxn.loopCoopFactor << '\n'
                        << oneRxn.length3Dto2D << '\n'
                        << oneRxn.area3Dto1D << '\n';
            restartFile << std::setprecision(20) << oneRxn.bindRadius << ' ' << oneRxn.assocAngles.theta1 << ' ' << oneRxn.assocAngles.theta2
                        << ' ' << oneRxn.assocAngles.phi1 << ' ' << oneRxn.assocAngles.phi2 << ' '
                        << oneRxn.assocAngles.omega << '\n';
            restartFile << oneRxn.norm1.x << ' ' << oneRxn.norm1.y << ' ' << oneRxn.norm1.z << '\n';
            restartFile << oneRxn.norm2.x << ' ' << oneRxn.norm2.y << ' ' << oneRxn.norm2.z << '\n';
            restartFile << oneRxn.excludeVolumeBound << '\n';
            restartFile << oneRxn.isCoupled;
            if (oneRxn.isCoupled)
                restartFile << ' ' << oneRxn.coupledRxn.absRxnIndex << ' ' << oneRxn.coupledRxn.relRxnIndex << ' ' << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.coupledRxn.rxnType) << ' ' << oneRxn.coupledRxn.label << ' ' << std::setprecision(20) << oneRxn.coupledRxn.probCoupled;
            restartFile << '\n';

            // integer reactants
            restartFile << oneRxn.intReactantList.size();
            for (const auto& oneReact : oneRxn.intReactantList)
                restartFile << ' ' << oneReact;
            restartFile << '\n';

            // integer products
            restartFile << oneRxn.intProductList.size();
            for (const auto& oneProd : oneRxn.intProductList)
                restartFile << ' ' << oneProd;
            restartFile << '\n';

            // reactant list
            restartFile << oneRxn.reactantListNew.size() << '\n';
            for (const auto& oneReact : oneRxn.reactantListNew) {
                restartFile << oneReact.molTypeIndex << '\n';
                restartFile << oneReact.ifaceName << ' ' << oneReact.absIfaceIndex << ' ' << oneReact.relIfaceIndex
                            << '\n';
                restartFile << oneReact.requiresState << ' ' << oneReact.requiresInteraction << '\n';
            }

            // product list
            restartFile << oneRxn.productListNew.size() << '\n';
            for (const auto& oneProd : oneRxn.productListNew) {
                restartFile << oneProd.molTypeIndex << '\n';
                restartFile << oneProd.ifaceName << ' ' << oneProd.absIfaceIndex << ' ' << oneProd.relIfaceIndex
                            << '\n';
                restartFile << oneProd.requiresState << ' ' << oneProd.requiresInteraction << '\n';
            }

            // rate list
            restartFile << oneRxn.rateList.size() << '\n';
            for (auto& oneRate : oneRxn.rateList) {
                restartFile << std::setprecision(20) << oneRate.rate << '\n';
                restartFile << oneRate.otherIfaceLists.size() << '\n';
                for (const auto& otherIfaceList : oneRate.otherIfaceLists) {
                    restartFile << otherIfaceList.size() << '\n';
                    for (const auto& anccIface : otherIfaceList) {
                        restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                                    << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                                    << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                    }
                }
            }
        }

        // back reactions
        for (const auto& oneRxn : backRxns) {
            restartFile << oneRxn.absRxnIndex << ' ' << oneRxn.relRxnIndex << '\n';
            restartFile << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType) << ' '
                        << oneRxn.isSymmetric << ' ' << oneRxn.isOnMem << ' ' << oneRxn.hasStateChange << '\n';
            restartFile << oneRxn.isObserved << ' ' << oneRxn.observeLabel << '\n';
            restartFile << oneRxn.conjForwardRxnIndex << '\n';
            restartFile << oneRxn.isCoupled;
            if (oneRxn.isCoupled)
                restartFile << ' ' << oneRxn.coupledRxn.absRxnIndex << ' ' << oneRxn.coupledRxn.relRxnIndex << ' ' << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.coupledRxn.rxnType) << ' ' << oneRxn.coupledRxn.label << ' ' << std::setprecision(20) << oneRxn.coupledRxn.probCoupled;
            restartFile << '\n';

            // integer reactants
            restartFile << oneRxn.intReactantList.size();
            for (const auto& oneReact : oneRxn.intReactantList)
                restartFile << ' ' << oneReact;
            restartFile << '\n';

            // integer products
            restartFile << oneRxn.intProductList.size();
            for (const auto& oneProd : oneRxn.intProductList)
                restartFile << ' ' << oneProd;
            restartFile << '\n';

            // reactant list
            restartFile << oneRxn.reactantListNew.size() << '\n';
            for (const auto& oneReact : oneRxn.reactantListNew) {
                restartFile << oneReact.molTypeIndex << '\n';
                restartFile << oneReact.ifaceName << ' ' << oneReact.absIfaceIndex << ' ' << oneReact.relIfaceIndex
                            << '\n';
                restartFile << oneReact.requiresState << ' ' << oneReact.requiresInteraction << '\n';
            }

            // product list
            restartFile << oneRxn.productListNew.size() << '\n';
            for (const auto& oneProd : oneRxn.productListNew) {
                restartFile << oneProd.molTypeIndex << '\n';
                restartFile << oneProd.ifaceName << ' ' << oneProd.absIfaceIndex << ' ' << oneProd.relIfaceIndex
                            << '\n';
                restartFile << oneProd.requiresState << ' ' << oneProd.requiresInteraction << '\n';
            }

            // rate list
            restartFile << oneRxn.rateList.size() << '\n';
            for (const auto& oneRate : oneRxn.rateList) {
                restartFile << std::setprecision(20) << oneRate.rate << '\n';
                restartFile << oneRate.otherIfaceLists.size() << '\n';
                for (const auto& oneList : oneRate.otherIfaceLists) {
                    restartFile << oneList.size() << '\n';
                    for (const auto& anccIface : oneList) {
                        restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                                    << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                                    << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                    }
                }
            }
        }

        // creation and destruction reactions
        for (auto& oneRxn : createDestructRxns) {
            restartFile << oneRxn.absRxnIndex << ' ' << oneRxn.relRxnIndex << '\n';
            restartFile << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType) << ' '
                        << oneRxn.isOnMem << '\n';
            restartFile << oneRxn.isObserved << ' ' << oneRxn.observeLabel << '\n';
            restartFile << oneRxn.creationRadius << '\n';

            // integer reactants
            restartFile << oneRxn.intReactantList.size();
            for (const auto& oneReact : oneRxn.intReactantList)
                restartFile << ' ' << oneReact;
            restartFile << '\n';

            // integer products
            restartFile << oneRxn.intProductList.size();
            for (const auto& oneProd : oneRxn.intProductList)
                restartFile << ' ' << oneProd;
            restartFile << std::endl;

            // reactant list
            restartFile << oneRxn.reactantMolList.size() << '\n';
            for (auto& oneReact : oneRxn.reactantMolList) {
                restartFile << oneReact.molTypeIndex << ' ' << oneReact.molName << ' ' << oneReact.interfaceList.size()
                            << '\n';
                for (auto& oneIface : oneReact.interfaceList) {
                    restartFile << oneIface.molTypeIndex << '\n';
                    restartFile << oneIface.ifaceName << ' ' << oneIface.absIfaceIndex << ' ' << oneIface.relIfaceIndex
                                << '\n';
                    restartFile << oneIface.requiresState << ' ' << oneIface.requiresInteraction << '\n';
                }
            }

            // product list
            restartFile << oneRxn.productMolList.size() << '\n';
            for (auto& oneProd : oneRxn.productMolList) {
                restartFile << oneProd.molTypeIndex << ' ' << oneProd.molName << ' ' << oneProd.interfaceList.size()
                            << '\n';
                for (auto& oneIface : oneProd.interfaceList) {
                    restartFile << oneIface.molTypeIndex << '\n';
                    restartFile << oneIface.ifaceName << ' ' << oneIface.absIfaceIndex << ' ' << oneIface.relIfaceIndex
                                << '\n';
                    restartFile << oneIface.requiresState << ' ' << oneIface.requiresInteraction << '\n';
                }
            }

            restartFile << oneRxn.rateList.size() << '\n';
            for (auto& oneRate : oneRxn.rateList) {
                restartFile << std::setprecision(20) << oneRate.rate << '\n';
                restartFile << oneRate.otherIfaceLists.size() << '\n';
                // if (oneRate.otherIfaceLists.size() != 0) {
                //     for (const auto& anccIface : oneRate.otherIfaceLists[0]) {
                //         restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                //                     << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                //                     << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                //     }
                // }
                //  restartFile << oneRate.otherIfaceLists.size() << '\n';
                for (const auto& otherIfaceList : oneRate.otherIfaceLists) {
                    restartFile << otherIfaceList.size() << '\n';
                    for (const auto& anccIface : otherIfaceList) {
                        restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                                    << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                                    << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                    }
                }
            }
        }

        // print out transmission reactions
        for (auto& oneRxn : transmissionRxns) {
            restartFile << oneRxn.absRxnIndex << ' ' << oneRxn.relRxnIndex << '\n';
            restartFile << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType) << ' '
                        << oneRxn.isOnMem << '\n';
            restartFile << oneRxn.isObserved << ' ' << oneRxn.observeLabel << '\n';

            // integer reactants
            restartFile << oneRxn.intReactantList.size();
            for (const auto& oneReact : oneRxn.intReactantList)
                restartFile << ' ' << oneReact;
            restartFile << '\n';

            // integer products
            restartFile << oneRxn.intProductList.size();
            for (const auto& oneProd : oneRxn.intProductList)
                restartFile << ' ' << oneProd;
            restartFile << std::endl;

            // reactant list
            restartFile << oneRxn.reactantMolList.size() << '\n';
            for (auto& oneReact : oneRxn.reactantMolList) {
                restartFile << oneReact.molTypeIndex << ' ' << oneReact.molName << ' ' << oneReact.interfaceList.size()
                            << '\n';
                for (auto& oneIface : oneReact.interfaceList) {
                    restartFile << oneIface.molTypeIndex << '\n';
                    restartFile << oneIface.ifaceName << ' ' << oneIface.absIfaceIndex << ' ' << oneIface.relIfaceIndex
                                << '\n';
                    restartFile << oneIface.requiresState << ' ' << oneIface.requiresInteraction << '\n';
                }
            }

            // product list
            restartFile << oneRxn.productMolList.size() << '\n';
            for (auto& oneProd : oneRxn.productMolList) {
                restartFile << oneProd.molTypeIndex << ' ' << oneProd.molName << ' ' << oneProd.interfaceList.size()
                            << '\n';
                for (auto& oneIface : oneProd.interfaceList) {
                    restartFile << oneIface.molTypeIndex << '\n';
                    restartFile << oneIface.ifaceName << ' ' << oneIface.absIfaceIndex << ' ' << oneIface.relIfaceIndex
                                << '\n';
                    restartFile << oneIface.requiresState << ' ' << oneIface.requiresInteraction << '\n';
                }
            }

            restartFile << oneRxn.rateList.size() << '\n';
            for (auto& oneRate : oneRxn.rateList) {
                restartFile << std::setprecision(20) << oneRate.rate << '\n';
                restartFile << oneRate.otherIfaceLists.size() << '\n';
                // if (oneRate.otherIfaceLists.size() != 0) {
                //     for (const auto& anccIface : oneRate.otherIfaceLists[0]) {
                //         restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                //                     << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                //                     << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                //     }
                // }
                //  restartFile << oneRate.otherIfaceLists.size() << '\n';
                for (const auto& otherIfaceList : oneRate.otherIfaceLists) {
                    restartFile << otherIfaceList.size() << '\n';
                    for (const auto& anccIface : otherIfaceList) {
                        restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                                    << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                                    << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                    }
                }
            }

            // The compartment code takes the crossing interface from
            // reactantListNew and the capture radius from bindRadius, and neither
            // can be rebuilt from the lists above; without them a restart read
            // reactantListNew[0] out of an empty vector in
            // initialize_paramters_for_implicitlipid_and_compartment_model().
            // productListNew goes with it, as in every other reaction record.
            // They follow the rest of the record behind a tag, so that a file
            // written before they existed is refused on read, not misparsed.
            restartFile << "bindRadius = " << oneRxn.bindRadius << '\n';

            // reactant list
            restartFile << oneRxn.reactantListNew.size() << '\n';
            for (const auto& oneReact : oneRxn.reactantListNew) {
                restartFile << oneReact.molTypeIndex << '\n';
                restartFile << oneReact.ifaceName << ' ' << oneReact.absIfaceIndex << ' ' << oneReact.relIfaceIndex
                            << '\n';
                restartFile << oneReact.requiresState << ' ' << oneReact.requiresInteraction << '\n';
            }

            // product list
            restartFile << oneRxn.productListNew.size() << '\n';
            for (const auto& oneProd : oneRxn.productListNew) {
                restartFile << oneProd.molTypeIndex << '\n';
                restartFile << oneProd.ifaceName << ' ' << oneProd.absIfaceIndex << ' ' << oneProd.relIfaceIndex
                            << '\n';
                restartFile << oneProd.requiresState << ' ' << oneProd.requiresInteraction << '\n';
            }
        }

    }

    // write Molecules
    {
        restartFile << "#All Molecules and coordinates \n";
        restartFile << moleculeList.size() << ' ' << Molecule::numberOfMolecules << '\n';
        for (auto& oneMol : moleculeList) {
            restartFile << oneMol.index << ' ' << oneMol.isEmpty << ' ' << oneMol.myComIndex << ' '
                        << oneMol.molTypeIndex << ' ' << oneMol.mySubVolIndex << '\n';
            restartFile << oneMol.mass << ' ' << oneMol.isLipid << ' ' << oneMol.isImplicitLipid 
                        << ' ' << oneMol.linksToSurface << ' ' << oneMol.isPromoter << ' ' << oneMol.isEmpty;
            // Set the first time a molecule is kept from crossing and never
            // cleared, so it is state, not a per-step flag: a molecule outside
            // the compartment reflects off it only while this is set.  Only a
            // compartment sets it, and only a compartment's file carries it.
            if (membraneObject.hasCompartment)
                restartFile << ' ' << oneMol.enforceCompartmentBC;
            restartFile << '\n';
            // center of mass
            restartFile << std::fixed << oneMol.comCoord.x << ' ' << oneMol.comCoord.y << ' ' << oneMol.comCoord.z
                        << '\n';

            // interface lists
            restartFile << oneMol.freelist.size();
            for (const auto& oneIface : oneMol.freelist)
                restartFile << ' ' << oneIface;
            restartFile << '\n';

            restartFile << oneMol.bndlist.size();
            for (const auto& oneIface : oneMol.bndlist)
                restartFile << ' ' << oneIface;
            restartFile << '\n';
            restartFile << oneMol.bndpartner.size();
            for (const auto& oneIface : oneMol.bndpartner)
                restartFile << ' ' << oneIface;
            restartFile << '\n';

            // interfaces
            restartFile << oneMol.interfaceList.size() << '\n';
            for (auto& oneIface : oneMol.interfaceList) {
                restartFile << oneIface.index << ' ' << oneIface.relIndex << ' ' << oneIface.molTypeIndex << ' '
                            << oneIface.stateIndex << ' ' << oneIface.stateIden << ' ' << oneIface.isBound << '\n';
                restartFile << std::fixed << oneIface.coord.x << ' ' << oneIface.coord.y << ' ' << oneIface.coord.z
                            << '\n';

                if (oneIface.isBound) {
                    restartFile << oneIface.interaction.partnerIndex << ' ' << oneIface.interaction.partnerIfaceIndex
                                << ' ' << oneIface.interaction.conjBackRxn << '\n';
                }
            }

            // Reweighting lists.  These are one std::vector<ReweightEntry> in
            // memory now, but the file still carries the six parallel arrays it
            // always did, in the same order and with the same per-line counts,
            // so restart files written by either build are interchangeable.
            restartFile << oneMol.prevReweight.size();
            for (const auto& oneEntry : oneMol.prevReweight)
                restartFile << ' ' << oneEntry.partner;
            restartFile << '\n';
            restartFile << oneMol.prevReweight.size();
            for (const auto& oneEntry : oneMol.prevReweight)
                restartFile << ' ' << oneEntry.myFace;
            restartFile << '\n';
            restartFile << oneMol.prevReweight.size();
            for (const auto& oneEntry : oneMol.prevReweight)
                restartFile << ' ' << oneEntry.partnerFace;
            restartFile << '\n';
            restartFile << oneMol.prevReweight.size();
            for (const auto& oneEntry : oneMol.prevReweight)
                restartFile << ' ' << oneEntry.norm;
            restartFile << '\n';
            restartFile << oneMol.prevReweight.size();
            for (const auto& oneEntry : oneMol.prevReweight)
                restartFile << ' ' << oneEntry.survProb;
            restartFile << '\n';
            restartFile << oneMol.prevReweight.size();
            for (const auto& oneEntry : oneMol.prevReweight)
                restartFile << ' ' << oneEntry.sep;
            restartFile << '\n';
        }

        restartFile << Molecule::emptyMolList.size();
        for (auto& index : Molecule::emptyMolList)
            restartFile << ' ' << index;
        restartFile << '\n';
    }

    // write Complexes
    {
        restartFile << "#All Complexes and their components \n";
        restartFile << complexList.size() << ' ' << Complex::numberOfComplexes << '\n';
        for (const auto& oneCom : complexList) {
            restartFile << oneCom.index << ' ' << oneCom.isEmpty << ' ' << oneCom.radius << ' ' << oneCom.mass << '\n';
            restartFile << oneCom.linksToSurface << ' ' << oneCom.iLipidIndex << ' ' << oneCom.OnSurface << '\n';
            restartFile << oneCom.onFiber << '\n';
            restartFile << std::fixed << oneCom.comCoord.x << ' ' << oneCom.comCoord.y << ' ' << oneCom.comCoord.z
                        << '\n';
            restartFile << std::fixed << oneCom.D.x << ' ' << oneCom.D.y << ' ' << oneCom.D.z << '\n';
            restartFile << std::fixed << oneCom.Dr.x << ' ' << oneCom.Dr.y << ' ' << oneCom.Dr.z << '\n';

            // member molecule lists
            restartFile << oneCom.memberList.size();
            for (const auto& memMol : oneCom.memberList)
                restartFile << ' ' << memMol;
            restartFile << '\n';
            restartFile << oneCom.numEachMol.size();
            for (const auto& memMol : oneCom.numEachMol)
                restartFile << ' ' << memMol;
            restartFile << '\n';

            restartFile << oneCom.lastNumberUpdateItrEachMol.size();
            for (const auto& memMol : oneCom.lastNumberUpdateItrEachMol)
                restartFile << ' ' << memMol;
            restartFile << '\n';
        }

        restartFile << Complex::emptyComList.size();
        for (auto& index : Complex::emptyComList)
            restartFile << ' ' << index;
        restartFile << '\n';
    }

    // Write observables
    {
        restartFile << "#Observables \n";
        restartFile << observablesList.size() << '\n';
        for (auto& observable : observablesList)
            restartFile << observable.first << ' ' << observable.second << '\n';
    }

    // Write counterArrays
    {
        restartFile << "#counterArrays.NLoops .nCancels\n";
        restartFile << counterArrays.nLoops << ' ' << counterArrays.nCancelOverlapPartner << ' ' << counterArrays.nCancelOverlapSystem << ' ' << counterArrays.nCancelDisplace2D << ' ' << counterArrays.nCancelDisplace3D << ' ' << counterArrays.nCancelDisplace3Dto2D << ' ' << counterArrays.nCancelSpanBox << ' ' << counterArrays.nAssocSuccess << ' ' << counterArrays.eventArraySize << '\n';
        //write events3D, 3Dto2D, and 2D
        for (auto& event : counterArrays.events3D)
            restartFile << event << ' ';
        restartFile << '\n';
        for (auto& event : counterArrays.events3Dto2D)
            restartFile << event << ' ';
        restartFile << '\n';
        for (auto& event : counterArrays.events2D)
            restartFile << event << ' ';
        restartFile << '\n';
        // Write bindPairList
        restartFile << "#counterArrays.bindPairList \n";
        restartFile << counterArrays.bindPairList.size() << '\n';
        for (auto& bindPair : counterArrays.bindPairList) {
            restartFile << bindPair.size() << '\n';
            for (auto& oneElem : bindPair)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
        }
        // for (auto& bindPair : counterArrays.bindPairListIL2D) {
        //     restartFile << bindPair.size() << '\n';
        //     for (auto& oneElem : bindPair)
        //         restartFile << ' ' << oneElem;
        //     restartFile << '\n';
        // }
        // for (auto& bindPair : counterArrays.bindPairListIL3D) {
        //     restartFile << bindPair.size() << '\n';
        //     for (auto& oneElem : bindPair)
        //         restartFile << ' ' << oneElem;
        //     restartFile << '\n';
        // }
    }

    // The implicit lipid's 2D binding table, and the protein counts its entries
    // are built from.  Each entry is computed the first time it is needed,
    // from the free-lipid count at that step, and kept for the rest of the
    // run; a restart that rebuilt it would use the count at the restart
    // instead and bind with different probabilities from then on.  The protein
    // counts are those at step 0, which a restart cannot recount from its own
    // molecules once any have been created, destroyed or changed state.
    // Written only for an implicit-lipid model, so no other model's file
    // changes, and last, so a build that predates it reads everything above
    // and never looks here.  Scientific, so every double reads back exactly.
    if (membraneObject.implicitLipid) {
        const std::ios_base::fmtflags savedFlags { restartFile.flags() };
        restartFile << std::scientific;
        restartFile << "#ImplicitLipid\n";
        restartFile << "numberOfProteinEachState =";
        for (int count : membraneObject.numberOfProteinEachState)
            restartFile << ' ' << count;
        restartFile << '\n';
        restartFile << "binding2DTable = " << membraneObject.IL2DbindingVec.size() << '\n';
        for (std::size_t entry { 0 }; entry < membraneObject.IL2DbindingVec.size(); ++entry)
            restartFile << membraneObject.ILTableIDs[3 * entry] << ' ' << membraneObject.ILTableIDs[3 * entry + 1] << ' '
                        << membraneObject.ILTableIDs[3 * entry + 2] << ' ' << membraneObject.IL2DbindingVec[entry] << '\n';
        restartFile.flags(savedFlags);
    }
}
