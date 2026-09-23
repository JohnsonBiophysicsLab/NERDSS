#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <string>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "json.hpp"

namespace {

/*! \brief Consume the next `key =` and abort unless it is the expected key.
 *
 * The reader used to skip to the next `=` without looking at what preceded it,
 * so the header was matched by position alone.  Every field added to the format
 * since therefore shifted an older file's values silently: a restart file
 * written before `accocDissocWrite` existed has its `checkPoint` read into
 * `assocDissocWrite`, its `scaleMaxDisplace` into `checkPoint`, and so on down
 * the block, and the run proceeds with the wrong output intervals and no
 * warning at all.
 * `sample_inputs/gagLatticeRemodeling/4-runRemodelingSim/restart1.dat` is such
 * a file: it predates `accocDissocWrite`, `RNGwrite`, `bondedComplexWrite` and
 * the two compartment fields of `implicitLipidsParams`.
 *
 * Checking the key does not make those files readable - that would need real
 * versioning of the format, and there is no reference output to validate such a
 * reader against - but it turns silent corruption into a diagnostic that names
 * the field which went missing.
 */
void expect_restart_key(std::ifstream& restartFile, const std::string& expected)
{
    // A short or malformed value on the PREVIOUS line leaves the stream in a
    // failed state, and every read after that is a no-op - which is the other
    // half of how an old file used to sail through unnoticed.  Report it here,
    // where the name of the field that was being looked for says where the file
    // stopped making sense.
    if (!restartFile) {
        std::cerr << "Cannot read this restart file: ran out of values before the field '"
                  << expected << "'. The field before it is missing values or malformed.\n"
                  << "This build's format may be newer than the one that wrote the file; "
                     "re-run from the input file instead." << std::endl;
        exit(1);
    }

    // Bounded. The ignore(max, '=') this replaced allocated nothing, so a file
    // corrupted such that no '=' ever follows used to run off the end harmlessly;
    // accumulating into a string would instead buffer the entire remainder, and
    // restart files for large systems run to hundreds of megabytes. No real key
    // is anywhere near this long, so hitting the cap is itself a malformed file.
    const std::size_t maxKeyLength { 256 };
    std::string key;
    char c;
    while (key.size() < maxKeyLength && restartFile.get(c)) {
        if (c == '=')
            break;
        key += c;
    }

    const std::size_t first { key.find_first_not_of(" \t\r\n") };
    const std::string found { first == std::string::npos
            ? std::string()
            : key.substr(first, key.find_last_not_of(" \t\r\n") - first + 1) };

    if (found == expected)
        return;

    // Inside a record, where there is no key, the text up to the next '=' is
    // many lines of data; its first line is enough to show what was there.
    std::cerr << "Cannot read this restart file: expected the field '" << expected
              << "' next, but found '" << found.substr(0, found.find('\n')) << "'.\n"
              << "The restart format has gained fields over time and is matched by "
                 "position, so a file written by an older build cannot be read by "
                 "this one. Re-run from the input file instead." << std::endl;
    exit(1);
}

} // namespace

/* JSON restart files; see write_restart.cpp for the layout.  Every field is
 * looked up by name, so a missing or mistyped one is reported by name, and a
 * field this branch added is simply absent from a file written without it.
 */
namespace {

using json = nlohmann::json;

// A null is a double that was NaN when written; see maybe_null() in
// write_restart.cpp.
double read_double(const json& j)
{
    return j.is_null() ? std::numeric_limits<double>::quiet_NaN() : j.get<double>();
}

void read_vec3(const json& j, Vec3D& v)
{
    v.x = read_double(j.at(0));
    v.y = read_double(j.at(1));
    v.z = read_double(j.at(2));
}

template <typename T>
std::vector<T> get_vector(const json& j, const std::string& key)
{
    if (j.contains(key))
        return j.at(key).get<std::vector<T>>();
    return {};
}

std::vector<double> get_double_vector(const json& j, const std::string& key)
{
    std::vector<double> values;
    if (j.contains(key)) {
        for (const auto& elem : j.at(key))
            values.push_back(read_double(elem));
    }
    return values;
}

void rxniface_from_json(const json& j, RxnIface& iface)
{
    iface.molTypeIndex = j.at("molTypeIndex");
    iface.ifaceName = j.at("ifaceName").get<std::string>();
    iface.absIfaceIndex = j.at("absIfaceIndex");
    iface.relIfaceIndex = j.at("relIfaceIndex");
    iface.requiresState = static_cast<char>(j.at("requiresState").get<int>());
    iface.requiresInteraction = j.at("requiresInteraction");
}

void iface_list_from_json(const json& j, std::vector<RxnIface>& ifaceList)
{
    for (const auto& ji : j) {
        RxnIface iface {};
        rxniface_from_json(ji, iface);
        ifaceList.push_back(iface);
    }
}

void state_change_iface_from_json(const json& j, std::pair<RxnIface, RxnIface>& stateChangeIface)
{
    rxniface_from_json(j.at(0), stateChangeIface.first);
    rxniface_from_json(j.at(1), stateChangeIface.second);
}

void rate_list_from_json(const json& j, std::vector<RxnBase::RateState>& rateList)
{
    for (const auto& jr : j) {
        RxnBase::RateState rate {};
        rate.rate = read_double(jr.at("rate"));
        if (jr.contains("otherIfaceLists")) {
            for (const auto& list : jr.at("otherIfaceLists")) {
                std::vector<RxnIface> ifaceList;
                iface_list_from_json(list, ifaceList);
                rate.otherIfaceLists.push_back(ifaceList);
            }
        }
        rateList.push_back(rate);
    }
}

void coupled_rxn_from_json(const json& j, RxnBase::CoupledRxn& coupled)
{
    coupled.absRxnIndex = j.at("absRxnIndex");
    coupled.relRxnIndex = j.at("relRxnIndex");
    coupled.rxnType = static_cast<ReactionType>(j.at("rxnType").get<int>());
    coupled.label = j.at("label").get<std::string>();
    coupled.probCoupled = j.at("probCoupled");
}

// CreateDestructRxn::CreateDestructMol and TransmissionRxn::TransmissionMol.
template <typename MolWithIfaces>
void mol_with_ifaces_from_json(const json& j, MolWithIfaces& mol)
{
    mol.molTypeIndex = j.at("molTypeIndex");
    mol.molName = j.at("molName").get<std::string>();
    // nerdss_development's json-restarts branch writes this list as
    // "interfaceList" and reads it as "interfaces", which fails an assertion
    // on any file with a creation, destruction or transmission reaction.  The
    // writer's name is the one in every file.
    iface_list_from_json(j.at("interfaceList"), mol.interfaceList);
}

void iface_state_from_json(const json& j, Interface::State& state)
{
    state.index = j.at("index");
    state.iden = static_cast<char>(j.at("iden").get<int>());
    state.rxnPartners = get_vector<unsigned>(j, "rxnPartners");
    state.myForwardRxns = get_vector<unsigned>(j, "myForwardRxns");
    state.myCreateDestructRxns = get_vector<unsigned>(j, "myCreateDestructRxns");
    if (j.contains("stateChangeRxns")) {
        for (const auto& pair : j.at("stateChangeRxns"))
            state.stateChangeRxns.emplace_back(pair.at(0).get<int>(), pair.at(1).get<int>());
    }
}

void iface_from_json(const json& j, Interface& iface)
{
    iface.index = j.at("index");
    iface.name = j.at("name").get<std::string>();
    read_vec3(j.at("coord"), iface.iCoord);
    for (const auto& js : j.at("states")) {
        Interface::State state {};
        iface_state_from_json(js, state);
        iface.stateList.push_back(state);
    }
}

void moltemplate_from_json(const json& j, MolTemplate& mt)
{
    mt.molTypeIndex = j.at("molTypeIndex");
    mt.molName = j.at("molName").get<std::string>();
    mt.copies = j.at("copies");
    mt.mass = j.at("mass");
    mt.radius = j.at("radius");
    mt.isLipid = j.at("isLipid");
    mt.isImplicitLipid = j.at("isImplicitLipid");
    mt.isRod = j.at("isRod");
    mt.isPoint = j.at("isPoint");
    mt.checkOverlap = j.at("checkOverlap");
    mt.countTransition = j.at("countTransition");
    mt.transitionMatrixSize = j.at("transitionMatrixSize");
    mt.outsideCompartment = j.at("outsideCompartment");
    mt.insideCompartment = j.at("insideCompartment");
    mt.crossesCompartment = j.at("crossesCompartment");
    mt.transmissionRxnIndex = j.at("transmissionRxnIndex");
    read_vec3(j.at("comCoord"), mt.comCoord);
    read_vec3(j.at("D"), mt.D);
    read_vec3(j.at("Dr"), mt.Dr);
    // Without this invCbrtDr stays zero, Complex::update_properties() sums it
    // to zero and divides by that, and every complex with a rotating member
    // comes out with Dr = inf and NaN coordinates after its first step.
    mt.cache_diffusion_derivatives();

    mt.rxnPartners = get_vector<int>(j, "rxnPartners");
    if (j.contains("bondList")) {
        for (const auto& bond : j.at("bondList"))
            mt.bondList.push_back(std::array<int, 2> { { bond.at(0).get<int>(), bond.at(1).get<int>() } });
    }
    for (const auto& ji : j.at("interfaces")) {
        Interface iface {};
        iface_from_json(ji, iface);
        mt.interfaceList.push_back(iface);
    }
    mt.ifacesWithStates = get_vector<int>(j, "ifacesWithStates");
    mt.monomerList = get_vector<int>(j, "monomerList");
    if (mt.countTransition) {
        mt.lifeTime = j.at("lifeTime").get<std::vector<std::vector<double>>>();
        mt.transitionMatrix = j.at("transitionMatrix").get<std::vector<std::vector<long long int>>>();
    }
}

void forward_rxn_from_json(const json& j, ForwardRxn& rxn)
{
    rxn.absRxnIndex = j.at("absRxnIndex");
    rxn.relRxnIndex = j.at("relRxnIndex");
    rxn.rxnLabel = j.at("rxnLabel").get<std::string>();
    rxn.rxnType = static_cast<ReactionType>(j.at("rxnType").get<int>());
    rxn.isSymmetric = j.at("isSymmetric");
    rxn.isOnMem = j.at("isOnMem");
    rxn.hasStateChange = j.at("hasStateChange");
    rxn.isObserved = j.at("isObserved");
    rxn.observeLabel = j.at("observeLabel").get<std::string>();
    rxn.productName = j.at("productName").get<std::string>();
    rxn.isReversible = j.at("isReversible");
    rxn.conjBackRxnIndex = j.at("conjBackRxnIndex");
    rxn.irrevRingClosure = j.at("irrevRingClosure");
    rxn.bindRadSameCom = j.at("bindRadSameCom");
    rxn.loopCoopFactor = j.at("loopCoopFactor");
    rxn.length3Dto2D = j.at("length3Dto2D");
    rxn.area3Dto1D = j.at("area3Dto1D");
    rxn.bindRadius = j.at("bindRadius");
    const json& jAngles = j.at("assocAngles");
    rxn.assocAngles.theta1 = read_double(jAngles.at("theta1"));
    rxn.assocAngles.theta2 = read_double(jAngles.at("theta2"));
    rxn.assocAngles.phi1 = read_double(jAngles.at("phi1"));
    rxn.assocAngles.phi2 = read_double(jAngles.at("phi2"));
    rxn.assocAngles.omega = read_double(jAngles.at("omega"));
    read_vec3(j.at("norm1"), rxn.norm1);
    read_vec3(j.at("norm2"), rxn.norm2);
    rxn.excludeVolumeBound = j.at("excludeVolumeBound");
    rxn.isCoupled = j.at("isCoupled");
    if (rxn.isCoupled)
        coupled_rxn_from_json(j.at("coupledRxn"), rxn.coupledRxn);
    rxn.intReactantList = get_vector<int>(j, "intReactantList");
    rxn.intProductList = get_vector<int>(j, "intProductList");
    iface_list_from_json(j.at("reactantList"), rxn.reactantListNew);
    iface_list_from_json(j.at("productList"), rxn.productListNew);
    rate_list_from_json(j.at("rateList"), rxn.rateList);
    if (j.contains("stateChangeIface"))
        state_change_iface_from_json(j.at("stateChangeIface"), rxn.stateChangeIface);
}

void back_rxn_from_json(const json& j, BackRxn& rxn)
{
    rxn.absRxnIndex = j.at("absRxnIndex");
    rxn.relRxnIndex = j.at("relRxnIndex");
    rxn.rxnType = static_cast<ReactionType>(j.at("rxnType").get<int>());
    rxn.isSymmetric = j.at("isSymmetric");
    rxn.isOnMem = j.at("isOnMem");
    rxn.hasStateChange = j.at("hasStateChange");
    rxn.isObserved = j.at("isObserved");
    rxn.observeLabel = j.at("observeLabel").get<std::string>();
    rxn.conjForwardRxnIndex = j.at("conjForwardRxnIndex");
    rxn.isCoupled = j.at("isCoupled");
    if (rxn.isCoupled)
        coupled_rxn_from_json(j.at("coupledRxn"), rxn.coupledRxn);
    rxn.intReactantList = get_vector<int>(j, "intReactantList");
    rxn.intProductList = get_vector<int>(j, "intProductList");
    iface_list_from_json(j.at("reactantList"), rxn.reactantListNew);
    iface_list_from_json(j.at("productList"), rxn.productListNew);
    rate_list_from_json(j.at("rateList"), rxn.rateList);
    // Absent from a file written by nerdss_development; see write_restart().
    if (j.contains("stateChangeIface"))
        state_change_iface_from_json(j.at("stateChangeIface"), rxn.stateChangeIface);
}

void createdestruct_rxn_from_json(const json& j, CreateDestructRxn& rxn)
{
    rxn.absRxnIndex = j.at("absRxnIndex");
    rxn.relRxnIndex = j.at("relRxnIndex");
    rxn.rxnType = static_cast<ReactionType>(j.at("rxnType").get<int>());
    rxn.isOnMem = j.at("isOnMem");
    rxn.isObserved = j.at("isObserved");
    rxn.observeLabel = j.at("observeLabel").get<std::string>();
    rxn.creationRadius = j.at("creationRadius");
    rxn.intReactantList = get_vector<int>(j, "intReactantList");
    rxn.intProductList = get_vector<int>(j, "intProductList");
    for (const auto& jm : j.at("reactantMolList")) {
        CreateDestructRxn::CreateDestructMol mol {};
        mol_with_ifaces_from_json(jm, mol);
        rxn.reactantMolList.push_back(mol);
    }
    for (const auto& jm : j.at("productMolList")) {
        CreateDestructRxn::CreateDestructMol mol {};
        mol_with_ifaces_from_json(jm, mol);
        rxn.productMolList.push_back(mol);
    }
    rate_list_from_json(j.at("rateList"), rxn.rateList);
}

void transmission_rxn_from_json(const json& j, TransmissionRxn& rxn)
{
    rxn.absRxnIndex = j.at("absRxnIndex");
    rxn.relRxnIndex = j.at("relRxnIndex");
    rxn.rxnType = static_cast<ReactionType>(j.at("rxnType").get<int>());
    rxn.isOnMem = j.at("isOnMem");
    rxn.isObserved = j.at("isObserved");
    rxn.observeLabel = j.at("observeLabel").get<std::string>();
    rxn.intReactantList = get_vector<int>(j, "intReactantList");
    rxn.intProductList = get_vector<int>(j, "intProductList");
    for (const auto& jm : j.at("reactantMolList")) {
        TransmissionRxn::TransmissionMol mol {};
        mol_with_ifaces_from_json(jm, mol);
        rxn.reactantMolList.push_back(mol);
    }
    for (const auto& jm : j.at("productMolList")) {
        TransmissionRxn::TransmissionMol mol {};
        mol_with_ifaces_from_json(jm, mol);
        rxn.productMolList.push_back(mol);
    }
    rate_list_from_json(j.at("rateList"), rxn.rateList);
    // Required, as the .dat reader requires its bindRadius tag: without these
    // a restart read reactantListNew[0] out of an empty vector in
    // initialize_paramters_for_implicitlipid_and_compartment_model().  A file
    // written by nerdss_development lacks them and is refused here by name.
    rxn.bindRadius = j.at("bindRadius");
    iface_list_from_json(j.at("reactantList"), rxn.reactantListNew);
    iface_list_from_json(j.at("productList"), rxn.productListNew);
}

void numerics_from_json(const json& j, NumericalSettings& numerics)
{
    numerics.integration.tableAbsoluteError = j.at("integrationAbsError");
    numerics.integration.tableRelativeError = j.at("integrationRelError");
    numerics.integration.fallbackError = j.at("integrationFallbackError");
    numerics.integration.tailCutoff = j.at("integrationTailCutoff");
    numerics.integration.normalizationAbsoluteError = j.at("normalizationAbsError");
    numerics.integration.normalizationRelativeError = j.at("normalizationRelError");
    numerics.tableLookup.reactionRate.absolute = j.at("tableRateAbsTolerance");
    numerics.tableLookup.reactionRate.relative = j.at("tableRateRelTolerance");
    numerics.tableLookup.diffusionCoefficient.absolute = j.at("tableDiffusionAbsTolerance");
    numerics.tableLookup.diffusionCoefficient.relative = j.at("tableDiffusionRelTolerance");
    numerics.classification.explicitLipidFlatDiffusion = j.at("explicitLipidFlatDiffusion");
    numerics.classification.implicitLipidFlatDiffusion = j.at("implicitLipidFlatDiffusion");
    numerics.associationAngles.sameAngle.absolute = j.at("associationSameAngleAbsTolerance");
    numerics.associationAngles.sameAngle.relative = j.at("associationSameAngleRelTolerance");
    numerics.associationAngles.rotationConvergenceTolerance = j.at("associationRotationTolerance");
    numerics.associationAngles.endpointSignTolerance = j.at("associationEndpointSignTolerance");
    numerics.vec3D.coordinateEqualityPrecision = j.at("vec3DCoordinatePrecision");
}

void read_json_restart(long long int& simItr, std::ifstream& restartFile, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<TransmissionRxn>& transmissionRxns,
    std::map<std::string, int>& observablesList, Membrane& membraneObject, copyCounters& counterArrays)
{
    json j;
    restartFile >> j;

    // parameters
    std::cout << "Reading parameters..." << std::endl;
    const json& jp = j.at("parameters");
    params.nItr = jp.at("nItr");
    simItr = jp.at("simItr");
    params.itrRestartFrom = simItr;
    params.timeRestartFrom = jp.at("currSimTime");
    std::cout << "Restarting simulation from iteration " << simItr << '\n';
    std::cout << "Current simulation time (s): " << params.timeRestartFrom << '\n';
    params.numMolTypes = jp.at("numMolTypes");
    params.numTotalSpecies = jp.at("numTotalSpecies");
    params.numTotalComplex = jp.at("numTotalComplex");
    params.numTotalUnits = jp.at("numTotalUnits");
    params.numLipids = jp.at("numLipids");
    params.timeStep = jp.at("timeStep");
    params.max2DRxns = jp.at("max2DRxns");
    params.overlapSepLimit = jp.at("overlapSepLimit");
    params.rMaxLimit = jp.at("rMaxLimit");
    params.timeWrite = jp.at("timeWrite");
    params.trajWrite = jp.at("trajWrite");
    params.restartWrite = jp.at("restartWrite");
    params.pdbWrite = jp.at("pdbWrite");
    params.assocDissocWrite = jp.at("assocDissocWrite");
    params.checkPoint = jp.at("checkPoint");
    params.scaleMaxDisplace = jp.at("scaleMaxDisplace");
    params.transitionWrite = jp.at("transitionWrite");
    params.clusterOverlapCheck = jp.at("clusterOverlapCheck");
    params.rngwrite = jp.at("rngwrite");
    params.bondedComplexWrite = jp.at("bondedComplexWrite");
    Parameters::lastUpdateTransition = jp.at("lastUpdateTransition").get<std::vector<long long int>>();

    // Always written by this build.  A file from a build without them keeps
    // the defaults, as a .dat file without the #NumericalSettings block does.
    if (jp.contains("numerics")) {
        numerics_from_json(jp.at("numerics"), params.numerics);
        try {
            params.numerics.validate();
        } catch (const std::invalid_argument& error) {
            throw std::runtime_error(std::string("invalid numerical settings: ") + error.what());
        }
    }

    // membrane
    std::cout << "Reading membrane..." << std::endl;
    const json& jm = jp.at("membrane");
    const std::vector<double> boxDimensions { jp.at("waterBox").get<std::vector<double>>() };
    if (boxDimensions.size() != 3)
        throw std::runtime_error("waterBox must hold three lengths");
    // Built by the constructor, as parse_input() builds it, so that xLeft and
    // xRight are set along with the volume: create_random_coords() places a
    // molecule created in a box at x = xLeft + (xRight - xLeft) * rand.
    // nerdss_mpi overwrites both with its rank's bounds in prepare.cpp.
    membraneObject.waterBox = Membrane::WaterBox(boxDimensions);
    membraneObject.implicitlipidIndex = jm.at("implicitlipidIndex");
    membraneObject.nSites = jm.at("nSites");
    membraneObject.nStates = jm.at("nStates");
    membraneObject.No_free_lipids = jm.at("No_free_lipids");
    membraneObject.No_protein = jm.at("No_protein");
    membraneObject.totalSA = jm.at("totalSA");
    membraneObject.numberOfFreeLipidsEachState = jp.at("numberOfFreeLipidsEachState").get<std::vector<int>>();

    const json& jil = jp.at("implicitLipidParams");
    membraneObject.implicitLipid = jil.at("implicitLipid");
    membraneObject.TwoD = jil.at("TwoD");
    membraneObject.sphereR = jil.at("sphereR");
    membraneObject.hasCompartment = jil.at("hasCompartment");
    membraneObject.compartmentR = jil.at("compartmentR");
    // The two flags map back one-to-one: the sphere flag is the geometry, the
    // box flag the waterBox provenance.
    const bool isBoxFlag { jil.at("isBox").get<bool>() };
    const bool isSphereFlag { jil.at("isSphere").get<bool>() };
    membraneObject.waterBoxGiven = isBoxFlag;
    membraneObject.shape = isSphereFlag ? BoundaryShape::Sphere
        : isBoxFlag                     ? BoundaryShape::Box
                                        : BoundaryShape::Unspecified;

    // A compartment file from a build that did not write them stops here,
    // rather than restarting with the sites' D and density at zero.
    if (jm.contains("compartmentSiteD")) {
        membraneObject.droplet.D = jm.at("compartmentSiteD");
        membraneObject.droplet.rho = jm.at("compartmentSiteRho");
    } else if (membraneObject.hasCompartment) {
        throw std::runtime_error("this compartment restart file does not record the compartment sites' "
                                 "diffusion constant and density (parameters.membrane.compartmentSiteD and "
                                 "compartmentSiteRho); re-run from the input file instead");
    }

    // molecule templates
    std::cout << "Reading molecule templates..." << std::endl;
    const json& jtemplates = j.at("molTemplates");
    MolTemplate::absToRelIface = jtemplates.at("absToRelIface").get<std::vector<int>>();
    for (const auto& jt : jtemplates.at("templates")) {
        MolTemplate oneTemp {};
        moltemplate_from_json(jt, oneTemp);
        molTemplateList.push_back(oneTemp);
    }
    MolTemplate::numMolTypes = molTemplateList.size();
    // Kept by the reactions that create and destroy molecules; recounted from
    // the molecules below only for a file that does not carry it.
    const bool hasNumEachMolType { jtemplates.contains("numEachMolType") };
    if (hasNumEachMolType)
        MolTemplate::numEachMolType = jtemplates.at("numEachMolType").get<std::vector<int>>();
    else
        MolTemplate::numEachMolType = std::vector<int>(MolTemplate::numMolTypes, 0);
    if (MolTemplate::numEachMolType.size() != MolTemplate::numMolTypes)
        throw std::runtime_error("molTemplates.numEachMolType does not have one count per template");
    // The states an add file's templates bring are numbered from here
    // (parse_molFile.cpp).  Every state has an index, so a file without it
    // gets the count of states.
    if (jtemplates.contains("totalNumOfStates")) {
        Interface::State::totalNumOfStates = jtemplates.at("totalNumOfStates");
    } else {
        int totalStates { 0 };
        for (const auto& oneTemp : molTemplateList)
            for (const auto& oneIface : oneTemp.interfaceList)
                totalStates += static_cast<int>(oneIface.stateList.size());
        Interface::State::totalNumOfStates = totalStates;
    }

    // reactions
    std::cout << "Reading reactions..." << std::endl;
    const json& jr = j.at("reactions");
    RxnBase::numberOfRxns = jr.at("numberOfRxns");
    RxnBase::totRxnSpecies = jr.at("totRxnSpecies");
    for (const auto& jrxn : jr.at("forward")) {
        ForwardRxn rxn;
        forward_rxn_from_json(jrxn, rxn);
        forwardRxns.push_back(rxn);
    }
    for (const auto& jrxn : jr.at("back")) {
        BackRxn rxn;
        back_rxn_from_json(jrxn, rxn);
        backRxns.push_back(rxn);
    }
    for (const auto& jrxn : jr.at("createDestruct")) {
        CreateDestructRxn rxn {};
        createdestruct_rxn_from_json(jrxn, rxn);
        createDestructRxns.push_back(rxn);
    }
    for (const auto& jrxn : jr.at("transmission")) {
        TransmissionRxn rxn {};
        transmission_rxn_from_json(jrxn, rxn);
        transmissionRxns.push_back(rxn);
    }

    // molecules
    std::cout << "Reading molecules..." << std::endl;
    const json& jmols = j.at("molecules");
    for (const auto& jmol : jmols.at("list")) {
        Molecule oneMol {};
        oneMol.index = jmol.at("index");
        oneMol.isEmpty = jmol.at("isEmpty");
        oneMol.myComIndex = jmol.at("myComIndex");
        oneMol.molTypeIndex = jmol.at("molTypeIndex");
        oneMol.mySubVolIndex = jmol.at("mySubVolIndex");
        oneMol.mass = jmol.at("mass");
        oneMol.isLipid = jmol.at("isLipid");
        oneMol.isImplicitLipid = jmol.at("isImplicitLipid");
        oneMol.linksToSurface = jmol.at("linksToSurface");
        oneMol.isPromoter = jmol.at("isPromoter");
        if (jmol.contains("enforceCompartmentBC")) {
            oneMol.enforceCompartmentBC = jmol.at("enforceCompartmentBC");
        } else if (membraneObject.hasCompartment) {
            throw std::runtime_error("this compartment restart file does not record enforceCompartmentBC for "
                                     "its molecules; re-run from the input file instead");
        }
        read_vec3(jmol.at("comCoord"), oneMol.comCoord);
        oneMol.freelist = get_vector<int>(jmol, "freelist");
        oneMol.bndlist = get_vector<int>(jmol, "bndlist");
        oneMol.bndpartner = get_vector<int>(jmol, "bndpartner");

        for (const auto& jiface : jmol.at("interfaceList")) {
            Molecule::Iface iface {};
            iface.index = jiface.at("index");
            iface.relIndex = jiface.at("relIndex");
            iface.molTypeIndex = jiface.at("molTypeIndex");
            iface.stateIndex = jiface.at("stateIndex");
            iface.stateIden = static_cast<char>(jiface.at("stateIden").get<int>());
            iface.isBound = jiface.at("isBound");
            read_vec3(jiface.at("coord"), iface.coord);
            if (iface.isBound) {
                const json& jinteraction = jiface.at("interaction");
                iface.interaction.partnerIndex = jinteraction.at("partnerIndex");
                iface.interaction.partnerIfaceIndex = jinteraction.at("partnerIfaceIndex");
                iface.interaction.conjBackRxn = jinteraction.at("conjBackRxn");
            }
            oneMol.interfaceList.push_back(iface);
        }

        // Reweighting lists: the six parallel arrays of the file are assembled
        // into the single prevReweight vector.
        const std::vector<int> prevlist { get_vector<int>(jmol, "prevlist") };
        const std::vector<int> prevmyface { get_vector<int>(jmol, "prevmyface") };
        const std::vector<int> prevpface { get_vector<int>(jmol, "prevpface") };
        const std::vector<double> prevnorm { get_double_vector(jmol, "prevnorm") };
        const std::vector<double> ps_prev { get_double_vector(jmol, "ps_prev") };
        const std::vector<double> prevsep { get_double_vector(jmol, "prevsep") };
        const std::size_t numEntries { prevlist.size() };
        if (prevmyface.size() != numEntries || prevpface.size() != numEntries || prevnorm.size() != numEntries
            || ps_prev.size() != numEntries || prevsep.size() != numEntries) {
            throw std::runtime_error("molecule " + std::to_string(oneMol.index)
                + ": the six reweighting arrays differ in length");
        }
        oneMol.prevReweight.resize(numEntries);
        for (std::size_t entry { 0 }; entry < numEntries; ++entry) {
            oneMol.prevReweight[entry].partner = prevlist[entry];
            oneMol.prevReweight[entry].myFace = prevmyface[entry];
            oneMol.prevReweight[entry].partnerFace = prevpface[entry];
            oneMol.prevReweight[entry].norm = prevnorm[entry];
            oneMol.prevReweight[entry].survProb = ps_prev[entry];
            oneMol.prevReweight[entry].sep = prevsep[entry];
        }

        // trajStatus is not in the file, and every other molecule can do
        // without it: the end of each timestep resets theirs to `none`, which
        // is what they all hold at a checkpoint.  That reset skips the implicit
        // lipid's one representative molecule, which holds `propagated` from
        // its first step on; left at `none`, a restart propagated it once more
        // and diverged.  A file written at step 0 is from before that first
        // propagation, where `none` is right.
        if (oneMol.isImplicitLipid && simItr > 0)
            oneMol.trajStatus = TrajStatus::propagated;

        if (!hasNumEachMolType && !oneMol.isEmpty) {
            if (oneMol.molTypeIndex < 0 || oneMol.molTypeIndex >= static_cast<int>(MolTemplate::numMolTypes))
                throw std::runtime_error("molecule " + std::to_string(oneMol.index) + " has no template");
            ++MolTemplate::numEachMolType[oneMol.molTypeIndex];
        }
        moleculeList.push_back(oneMol);
    }
    // The list keeps a slot for every molecule ever created; the count is of
    // those not destroyed.
    if (jmols.contains("numberOfMolecules"))
        Molecule::numberOfMolecules = jmols.at("numberOfMolecules");
    else
        Molecule::numberOfMolecules = static_cast<int>(moleculeList.size());
    Molecule::emptyMolList = get_vector<int>(jmols, "emptyMolList");

    // complexes
    std::cout << "Reading complexes..." << std::endl;
    const json& jcoms = j.at("complexes");
    Complex::numberOfComplexes = jcoms.at("numberOfComplexes");
    for (const auto& jc : jcoms.at("list")) {
        Complex oneCom {};
        oneCom.index = jc.at("index");
        oneCom.isEmpty = jc.at("isEmpty");
        oneCom.radius = jc.at("radius");
        oneCom.mass = jc.at("mass");
        oneCom.linksToSurface = jc.at("linksToSurface");
        oneCom.iLipidIndex = jc.at("iLipidIndex");
        oneCom.OnSurface = jc.at("OnSurface");
        oneCom.onFiber = jc.at("onFiber");
        read_vec3(jc.at("comCoord"), oneCom.comCoord);
        read_vec3(jc.at("D"), oneCom.D);
        read_vec3(jc.at("Dr"), oneCom.Dr);
        oneCom.memberList = get_vector<int>(jc, "memberList");
        oneCom.numEachMol = get_vector<int>(jc, "numEachMol");
        oneCom.lastNumberUpdateItrEachMol = get_vector<long long int>(jc, "lastNumberUpdateItrEachMol");
        complexList.push_back(oneCom);
    }
    Complex::emptyComList = get_vector<int>(jcoms, "emptyComList");

    // observables
    std::cout << "Reading observables..." << std::endl;
    observablesList.clear();
    for (const auto& observable : j.at("observables"))
        observablesList.emplace(observable.at("name").get<std::string>(), observable.at("value").get<int>());

    // counter arrays
    std::cout << "Reading counter arrays..." << std::endl;
    const json& jcounters = j.at("counterArrays");
    counterArrays.nLoops = jcounters.at("nLoops");
    counterArrays.nCancelOverlapPartner = jcounters.at("nCancelOverlapPartner");
    counterArrays.nCancelOverlapSystem = jcounters.at("nCancelOverlapSystem");
    counterArrays.nCancelDisplace2D = jcounters.at("nCancelDisplace2D");
    counterArrays.nCancelDisplace3D = jcounters.at("nCancelDisplace3D");
    counterArrays.nCancelDisplace3Dto2D = jcounters.at("nCancelDisplace3Dto2D");
    counterArrays.nCancelSpanBox = jcounters.at("nCancelSpanBox");
    counterArrays.nAssocSuccess = jcounters.at("nAssocSuccess");
    counterArrays.eventArraySize = jcounters.at("eventArraySize");
    counterArrays.events3D = jcounters.at("events3D").get<std::vector<int>>();
    counterArrays.events3Dto2D = jcounters.at("events3Dto2D").get<std::vector<int>>();
    counterArrays.events2D = jcounters.at("events2D").get<std::vector<int>>();
    counterArrays.bindPairList = jcounters.at("bindPairList").get<std::vector<std::vector<int>>>();

    // The implicit lipid's 2D binding table and the protein counts it is built
    // from; see write_restart().  A file from before they were written still
    // restarts, but both are then rebuilt from the state at the restart, so the
    // run cannot continue exactly.
    if (j.contains("implicitLipid")) {
        const json& jlipid = j.at("implicitLipid");
        membraneObject.numberOfProteinEachState = jlipid.at("numberOfProteinEachState").get<std::vector<int>>();
        membraneObject.ILTableIDs.clear();
        membraneObject.IL2DbindingVec.clear();
        for (const auto& row : jlipid.at("binding2DTable")) {
            membraneObject.ILTableIDs.push_back(read_double(row.at("ka")));
            membraneObject.ILTableIDs.push_back(read_double(row.at("Dtot")));
            membraneObject.ILTableIDs.push_back(read_double(row.at("kb")));
            membraneObject.IL2DbindingVec.push_back(read_double(row.at("probability")));
        }
    } else if (membraneObject.implicitLipid) {
        std::cout << "This restart file has no implicitLipid section, so it predates saving the implicit lipid's 2D "
                     "binding table and protein counts. They are rebuilt from the state at the restart, and the run "
                     "will not continue exactly as the one that wrote this file would have."
                  << std::endl;
    }
}

/*! \brief Whether the restart file is JSON, by its first non-blank byte.
 *
 * A JSON restart file is one object and starts with '{'; a .dat file starts
 * with its "#Parameters" line.  The name is not consulted: nerdss_mpi appends
 * the rank after the extension (restart.dat0), so an extension test would
 * send every parallel restart down the wrong path.
 */
bool restart_file_is_json(std::istream& restartFile)
{
    const std::streampos start { restartFile.tellg() };
    char c { '\0' };
    while (restartFile.get(c) && std::isspace(static_cast<unsigned char>(c))) { }
    const bool isJson { restartFile.good() && c == '{' };
    restartFile.clear();
    restartFile.seekg(start);
    return isJson;
}

} // namespace

void read_restart(long long int& simItr, std::ifstream& restartFile, Parameters& params, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<TransmissionRxn>& transmissionRxns,
    std::map<std::string, int>& observablesList, Membrane& membraneObject, copyCounters& counterArrays)
{
    // A restart keeps writing the format it read, so a run that writes .dat
    // files goes on doing so across restarts without carrying the flag in the
    // .dat format, whose layout must not change.  An add file parsed after
    // this can set legacyRestartFormat either way.
    if (!restart_file_is_json(restartFile)) {
        std::cout << "The restart file is in the legacy .dat format." << std::endl;
        LEGACY_read_restart(simItr, restartFile, params, simulVolume, moleculeList, complexList, molTemplateList,
            forwardRxns, backRxns, createDestructRxns, transmissionRxns, observablesList, membraneObject,
            counterArrays);
        params.legacyRestartFormat = true;
        return;
    }

    std::cout << "The restart file is in the JSON format." << std::endl;
    try {
        read_json_restart(simItr, restartFile, params, moleculeList, complexList, molTemplateList, forwardRxns,
            backRxns, createDestructRxns, transmissionRxns, observablesList, membraneObject, counterArrays);
    } catch (const std::exception& e) {
        // json::exception names the key or the type that was wrong.
        std::cerr << "Cannot read this JSON restart file: " << e.what() << std::endl;
        exit(1);
    }
    params.legacyRestartFormat = false;
    std::cout << "Finished reading restart file." << std::endl;
}

/* The .dat format: positional, one value after another, matched by the order
 * the writer used.  Files in it come from earlier builds, or from a run with
 * legacyRestartFormat set; read_restart() sends each file to the reader for
 * its format.
 */

void LEGACY_read_restart(long long int& simItr, std::ifstream& restartFile, Parameters& params, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<TransmissionRxn>& transmissionRxns,
    std::map<std::string, int>& observablesList, Membrane& membraneObject, copyCounters& counterArrays)
{
    // TRACE();
    try {
        // Read parameters
        std::cout << "READ IN PARMATERS from restart file" << std::endl;
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        {
            expect_restart_key(restartFile, "numItr");
            restartFile >> params.nItr;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "currItr");
            restartFile >> simItr;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Restarting simulation from iteration " << simItr << '\n';
            params.itrRestartFrom = simItr;

            expect_restart_key(restartFile, "currSimTime (s)");
            restartFile >> params.timeRestartFrom;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Current simulation time (s): " << params.timeRestartFrom << '\n';

            expect_restart_key(restartFile, "numMolTypes");
            restartFile >> params.numMolTypes;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "numTotalSpecies");
            restartFile >> params.numTotalSpecies;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "numComplexs");
            restartFile >> params.numTotalComplex;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "numTotalUnits");
            restartFile >> params.numTotalUnits;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "numLipids");
            restartFile >> params.numLipids;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "timestep");
            restartFile >> params.timeStep;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "max2DRxns");
            restartFile >> params.max2DRxns;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "simulDimensions");
            std::vector<double> boxDimensions(3);
            restartFile >> boxDimensions[0] >> boxDimensions[1] >> boxDimensions[2];
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // Built by the constructor, as parse_input() builds it, so that
            // xLeft and xRight are set along with the volume.  Reading x, y and
            // z alone left both at zero, and create_random_coords() places a
            // molecule created in a box at x = xLeft + (xRight - xLeft) * rand:
            // after a restart, every molecule a creation reaction made landed on
            // the plane x = 0.  nerdss_mpi overwrites both with its rank's
            // bounds in prepare.cpp.
            membraneObject.waterBox = Membrane::WaterBox(boxDimensions);

            expect_restart_key(restartFile, "membrane");
            restartFile >> membraneObject.implicitlipidIndex >> membraneObject.nSites >> membraneObject.nStates >> membraneObject.No_free_lipids >> membraneObject.No_protein >> membraneObject.totalSA;

            expect_restart_key(restartFile, "implicitLipidStates");
            for (int i = 0; i < membraneObject.nStates; i++) {
                membraneObject.numberOfFreeLipidsEachState.emplace_back(0);
                restartFile >> membraneObject.numberOfFreeLipidsEachState[i];
            }

            expect_restart_key(restartFile, "implicitLipidsParams");
            // The two flags this format carries map back one-to-one: the
            // sphere bit is the geometry, the box bit is the waterBox
            // provenance.  Nothing is folded away, so a file naming both round
            // trips as `1 1` exactly as it did before.
            //
            // NOTE: this line wants seven fields and restart files under
            // sample_inputs/ carry five, which sets failbit here and silently
            // no-ops every read below it.  That predates this change and is
            // not addressed here.
            bool isBoxFlag { false };
            bool isSphereFlag { false };
            restartFile >> membraneObject.implicitLipid >> membraneObject.TwoD >> isBoxFlag >> isSphereFlag >> membraneObject.sphereR >> membraneObject.hasCompartment >> membraneObject.compartmentR;
            membraneObject.waterBoxGiven = isBoxFlag;
            membraneObject.shape = isSphereFlag ? BoundaryShape::Sphere
                : isBoxFlag                     ? BoundaryShape::Box
                                                : BoundaryShape::Unspecified;

            // Present only for a compartment; see write_restart().  A
            // compartment file from before they were written stops here, rather
            // than restarting with the sites' D and density at zero.
            if (membraneObject.hasCompartment) {
                expect_restart_key(restartFile, "compartmentSiteD");
                restartFile >> membraneObject.droplet.D;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                expect_restart_key(restartFile, "compartmentSiteRho");
                restartFile >> membraneObject.droplet.rho;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }

            expect_restart_key(restartFile, "ifaceOverlapSepLimit");
            restartFile >> params.overlapSepLimit;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "rMaxLimit");
            restartFile >> params.rMaxLimit;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "timeWrite");
            restartFile >> params.timeWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "trajWrite");
            restartFile >> params.trajWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "restartWrite");
            restartFile >> params.restartWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "pdbWrite");
            restartFile >> params.pdbWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "accocDissocWrite");
            restartFile >> params.assocDissocWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "checkPoint");
            restartFile >> params.checkPoint;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "scaleMaxDisplace");
            restartFile >> params.scaleMaxDisplace;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "transitionWrite");
            restartFile >> params.transitionWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "clusterOverlapCheck");
            restartFile >> params.clusterOverlapCheck;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "RNGwrite");
            restartFile >> params.rngwrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            expect_restart_key(restartFile, "bondedComplexWrite");
            restartFile >> params.bondedComplexWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            unsigned long lastUpdateTransitionSize { 0 };
            restartFile >> lastUpdateTransitionSize;
            for (unsigned itr { 0 }; itr < lastUpdateTransitionSize; ++itr) {
                int index { 0 };
                restartFile >> index;
                Parameters::lastUpdateTransition.push_back(index);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        std::cout << "restart write, pdbWrite: " << params.restartWrite << ' ' << params.pdbWrite << std::endl;
        /*	std::cout<<"READ IN SUB volume PARTITIONING from restart file"<<std::endl;
        // Read Simulation Volume
        {
            restartFile >> simulVolume.numSubCells.x >> simulVolume.numSubCells.y >> simulVolume.numSubCells.z
                >> simulVolume.numSubCells.tot;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile >> simulVolume.subCellSize.x >> simulVolume.subCellSize.y >> simulVolume.subCellSize.z;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            for (unsigned subCellItr { 0 }; subCellItr < simulVolume.numSubCells.tot; ++subCellItr) {
                SimulVolume::SubVolume subCell {};
                restartFile >> subCell.absIndex >> subCell.xIndex >> subCell.yIndex >> subCell.zIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned neighborListSize { 0 };
                restartFile >> neighborListSize;
                for (unsigned itr { 0 }; itr < neighborListSize; ++itr) {
                    unsigned neighbor { 0 };
                    restartFile >> neighbor;
                    subCell.neighborList.push_back(neighbor);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned memberListSize { 0 };
                restartFile >> memberListSize;
                for (unsigned itr { 0 }; itr < memberListSize; ++itr) {
                    unsigned member { 0 };
                    restartFile >> member;
                    subCell.memberMolList.push_back(member);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                simulVolume.subCellList.push_back(subCell);
            }
        }
	*/
        std::string nextSection;
        std::getline(restartFile, nextSection);
        if (nextSection.find("#NumericalSettings") == 0) {
            const bool isVersion1 { nextSection.find("version = 1") != std::string::npos };
            const bool isVersion2 { nextSection.find("version = 2") != std::string::npos };
            const bool isVersion3 { nextSection.find("version = 3") != std::string::npos };
            if (!isVersion1 && !isVersion2 && !isVersion3)
                throw std::string("Unsupported numerical-settings restart section: ") + nextSection;

            const auto readNumericalValue = [&restartFile](double& value) {
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
                restartFile >> value;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            };
            const auto readNumericalInteger = [&restartFile](int& value) {
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
                restartFile >> value;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            };

            readNumericalValue(params.numerics.integration.tableAbsoluteError);
            readNumericalValue(params.numerics.integration.tableRelativeError);
            readNumericalValue(params.numerics.integration.fallbackError);
            readNumericalValue(params.numerics.integration.tailCutoff);
            readNumericalValue(params.numerics.integration.normalizationAbsoluteError);
            readNumericalValue(params.numerics.integration.normalizationRelativeError);
            readNumericalValue(params.numerics.tableLookup.reactionRate.absolute);
            readNumericalValue(params.numerics.tableLookup.reactionRate.relative);
            readNumericalValue(params.numerics.tableLookup.diffusionCoefficient.absolute);
            readNumericalValue(params.numerics.tableLookup.diffusionCoefficient.relative);
            readNumericalValue(params.numerics.classification.explicitLipidFlatDiffusion);
            readNumericalValue(params.numerics.classification.implicitLipidFlatDiffusion);
            if (isVersion2 || isVersion3) {
                readNumericalValue(params.numerics.associationAngles.sameAngle.absolute);
                readNumericalValue(params.numerics.associationAngles.sameAngle.relative);
                readNumericalValue(params.numerics.associationAngles.rotationConvergenceTolerance);
                readNumericalValue(params.numerics.associationAngles.endpointSignTolerance);
            }
            if (isVersion3)
                readNumericalInteger(params.numerics.vec3D.coordinateEqualityPrecision);

            std::getline(restartFile, nextSection);
            if (nextSection.find("#EndNumericalSettings") != 0)
                throw std::string("Expected #EndNumericalSettings, found: ") + nextSection;
            try {
                params.numerics.validate();
            } catch (const std::invalid_argument& error) {
                throw std::string("Invalid numerical settings in restart file: ") + error.what();
            }
            std::getline(restartFile, nextSection); // #MolTemplates
        }
        if (nextSection.find("#MolTemplates") != 0)
            throw std::string("Expected #MolTemplates section, found: ") + nextSection;

        std::cout << "READ IN MOL TEMPLATE from restart file" << std::endl;
        // Read MolTemplates
        {
            restartFile >> MolTemplate::numMolTypes;
            for (unsigned itr { 0 }; itr < MolTemplate::numMolTypes; ++itr) {
                unsigned num { 0 };
                restartFile >> num;
                MolTemplate::numEachMolType.push_back(num);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << " Num moltypes: " << MolTemplate::numMolTypes << '\n';
            unsigned absToRelIfaceSize { 0 };
            restartFile >> absToRelIfaceSize;
            for (unsigned itr { 0 }; itr < absToRelIfaceSize; ++itr) {
                unsigned iface { 0 };
                restartFile >> iface;
                MolTemplate::absToRelIface.push_back(iface);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile >> Interface::State::totalNumOfStates;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (unsigned itr { 0 }; itr < MolTemplate::numMolTypes; ++itr) {
                MolTemplate oneTemp {};
                restartFile >> oneTemp.molTypeIndex >> oneTemp.molName;
                std::cout << " protein index, name: " << oneTemp.molTypeIndex << ' ' << oneTemp.molName << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.copies >> oneTemp.mass >> oneTemp.radius;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.isLipid >> oneTemp.isImplicitLipid >> oneTemp.isRod >> oneTemp.isPoint >> oneTemp.checkOverlap >> oneTemp.countTransition >> oneTemp.transitionMatrixSize 
                >> oneTemp.outsideCompartment >> oneTemp.insideCompartment >> oneTemp.crossesCompartment >> oneTemp.transmissionRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.comCoord.x >> oneTemp.comCoord.y >> oneTemp.comCoord.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.D.x >> oneTemp.D.y >> oneTemp.D.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.Dr.x >> oneTemp.Dr.y >> oneTemp.Dr.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                // Without this invCbrtDr stays zero, Complex::update_properties()
                // sums it to zero and divides by that, and every complex with a
                // rotating member comes out with Dr = inf and NaN coordinates
                // after its first step.
                oneTemp.cache_diffusion_derivatives();

                // reaction partners
                {
                    unsigned partSize { 0 };
                    restartFile >> partSize;

                    for (unsigned partItr { 0 }; partItr < partSize; ++partItr) {
                        int partner { 0 };
                        restartFile >> partner;
                        oneTemp.rxnPartners.push_back(partner);
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // optional bonds
                {
                    unsigned bondSize { 0 };
                    restartFile >> bondSize;

                    for (unsigned bondItr { 0 }; bondItr < bondSize; ++bondItr) {
                        int iface1 { 0 };
                        int iface2 { 0 };
                        restartFile >> iface1 >> iface2;
                        oneTemp.bondList.push_back(std::array<int, 2> { { iface1, iface2 } });
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // write interfaces
                unsigned oneTempIfaceSize { 0 };
                restartFile >> oneTempIfaceSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                for (unsigned ifaceItr { 0 }; ifaceItr < oneTempIfaceSize; ++ifaceItr) {
                    Interface tmpIface {};
                    restartFile >> tmpIface.index >> tmpIface.name;
                    std::cout << " iface index, name: " << tmpIface.index << ' ' << tmpIface.name << '\n';
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> tmpIface.iCoord.x >> tmpIface.iCoord.y >> tmpIface.iCoord.z;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned stateListSize { 0 };
                    restartFile >> stateListSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned stateItr { 0 }; stateItr < stateListSize; ++stateItr) {
                        Interface::State tmpState {};
                        restartFile >> tmpState.index >> tmpState.iden;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        // partner list
                        unsigned rxnPartnersSize { 0 };
                        restartFile >> rxnPartnersSize;
                        for (unsigned partItr { 0 }; partItr < rxnPartnersSize; ++partItr) {
                            unsigned partner { 0 };
                            restartFile >> partner;
                            tmpState.rxnPartners.push_back(partner);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        // reaction lists
                        unsigned myForwardRxnsSize { 0 };
                        restartFile >> myForwardRxnsSize;
                        for (unsigned rxnItr { 0 }; rxnItr < myForwardRxnsSize; ++rxnItr) {
                            int rxn { 0 };
                            restartFile >> rxn;
                            tmpState.myForwardRxns.push_back(rxn);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        unsigned myCreateDestructRxnsSize { 0 };
                        restartFile >> myCreateDestructRxnsSize;
                        for (unsigned rxnItr { 0 }; rxnItr < myCreateDestructRxnsSize; ++rxnItr) {
                            int rxn { 0 };
                            restartFile >> rxn;
                            tmpState.myCreateDestructRxns.push_back(rxn);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        unsigned stateChangeRxnsSize { 0 };
                        restartFile >> stateChangeRxnsSize;
                        for (unsigned rxnItr { 0 }; rxnItr < stateChangeRxnsSize; ++rxnItr) {
                            unsigned elem1 { 0 };
                            unsigned elem2 { 0 };
                            restartFile >> elem1 >> elem2;
                            tmpState.stateChangeRxns.emplace_back(elem1, elem2);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        // done reading state, add to state list
                        tmpIface.stateList.emplace_back(tmpState);
                    }
                    // done reading interface, add to interface list
                    oneTemp.interfaceList.emplace_back(tmpIface);
                }

                // read ifaceWithStates
                unsigned ifacesWithStatesSize { 0 };
                restartFile >> ifacesWithStatesSize;
                for (unsigned rxnItr { 0 }; rxnItr < ifacesWithStatesSize; ++rxnItr) {
                    unsigned elem { 0 };
                    restartFile >> elem;
                    oneTemp.ifacesWithStates.push_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // read monomerList
                unsigned monomerListSize { 0 };
                restartFile >> monomerListSize;
                for (unsigned rxnItr { 0 }; rxnItr < monomerListSize; ++rxnItr) {
                    unsigned elem { 0 };
                    restartFile >> elem;
                    oneTemp.monomerList.push_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // read lifetime
                if(oneTemp.countTransition == true) {
                    unsigned lifeTimeSize {0};
                    oneTemp.lifeTime.resize(oneTemp.transitionMatrixSize);
                    for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne) {
                        restartFile >> lifeTimeSize;
                        for (unsigned rxnItr { 0 }; rxnItr < lifeTimeSize; ++rxnItr) {
                            double elem { 0.0 };
                            restartFile >> elem;
                            oneTemp.lifeTime[indexOne].push_back(elem);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // read transition matrix
                if(oneTemp.countTransition == true) {
                    oneTemp.transitionMatrix.resize(oneTemp.transitionMatrixSize);
                    for (int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne) {
                        oneTemp.transitionMatrix[indexOne].resize(oneTemp.transitionMatrixSize);
                    }

                    for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne){
                        for (int indexTwo = 0; indexTwo < oneTemp.transitionMatrixSize; ++indexTwo){
                            restartFile >> oneTemp.transitionMatrix[indexOne][indexTwo];
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // done reading template, add to template list
                molTemplateList.emplace_back(oneTemp);
            }
        }
        std::cout << "READ IN REACTIONs from restart file" << std::endl;
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // write Reactions
        {
            unsigned forwardRxnsSize { 0 };
            unsigned backRxnsSize { 0 };
            unsigned createDestructRxnsSize { 0 };
            unsigned transmissionRxnSize { 0 };
            restartFile >> RxnBase::numberOfRxns >> forwardRxnsSize >> backRxnsSize >> createDestructRxnsSize >> transmissionRxnSize
                >> RxnBase::totRxnSpecies;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << " Num rxns: " << RxnBase::numberOfRxns << " forwardRxns: " << forwardRxnsSize << '\n';
            // forward reactions
            for (unsigned rxnItr { 0 }; rxnItr < forwardRxnsSize; ++rxnItr) {
                ForwardRxn tmpRxn;
                restartFile >> tmpRxn.absRxnIndex >> tmpRxn.relRxnIndex >> tmpRxn.rxnLabel;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::cout << " itr: " << rxnItr << " abs index, relindex: " << tmpRxn.absRxnIndex << ' ' << tmpRxn.relRxnIndex << '\n';
                int rxnType { -1 };
                restartFile >> rxnType >> tmpRxn.isSymmetric >> tmpRxn.isOnMem >> tmpRxn.hasStateChange;
                std::cout << "Rxntype: " << rxnType << '\n';
                if (rxnType != -1) {
                    tmpRxn.rxnType = static_cast<ReactionType>(rxnType); // turn the int rxnType into ReactionType
                } else {
                    std::cerr << "ERROR: Cannot parse reaction type for reaction " << tmpRxn.absRxnIndex
                              << ". Exiting.\n";
                    exit(1);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isObserved;
                if (tmpRxn.isObserved)
                    restartFile >> tmpRxn.observeLabel;
                restartFile >> tmpRxn.productName;
                std::cout << "Product name: " << tmpRxn.productName << std::endl;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isReversible >> tmpRxn.conjBackRxnIndex >>
                    tmpRxn.irrevRingClosure >> tmpRxn.bindRadSameCom >>
                    tmpRxn.loopCoopFactor >> tmpRxn.length3Dto2D >> tmpRxn.area3Dto1D;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.bindRadius;
                std::string th1, th2, ph1, ph2, omega;
                restartFile >> th1 >> th2 >> ph1 >> ph2 >> omega;

                tmpRxn.assocAngles.theta1 = std::stod(th1);
                tmpRxn.assocAngles.theta2 = std::stod(th2);
                tmpRxn.assocAngles.phi1 = std::stod(ph1);
                tmpRxn.assocAngles.phi2 = std::stod(ph2);
                tmpRxn.assocAngles.omega = std::stod(omega);
                if (th1 == "nan")
                    tmpRxn.assocAngles.theta1 = std::numeric_limits<double>::quiet_NaN();
                if (th2 == "nan")
                    tmpRxn.assocAngles.theta2 = std::numeric_limits<double>::quiet_NaN();
                if (ph1 == "nan")
                    tmpRxn.assocAngles.phi1 = std::numeric_limits<double>::quiet_NaN();
                if (ph2 == "nan")
                    tmpRxn.assocAngles.phi2 = std::numeric_limits<double>::quiet_NaN();
                if (omega == "nan")
                    tmpRxn.assocAngles.omega = std::numeric_limits<double>::quiet_NaN();
                //>> tmpRxn.assocAngles.theta1 >> tmpRxn.assocAngles.theta2
                //>> tmpRxn.assocAngles.phi1 >> tmpRxn.assocAngles.phi2 >> tmpRxn.assocAngles.omega;
                std::cout << "RXN angles " << tmpRxn.bindRadius << ' ' << tmpRxn.assocAngles.theta1 << ' ' << tmpRxn.assocAngles.theta2
                          << ' ' << tmpRxn.assocAngles.phi1 << ' ' << tmpRxn.assocAngles.phi2 << ' ' << tmpRxn.assocAngles.omega << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.norm1.x >> tmpRxn.norm1.y >> tmpRxn.norm1.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.norm2.x >> tmpRxn.norm2.y >> tmpRxn.norm2.z;
                std::cout << " norm 2: " << tmpRxn.norm2.x << ' ' << tmpRxn.norm2.y << ' ' << tmpRxn.norm2.z << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.excludeVolumeBound;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isCoupled;
                std::cout << " reaction is coupled? " << tmpRxn.isCoupled << std::endl;
                if (tmpRxn.isCoupled) {
                    rxnType = -1;
                    std::cout << "did not enter iscoupled loop " << '\n';
                    restartFile >> tmpRxn.coupledRxn.absRxnIndex >> tmpRxn.coupledRxn.relRxnIndex >> rxnType >> tmpRxn.coupledRxn.label >> tmpRxn.coupledRxn.probCoupled;
                    if (rxnType != -1) {
                        tmpRxn.coupledRxn.rxnType = static_cast<ReactionType>(rxnType);
                    } else {
                        std::cerr << "ERROR: Cannot parse reaction type for reaction coupled to reaction "
                                  << tmpRxn.absRxnIndex << ". Exiting.\n";
                        exit(1);
                    }
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer reactants
                unsigned intReactantListSize { 0 };
                restartFile >> intReactantListSize;
                std::cout << "Nreactant first round " << intReactantListSize << '\n';
                for (unsigned itr { 0 }; itr < intReactantListSize; ++itr) {
                    int reactant { -1 };
                    restartFile >> reactant;
                    tmpRxn.intReactantList.push_back(reactant);
                    std::cout << " reactant: " << reactant << '\n';
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer products
                unsigned intProductListSize { 0 };
                restartFile >> intProductListSize;
                for (unsigned itr { 0 }; itr < intProductListSize; ++itr) {
                    int product { -1 };
                    restartFile >> product;
                    tmpRxn.intProductList.push_back(product);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantListNewSize { 0 };
                restartFile >> reactantListNewSize;
                std::cout << "N reactants: " << reactantListNewSize << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantListNewSize; ++itr) {
                    RxnIface oneReact {};
                    restartFile >> oneReact.molTypeIndex;
                    std::cout << " molTypeIndex: " << oneReact.molTypeIndex << '\n';
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.ifaceName >> oneReact.absIfaceIndex >> oneReact.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.requiresState >> oneReact.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.reactantListNew.emplace_back(oneReact);
                    std::cout << " requiresState, requiresInteraction: " << oneReact.requiresState << ' ' << oneReact.requiresInteraction << '\n';
                }

                // product list
                unsigned productListNewSize { 0 };
                restartFile >> productListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productListNewSize; ++itr) {
                    RxnIface oneProd {};
                    restartFile >> oneProd.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.ifaceName >> oneProd.absIfaceIndex >> oneProd.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.requiresState >> oneProd.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.productListNew.emplace_back(oneProd);
                }

                // rate list
                unsigned rateListSize { 0 };
                restartFile >> rateListSize;
                std::cout << "Nrates: " << rateListSize << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < rateListSize; ++itr) {
                    RxnBase::RateState oneRate {};
                    restartFile >> oneRate.rate;
                    std::cout << " rate: " << oneRate.rate << '\n';
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned otherIfaceListsSize { 0 };
                    unsigned ifaceItr { 0 };
                    restartFile >> otherIfaceListsSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned listItr { 0 }; listItr < otherIfaceListsSize; ++listItr) {
                        unsigned oneListSize { 0 };
                        std::vector<RxnIface> tmpIfaceVec {};
                        restartFile >> oneListSize;
                        std::cout << "onelistsize: " << oneListSize << '\n';
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        for (unsigned anccIfaceItr { 0 }; anccIfaceItr < oneListSize; ++anccIfaceItr) {
                            RxnIface otherIface {};
                            restartFile >> otherIface.molTypeIndex >> otherIface.ifaceName >> otherIface.absIfaceIndex
                                >> otherIface.relIfaceIndex >> otherIface.requiresState
                                >> otherIface.requiresInteraction;
                            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            tmpIfaceVec.emplace_back(otherIface);
                        }
                        oneRate.otherIfaceLists.push_back(tmpIfaceVec);
                    }
                    tmpRxn.rateList.emplace_back(oneRate);
                }
                tmpRxn.display();
                forwardRxns.emplace_back(tmpRxn);
            }
            std::cout << " Done with forward reactions " << '\n';
            // backRxns
            for (unsigned rxnItr { 0 }; rxnItr < backRxnsSize; ++rxnItr) {
                BackRxn tmpRxn;
                restartFile >> tmpRxn.absRxnIndex >> tmpRxn.relRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                int rxnType { -1 };
                restartFile >> rxnType >> tmpRxn.isSymmetric >> tmpRxn.isOnMem >> tmpRxn.hasStateChange;
                tmpRxn.rxnType = static_cast<ReactionType>(rxnType); // turn the int rxnType into ReactionType
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.isObserved;
                if (tmpRxn.isObserved)
                    restartFile >> tmpRxn.observeLabel;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.conjForwardRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.isCoupled;
                std::cout << " reaction is coupled? " << tmpRxn.isCoupled << std::endl;
                if (tmpRxn.isCoupled) {
                    rxnType = -1;
                    std::cout << "did not enter iscoupled loop " << '\n';
                    restartFile >> tmpRxn.coupledRxn.absRxnIndex >> tmpRxn.coupledRxn.relRxnIndex >> rxnType >> tmpRxn.coupledRxn.label >> tmpRxn.coupledRxn.probCoupled;
                    if (rxnType != -1) {
                        tmpRxn.coupledRxn.rxnType = static_cast<ReactionType>(rxnType);
                    } else {
                        std::cerr << "ERROR: Cannot parse reaction type for reaction coupled to reaction "
                                  << tmpRxn.absRxnIndex << ". Exiting.\n";
                        exit(1);
                    }
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer reactants
                unsigned intReactantListSize { 0 };
                restartFile >> intReactantListSize;
                for (unsigned itr { 0 }; itr < intReactantListSize; ++itr) {
                    int reactant { -1 };
                    restartFile >> reactant;
                    tmpRxn.intReactantList.push_back(reactant);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer products
                unsigned intProductListSize { 0 };
                restartFile >> intProductListSize;
                for (unsigned itr { 0 }; itr < intProductListSize; ++itr) {
                    int product { -1 };
                    restartFile >> product;
                    tmpRxn.intProductList.push_back(product);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantListNewSize { 0 };
                restartFile >> reactantListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantListNewSize; ++itr) {
                    RxnIface oneReact {};
                    restartFile >> oneReact.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.ifaceName >> oneReact.absIfaceIndex >> oneReact.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.requiresState >> oneReact.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.reactantListNew.emplace_back(oneReact);
                }

                // product list
                unsigned productListNewSize { 0 };
                restartFile >> productListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productListNewSize; ++itr) {
                    RxnIface oneProd {};
                    restartFile >> oneProd.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.ifaceName >> oneProd.absIfaceIndex >> oneProd.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.requiresState >> oneProd.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.productListNew.emplace_back(oneProd);
                }

                // rate list
                unsigned rateListSize { 0 };
                restartFile >> rateListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < rateListSize; ++itr) {
                    RxnBase::RateState oneRate {};
                    restartFile >> oneRate.rate;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned otherIfaceListSize { 0 };
                    restartFile >> otherIfaceListSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned listItr { 0 }; listItr < otherIfaceListSize; ++listItr) {
                        unsigned oneListSize { 0 };
                        restartFile >> oneListSize;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        std::vector<RxnIface> tmpIfaceVec;
                        for (unsigned anccIfaceItr { 0 }; anccIfaceItr < oneListSize; ++anccIfaceItr) {
                            RxnIface otherIface {};
                            restartFile >> otherIface.molTypeIndex >> otherIface.ifaceName >> otherIface.absIfaceIndex
                                >> otherIface.relIfaceIndex >> otherIface.requiresState
                                >> otherIface.requiresInteraction;
                            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            tmpIfaceVec.emplace_back(otherIface);
                        }
                        oneRate.otherIfaceLists.push_back(tmpIfaceVec);
                    }
                    tmpRxn.rateList.emplace_back(oneRate);
                }
                backRxns.emplace_back(tmpRxn);
            }
            std::cout << "Done with back reactions " << '\n';
            // creation and destruction reactions
            std::cout << "Now creation and destruction " << '\n';
            for (unsigned rxnItr { 0 }; rxnItr < createDestructRxnsSize; ++rxnItr) {
                CreateDestructRxn tmpRxn {};
                restartFile >> tmpRxn.absRxnIndex >> tmpRxn.relRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                int rxnType { -1 };
                restartFile >> rxnType >> tmpRxn.isOnMem;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isObserved; //>>tmpRxn.observeLabel;
                if (tmpRxn.isObserved)
                    restartFile >> tmpRxn.observeLabel;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                tmpRxn.rxnType = static_cast<ReactionType>(rxnType);
                restartFile >> tmpRxn.creationRadius;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer reactants
                unsigned intReactantListSize { 0 };
                restartFile >> intReactantListSize;
                for (unsigned itr { 0 }; itr < intReactantListSize; ++itr) {
                    int reactant { -1 };
                    restartFile >> reactant;
                    tmpRxn.intReactantList.push_back(reactant);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer products
                unsigned intProductListSize { 0 };
                restartFile >> intProductListSize;
                for (unsigned itr { 0 }; itr < intProductListSize; ++itr) {
                    int product { -1 };
                    restartFile >> product;
                    tmpRxn.intProductList.push_back(product);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantMolListSize { 0 };
                restartFile >> reactantMolListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantMolListSize; ++itr) {
                    CreateDestructRxn::CreateDestructMol oneMol {};
                    unsigned interfaceListSize { 0 };
                    restartFile >> oneMol.molTypeIndex >> oneMol.molName >> interfaceListSize;
                    std::cout << " Destroy molecule: " << oneMol.molName << std::endl;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned ifaceItr { 0 }; ifaceItr < interfaceListSize; ++ifaceItr) {
                        RxnIface tmpIface {};
                        restartFile >> tmpIface.molTypeIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.ifaceName >> tmpIface.absIfaceIndex >> tmpIface.relIfaceIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        oneMol.interfaceList.emplace_back(tmpIface);
                    }
                    tmpRxn.reactantMolList.emplace_back(oneMol);
                }

                // product list
                unsigned productMolListSize { 0 };
                restartFile >> productMolListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productMolListSize; ++itr) {
                    CreateDestructRxn::CreateDestructMol oneMol {};
                    unsigned interfaceListSize { 0 };
                    restartFile >> oneMol.molTypeIndex >> oneMol.molName >> interfaceListSize;
                    std::cout << " CREATION OF MOL: " << oneMol.molName << std::endl;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned ifaceItr { 0 }; ifaceItr < interfaceListSize; ++ifaceItr) {
                        RxnIface tmpIface {};
                        restartFile >> tmpIface.molTypeIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.ifaceName >> tmpIface.absIfaceIndex >> tmpIface.relIfaceIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        oneMol.interfaceList.emplace_back(tmpIface);
                    }
                    tmpRxn.productMolList.emplace_back(oneMol);
                }

                unsigned rateListSize { 0 };
                restartFile >> rateListSize;
                for (unsigned itr { 0 }; itr < rateListSize; ++itr) {
                    RxnBase::RateState tmpRate {};
                    restartFile >> tmpRate.rate;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned otherIfaceListsSize { 0 };
                    restartFile >> otherIfaceListsSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    // for (unsigned ifaceItr { 0 }; ifaceItr < otherIfaceListSize; ++ifaceItr) {
                    //     RxnIface tmpIface {};
                    //     restartFile >> tmpIface.molTypeIndex >> tmpIface.ifaceName >> tmpIface.absIfaceIndex
                    //         >> tmpIface.relIfaceIndex >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                    //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    //     tmpRate.otherIfaceLists.emplace_back(std::vector<RxnIface> { tmpIface });
                    // }

                    // restartFile >> otherIfaceListsSize;
                    // restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned listItr { 0 }; listItr < otherIfaceListsSize; ++listItr) {
                        unsigned oneListSize { 0 };
                        std::vector<RxnIface> tmpIfaceVec {};
                        restartFile >> oneListSize;
                        std::cout << "onelistsize: " << oneListSize << '\n';
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        for (unsigned anccIfaceItr { 0 }; anccIfaceItr < oneListSize; ++anccIfaceItr) {
                            RxnIface otherIface {};
                            restartFile >> otherIface.molTypeIndex >> otherIface.ifaceName >> otherIface.absIfaceIndex
                                >> otherIface.relIfaceIndex >> otherIface.requiresState
                                >> otherIface.requiresInteraction;
                            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            tmpIfaceVec.emplace_back(otherIface);
                        }
                        tmpRate.otherIfaceLists.push_back(tmpIfaceVec);
                    }
                    //copied new version up to here.
                    tmpRxn.rateList.emplace_back(tmpRate);
                }
                createDestructRxns.emplace_back(tmpRxn);
            }

            // tranmission reaction 
            std::cout << "Now transmission reaction " << '\n';
            for (unsigned rxnItr { 0 }; rxnItr < transmissionRxnSize; ++rxnItr) {
                TransmissionRxn tmpRxn {};
                restartFile >> tmpRxn.absRxnIndex >> tmpRxn.relRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                int rxnType { -1 };
                restartFile >> rxnType >> tmpRxn.isOnMem;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isObserved; //>>tmpRxn.observeLabel;
                if (tmpRxn.isObserved)
                    restartFile >> tmpRxn.observeLabel;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                tmpRxn.rxnType = static_cast<ReactionType>(rxnType);

                // integer reactants
                unsigned intReactantListSize { 0 };
                restartFile >> intReactantListSize;
                for (unsigned itr { 0 }; itr < intReactantListSize; ++itr) {
                    int reactant { -1 };
                    restartFile >> reactant;
                    tmpRxn.intReactantList.push_back(reactant);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer products
                unsigned intProductListSize { 0 };
                restartFile >> intProductListSize;
                for (unsigned itr { 0 }; itr < intProductListSize; ++itr) {
                    int product { -1 };
                    restartFile >> product;
                    tmpRxn.intProductList.push_back(product);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantMolListSize { 0 };
                restartFile >> reactantMolListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantMolListSize; ++itr) {
                    TransmissionRxn::TransmissionMol oneMol {};
                    unsigned interfaceListSize { 0 };
                    restartFile >> oneMol.molTypeIndex >> oneMol.molName >> interfaceListSize;
                    std::cout << " Destroy molecule: " << oneMol.molName << std::endl;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned ifaceItr { 0 }; ifaceItr < interfaceListSize; ++ifaceItr) {
                        RxnIface tmpIface {};
                        restartFile >> tmpIface.molTypeIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.ifaceName >> tmpIface.absIfaceIndex >> tmpIface.relIfaceIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        oneMol.interfaceList.emplace_back(tmpIface);
                    }
                    tmpRxn.reactantMolList.emplace_back(oneMol);
                }

                // product list
                unsigned productMolListSize { 0 };
                restartFile >> productMolListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productMolListSize; ++itr) {
                    TransmissionRxn::TransmissionMol oneMol {};
                    unsigned interfaceListSize { 0 };
                    restartFile >> oneMol.molTypeIndex >> oneMol.molName >> interfaceListSize;
                    std::cout << " CREATION OF MOL: " << oneMol.molName << std::endl;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned ifaceItr { 0 }; ifaceItr < interfaceListSize; ++ifaceItr) {
                        RxnIface tmpIface {};
                        restartFile >> tmpIface.molTypeIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.ifaceName >> tmpIface.absIfaceIndex >> tmpIface.relIfaceIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        oneMol.interfaceList.emplace_back(tmpIface);
                    }
                    tmpRxn.productMolList.emplace_back(oneMol);
                }

                unsigned rateListSize { 0 };
                restartFile >> rateListSize;
                for (unsigned itr { 0 }; itr < rateListSize; ++itr) {
                    RxnBase::RateState tmpRate {};
                    restartFile >> tmpRate.rate;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned otherIfaceListsSize { 0 };
                    restartFile >> otherIfaceListsSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    // for (unsigned ifaceItr { 0 }; ifaceItr < otherIfaceListSize; ++ifaceItr) {
                    //     RxnIface tmpIface {};
                    //     restartFile >> tmpIface.molTypeIndex >> tmpIface.ifaceName >> tmpIface.absIfaceIndex
                    //         >> tmpIface.relIfaceIndex >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                    //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    //     tmpRate.otherIfaceLists.emplace_back(std::vector<RxnIface> { tmpIface });
                    // }

                    // restartFile >> otherIfaceListsSize;
                    // restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned listItr { 0 }; listItr < otherIfaceListsSize; ++listItr) {
                        unsigned oneListSize { 0 };
                        std::vector<RxnIface> tmpIfaceVec {};
                        restartFile >> oneListSize;
                        std::cout << "onelistsize: " << oneListSize << '\n';
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        for (unsigned anccIfaceItr { 0 }; anccIfaceItr < oneListSize; ++anccIfaceItr) {
                            RxnIface otherIface {};
                            restartFile >> otherIface.molTypeIndex >> otherIface.ifaceName >> otherIface.absIfaceIndex
                                >> otherIface.relIfaceIndex >> otherIface.requiresState
                                >> otherIface.requiresInteraction;
                            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            tmpIfaceVec.emplace_back(otherIface);
                        }
                        tmpRate.otherIfaceLists.push_back(tmpIfaceVec);
                    }
                    //copied new version up to here.
                    tmpRxn.rateList.emplace_back(tmpRate);
                }

                // Appended to the record; see write_restart().  The tag is what
                // refuses a file written before these fields existed.
                expect_restart_key(restartFile, "bindRadius");
                restartFile >> tmpRxn.bindRadius;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantListNewSize { 0 };
                restartFile >> reactantListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantListNewSize; ++itr) {
                    RxnIface oneReact {};
                    restartFile >> oneReact.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.ifaceName >> oneReact.absIfaceIndex >> oneReact.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.requiresState >> oneReact.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.reactantListNew.emplace_back(oneReact);
                }

                // product list
                unsigned productListNewSize { 0 };
                restartFile >> productListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productListNewSize; ++itr) {
                    RxnIface oneProd {};
                    restartFile >> oneProd.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.ifaceName >> oneProd.absIfaceIndex >> oneProd.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.requiresState >> oneProd.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.productListNew.emplace_back(oneProd);
                }
                transmissionRxns.emplace_back(tmpRxn);
            }
        }

        std::cout << "Now read in coordinates " << std::endl;
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // write Molecules
        {
            int molListSize { 0 };

            restartFile >> molListSize >> Molecule::numberOfMolecules;
            std::cout << "Mol list size and molecule.numberofMolecules: " << molListSize << ' ' << Molecule::numberOfMolecules << std::endl;
            for (unsigned molItr { 0 }; molItr < molListSize; ++molItr) {
                Molecule tmpMol {};
                restartFile >> tmpMol.index >> tmpMol.isEmpty >> tmpMol.myComIndex >> tmpMol.molTypeIndex
                    >> tmpMol.mySubVolIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpMol.mass >> tmpMol.isLipid >> tmpMol.isImplicitLipid 
                            >> tmpMol.linksToSurface >> tmpMol.isPromoter >> tmpMol.isEmpty;
                if (membraneObject.hasCompartment) // see write_restart()
                    restartFile >> tmpMol.enforceCompartmentBC;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                // center of mass
                restartFile >> std::fixed >> tmpMol.comCoord.x >> tmpMol.comCoord.y >> tmpMol.comCoord.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // interface lists
                unsigned freeListSize { 0 };
                restartFile >> freeListSize;
                for (unsigned itr { 0 }; itr < freeListSize; ++itr) {
                    int freeIface { 0 };
                    restartFile >> freeIface;
                    tmpMol.freelist.emplace_back(freeIface);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned bndListSize { 0 };
                restartFile >> bndListSize;
                for (unsigned itr { 0 }; itr < bndListSize; ++itr) {
                    int bndIface { 0 };
                    restartFile >> bndIface;
                    tmpMol.bndlist.emplace_back(bndIface);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned bndPartnerListSize { 0 };
                restartFile >> bndPartnerListSize;
                for (unsigned itr { 0 }; itr < bndPartnerListSize; ++itr) {
                    int bndPartner { 0 };
                    restartFile >> bndPartner;
                    tmpMol.bndpartner.emplace_back(bndPartner);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // interfaces
                unsigned interfaceListSize { 0 };
                restartFile >> interfaceListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                for (unsigned itr { 0 }; itr < interfaceListSize; ++itr) {
                    Molecule::Iface tmpIface {};
                    restartFile >> tmpIface.index >> tmpIface.relIndex >> tmpIface.molTypeIndex >> tmpIface.stateIndex
                        >> tmpIface.stateIden >> tmpIface.isBound;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> std::fixed >> tmpIface.coord.x >> tmpIface.coord.y >> tmpIface.coord.z;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    if (tmpIface.isBound) {
                        restartFile >> tmpIface.interaction.partnerIndex >> tmpIface.interaction.partnerIfaceIndex
                            >> tmpIface.interaction.conjBackRxn;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    tmpMol.interfaceList.emplace_back(tmpIface);
                }

                // Reweighting lists.  The file carries six parallel arrays, as
                // it always has; they are read in the same order and assembled
                // into the single prevReweight vector.  Each array carries its
                // own count, and a file written by any build gives all six the
                // same count, so the first one sizes the vector and the rest
                // fill fields of entries that already exist.
                unsigned listSize { 0 };
                restartFile >> listSize;
                tmpMol.prevReweight.resize(listSize);
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int elem { 0 };
                    restartFile >> elem;
                    tmpMol.prevReweight[itr].partner = elem;
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int elem { 0 };
                    restartFile >> elem;
                    // Always consume, store only what fits: a malformed file
                    // whose arrays disagree in length must not leave numbers in
                    // the stream for the next record to misread.
                    if (itr < tmpMol.prevReweight.size())
                        tmpMol.prevReweight[itr].myFace = elem;
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int elem { 0 };
                    restartFile >> elem;
                    // Always consume, store only what fits: a malformed file
                    // whose arrays disagree in length must not leave numbers in
                    // the stream for the next record to misread.
                    if (itr < tmpMol.prevReweight.size())
                        tmpMol.prevReweight[itr].partnerFace = elem;
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    double elem { 0 };
                    restartFile >> elem;
                    // Always consume, store only what fits: a malformed file
                    // whose arrays disagree in length must not leave numbers in
                    // the stream for the next record to misread.
                    if (itr < tmpMol.prevReweight.size())
                        tmpMol.prevReweight[itr].norm = elem;
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    double elem { 0 };
                    restartFile >> elem;
                    // Always consume, store only what fits: a malformed file
                    // whose arrays disagree in length must not leave numbers in
                    // the stream for the next record to misread.
                    if (itr < tmpMol.prevReweight.size())
                        tmpMol.prevReweight[itr].survProb = elem;
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    double elem { 0 };
                    restartFile >> elem;
                    // Always consume, store only what fits: a malformed file
                    // whose arrays disagree in length must not leave numbers in
                    // the stream for the next record to misread.
                    if (itr < tmpMol.prevReweight.size())
                        tmpMol.prevReweight[itr].sep = elem;
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // trajStatus is not in the file, and every other molecule can do
                // without it: the end of each timestep resets theirs to `none`,
                // which is what they all hold at a checkpoint.  That reset skips
                // the implicit lipid's one representative molecule.  The overlap
                // loop propagates it in the first timestep, and from then on it
                // holds `propagated`, so the loop never propagates it again.
                // Left at `none`, a restart propagated it once more at its first
                // step, drawing random numbers the uninterrupted run never drew,
                // and every step after that diverged.  A file written at step 0
                // is from before that first propagation, where `none` is right.
                if (tmpMol.isImplicitLipid && simItr > 0)
                    tmpMol.trajStatus = TrajStatus::propagated;
                moleculeList.emplace_back(tmpMol);
                //std::cout <<"read in : "<<tmpMol.index<<" first interface z crd: "<<tmpMol.comCoord.z<<std::endl;
            }

            unsigned long emptyMolListSize { 0 };
            restartFile >> emptyMolListSize;
            for (unsigned itr { 0 }; itr < emptyMolListSize; ++itr) {
                int index { 0 };
                restartFile >> index;
                Molecule::emptyMolList.push_back(index);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "N empty molecules: " << emptyMolListSize << std::endl;
        }
        std::cout << "Now read in complexes from RESTART " << '\n';
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // read Complexes
        {
            int comListSize { 0 };
            restartFile >> comListSize >> Complex::numberOfComplexes;
            std::cout << " Ncomplexes including empties: " << comListSize << " N actual complexes: " << Complex::numberOfComplexes << std::endl;
            for (unsigned comItr { 0 }; comItr < comListSize; ++comItr) {
                Complex tmpCom {};
                restartFile >> tmpCom.index >> tmpCom.isEmpty >> tmpCom.radius >> tmpCom.mass;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.linksToSurface >> tmpCom.iLipidIndex >> tmpCom.OnSurface;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.onFiber;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.comCoord.x >> tmpCom.comCoord.y >> tmpCom.comCoord.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.D.x >> tmpCom.D.y >> tmpCom.D.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.Dr.x >> tmpCom.Dr.y >> tmpCom.Dr.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // member molecule lists
                unsigned listSize { 0 };
                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int memMol { -1 };
                    restartFile >> memMol;
                    tmpCom.memberList.emplace_back(memMol);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // numEachMol list
                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int memMol { -1 };
                    restartFile >> memMol;
                    tmpCom.numEachMol.emplace_back(memMol);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int memMol { -1 };
                    restartFile >> memMol;
                    tmpCom.lastNumberUpdateItrEachMol.emplace_back(memMol);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                complexList.emplace_back(tmpCom);
            }

            unsigned long emptyComListSize { 0 };
            restartFile >> emptyComListSize;
            std::cout << "N empty complexes " << emptyComListSize << '\t';
            for (unsigned itr { 0 }; itr < emptyComListSize; ++itr) {
                int index { 0 };
                restartFile >> index;
                Complex::emptyComList.push_back(index);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } //done reading complexes
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // read observables
        {
            int numObs { 0 };
            restartFile >> numObs;
            std::cout << "N observables " << numObs << '\t';
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (unsigned itr { 0 }; itr < numObs; ++itr) {
                std::string obsName;
                int obsCount;
                restartFile >> obsName >> obsCount;
                observablesList.emplace(obsName, obsCount);
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // read counterArrays
        {
            restartFile >> counterArrays.nLoops >> counterArrays.nCancelOverlapPartner >> counterArrays.nCancelOverlapSystem >> counterArrays.nCancelDisplace2D >> counterArrays.nCancelDisplace3D >> counterArrays.nCancelDisplace3Dto2D >> counterArrays.nCancelSpanBox >> counterArrays.nAssocSuccess >> counterArrays.eventArraySize;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //read events3D, 3Dto2D, and 2D
            for (int i = 0; i < counterArrays.eventArraySize; i++)
                restartFile >> counterArrays.events3D[i];
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (int i = 0; i < counterArrays.eventArraySize; i++)
                restartFile >> counterArrays.events3Dto2D[i];
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (int i = 0; i < counterArrays.eventArraySize; i++)
                restartFile >> counterArrays.events2D[i];
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            int numSpecies { 0 };
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            restartFile >> numSpecies;
            std::cout << "N species " << numSpecies << '\t';
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (unsigned itr { 0 }; itr < numSpecies; ++itr) {
                unsigned long listSize { 0 };
                restartFile >> listSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                counterArrays.bindPairList.emplace_back();
                std::cout << "Specie " << itr << " N bindPairs " << listSize << '\t';
                for (unsigned itr2 { 0 }; itr2 < listSize; ++itr2) {
                    int index { 0 };
                    restartFile >> index;
                    counterArrays.bindPairList[itr].push_back(index);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            // for (unsigned itr { 0 }; itr < numSpecies; ++itr) {
            //     unsigned long listSize { 0 };
            //     restartFile >> listSize;
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //     counterArrays.bindPairListIL2D.emplace_back();
            //     std::cout << "Specie " << itr << " N bindPairs " << listSize << '\t';
            //     for (unsigned itr2 { 0 }; itr2 < listSize; ++itr2) {
            //         int index { 0 };
            //         restartFile >> index;
            //         counterArrays.bindPairListIL2D[itr].push_back(index);
            //     }
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // }
            // for (unsigned itr { 0 }; itr < numSpecies; ++itr) {
            //     unsigned long listSize { 0 };
            //     restartFile >> listSize;
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //     counterArrays.bindPairListIL3D.emplace_back();
            //     std::cout << "Specie " << itr << " N bindPairs " << listSize << '\t';
            //     for (unsigned itr2 { 0 }; itr2 < listSize; ++itr2) {
            //         int index { 0 };
            //         restartFile >> index;
            //         counterArrays.bindPairListIL3D[itr].push_back(index);
            //     }
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // }
        }

        // The implicit lipid's 2D binding table and the protein counts it is
        // built from; see write_restart().  A file from before they were
        // written ends above.  That file still restarts, as it always did, but
        // both are then rebuilt from the state at the restart, so the run
        // cannot continue exactly.
        std::string implicitLipidSection;
        if (std::getline(restartFile, implicitLipidSection) && implicitLipidSection.find("#ImplicitLipid") == 0) {
            expect_restart_key(restartFile, "numberOfProteinEachState");
            membraneObject.numberOfProteinEachState.assign(membraneObject.nStates, 0);
            for (int& count : membraneObject.numberOfProteinEachState)
                restartFile >> count;

            expect_restart_key(restartFile, "binding2DTable");
            std::size_t tableSize { 0 };
            restartFile >> tableSize;
            for (std::size_t entry { 0 }; entry < tableSize; ++entry) {
                double ka { 0 };
                double Dtot { 0 };
                double kb { 0 };
                double value { 0 };
                restartFile >> ka >> Dtot >> kb >> value;
                membraneObject.ILTableIDs.push_back(ka);
                membraneObject.ILTableIDs.push_back(Dtot);
                membraneObject.ILTableIDs.push_back(kb);
                membraneObject.IL2DbindingVec.push_back(value);
            }
            if (!restartFile)
                throw std::string("Cannot read this restart file: its #ImplicitLipid section is truncated or malformed.");
        } else if (membraneObject.implicitLipid && restartFile.eof()) {
            std::cout << "This restart file has no #ImplicitLipid section, so it predates saving the implicit lipid's 2D "
                         "binding table and protein counts. They are rebuilt from the state at the restart, and the run "
                         "will not continue exactly as the one that wrote this file would have."
                      << std::endl;
        }
    } catch (const std::string& msg) {
        std::cerr << msg << '\n';
        exit(1);
    } catch (const std::length_error& e) {
        std::cerr << "Error in reading template vectors for " << e.what() << '\n';
        exit(1);
    }
}
