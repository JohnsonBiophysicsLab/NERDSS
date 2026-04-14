#include "io/io.hpp"
#include "json.hpp"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <tuple>

#include "reactions/association/association.hpp"

using json = nlohmann::json;

// forward declaration
Quat molecule_orientation(const MolTemplate& oneTemplate, const Molecule& targMol);

void write_bonded_complex_json(const std::string filename,
    std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns,
    const Membrane& membraneObject)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Unable to open JSON file for writing: " << filename << std::endl;
        return;
    }

    json result = json::array();

    // collect bonded complexes
    std::vector<size_t> bonded_complex_indices;
    for (size_t k = 0; k < complexList.size(); k++) {
        if (!complexList[k].isEmpty && complexList[k].memberList.size() > 1) {
            bonded_complex_indices.push_back(k);
        }
    }

    double xOffset = membraneObject.waterBox.x / 2.0;
    double yOffset = membraneObject.waterBox.y / 2.0;
    double zOffset = membraneObject.waterBox.z / 2.0;

    for (const size_t k : bonded_complex_indices)
    {
        const auto& memberList = complexList[k].memberList;

        // for quickly looking up a molecule's index in memberList using its index in moleculeList.
        // We will populate this in the first loop where we get the mol names
        // it is used for bonds later
        std::vector<int> localIndLookup(moleculeList.size(), -1);

        // the json object we will populate with data
        json complex_json;

        // molecule names
        json names = json::array();
        for (size_t i = 0; i < memberList.size(); i++) {
            int molInd = memberList[i];
            names.push_back(molTemplateList[moleculeList[molInd].molTypeIndex].molName);
            localIndLookup[molInd] = i;
        }
        complex_json["names"] = names;

        // center of mass coordinates
        json coords = json::array();
        for (size_t i = 0; i < memberList.size(); i++) {
            int molInd = memberList[i];
            coords.push_back({
                moleculeList[molInd].comCoord.x + xOffset,
                moleculeList[molInd].comCoord.y + yOffset,
                moleculeList[molInd].comCoord.z + zOffset
            });
        }
        complex_json["coords"] = coords;

        // orientation, given as the quaternion rotation about the CoM required
        // to turn the .mol template into the current orientation
        json rotations = json::array();
        for (size_t i = 0; i < memberList.size(); i++) {
            int molInd = memberList[i];
            Quat orientation =
                molecule_orientation(
                    molTemplateList[moleculeList[molInd].molTypeIndex],
                    moleculeList[molInd]
                ).inverse(); // inverse because molecule_orientation() gives rotation FROM instance TO template

            rotations.push_back({
                orientation.w,
                orientation.x,
                orientation.y,
                orientation.z
            });
        }
        complex_json["rotations"] = rotations;

        // bond types
        json bond_types = json::object();
        std::vector<std::string> bondNames;

        // track already found pairs (localindex1,localindex2,ifaceindex1,ifaceindex2)
        // we will then loop over these to output bond instances after the type loop
        std::vector<std::tuple<int,int,int,int>> boundPairs;

        for (size_t i = 0; i < memberList.size(); i++)
        {
            int molInd = memberList[i];

            // loop over molecule interfaces to find bonds
            for (size_t j = 0; j < moleculeList[molInd].interfaceList.size(); j++)
            {
                if (!moleculeList[molInd].interfaceList[j].isBound)
                    continue;

                int partner_molInd = moleculeList[molInd].interfaceList[j].interaction.partnerIndex;
                int partner_i = localIndLookup[partner_molInd]; // index of partner in memberList

                int backRxnInd = moleculeList[molInd].interfaceList[j].interaction.conjBackRxn;
                size_t rxnInd = backRxns[backRxnInd].conjForwardRxnIndex;

                std::string bName = forwardRxns[rxnInd].productName;

                if (std::find(bondNames.begin(), bondNames.end(), bName) == bondNames.end())
                {
                    // haven't seen this bond type yet, add it to the list
                    bondNames.push_back(bName);

                    auto& rxn = forwardRxns[rxnInd];

                    auto maybe_null = [](double val) -> json {
                        return std::isnan(val) ? json(nullptr) : json(val);
                    };

                    bond_types[bName] = {
                        {"n1", {rxn.norm1.x, rxn.norm1.y, rxn.norm1.z}},
                        {"n2", {rxn.norm2.x, rxn.norm2.y, rxn.norm2.z}},
                        {"sigma", rxn.bindRadius},
                        {"theta1", maybe_null(rxn.assocAngles.theta1)},
                        {"theta2", maybe_null(rxn.assocAngles.theta2)},
                        {"phi1",   maybe_null(rxn.assocAngles.phi1)},
                        {"phi2",   maybe_null(rxn.assocAngles.phi2)},
                        {"omega",  maybe_null(rxn.assocAngles.omega)}
                    };

                    boundPairs.emplace_back(
                        i, partner_i, j,
                        moleculeList[molInd].interfaceList[j].interaction.partnerIfaceIndex
                    );
                }
                else
                {
                    // we've already seen that bond *type*, but have we seen that *bond instance*?
                    bool alreadyFound = false;
                    for (const auto& pair : boundPairs) {
                        if (std::get<0>(pair) == partner_i && std::get<1>(pair) == i)
                        {
                            alreadyFound = true;
                            break;
                        }
                    }
                    if (!alreadyFound) {
                        boundPairs.emplace_back(
                            i, partner_i, j,
                            moleculeList[molInd].interfaceList[j].interaction.partnerIfaceIndex
                        );
                    }
                }
            }
        }

        complex_json["bond_types"] = bond_types;

        // bond instances (which refer to the above types)
        json bonds = json::array();

        for (const auto& pair : boundPairs)
        {
            int ind1 = std::get<0>(pair);
            int ind2 = std::get<1>(pair);
            int ifaceInd1 = std::get<2>(pair);
            int ifaceInd2 = std::get<3>(pair);

            std::string site1Name =
                molTemplateList[
                    moleculeList[memberList[ind1]].molTypeIndex
                ].interfaceList[ifaceInd1].name;

            std::string site2Name =
                molTemplateList[
                    moleculeList[memberList[ind2]].molTypeIndex
                ].interfaceList[ifaceInd2].name;

            int backRxnInd =
                moleculeList[memberList[ind1]]
                .interfaceList[ifaceInd1]
                .interaction.conjBackRxn;

            std::string bName =
                forwardRxns[
                    backRxns[backRxnInd].conjForwardRxnIndex
                ].productName;

            bonds.push_back({
                {"molindex1", ind1},
                {"molindex2", ind2},
                {"site1", site1Name},
                {"site2", site2Name},
                {"type", bName}
            });
        }

        complex_json["bonds"] = bonds;

        result.push_back(complex_json);
    }

    // pretty print (2-space indent)
    out << result.dump(2) << std::endl;
    out.close();
}

// this is a slightly modified version of orient_crds_to_template()
// This one does not make weird state modifications to the molecule
Quat molecule_orientation(const MolTemplate& oneTemplate, const Molecule& targMol)
{
    // TRACE();
    Quat firstRot {};
    Quat secondRot {};
    // make a temporary copy of the interface coordinates for us to modify
    std::vector<Coord> tempIfaceCoords;
    for (auto& iface : targMol.interfaceList)
    {
        tempIfaceCoords.push_back(Coord { iface.coord } );
    }
    { // first rotation
        // determine the interface-center of mass vector (v0, v1), the rotation vector (u), and the angle to rotate
        // (angle)
        // arbitrarily the center of mass to first interface of the MolTemplate
        Vector vec1 { oneTemplate.interfaceList[0].iCoord - oneTemplate.comCoord };
        // the center of mass to first interface of the target Molecule
        Vector vec2 { targMol.interfaceList[0].coord - targMol.comCoord };
        // the rotation vector formed by v0 cross v1
        vec1.calc_magnitude();
        vec2.calc_magnitude();
        Vector rotAxis { vec1.cross(vec2) };
        double angle { vec1.dot_theta(vec2) };
        // TODO: the above vectors are correct
        // determine if we need to flip the sign
        //        if (requiresSignFlip(rotAxis, vec1, vec2))
        //            angle = -angle;
        double sa = sin(angle / 2);
        firstRot = Quat { cos(angle / 2), sa * rotAxis.x, sa * rotAxis.y,
            sa * rotAxis.z };
        firstRot = firstRot.unit();
        {
            Vector tmpVec { vec2 };
            firstRot.rotate(tmpVec);
            if (std::abs(tmpVec.x - oneTemplate.interfaceList[0].iCoord.x) > 1E-8 || std::abs(tmpVec.y - oneTemplate.interfaceList[0].iCoord.y) > 1E-8 || std::abs(tmpVec.z - oneTemplate.interfaceList[0].iCoord.z) > 1E-8) {
                angle = -angle;
                sa = sin(angle / 2);
                firstRot = Quat { cos(angle / 2), sa * rotAxis.x, sa * rotAxis.y,
                    sa * rotAxis.z };
                firstRot = firstRot.unit();
            }
        }
        // generate intermediate interface coordinates after first rotation
        for (size_t i = 0; i < targMol.interfaceList.size(); i++)
        {
            Vector tmpVec { targMol.interfaceList[i].coord - targMol.comCoord };
            firstRot.rotate(tmpVec);
            tempIfaceCoords[i] = tmpVec;
        }
    }

    if (targMol.interfaceList.size() > 1) {
        // if the protein has more than one interface, use a second one to make sure all the interfaces line up
        
        size_t ifaceIndex { 1 };
        
        Vector ifaceVec1 { oneTemplate.interfaceList[0].iCoord - oneTemplate.comCoord };
        ifaceVec1.calc_magnitude();

        // First check to make sure they're not in a line. If so, use the last interface

        bool allcollinear{true};
        for (size_t tmpIndex=0; tmpIndex < oneTemplate.interfaceList.size(); tmpIndex++){
            // check all other interfaces to find one that is not collinear to the first one
            Vector ifaceVec2 {oneTemplate.interfaceList[tmpIndex].iCoord - oneTemplate.comCoord };
            ifaceVec2.calc_magnitude();
            double ang1 { ifaceVec1.dot_theta(ifaceVec2) };
            if (ang1 == 0 || ang1 == M_PI) {
                // check next
            } else {
                ifaceIndex = tmpIndex;
                allcollinear = false;
                break;
            }
        }
        if (allcollinear && !oneTemplate.isRod){
            // all interfaces are collinear, just use the first rotation
            // By Mankun Sang (msang2@jh.edu): I don't know what isRod does here,
            //   and I just preserved the logic. If this is a rod, just do what's remaining. 
            return firstRot;
        } else {
            Vector v0 { oneTemplate.interfaceList[ifaceIndex].iCoord - oneTemplate.comCoord };
            Vector v1 { tempIfaceCoords[ifaceIndex] - targMol.tmpComCoord };
            Vector rotAxis { tempIfaceCoords[0] - targMol.tmpComCoord };
            v0.calc_magnitude();
            v1.calc_magnitude();
            rotAxis.normalize();

            // project the current and desired iface-com vectors
            // onto a plane of which the first iface-com vector is normal to
            Vector projVec0(v0.vector_projection(rotAxis));
            Vector projVec1(v1.vector_projection(rotAxis));
            projVec0.calc_magnitude();
            projVec1.calc_magnitude();
            double angle { projVec0.dot_theta(projVec1) };

            if (requiresSignFlip(rotAxis, projVec0, projVec1))
                angle = -angle;

            double sa { std::sin(angle / 2) };
            secondRot = Quat(cos(angle / 2), sin(angle / 2) * rotAxis.x, sin(angle / 2) * rotAxis.y, sin(angle / 2) * rotAxis.z);
            secondRot = secondRot.unit();
        }
    }

    // now use the inverse quat product to rotate the norm to
    // its actual position
    return (targMol.interfaceList.size() > 1) ? (secondRot * firstRot) : firstRot;
}