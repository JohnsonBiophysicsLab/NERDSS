/*! \file restart_convert.cpp
 * \brief Reads a restart file of either format and writes it back in both.
 *
 *     restart_convert <restart file> <out.json> <out.dat>
 *
 * check.sh builds this against the objects of a `make serial` build, so the
 * reader and the two writers are exactly those of the binary under test.  The
 * three globals are defined by nerdss's own main() and are unused here.
 */
#include "io/io.hpp"
#include "reactions/association/association.hpp"

#include <cstdio>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <map>
#include <string>
#include <vector>

gsl_rng* r;
unsigned long totMatches = 0;
long long randNum = 0;

int main(int argc, char* argv[])
{
    if (argc != 4) {
        std::fprintf(stderr, "usage: restart_convert <restart file> <out.json> <out.dat>\n");
        return 2;
    }

    long long int simItr { 0 };
    Parameters params;
    SimulVolume simulVolume;
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<MolTemplate> molTemplateList;
    std::vector<ForwardRxn> forwardRxns;
    std::vector<BackRxn> backRxns;
    std::vector<CreateDestructRxn> createDestructRxns;
    std::vector<TransmissionRxn> transmissionRxns;
    std::map<std::string, int> observablesList;
    Membrane membraneObject;
    copyCounters counterArrays;
    // The .dat reader fills the event arrays in place, so they must exist.
    init_association_events(counterArrays);

    std::ifstream in { argv[1] };
    if (!in) {
        std::fprintf(stderr, "restart_convert: cannot open %s\n", argv[1]);
        return 1;
    }
    // Exits on a malformed file.
    read_restart(simItr, in, params, simulVolume, moleculeList, complexList, molTemplateList, forwardRxns, backRxns,
        createDestructRxns, transmissionRxns, observablesList, membraneObject, counterArrays);
    in.close();

    std::ofstream outJson { argv[2] };
    write_json_restart(simItr, outJson, params, simulVolume, moleculeList, complexList, molTemplateList, forwardRxns,
        backRxns, createDestructRxns, transmissionRxns, observablesList, membraneObject, counterArrays);
    const bool jsonOk { static_cast<bool>(outJson) };
    outJson.close();

    std::ofstream outDat { argv[3] };
    LEGACY_write_restart(simItr, outDat, params, simulVolume, moleculeList, complexList, molTemplateList,
        forwardRxns, backRxns, createDestructRxns, transmissionRxns, observablesList, membraneObject, counterArrays);
    const bool datOk { static_cast<bool>(outDat) };
    outDat.close();

    if (!jsonOk || !datOk) {
        std::fprintf(stderr, "restart_convert: could not write %s\n", jsonOk ? argv[3] : argv[2]);
        return 1;
    }
    return 0;
}
