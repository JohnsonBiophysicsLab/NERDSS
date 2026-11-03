#include "math/rand_gsl.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace {

enum class ReadResult {
    success,
    notFound,
    failed,
};

void report_io_error(const char* message, const std::string& filename, int errorNumber)
{
    std::cerr << message << " '" << filename << "'";
    if (errorNumber != 0) {
        std::cerr << ": " << std::strerror(errorNumber);
    }
    std::cerr << '\n';
}

[[noreturn]] void fail_rng_state_write(const std::string& filename, const char* message, int errorNumber)
{
    report_io_error(message, filename, errorNumber);
    std::exit(EXIT_FAILURE);
}

void write_rng_state_file(const std::string& filename)
{
    if (r == nullptr) {
        fail_rng_state_write(filename, "ERROR: RNG is not initialized; cannot write", 0);
    }

    errno = 0;
    FILE* stateOut = std::fopen(filename.c_str(), "wb");
    if (stateOut == nullptr) {
        fail_rng_state_write(filename, "ERROR: Could not open RNG state file for writing", errno);
    }

    const std::size_t stateSize = gsl_rng_size(r);
    errno = 0;
    const bool writeSucceeded = std::fwrite(gsl_rng_state(r), stateSize, 1, stateOut) == 1;
    int errorNumber = writeSucceeded ? 0 : errno;

    errno = 0;
    const bool closeSucceeded = std::fclose(stateOut) == 0;
    if (!closeSucceeded && errorNumber == 0) {
        errorNumber = errno;
    }

    if (!writeSucceeded || !closeSucceeded) {
        fail_rng_state_write(filename, "ERROR: Could not write RNG state file", errorNumber);
    }
}

ReadResult read_rng_state_file(const std::string& filename)
{
    if (r == nullptr) {
        report_io_error("ERROR: RNG is not initialized; cannot read", filename, 0);
        return ReadResult::failed;
    }

    errno = 0;
    FILE* stateIn = std::fopen(filename.c_str(), "rb");
    if (stateIn == nullptr) {
        const int errorNumber = errno;
        if (errorNumber == ENOENT) {
            return ReadResult::notFound;
        }
        report_io_error("ERROR: Could not open RNG state file for reading", filename, errorNumber);
        return ReadResult::failed;
    }

    const std::size_t stateSize = gsl_rng_size(r);
    std::vector<unsigned char> savedState(stateSize);

    errno = 0;
    const bool readSucceeded = std::fread(savedState.data(), stateSize, 1, stateIn) == 1;
    int errorNumber = readSucceeded ? 0 : errno;

    errno = 0;
    const bool closeSucceeded = std::fclose(stateIn) == 0;
    if (!closeSucceeded && errorNumber == 0) {
        errorNumber = errno;
    }

    if (!readSucceeded || !closeSucceeded) {
        report_io_error("ERROR: Could not read complete RNG state from", filename, errorNumber);
        return ReadResult::failed;
    }

    std::memcpy(gsl_rng_state(r), savedState.data(), stateSize);
    return ReadResult::success;
}

double inverse_rng_range(const gsl_rng* generator)
{
    static const gsl_rng_type* cachedType = nullptr;
    static double cachedInverseRange = 0.0;

    if (generator->type != cachedType) {
        cachedType = generator->type;
        cachedInverseRange = 1.0
            / (static_cast<double>(cachedType->max) - static_cast<double>(cachedType->min) + 1.0);
    }

    return cachedInverseRange;
}

void report_rng_state_read_failure()
{
    std::cerr << "RNG state was not restored; continuing with the current initialized RNG state.\n";
}

} // namespace

double rand_gsl()
{
    return gsl_rng_uniform(r);
}

double rand_gsl64()
{
    const double high = gsl_rng_uniform(r);
    const double low = gsl_rng_uniform(r);
    const double result = high + inverse_rng_range(r) * low;

    return result < 1.0 ? result : std::nextafter(1.0, 0.0);
}

void write_rng_state()
{
    write_rng_state_file("DATA/rng_state");
}

void write_rng_state(int rank)
{
    write_rng_state_file("DATA/rng_state_" + std::to_string(rank));
}

void write_rng_state_simItr(int simItr)
{
    write_rng_state_file("RESTARTS/rng_state" + std::to_string(simItr));
}

void read_rng_state()
{
    ReadResult result = read_rng_state_file("DATA/rng_state");
    if (result == ReadResult::notFound) {
        result = read_rng_state_file("rng_state");
    }

    if (result != ReadResult::success) {
        report_rng_state_read_failure();
    }
}

void read_rng_state(int rank)
{
    const std::string rankedFilename = "rng_state_" + std::to_string(rank);
    ReadResult result = read_rng_state_file("DATA/" + rankedFilename);
    if (result == ReadResult::notFound) {
        result = read_rng_state_file(rankedFilename);
    }

    if (result != ReadResult::success) {
        report_rng_state_read_failure();
    }
}

// ---- every rank's generator, for a restart at the same rank count ----
//
// nerdss_mpi gives each rank its own stream (rank 0 the seed, rank r the
// seed + r), and a checkpoint used to save rank 0's state alone under the
// serial build's names, so a restart at np > 1 could continue rank 0's
// stream at best and re-seeded the others.  write_gathered_restart() now has
// every rank write its state as well, under a name that carries the rank and
// the rank count: a state is only meaningful for the same rank of the same
// decomposition, and a restart at another count falls back to seeding.
//
// Names: DATA/rng_state_<rank>_np<nprocs> beside DATA/restart.json, and
// RESTARTS/rng_state<itr>_<rank>_np<nprocs> beside RESTARTS/restart<itr>.json.
// read_rng_state() below looks in DATA/, in the working directory, and in
// the directory the restart file was read from (so a checkpoint restarts from
// its own RESTARTS/ with nothing copied or renamed); at np=1 it also accepts
// the names an np=1 run and the serial build write, so an np=1 continuation
// of either is still exact.

namespace {
std::string rng_state_rank_tag(int rank, int nprocs)
{
    return "_" + std::to_string(rank) + "_np" + std::to_string(nprocs);
}
} // namespace

void write_rng_state(int rank, int nprocs)
{
    write_rng_state_file("DATA/rng_state" + rng_state_rank_tag(rank, nprocs));
}

void write_rng_state_simItr(long long simItr, int rank, int nprocs)
{
    write_rng_state_file("RESTARTS/rng_state" + std::to_string(simItr) + rng_state_rank_tag(rank, nprocs));
}

bool read_rng_state(int rank, int nprocs, const std::string& checkpointDir, long long simItr,
    std::string& fileRead)
{
    const std::string tagged = "rng_state" + rng_state_rank_tag(rank, nprocs);
    const std::string dir = checkpointDir.empty() ? std::string(".") : checkpointDir;
    std::vector<std::string> candidates {
        "DATA/" + tagged,
        tagged,
        dir + "/rng_state" + std::to_string(simItr) + rng_state_rank_tag(rank, nprocs),
    };
    if (nprocs == 1) {
        // what an np=1 run wrote before the tag existed, and the serial names
        candidates.push_back("DATA/rng_state_0");
        candidates.push_back("rng_state_0");
        candidates.push_back("DATA/rng_state");
        candidates.push_back("rng_state");
        candidates.push_back(dir + "/rng_state" + std::to_string(simItr));
    }
    for (const auto& name : candidates) {
        const ReadResult result = read_rng_state_file(name);
        if (result == ReadResult::success) {
            fileRead = name;
            return true;
        }
        // a file that exists but cannot be read has been reported; keep looking
    }
    fileRead.clear();
    return false;
}

std::string rng_state_digest()
{
    // FNV-1a over the generator's state bytes: two generators with the same
    // digest are (to the hash) the same generator, and will draw the same
    // sequence from here on.
    const unsigned char* bytes = static_cast<const unsigned char*>(gsl_rng_state(r));
    const std::size_t n = gsl_rng_size(r);
    unsigned long long h = 14695981039346656037ULL;
    for (std::size_t i = 0; i < n; i++) {
        h ^= bytes[i];
        h *= 1099511628211ULL;
    }
    char text[24];
    std::snprintf(text, sizeof text, "%016llx", h);
    return text;
}

void rng_probe_line(int rank, long long simItr)
{
    // NERDSS_RNG_PROBE in the environment (mpirun -x NERDSS_RNG_PROBE=1)
    // prints the digest at the top of every step to stderr, so that a
    // restarted run can be compared with the uninterrupted one at the same
    // step; see benchmarks/nerdss_optimized/mpi_validation/restart_np_check.sh.
    static const bool on = std::getenv("NERDSS_RNG_PROBE") != nullptr;
    if (!on) return;
    std::fprintf(stderr, "RNG STATE rank %d itr=%lld %s\n", rank, simItr, rng_state_digest().c_str());
}

double GaussV()
{
    // The Marsaglia polar method this used to implement draws ~2.55 uniforms on
    // average, takes a log, a sqrt and a division, and branches unpredictably
    // because it rejects samples outside the unit disc.  GSL's ziggurat sampler
    // answers from a table lookup about 98% of the time and only evaluates a
    // logarithm in the tail (issue #12).
#ifdef NERDSS_LEGACY_STREAM
    // master's Marsaglia polar sampler, verbatim
    double R { 2.0 };
    double V1 {};
    while (R >= 1.0) {
        V1 = 2.0 * gsl_rng_uniform(r) - 1.0;
        double V2 = 2.0 * gsl_rng_uniform(r) - 1.0;
        R = (V1 * V1) + (V2 * V2);
    }
    return (V1 * sqrt(-2.0 * log(R) / R));
#else
    return gsl_ran_gaussian_ziggurat(r, 1.0);
#endif
}
