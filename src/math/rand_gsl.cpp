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

void srand_gsl(int num)
{
    gsl_rng_set(r, num);
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

double GaussV()
{
    // The Marsaglia polar method this used to implement draws ~2.55 uniforms on
    // average, takes a log, a sqrt and a division, and branches unpredictably
    // because it rejects samples outside the unit disc.  GSL's ziggurat sampler
    // answers from a table lookup about 98% of the time and only evaluates a
    // logarithm in the tail (issue #12).
    return gsl_ran_gaussian_ziggurat(r, 1.0);
}
