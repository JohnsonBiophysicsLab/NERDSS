#pragma once

#include "gsl/gsl_rng.h"
#include <gsl/gsl_randist.h>

#include <string>

extern gsl_rng* r;

/*!
 * \brief Uses the previously initialized random number generator to return a random number
 * \param[out] double Uniformly distributed random double.
 */
double rand_gsl();
double rand_gsl64();

/*!
 * \brief Wrapper for the internal GSL RNG state read function.
 *
 * Reads a previously written binary file with the current status of the RNG, so restarting will give the same random
 * numbers as a continuous run would.
 */
void read_rng_state();
void read_rng_state(int rank);

/*!
 * \brief Wrapper for the internal GSL RNG state write function.
 *
 * Writes the current state of the RNG to a binary file.
 */
void write_rng_state();
void write_rng_state(int rank);

/*!
 * \brief Wrapper for the internal GSL RNG state write function.
 *
 * Writes the current state of the RNG to a binary file for check point.
 */
void write_rng_state_simItr(int simItr);

/*!
 * \brief Every rank's generator, named by rank and rank count, for a restart
 * at the same count (nerdss_mpi); see rand_gsl.cpp.
 */
void write_rng_state(int rank, int nprocs);
void write_rng_state_simItr(long long simItr, int rank, int nprocs);

/*!
 * \brief Restore this rank's generator from the state a run at the same rank
 * count saved: DATA/, the working directory, then the directory the restart
 * file came from (checkpointDir, for RESTARTS/rng_state<simItr>_<rank>_np<N>).
 * Returns false, with the generator untouched, when no such state exists;
 * fileRead names the file on success.
 */
bool read_rng_state(int rank, int nprocs, const std::string& checkpointDir, long long simItr,
    std::string& fileRead);

/*!
 * \brief A hash of the generator's state bytes, for comparing two runs' generators.
 */
std::string rng_state_digest();

/*!
 * \brief Prints the digest to stderr when NERDSS_RNG_PROBE is set; a no-op otherwise.
 */
void rng_probe_line(int rank, long long simItr);

/*!
 * \brief Returns a standard-normal random number (mean 0, variance 1).
 *
 * Backed by GSL's ziggurat sampler.  This replaced a hand-written Marsaglia
 * polar (rejection Box-Muller) implementation, which was slower and consumed a
 * variable number of uniform variates per call.
 *
 * \param[out] double Gaussian-distributed random double.
 */
double GaussV();
