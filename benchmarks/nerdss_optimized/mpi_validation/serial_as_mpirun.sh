#!/usr/bin/env bash
# Stand-in for "mpirun --oversubscribe" so run_mpi_suite.sh can drive the
# serial build: drops "-np N <binary>", runs $SERIAL_BIN with the rest, and
# copies the serial outputs to where the MPI merge would put them.
set -u
shift 2   # -np N
shift     # nerdss_mpi path
status=0
"$SERIAL_BIN" "$@" || status=$?
mkdir -p mergeOUT
cp DATA/copy_numbers_time.dat DATA/histogram_complexes_time.dat mergeOUT/ 2>/dev/null
exit $status
