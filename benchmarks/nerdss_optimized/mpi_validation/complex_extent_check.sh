#!/usr/bin/env bash
# Checks the report for a complex wider than the stripe carries
# (check_complex_extent() in src/mpi/check_structure_invariant.cpp).
#
# 1. models/hex_chain, a linear polymer on an 8-column box, at np=3: a chain
#    of six or more molecules lying along x reaches past both ends of the
#    middle rank's three columns.  The -O3 build has to print one WARNING line
#    per rank the first time and a count at the end of the run; a checker
#    build (CHECK_STRUCTURE) prints every such complex under the STRUCTURE
#    VIOLATION label instead, and the run goes on.  The probe reads which
#    build it was given from the output.
# 2. hetTrimer at np=3 and np=4 (trimers of nanometre-sized molecules on
#    40 nm columns): nothing may be printed, and the run must end on its
#    usual final row.
#
# Usage: complex_extent_check.sh <nerdss_mpi> [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
MPI_BIN=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
WORK=${2:-$(mktemp -d)}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
SEED=${SEED:-2}
fail=0
mkdir -p "$WORK"

# --- 1
d=$WORK/hex_chain
rm -rf "$d"; mkdir -p "$d"; cp "$SCRIPT_DIR"/models/hex_chain/* "$d/"
(cd "$d" && $MPIRUN -np 3 "$MPI_BIN" -f parms.inp -s "$SEED" > stdout.log 2> stderr.log < /dev/null)
status=$?
warnings=$(grep -c '^WARNING rank [0-9]* itr=[0-9]*: complex id=' "$d/stderr.log" || true)
counts=$(grep -c '^WARNING rank [0-9]*: a complex wider than the stripe carries on' "$d/stderr.log" || true)
violations=$(grep -c 'STRUCTURE VIOLATION .* at=complex extent' "$d/stderr.log" || true)
first=$(grep 'complex id=' "$d/stderr.log" | grep -o 'itr=[0-9]*' | head -n 1)
rows=$(grep -c '' "$d/mergeOUT/copy_numbers_time.dat" 2>/dev/null || echo 0)
if [[ $status == 0 && $rows -ge 13 && ( ( $warnings -ge 1 && $warnings == "$counts" ) || $violations -ge 1 ) ]]; then
    if [[ $violations -ge 1 ]]; then
        echo "1 ok: hex_chain np=3 (checker build) reported $violations complex copies past both ends of the middle rank, first at $first"
    else
        echo "1 ok: hex_chain np=3 reported the first such complex at $first on $warnings rank(s), with a count at the end"
    fi
    grep -m 2 'complex id=' "$d/stderr.log" | cut -c1-160 | sed 's/^/    /'
    grep 'wider than the stripe' "$d/stderr.log" | sed 's/^/    /'
else
    echo "1 FAIL: hex_chain np=3 exit=$status rows=$rows warnings=$warnings counts=$counts violations=$violations"
    fail=1
fi

# --- 2
HET=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/hetTrimer
for np in 3 4; do
    d=$WORK/hetTrimer_np$np
    rm -rf "$d"; mkdir -p "$d"; cp "$HET"/*.mol "$d/"
    awk '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 3000"; next } { print }' \
        "$HET/parm_autodiff_hetTri.inp" > "$d/parm.inp"
    (cd "$d" && $MPIRUN -np "$np" "$MPI_BIN" -f parm.inp -s "$SEED" > stdout.log 2> stderr.log < /dev/null)
    status=$?
    noise=$(grep -c 'complex id=\|wider than the stripe\|complex extent' "$d/stderr.log" || true)
    last=$(tail -n 1 "$d/mergeOUT/copy_numbers_time.dat" 2>/dev/null | cut -d, -f1)
    if [[ $status == 0 && $noise == 0 && $last == 0.0008997 ]]; then
        echo "2 ok: hetTrimer np=$np silent, final row at $last"
    else
        echo "2 FAIL: hetTrimer np=$np exit=$status reports=$noise last row '$last'"
        fail=1
    fi
done

[[ $fail == 0 ]] && echo "complex_extent_check: ok"
exit $fail
