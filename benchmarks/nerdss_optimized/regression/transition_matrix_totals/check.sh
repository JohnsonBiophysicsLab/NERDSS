#!/usr/bin/env bash
# Run the reversible A + R model with transition counting switched on and check
# each transition matrix against the one total it must have.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# A matrix is filled lazily.  When an event changes how many R a complex holds,
# every complex holding R is credited with the steps since the previous update
# for R: the complexes that took part get the stays before this step plus one
# transition, every other complex gets the whole interval as stays.  Either way
# each complex holding R gains one entry per step (a complex that splits into two
# which both hold R gains two, but that cannot happen here).  In this model every
# molecule has a single site, so a complex holds at most one R and the 1000
# copies of R always sit in exactly 1000 complexes.  The matrix for R must
# therefore total 1000 x lastUpdateTransition[R] at every write, whatever the
# trajectory, and likewise for A.
#
# A build that credits a destroyed complex breaks the total: associate_box()
# destroys the second reactant's complex and then counted its old contents a
# second time as an unchanged complex.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../../.." && pwd)
MODEL_DIR=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
SEED=${SEED:-20260810}
NITR=${NITR:-20000}
COPIES=1000

mkdir -p "$WORK"
cd "$WORK"

# The stock model does not count transitions; switch it on for both types.  A
# complex holds at most one of each, so a 2 x 2 matrix is ample (the start-up
# guardrail still warns, since it compares the size against the copy number).
for mol in A R; do
    awk '{ print } /^Name[[:space:]]*=/ { print "countTransition = true"; print "transitionMatrixSize = 2" }' \
        "$MODEL_DIR/$mol.mol" > "$mol.mol"
done
awk -v n="$NITR" '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next } { print }' \
    "$MODEL_DIR/parms3d.inp" > parms3d.inp

echo "model:  $MODEL_DIR (countTransition switched on, nItr = $NITR)"
echo "binary: $BINARY"
echo "seed:   $SEED"
echo "work:   $WORK"
echo

status=0
"$BINARY" -f parms3d.inp -s "$SEED" > stdout.log 2> stderr.log || status=$?
if [[ $status -ne 0 ]]; then
    echo "FAIL: simulation exited $status (see $WORK/stdout.log)" >&2
    exit 1
fi

# lastUpdateTransition is written on the line after bondedComplexWrite in the
# restart file, as its size followed by one step number per molecule type.  Both
# types count transitions here, so the j-th matrix in the transition file (they
# are written in template order) pairs with the j-th step number.
last_update=$(grep -a -A1 '^bondedComplexWrite' DATA/restart.dat | tail -1)

# Records are keyed by their position, not their time: the final write can
# print the same time as the periodic write before it.
awk -v copies="$COPIES" -v last="$last_update" '
    /^time: /                  { when[++nRecords] = $2; next }
    /^transition matrix/       { inMatrix = 1; next }
    /^lifetime/                { inMatrix = 0; next }
    inMatrix && NF == 1        { name = $1; if (!(name in index_of)) { index_of[name] = ++nTypes; names[nTypes] = name } next }
    inMatrix && NF > 1         { for (i = 1; i <= NF; i++) total[nRecords, name] += $i; next }
    END {
        if (nRecords == 0 || nTypes == 0) { print "FAIL: no transition matrices in DATA/transition_matrix_time.dat" > "/dev/stderr"; exit 1 }
        if (split(last, lut, " ") != nTypes + 1 || lut[1] != nTypes) {
            printf "FAIL: expected %d lastUpdateTransition entries in DATA/restart.dat, read \"%s\"\n",
                nTypes, last > "/dev/stderr"
            exit 1
        }
        bad = 0; nonzero = 0
        # Every record must be a whole number of steps for each of the copies.
        for (k = 1; k <= nRecords; k++)
            for (j = 1; j <= nTypes; j++)
                if (total[k, names[j]] % copies != 0) {
                    printf "FAIL: %s matrix totals %d at t = %s, not a multiple of %d\n",
                        names[j], total[k, names[j]], when[k], copies > "/dev/stderr"
                    bad = 1
                }
        # The final record must equal copies x lastUpdateTransition exactly.
        for (j = 1; j <= nTypes; j++) {
            got = total[nRecords, names[j]]
            want = copies * lut[j + 1]
            if (got > 0) nonzero = 1
            printf "%s: matrix total %d, expected %d x %d = %d  %s\n",
                names[j], got, copies, lut[j + 1], want, (got == want) ? "ok" : "MISMATCH"
            if (got != want) bad = 1
        }
        if (!nonzero) { print "FAIL: every matrix is empty, so nothing was checked" > "/dev/stderr"; bad = 1 }
        if (bad) exit 1
        printf "\nPASS: %d records, every matrix totals %d entries per step\n", nRecords, copies
    }
' DATA/transition_matrix_time.dat
