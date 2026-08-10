#!/usr/bin/env bash
# Run the reversible unimolecular state-change model and check both state changes
# against the equilibrium fraction they must reach.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# For U <-> P with forward rate kf and back rate kb, the fraction in P at
# equilibrium is kf/(kf+kb), independently of anything else in the model.  The
# input uses kf = 1000 s-1 and kb = 3000 s-1, so the expected fraction is 0.25.
#
# The check is a real test rather than a print: it fails if either species misses
# the prediction, and it fails if A(tag~P) is nonzero, since that reaction has
# rate 0 and exists only to offset the forward and back index spaces.
#
# A build whose reverse path does not fire makes P absorbing, so it lands near
# 1.00 instead of 0.25 and fails loudly.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../../.." && pwd)
MODEL_DIR=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
SEED=${SEED:-20260810}

# Predicted equilibrium fraction in the P state, and how far off we tolerate.
# 0.06 is about 3.5 binomial sigma on 200 copies (sigma = 0.0306), loose enough
# not to flake and far tighter than the ~1.0 a broken reverse path produces.
EXPECTED=0.25
TOLERANCE=0.06
COPIES=200

mkdir -p "$WORK"
cp "$MODEL_DIR"/*.mol "$MODEL_DIR"/uni_state_rev.inp "$WORK/"
cd "$WORK"

echo "model:  $MODEL_DIR"
echo "binary: $BINARY"
echo "seed:   $SEED"
echo "work:   $WORK"
echo

status=0
"$BINARY" -f uni_state_rev.inp -s "$SEED" > stdout.log 2> stderr.log || status=$?
if [[ $status -ne 0 ]]; then
    echo "FAIL: simulation exited $status (see $WORK/stdout.log)" >&2
    exit 1
fi

# Time-average the second half of the trajectory, which is past ~4 relaxation
# times, so the comparison is against a mean rather than one noisy frame.
# Columns: 1 Time, 4 A(ser~P), 6 A(tag~P), 9 B(thr~P).
awk -F, -v exp_frac="$EXPECTED" -v tol="$TOLERANCE" -v copies="$COPIES" '
    NR == 1 { next }
    $1 > 0.001 { n++; serP += $4; thrP += $9; if ($6 > tagPmax) tagPmax = $6 }
    END {
        if (n == 0) { print "FAIL: no equilibrated frames in DATA/copy_numbers_time.dat" > "/dev/stderr"; exit 1 }
        ser = serP / n / copies
        thr = thrP / n / copies
        printf "frames averaged:        %d\n", n
        printf "A(ser~P) fraction:      %.4f  (expected %.4f)\n", ser, exp_frac
        printf "B(thr~P) fraction:      %.4f  (expected %.4f)\n", thr, exp_frac
        printf "A(tag~P) peak:          %d  (expected 0, rate is 0)\n", tagPmax
        bad = 0
        d = ser - exp_frac; if (d < 0) d = -d
        if (d > tol) { printf "FAIL: A(ser~P) off by %.4f, tolerance %.4f\n", d, tol > "/dev/stderr"; bad = 1 }
        d = thr - exp_frac; if (d < 0) d = -d
        if (d > tol) { printf "FAIL: B(thr~P) off by %.4f, tolerance %.4f\n", d, tol > "/dev/stderr"; bad = 1 }
        if (tagPmax > 0) { printf "FAIL: A(tag~P) reached %d but its reaction has rate 0\n", tagPmax > "/dev/stderr"; bad = 1 }
        if (bad) exit 1
        print "\nPASS: both reversible state changes reached kf/(kf+kb)"
    }
' DATA/copy_numbers_time.dat
