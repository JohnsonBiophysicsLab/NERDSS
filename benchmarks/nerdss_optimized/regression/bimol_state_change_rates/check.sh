#!/usr/bin/env bash
# Run the reversible bimolecular state-change model and check its equilibrium
# against the ratio of the two rates.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# K(k) + A(a~X) <-> K(k) + A(a~Y) runs X -> Y at ka = 10 and Y -> X at kb = 50
# nm^3/us.  K is unchanged by either direction, so detailed balance puts the
# equilibrium at Y/X = ka/kb = 0.2, a fraction in Y of 1/6, whatever [K] and the
# diffusion constants are.  A build that runs Y -> X at the forward rate settles
# at 0.5 instead, and one whose back direction never fires drains X entirely.
#
# The check also fails if the copy counters stop adding up: K(k) and K(kd) must
# stay at 100 and A(a~X) + A(a~Y) at 200, which a state change counted under the
# relative state index broke.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
SEED=${SEED:-20260810}

# Predicted fraction in Y and how far off one run may land.  Across 6 seeds the
# per-run average from 2 ms on has a standard deviation of 0.005 (see README.md),
# so 0.04 is 8 sigma wide and still far from the 0.5 of a build that uses the
# forward rate both ways.
EXPECTED=0.1667
TOLERANCE=0.04
COPIES=200
K_COPIES=100

mkdir -p "$WORK"
cp "$SCRIPT_DIR"/*.mol "$SCRIPT_DIR"/parms.inp "$WORK/"
cd "$WORK"

echo "binary: $BINARY"
echo "seed:   $SEED"
echo "work:   $WORK"
echo

status=0
"$BINARY" -f parms.inp -s "$SEED" > stdout.log 2> stderr.log || status=$?
if [[ $status -ne 0 ]]; then
    echo "FAIL: simulation exited $status (see $WORK/stdout.log)" >&2
    exit 1
fi

# Columns: 1 Time, 2 K(k), 3 K(kd), 4 A(a~X), 5 A(a~Y), 6 K(kd!1).K(kd!1).
awk -F, -v exp_frac="$EXPECTED" -v tol="$TOLERANCE" -v copies="$COPIES" -v kcopies="$K_COPIES" '
    NR == 1 { next }
    $2 != kcopies || $3 != kcopies || $4 + $5 != copies { badCount++ }
    $1 >= 0.002 { n++; y += $5 }
    END {
        if (n == 0) { print "FAIL: no equilibrated frames in DATA/copy_numbers_time.dat" > "/dev/stderr"; exit 1 }
        frac = y / n / copies
        printf "frames averaged:        %d\n", n
        printf "A(a~Y) fraction:        %.4f  (expected %.4f)\n", frac, exp_frac
        printf "Y/X:                    %.4f  (expected %.4f)\n", frac / (1 - frac), exp_frac / (1 - exp_frac)
        bad = 0
        d = frac - exp_frac; if (d < 0) d = -d
        if (d > tol) { printf "FAIL: A(a~Y) fraction off by %.4f, tolerance %.4f\n", d, tol > "/dev/stderr"; bad = 1 }
        if (badCount > 0) { printf "FAIL: %d frames where the copy counters do not add up\n", badCount > "/dev/stderr"; bad = 1 }
        if (bad) exit 1
        print "\nPASS: the state change settled on the ratio of its two rates"
    }
' DATA/copy_numbers_time.dat
