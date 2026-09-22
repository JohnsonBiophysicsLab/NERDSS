#!/usr/bin/env bash
# Create 36000 B in a box already holding 1000 A and check that no B was left
# within a bindRadius of an A.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# moleculeOverlaps() rejects a place whose interfaces are within their
# reaction's bindRadius of another molecule's, and
# create_molecule_and_complex_from_rxn() draws again while it does.  Both
# templates put their one interface on the center and nothing diffuses, so the
# final coordinates are exactly the places that test accepted, and a brute-force
# count of A-B pairs closer than sigma over them is the whole verdict: it has to
# be zero.
#
# A build whose test walks only the SubBox the place falls in leaves thousands
# of them.  SubBoxes here are 5 nm on a side, sigma is 5 nm, and 523.6 nm^3 of
# exclusion sphere against a 125 nm^3 SubBox is why: 77% of each sphere lies in
# the 26 SubBoxes around it.
#
# Such a build is also slow -- the pairs it leaves behind are offered to the
# reaction machinery on every step of the run -- so give it minutes, not the
# twenty seconds a build that places the molecules properly takes.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
SEED=${SEED:-20260810}
SIGMA=${SIGMA:-5.0}

mkdir -p "$WORK"
cp "$SCRIPT_DIR"/A.mol "$SCRIPT_DIR"/B.mol "$SCRIPT_DIR"/parms.inp "$WORK"
cd "$WORK"

echo "model:  $SCRIPT_DIR"
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

# The grid the model asks for is worth reporting: the claim above is about the
# size of a SubBox against sigma, and this is where that shows.
grep -E 'Rmaxlimit|Dimensions|Sub-volume size' stdout.log || true
echo

# final_coords.xyz writes each molecule twice, its center and then its one
# interface, which here are the same point.  The A are binned on a grid of sigma
# so each B looks at 27 cells rather than at all 1000 A.
status=0
awk -v sigma="$SIGMA" '
    function cellOf(value,    quotient) {
        quotient = int(value / sigma)
        # int() truncates toward zero, which would make the cell holding 0 two
        # cells wide and put two cells one index apart 2 * sigma from each other.
        if (value < 0 && quotient * sigma != value) quotient--
        return quotient
    }
    NR <= 2 { next }
    # Every second row is the interface of the molecule named on the row before
    # it, at the same coordinates, so half the rows are skipped.
    (++row) % 2 == 1 {
        if ($1 == "A") {
            ax[++nA] = $2; ay[nA] = $3; az[nA] = $4
            key = cellOf($2) "," cellOf($3) "," cellOf($4)
            cell[key] = cell[key] " " nA
        } else {
            bx[++nB] = $2; by[nB] = $3; bz[nB] = $4
        }
    }
    END {
        pairs = 0
        closest = 1e9
        for (b = 1; b <= nB; b++) {
            bi = cellOf(bx[b]); bj = cellOf(by[b]); bk = cellOf(bz[b])
            for (di = -1; di <= 1; di++)
            for (dj = -1; dj <= 1; dj++)
            for (dk = -1; dk <= 1; dk++) {
                key = (bi + di) "," (bj + dj) "," (bk + dk)
                if (!(key in cell)) continue
                n = split(cell[key], member, " ")
                for (m = 1; m <= n; m++) {
                    a = member[m]
                    d = sqrt((ax[a] - bx[b])^2 + (ay[a] - by[b])^2 + (az[a] - bz[b])^2)
                    if (d < sigma) pairs++
                    if (d < closest) closest = d
                }
            }
        }
        printf "%d A, %d B, %d A-B pairs within %.3f nm, closest %.6f nm\n", nA, nB, pairs, sigma, closest
        exit (pairs > 0) ? 1 : 0
    }
' DATA/final_coords.xyz > pairs.txt || status=$?
cat pairs.txt

if [[ $status -ne 0 ]]; then
    echo "FAIL: a created B was left within $SIGMA nm of an A" >&2
    exit 1
fi
echo "PASS"
