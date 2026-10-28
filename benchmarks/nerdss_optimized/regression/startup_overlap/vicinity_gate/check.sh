#!/usr/bin/env bash
# Start the P/Q model at several seeds and check that no P and Q interface pair
# is left within sigma of each other in DATA/initial_crds.xyz.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# generate_coordinates() places every molecule at random and then re-places any
# molecule whose interface lies within a reaction's sigma of a partner's.  It
# only looks at a pair whose centers pass are_molecules_in_vicinity(), which
# used to mean closer than r1 + r2 + 2 nm, 3 nm here.  With sigma = 5 an
# interface pair can be inside sigma while the centers are anywhere up to 6 nm
# apart, and every such pair was left where it was drawn.  The gate is now
# r1 + r2 + the largest sigma between the two types.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
SEEDS=${SEEDS:-"1 2 3 4 5 6 7 8 9 10"}
SIGMA=5.0
# Each template's interface is 0.5 nm from its center.
OLD_GATE=3.0

echo "model:  $SCRIPT_DIR"
echo "binary: $BINARY"
echo "work:   $WORK"
echo

failed=0
for seed in $SEEDS; do
    run=$WORK/seed_$seed
    rm -rf "$run"
    mkdir -p "$run"
    cp "$SCRIPT_DIR"/P.mol "$SCRIPT_DIR"/Q.mol "$SCRIPT_DIR"/parms.inp "$run"

    status=0
    (cd "$run" && "$BINARY" -f parms.inp -s "$seed" > stdout.log 2> stderr.log) || status=$?
    if [[ $status -ne 0 ]]; then
        echo "FAIL: seed $seed: simulation exited $status (see $run/stdout.log)" >&2
        failed=1
        continue
    fi

    # initial_crds.xyz has two header lines, then each molecule's center and
    # then its one interface.
    status=0
    awk -v seed="$seed" -v sigma="$SIGMA" -v oldGate="$OLD_GATE" '
        NR <= 2 { next }
        (++row) % 2 == 1 { name = $1; cx = $2; cy = $3; cz = $4; next }
        {
            n = ++count[name]
            comX[name, n] = cx; comY[name, n] = cy; comZ[name, n] = cz
            ifX[name, n] = $2; ifY[name, n] = $3; ifZ[name, n] = $4
        }
        END {
            pairs = 0; beyond = 0; closest = 1e9
            for (p = 1; p <= count["P"]; p++)
                for (q = 1; q <= count["Q"]; q++) {
                    d = sqrt((ifX["P", p] - ifX["Q", q])^2 + (ifY["P", p] - ifY["Q", q])^2 + (ifZ["P", p] - ifZ["Q", q])^2)
                    if (d < closest) closest = d
                    if (d >= sigma) continue
                    pairs++
                    com = sqrt((comX["P", p] - comX["Q", q])^2 + (comY["P", p] - comY["Q", q])^2 + (comZ["P", p] - comZ["Q", q])^2)
                    if (com >= oldGate) beyond++
                }
            printf "seed %s: %d P, %d Q, %d P-Q interface pairs within %.1f nm (%d with centers >= %.1f nm apart), closest %.3f nm\n",
                seed, count["P"], count["Q"], pairs, sigma, beyond, oldGate, closest
            exit (pairs > 0) ? 1 : 0
        }
    ' "$run"/DATA/initial_crds.xyz || status=$?
    [[ $status -eq 0 ]] || failed=1
done

echo
if [[ $failed -ne 0 ]]; then
    echo "FAIL: a P and a Q start within sigma of each other" >&2
    exit 1
fi
echo "PASS"
