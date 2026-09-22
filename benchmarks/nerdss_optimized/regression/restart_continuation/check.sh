#!/usr/bin/env bash
# Check that a restart continues the run it was written by, bit for bit.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# For each case in continuation_cases.tsv: run the model for 2N steps (or the
# case's own total) with checkPoint = N and RNGwrite = true; then, in a clean
# directory, restart from RESTARTS/restart<N>.dat with RESTARTS/rng_state<N>
# copied in as rng_state, so read_rng_state() resumes the same stream.  The
# case passes when
#   * the two final DATA/restart.dat files are identical apart from line 2
#     (numItr) and line 4 (currSimTime, which the restart reaches through a
#     different but equivalent expression), and
#   * every record the restart wrote to DATA/*_time.dat equals the
#     uninterrupted run's record for the same time (compare_series.py).
# Exit status 0 when every case passes.
#
# The restart has to run in a clean directory: with DATA/trajectory.xyz
# present, nerdss checks the trajectory's length against the restart file and
# exits.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../../.." && pwd)

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
mkdir -p "$WORK"

echo "binary: $BINARY"
echo "work:   $WORK"
echo

failures=0
while IFS=$'\t' read -r id dir parm n seed total; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    total=${total:-$((2 * n))}
    case_dir=$ROOT_DIR/$dir
    full=$WORK/$id/full
    cont=$WORK/$id/cont
    rm -rf "$WORK/$id"
    mkdir -p "$full" "$cont"

    cp "$case_dir"/*.mol "$full/"
    # Keys are case-insensitive to the parser, so match them that way.
    awk -v total="$total" -v n="$n" '
        { key = tolower($0) }
        key ~ /^[[:space:]]*nitr[[:space:]]*=/                  { print "    nItr = " total; next }
        key ~ /^[[:space:]]*(checkpoint|rngwrite)[[:space:]]*=/ { next }
        key ~ /^[[:space:]]*end[[:space:]]+parameters/ { print "    checkPoint = " n; print "    RNGwrite = true" }
        { print }
    ' "$case_dir/$parm" > "$full/parms.inp"

    if ! (cd "$full" && "$BINARY" -f parms.inp -s "$seed" > stdout.log 2> stderr.log); then
        printf '%-22s FAIL: the uninterrupted run exited non-zero\n' "$id"
        failures=$((failures + 1))
        continue
    fi
    cp "$full/RESTARTS/restart$n.dat" "$cont/restart.dat"
    cp "$full/RESTARTS/rng_state$n" "$cont/rng_state"
    # The seed is irrelevant here: read_rng_state() replaces the stream.
    if ! (cd "$cont" && "$BINARY" -r restart.dat -s 1 > stdout.log 2> stderr.log); then
        printf '%-22s FAIL: the restart exited non-zero\n' "$id"
        failures=$((failures + 1))
        continue
    fi

    # restart.dat holds NUL bytes (requiresState is written raw), so compare
    # it with sed and cmp, which pass them through; awk does not.
    if ! cmp -s <(sed '2d;4d' "$full/DATA/restart.dat") <(sed '2d;4d' "$cont/DATA/restart.dat"); then
        lines=$(diff -a <(sed '2d;4d' "$full/DATA/restart.dat") <(sed '2d;4d' "$cont/DATA/restart.dat") | grep -c '^<')
        printf '%-22s FAIL: final restart.dat differs in %s lines\n' "$id" "$lines"
        failures=$((failures + 1))
    elif ! series=$(python3 "$SCRIPT_DIR/compare_series.py" "$full/DATA" "$cont/DATA"); then
        printf '%-22s FAIL: DATA time series differ\n' "$id"
        echo "$series" | grep -v ' ok ' | sed 's/^/    /'
        failures=$((failures + 1))
    else
        printf '%-22s ok (checkpoint %s of %s steps, seed %s)\n' "$id" "$n" "$total" "$seed"
    fi
done < "$SCRIPT_DIR/continuation_cases.tsv"

echo
if [[ $failures -eq 0 ]]; then
    echo "every case continued exactly"
else
    echo "$failures case(s) did not continue exactly"
    exit 1
fi
