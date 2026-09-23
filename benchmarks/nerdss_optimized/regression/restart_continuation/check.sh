#!/usr/bin/env bash
# Check that a restart continues the run it was written by, bit for bit.
#
# Usage: [RESTART_FORMAT=json|dat] check.sh <path-to-nerdss> [work-dir]
#
# For each case in continuation_cases.tsv: run the model for 2N steps (or the
# case's own total) with checkPoint = N and RNGwrite = true; then, in a clean
# directory, restart from RESTARTS/restart<N>.<format> with
# RESTARTS/rng_state<N> copied in as rng_state, so read_rng_state() resumes
# the same stream.  The case passes when
#   * the two final DATA/restart.<format> files agree on every value but
#     currSimTime, which the restart reaches through a different but
#     equivalent expression (compare_restart_json.py for JSON; for .dat, cmp
#     with lines 2 and 4, numItr and currSimTime, left out), and
#   * every record the restart wrote to DATA/*_time.dat equals the
#     uninterrupted run's record for the same time (compare_series.py).
# Exit status 0 when every case passes.
#
# RESTART_FORMAT selects the restart file format: json (the default) or dat,
# which adds `legacyRestartFormat = true` to every case so that the checkpoint
# and the final files are written in the .dat format of earlier builds, and
# the restart, having read a .dat file, keeps writing it.
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
FORMAT=${RESTART_FORMAT:-json}
case $FORMAT in
    json|dat) ;;
    *) echo "RESTART_FORMAT must be json or dat, not '$FORMAT'" >&2; exit 2 ;;
esac
mkdir -p "$WORK"

echo "binary: $BINARY"
echo "format: $FORMAT"
echo "work:   $WORK"
echo

# Sets `report` and returns 0 when the two final restart files agree.
final_files_agree() { # <full run's file> <restart's file>
    if [[ $FORMAT == json ]]; then
        report=$(python3 "$SCRIPT_DIR/compare_restart_json.py" "$1" "$2" --ignore currSimTime)
    else
        # restart.dat holds NUL bytes (requiresState is written raw), so compare
        # it with sed and cmp, which pass them through; awk does not.
        if cmp -s <(sed '2d;4d' "$1") <(sed '2d;4d' "$2"); then
            report="identical (ignoring lines 2 and 4)"
        else
            report="differs in $(diff -a <(sed '2d;4d' "$1") <(sed '2d;4d' "$2") | grep -c '^<') lines"
            return 1
        fi
    fi
}

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
    awk -v total="$total" -v n="$n" -v legacy="$([[ $FORMAT == dat ]] && echo 1 || echo 0)" '
        { key = tolower($0) }
        key ~ /^[[:space:]]*nitr[[:space:]]*=/                                      { print "    nItr = " total; next }
        key ~ /^[[:space:]]*(checkpoint|rngwrite|legacyrestartformat)[[:space:]]*=/ { next }
        key ~ /^[[:space:]]*end[[:space:]]+parameters/ {
            print "    checkPoint = " n
            print "    RNGwrite = true"
            if (legacy) print "    legacyRestartFormat = true"
        }
        { print }
    ' "$case_dir/$parm" > "$full/parms.inp"

    if ! (cd "$full" && "$BINARY" -f parms.inp -s "$seed" > stdout.log 2> stderr.log); then
        printf '%-22s FAIL: the uninterrupted run exited non-zero\n' "$id"
        failures=$((failures + 1))
        continue
    fi
    if [[ ! -s "$full/RESTARTS/restart$n.$FORMAT" ]]; then
        printf '%-22s FAIL: no RESTARTS/restart%s.%s was written\n' "$id" "$n" "$FORMAT"
        failures=$((failures + 1))
        continue
    fi
    cp "$full/RESTARTS/restart$n.$FORMAT" "$cont/restart.$FORMAT"
    cp "$full/RESTARTS/rng_state$n" "$cont/rng_state"
    # The seed is irrelevant here: read_rng_state() replaces the stream.
    if ! (cd "$cont" && "$BINARY" -r "restart.$FORMAT" -s 1 > stdout.log 2> stderr.log); then
        printf '%-22s FAIL: the restart exited non-zero\n' "$id"
        failures=$((failures + 1))
        continue
    fi
    if [[ ! -s "$cont/DATA/restart.$FORMAT" ]]; then
        printf '%-22s FAIL: the restart did not write DATA/restart.%s\n' "$id" "$FORMAT"
        failures=$((failures + 1))
        continue
    fi

    if ! final_files_agree "$full/DATA/restart.$FORMAT" "$cont/DATA/restart.$FORMAT"; then
        printf '%-22s FAIL: final restart.%s differs\n' "$id" "$FORMAT"
        echo "$report" | sed 's/^/    /'
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
