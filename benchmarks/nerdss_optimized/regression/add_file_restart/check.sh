#!/usr/bin/env bash
# Check that a restart with an add file (`nerdss -r DATA/restart.json -a
# add.inp`) enters the add file's molecules and reactions and runs
# reproducibly, and, given a reference build, that it produces the same bytes.
#
# Usage: [REFERENCE_BINARY=<path>] check.sh <path-to-nerdss> [work-dir]
#
# For each case in cases.tsv the model is run twice in one directory:
#
#   leg 1   nerdss -f parms.inp -s <seed1>              nItr = N, restartWrite = N
#   leg 2   nerdss -r DATA/restart.json -a add.inp -s <seed2>   add.inp's nItr = 2N
#
# and every file under DATA/ and RESTARTS/ is hashed into manifest.sha256.  The
# case passes when
#   * both legs exit 0 and leg 2 reports the add file on stdout,
#   * the header of DATA/copy_numbers_time.dat written by leg 2 carries every
#     species the case's `adds` column names, and the header leg 1 wrote
#     carries none of them: the add file's reactions were entered,
#   * the final DATA/restart.json has simItr = 2N - 1: leg 2 ran to the add
#     file's nItr,
#   * a second pair of legs with the same seeds gives an identical manifest,
#   * with REFERENCE_BINARY set, that build's manifest is identical as well.
# Exit status 0 when every case passes.
#
# Leg 2 MUST be given its own seed.  No DATA/rng_state is written unless
# RNGwrite = true, so the restart does not restore the RNG and seeds from -s,
# or from the clock when -s is absent: two unseeded restarts differ in
# final_coords.xyz, trajectory.xyz and restart.json.  The seed is what makes
# the two-run comparison, and the comparison between builds, meaningful.
#
# Leg 2 runs in leg 1's directory on purpose: that is how an add-file restart
# is used, and the trajectory-length check reads ./trajectory.xyz, not
# DATA/trajectory.xyz, so it finds no file and only warns.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../../.." && pwd)

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
REFERENCE_BINARY=${REFERENCE_BINARY:-}
if [[ -n $REFERENCE_BINARY ]]; then
    REFERENCE_BINARY=$(cd "$(dirname "$REFERENCE_BINARY")" && pwd)/$(basename "$REFERENCE_BINARY")
fi
mkdir -p "$WORK"

echo "binary:    $BINARY"
[[ -n $REFERENCE_BINARY ]] && echo "reference: $REFERENCE_BINARY"
echo "work:      $WORK"
echo

if command -v shasum >/dev/null 2>&1; then
    sha_files() { xargs shasum -a 256; }
else
    sha_files() { xargs sha256sum; }
fi

# Runs both legs of a case in <run-dir>.  Sets `why` and returns non-zero when
# a leg fails.  Leaves header1 (leg 1's copy_numbers header) and
# manifest.sha256 in the directory.
run_case() { # <binary> <run-dir> <case-dir> <N> <seed1> <seed2>
    local binary=$1 run=$2 case_dir=$3 n=$4 seed1=$5 seed2=$6
    rm -rf "$run"
    mkdir -p "$run"
    cp "$case_dir"/*.mol "$run/"
    # Keys are case-insensitive to the parser, so match them that way.
    awk -v n="$n" '
        { key = tolower($0) }
        key ~ /^[[:space:]]*nitr[[:space:]]*=/         { print "    nItr = " n; next }
        key ~ /^[[:space:]]*restartwrite[[:space:]]*=/ { print "    restartWrite = " n; next }
        { print }
    ' "$case_dir/parms.inp" > "$run/parms.inp"
    awk -v n="$((2 * n))" '
        tolower($0) ~ /^[[:space:]]*nitr[[:space:]]*=/ { print "    nItr = " n; next }
        { print }
    ' "$case_dir/add.inp" > "$run/add.inp"

    if ! (cd "$run" && "$binary" -f parms.inp -s "$seed1" > stdout1.log 2> stderr1.log); then
        why="leg 1 (the run that writes the restart file) exited non-zero"; return 1
    fi
    if [[ ! -s "$run/DATA/restart.json" ]]; then
        why="leg 1 did not write DATA/restart.json"; return 1
    fi
    head -1 "$run/DATA/copy_numbers_time.dat" > "$run/header1"
    if ! (cd "$run" && "$binary" -r DATA/restart.json -a add.inp -s "$seed2" > stdout2.log 2> stderr2.log); then
        why="leg 2 (the restart with the add file) exited non-zero"; return 1
    fi
    if ! grep -q 'restart simulation with add file' "$run/stdout2.log"; then
        why="leg 2 did not report the add file"; return 1
    fi
    (cd "$run" && find DATA RESTARTS -type f | LC_ALL=C sort | sha_files > manifest.sha256)
}

# Sets `why` and returns non-zero when the header lacks (all) or has (none)
# one of the comma-separated names.  The header is a CSV line whose first
# field is the time column.
header_has() { # <header line> <comma-separated names> <all|none>
    local header="$1," names=$2 want=$3 name
    IFS=',' read -ra list <<< "$names"
    for name in "${list[@]}"; do
        if [[ $header == *",$name,"* ]]; then
            [[ $want == none ]] && { why="$name was already in leg 1's copy_numbers header"; return 1; }
        else
            [[ $want == all ]] && { why="$name is missing from leg 2's copy_numbers header"; return 1; }
        fi
    done
    return 0
}

final_sim_itr() { # <restart.json>
    python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["parameters"]["simItr"])' "$1"
}

# Lists the files whose hashes differ between two manifests.
manifests_differ() { # <manifest a> <manifest b>
    diff <(awk '{print $2, $1}' "$1") <(awk '{print $2, $1}' "$2") \
        | awk '$1 == "<" || $1 == ">" {print $2}' | sort -u | sed 's/^/    /'
}

failures=0
while IFS=$'\t' read -r id dir n seed1 seed2 adds; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    case_dir=$ROOT_DIR/$dir
    first=$WORK/$id/first
    again=$WORK/$id/again
    ref=$WORK/$id/reference
    rm -rf "$WORK/$id"

    if ! run_case "$BINARY" "$first" "$case_dir" "$n" "$seed1" "$seed2"; then
        printf '%-12s FAIL: %s\n' "$id" "$why"; failures=$((failures + 1)); continue
    fi
    if ! header_has "$(cat "$first/header1")" "$adds" none \
        || ! header_has "$(head -1 "$first/DATA/copy_numbers_time.dat")" "$adds" all; then
        printf '%-12s FAIL: %s\n' "$id" "$why"; failures=$((failures + 1)); continue
    fi
    last=$(final_sim_itr "$first/DATA/restart.json")
    if [[ $last != $((2 * n - 1)) ]]; then
        printf '%-12s FAIL: final restart.json has simItr %s, not %s\n' "$id" "$last" "$((2 * n - 1))"
        failures=$((failures + 1)); continue
    fi
    if ! run_case "$BINARY" "$again" "$case_dir" "$n" "$seed1" "$seed2"; then
        printf '%-12s FAIL: second run: %s\n' "$id" "$why"; failures=$((failures + 1)); continue
    fi
    if ! cmp -s "$first/manifest.sha256" "$again/manifest.sha256"; then
        printf '%-12s FAIL: two runs of %s differ in\n' "$id" "$BINARY"
        manifests_differ "$first/manifest.sha256" "$again/manifest.sha256"
        failures=$((failures + 1)); continue
    fi
    if [[ -n $REFERENCE_BINARY ]]; then
        if ! run_case "$REFERENCE_BINARY" "$ref" "$case_dir" "$n" "$seed1" "$seed2"; then
            printf '%-12s FAIL: reference: %s\n' "$id" "$why"; failures=$((failures + 1)); continue
        fi
        if ! cmp -s "$first/manifest.sha256" "$ref/manifest.sha256"; then
            printf '%-12s FAIL: differs from the reference build in\n' "$id"
            manifests_differ "$first/manifest.sha256" "$ref/manifest.sha256"
            failures=$((failures + 1)); continue
        fi
    fi
    printf '%-12s ok (%s files, %s + %s steps, seeds %s / %s, adds %s%s)\n' "$id" \
        "$(wc -l < "$first/manifest.sha256" | tr -d ' ')" "$n" "$n" "$seed1" "$seed2" "$adds" \
        "$([[ -n $REFERENCE_BINARY ]] && echo ', matches the reference')"
done < "$SCRIPT_DIR/cases.tsv"

echo
if [[ $failures -eq 0 ]]; then
    echo "every add-file restart passed"
else
    echo "$failures case(s) failed"
    exit 1
fi
