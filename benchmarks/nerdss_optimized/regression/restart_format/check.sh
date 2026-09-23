#!/usr/bin/env bash
# Check that the JSON restart format and the .dat format it replaced carry the
# same state, and that each reads back everything it writes.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# restart_convert.cpp is built against <repo>/obj/serial, the objects the
# binary was linked from (build with `make serial` first), so the reader and
# writers under test are the binary's own.  For each case in cases.tsv the
# model is run, and then
#
#     J1 = DATA/restart.json, written at the end of the run
#     restart_convert J1  ->  J2 (JSON) and D1 (.dat)
#     restart_convert D1  ->  J3 (JSON) and D2 (.dat)
#
# must give J1 == J2 and D1 == D2 byte for byte, and J1 == J3 on every value
# but stateChangeIface, which only the JSON format carries, to the 1e-20 that
# the .dat format keeps (it prints most doubles in fixed notation with twenty
# decimal places, so a template interface coordinate of 4.7e-15 comes back
# from it with six significant digits):
#   J1 == J2  the JSON reader consumes everything the JSON writer wrote;
#   D1 == D2  the .dat reader consumes everything the .dat writer wrote;
#   J1 == J3  the state a .dat file holds is the state the JSON file holds.
#
# With REFERENCE_BINARY set to a build that still writes DATA/restart.dat, the
# same model is also run with it and its final restart.dat must equal D1 byte
# for byte: the JSON file, converted, reproduces the old format exactly.
#
# Exit status 0 when every case passes.  Needs bash, awk, python3 and a C++
# compiler (CXX, default g++).
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../../.." && pwd)
COMPARE=$SCRIPT_DIR/../restart_continuation/compare_restart_json.py

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
OBJ_DIR=${OBJ_DIR:-$ROOT_DIR/obj/serial}
REFERENCE_BINARY=${REFERENCE_BINARY:-}
mkdir -p "$WORK"

echo "binary:  $BINARY"
echo "objects: $OBJ_DIR"
[[ -n $REFERENCE_BINARY ]] && echo "reference: $REFERENCE_BINARY"
echo "work:    $WORK"
echo

CONVERT=$WORK/restart_convert
if ! "${CXX:-g++}" -O2 -std=c++0x $(gsl-config --cflags) -I"$ROOT_DIR/include" \
        -o "$CONVERT" "$SCRIPT_DIR/restart_convert.cpp" "$OBJ_DIR"/*/*.o $(gsl-config --libs) \
        > "$WORK/build.log" 2>&1; then
    echo "could not build restart_convert (see $WORK/build.log)" >&2
    exit 1
fi

run_model() { # <binary> <run-dir> <case-dir> <parm> <nItr> <seed>
    local binary=$1 run_dir=$2 case_dir=$3 parm=$4 n_itr=$5 seed=$6
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    cp "$case_dir"/*.mol "$run_dir/"
    awk -v n="$n_itr" '
        tolower($0) ~ /^[[:space:]]*nitr[[:space:]]*=/ { print "    nItr = " n; next }
        { print }
    ' "$case_dir/$parm" > "$run_dir/$parm"
    (cd "$run_dir" && "$binary" -f "$parm" -s "$seed" > stdout.log 2> stderr.log)
}

failures=0
while IFS=$'\t' read -r id dir parm n_itr seed; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    case_dir=$ROOT_DIR/$dir
    run=$WORK/$id/run
    conv=$WORK/$id/convert
    rm -rf "$WORK/$id"
    mkdir -p "$conv"

    if ! run_model "$BINARY" "$run" "$case_dir" "$parm" "$n_itr" "$seed"; then
        printf '%-20s FAIL: the run exited non-zero\n' "$id"
        failures=$((failures + 1))
        continue
    fi
    if [[ ! -s $run/DATA/restart.json ]]; then
        printf '%-20s FAIL: the run wrote no DATA/restart.json\n' "$id"
        failures=$((failures + 1))
        continue
    fi

    j1=$run/DATA/restart.json
    if ! "$CONVERT" "$j1" "$conv/J2.json" "$conv/D1.dat" > "$conv/convert1.log" 2>&1 \
       || ! "$CONVERT" "$conv/D1.dat" "$conv/J3.json" "$conv/D2.dat" > "$conv/convert2.log" 2>&1; then
        printf '%-20s FAIL: restart_convert exited non-zero (see %s)\n' "$id" "$conv"
        failures=$((failures + 1))
        continue
    fi

    problems=()
    cmp -s "$j1" "$conv/J2.json" || problems+=("JSON -> JSON changed the file")
    cmp -s "$conv/D1.dat" "$conv/D2.dat" || problems+=(".dat -> .dat changed the file")
    if ! report=$(python3 "$COMPARE" "$j1" "$conv/J3.json" --ignore stateChangeIface --abs-tol 1e-19); then
        problems+=("JSON -> .dat -> JSON lost state:")
        while IFS= read -r line; do problems+=("    $line"); done <<< "$report"
    fi
    if [[ -n $REFERENCE_BINARY ]]; then
        ref=$WORK/$id/reference
        if ! run_model "$REFERENCE_BINARY" "$ref" "$case_dir" "$parm" "$n_itr" "$seed"; then
            problems+=("the reference run exited non-zero")
        elif [[ ! -s $ref/DATA/restart.dat ]]; then
            problems+=("the reference run wrote no DATA/restart.dat")
        elif ! cmp -s "$ref/DATA/restart.dat" "$conv/D1.dat"; then
            lines=$(diff -a "$ref/DATA/restart.dat" "$conv/D1.dat" | grep -c '^<')
            problems+=("JSON converted to .dat differs from the reference build's restart.dat in $lines lines")
        fi
    fi

    if [[ ${#problems[@]} -eq 0 ]]; then
        printf '%-20s ok (%s steps, seed %s)\n' "$id" "$n_itr" "$seed"
    else
        printf '%-20s FAIL\n' "$id"
        printf '    %s\n' "${problems[@]}"
        failures=$((failures + 1))
    fi
done < "$SCRIPT_DIR/cases.tsv"

echo
if [[ $failures -eq 0 ]]; then
    echo "every case round-tripped"
else
    echo "$failures case(s) failed"
    exit 1
fi
