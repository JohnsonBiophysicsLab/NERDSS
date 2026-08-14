#!/usr/bin/env bash
# Run every case in a case table with one nerdss build and record, per case:
#   * wall-clock seconds for each repetition
#   * a content hash of every output file
#
# Usage: run_suite.sh <path-to-nerdss> <label> [repetitions] [result-root]
#
# The hash manifest is what the bitwise comparison uses.  Line 1 of every PDB
# frame carries the wall-clock creation time, so PDB files are hashed from
# line 2 onward while every other output is hashed raw.  Nothing else in the
# outputs depends on the clock, which `compare_suites.sh` relies on.
#
# CASES_FILE selects the table; it defaults to cases.tsv.  A table may carry an
# optional fifth column holding one parameter line to inject into the model's
# `parameters` block, which is how coverage_cases.tsv reaches code paths that
# are off by default:
#
#   CASES_FILE=coverage_cases.tsv ./run_suite.sh ../../bin/nerdss cov 1
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../.." && pwd)

BINARY=${1:?usage: run_suite.sh <path-to-nerdss> <label> [repetitions] [result-root]}
LABEL=${2:?missing label}
REPS=${3:-1}
RESULT_ROOT=${4:-$SCRIPT_DIR/results}

BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
SEED=${SUITE_SEED:-20260810}
INPUT_ROOT=$ROOT_DIR/sample_inputs

# Relative names are resolved against this directory, so CASES_FILE can be given
# as a bare table name.
CASES_FILE=${CASES_FILE:-cases.tsv}
[[ $CASES_FILE == /* ]] || CASES_FILE=$SCRIPT_DIR/$CASES_FILE
if [[ ! -f $CASES_FILE ]]; then
    echo "no such case table: $CASES_FILE" >&2
    exit 1
fi

OUT_DIR=$RESULT_ROOT/$LABEL
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

if command -v shasum >/dev/null 2>&1; then
    sha_file() { shasum -a 256 "$1" | cut -d' ' -f1; }
    sha_stdin() { shasum -a 256 | cut -d' ' -f1; }
elif command -v sha256sum >/dev/null 2>&1; then
    sha_file() { sha256sum "$1" | cut -d' ' -f1; }
    sha_stdin() { sha256sum | cut -d' ' -f1; }
else
    echo "need shasum or sha256sum" >&2
    exit 1
fi

if command -v perl >/dev/null 2>&1; then
    now() { perl -MTime::HiRes=time -e 'printf "%.6f\n", time'; }
else
    now() { python3 -c 'import time; print("%.6f" % time.time())'; }
fi

# macOS ships neither GNU timeout nor gtimeout, and at least one sample input
# makes the master build loop forever ("CELL PAIR MAX EXCEEDED"), so the suite
# needs its own watchdog.  A separate killer process is used instead of polling
# so that `wait` still returns the instant the simulation exits and the measured
# wall time stays exact.  A killed run reports status 137.
CASE_TIMEOUT=${CASE_TIMEOUT:-900}
run_with_timeout() {
    local limit=$1
    shift
    "$@" &
    local child=$!

    # The watchdog polls so that it retires by itself shortly after the
    # simulation exits.  Killing it from here would leave its `sleep` orphaned,
    # and that orphan would later signal whatever process inherited the recycled
    # pid.  The poll interval does not affect the measurement: the `wait` below
    # returns the moment the simulation exits.
    (
        local waited=0
        while [[ $waited -lt $limit ]] && kill -0 "$child" 2>/dev/null; do
            sleep 1
            waited=$((waited + 1))
        done
        if kill -0 "$child" 2>/dev/null; then
            kill -9 "$child" 2>/dev/null
        fi
    ) &
    local watchdog=$!

    local status=0
    wait "$child" || status=$?
    wait "$watchdog" 2>/dev/null || true
    return $status
}

{
    printf 'label\t%s\n' "$LABEL"
    printf 'binary\t%s\n' "$BINARY"
    printf 'binary_sha256\t%s\n' "$(sha_file "$BINARY")"
    printf 'seed\t%s\n' "$SEED"
    printf 'repetitions\t%s\n' "$REPS"
    printf 'git_branch\t%s\n' "$(git -C "$ROOT_DIR" rev-parse --abbrev-ref HEAD)"
    printf 'git_commit\t%s\n' "$(git -C "$ROOT_DIR" rev-parse --short HEAD)"
    printf 'host\t%s\n' "$(uname -srm)"
    if [[ "$(uname)" == Darwin ]]; then
        printf 'cpu\t%s\n' "$(sysctl -n machdep.cpu.brand_string)"
        printf 'physical_cpus\t%s\n' "$(sysctl -n hw.physicalcpu)"
    elif command -v lscpu >/dev/null 2>&1; then
        printf 'cpu\t%s\n' "$(lscpu | sed -n 's/^Model name: *//p' | head -1)"
    fi
    printf 'compiler\t%s\n' "$("${CXX:-g++}" --version | head -1)"
} > "$OUT_DIR/metadata.tsv"

printf 'case\trepetition\tnItr\tseconds\texit_status\n' > "$OUT_DIR/timings.tsv"
: > "$OUT_DIR/manifest.sha256"

while IFS=$'\t' read -r id dir parm n_itr extra; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    extra=${extra:-}

    case_dir=$INPUT_ROOT/$dir
    if [[ ! -f "$case_dir/$parm" ]]; then
        echo "SKIP $id: missing $case_dir/$parm" >&2
        continue
    fi

    for rep in $(seq 1 "$REPS"); do
        run_dir=$OUT_DIR/runs/$id/rep$rep
        rm -rf "$run_dir"
        mkdir -p "$run_dir"

        # Copy only the model files the simulation reads; the sample directories
        # also hold notebooks, figures and reference output we must not pick up.
        find "$case_dir" -maxdepth 1 -name '*.mol' -exec cp {} "$run_dir/" \;

        # Force the iteration count so every build runs exactly the same length.
        # SUITE_NITR overrides the table for calibration runs.  A fifth column in
        # the table adds one parameter line, which is how a case switches on a
        # feature the stock model leaves off.
        awk -v n="${SUITE_NITR:-$n_itr}" -v extra="$extra" '
            /^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next }
            /^[[:space:]]*end parameters/    { if (extra != "") print "    " extra; print; next }
            { print }
        ' "$case_dir/$parm" > "$run_dir/$parm"

        start=$(now)
        status=0
        (
            cd "$run_dir"
            run_with_timeout "$CASE_TIMEOUT" "$BINARY" -f "$parm" -s "$SEED" \
                > stdout.log 2> stderr.log
        ) || status=$?
        finish=$(now)

        seconds=$(awk -v a="$start" -v b="$finish" 'BEGIN { printf "%.3f", b - a }')
        printf '%s\t%s\t%s\t%s\t%s\n' "$id" "$rep" "${SUITE_NITR:-$n_itr}" "$seconds" "$status" \
            >> "$OUT_DIR/timings.tsv"

        if [[ $status -ne 0 ]]; then
            echo "FAIL $id rep$rep exit=$status" >&2
        fi

        # Hash outputs from the first repetition only: repetitions exist to
        # average timing noise and every one of them produces the same bytes.
        if [[ $rep -eq 1 ]]; then
            (
                cd "$run_dir"
                find DATA RESTARTS -type f 2>/dev/null | LC_ALL=C sort |
                    while IFS= read -r f; do
                        printf '%s  %s/%s\n' "$(sha_file "$f")" "$id" "$f"
                    done
                find PDB -type f 2>/dev/null | LC_ALL=C sort |
                    while IFS= read -r f; do
                        printf '%s  %s/%s\n' "$(tail -n +2 "$f" | sha_stdin)" "$id" "$f"
                    done
            ) >> "$OUT_DIR/manifest.sha256"
        fi

        printf '%-22s rep%-3s %8ss exit=%s\n' "$id" "$rep" "$seconds" "$status"
    done
done < "$CASES_FILE"

echo
echo "cases:   $CASES_FILE"
echo "results: $OUT_DIR"
