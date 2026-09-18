#!/usr/bin/env bash
# Regression test: no bond is written between interfaces that are far apart.
#
# Runs the 75-molecule actin model in this directory at the seeds that used to
# expose the defect described in README.md, and checks every frame of
# DATA/COMPLEXES with check_bond_site_separation.py.
#
# Usage: run_test.sh [path-to-nerdss] [work-dir]
#
# Defaults to ../../bin/nerdss and a temporary work directory that is removed on
# success and kept on failure, so a failing run can be inspected.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BINARY=${1:-$SCRIPT_DIR/../../bin/nerdss}
WORK_DIR=${2:-}

# The five seeds below are the ones in 1..60 whose unfixed run first shows a
# far-apart bond before iteration 40000 (parms.inp's nItr), measured on
# nerdss-optimized at 8025fd7.  The defect needs two associations in one timestep
# whose complexes overlap, so most seeds never reach it and a single seed is not
# enough to call this a regression test.  The seeds are only meaningful while the
# trajectories they produce stay the same; ../loop_closure_separation_check.cpp
# pins the criterion itself without depending on any trajectory.
SEEDS=${SEEDS:-"15 21 32 46 54"}
TOLERANCE=${TOLERANCE:-1.15}

if [[ ! -x $BINARY ]]; then
    echo "no nerdss binary at $BINARY -- build one with 'make serial'" >&2
    exit 2
fi
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")

if [[ -z $WORK_DIR ]]; then
    WORK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/bond_site_separation.XXXXXX")
    CLEAN_ON_SUCCESS=1
else
    mkdir -p "$WORK_DIR"
    CLEAN_ON_SUCCESS=0
fi

failures=0
for seed in $SEEDS; do
    run_dir=$WORK_DIR/seed$seed
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    cp "$SCRIPT_DIR/A.mol" "$SCRIPT_DIR/parms.inp" "$run_dir/"

    if ! ( cd "$run_dir" && "$BINARY" -f parms.inp -s "$seed" > run.log 2>&1 ); then
        echo "FAIL seed $seed: nerdss exited nonzero, see $run_dir/run.log" >&2
        failures=$((failures + 1))
        continue
    fi

    # Not a pipeline: the checker's exit status is the verdict, and a pipe would
    # replace it with sed's.
    report=$(python3 "$SCRIPT_DIR/check_bond_site_separation.py" "$run_dir" \
        --tolerance "$TOLERANCE" 2>&1)
    status=$?
    sed "s/^/seed $seed: /" <<< "$report"
    (( status )) && failures=$((failures + 1))
done

if (( failures )); then
    echo "$failures of $(wc -w <<< "$SEEDS") seeds failed; runs kept in $WORK_DIR" >&2
    exit 1
fi

echo "all seeds passed"
(( CLEAN_ON_SUCCESS )) && rm -rf "$WORK_DIR"
exit 0
