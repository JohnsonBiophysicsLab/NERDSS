#!/usr/bin/env bash
# Statistical equivalence check for the changes that alter the random stream
# (issues #10 and #12), where bitwise comparison is meaningless by construction.
#
# Runs the same model under two builds with several independent seeds, then
# compares the seed-averaged species copy numbers.  If the two builds agree
# within the seed-to-seed scatter, the RNG changes have not moved the physics.
#
# Usage: statistical_check.sh <baseline-binary> <candidate-binary> [seeds] [result-root]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=${ROOT_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)}

BASE_BIN=${1:?usage: statistical_check.sh <baseline-binary> <candidate-binary> [seeds] [result-root]}
CAND_BIN=${2:?missing candidate binary}
SEEDS=${3:-6}
RESULT_ROOT=${4:-$SCRIPT_DIR/results}

BASE_BIN=$(cd "$(dirname "$BASE_BIN")" && pwd)/$(basename "$BASE_BIN")
CAND_BIN=$(cd "$(dirname "$CAND_BIN")" && pwd)/$(basename "$CAND_BIN")

OUT_DIR=$RESULT_ROOT/statistical
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

# Models chosen because they reach a steady state inside a short run and report
# species copy numbers that a rate change would visibly move.
# id | input dir under sample_inputs/ | parm file | nItr
STAT_CASES=(
    "implicit_lipid|VALIDATE_SUITE/implicit_lipid|parms.inp|100000"
    "michaelis_menten|VALIDATE_SUITE/michaelis_menten|michaelis.inp|150000"
    "rev_3D|VALIDATE_SUITE/bimolecular_reversible/rev_3D|parms3d.inp|60000"
)

run_one() {
    local variant=$1 binary=$2 id=$3 dir=$4 parm=$5 n_itr=$6 seed=$7
    local run_dir=$OUT_DIR/$id/$variant/seed$seed
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    find "$ROOT_DIR/sample_inputs/$dir" -maxdepth 1 -name '*.mol' -exec cp {} "$run_dir/" \;
    awk -v n="$n_itr" '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next } { print }' \
        "$ROOT_DIR/sample_inputs/$dir/$parm" > "$run_dir/$parm"
    (
        cd "$run_dir"
        "$binary" -f "$parm" -s "$seed" > stdout.log 2> stderr.log
    )
}

for entry in "${STAT_CASES[@]}"; do
    IFS='|' read -r id dir parm n_itr <<< "$entry"
    echo "=== $id  (nItr=$n_itr, seeds=$SEEDS)"
    for seed in $(seq 1 "$SEEDS"); do
        actual_seed=$((20260810 + seed))
        run_one baseline "$BASE_BIN" "$id" "$dir" "$parm" "$n_itr" "$actual_seed"
        run_one candidate "$CAND_BIN" "$id" "$dir" "$parm" "$n_itr" "$actual_seed"
        printf '  seed %s done\n' "$actual_seed"
    done
    printf '%s\t%s\n' "$id" "$n_itr" >> "$OUT_DIR/cases_run.tsv"
done

python3 "$SCRIPT_DIR/compare_observables.py" "$OUT_DIR" | tee "$OUT_DIR/summary.txt"
