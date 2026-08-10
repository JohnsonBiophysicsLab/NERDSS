#!/usr/bin/env bash
# Timing pass that interleaves the builds instead of running each build's whole
# suite in turn.
#
# run_suite.sh measures one build at a time, which is fine for producing the
# output hashes but leaves the timings exposed to anything else happening on the
# machine: a desktop that becomes busy halfway through penalises whichever build
# happened to be running then.  This script runs case A build 1, case A build 2,
# case A build 3, then repeats, so every build meets the same conditions, and
# reports the median over repetitions.
#
# Usage: interleaved_timing.sh <reps> <label>=<binary> [<label>=<binary> ...]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../.." && pwd)

REPS=${1:?usage: interleaved_timing.sh <reps> <label>=<binary> ...}
shift
[[ $# -ge 2 ]] || { echo "need at least two builds to compare" >&2; exit 1; }

LABELS=()
BINARIES=()
for spec in "$@"; do
    LABELS+=("${spec%%=*}")
    path=${spec#*=}
    BINARIES+=("$(cd "$(dirname "$path")" && pwd)/$(basename "$path")")
done

SEED=${SUITE_SEED:-20260810}
WORK_DIR=${WORK_DIR:-$SCRIPT_DIR/results/interleaved}
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"

if command -v perl >/dev/null 2>&1; then
    now() { perl -MTime::HiRes=time -e 'printf "%.6f\n", time'; }
else
    now() { python3 -c 'import time; print("%.6f" % time.time())'; }
fi

printf 'case\tnItr\tvariant\trepetition\tseconds\n' > "$WORK_DIR/timings.tsv"

while IFS=$'\t' read -r id dir parm n_itr; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    case_dir=$ROOT_DIR/sample_inputs/$dir
    [[ -f "$case_dir/$parm" ]] || continue

    for rep in $(seq 1 "$REPS"); do
        for i in "${!LABELS[@]}"; do
            run_dir=$WORK_DIR/$id/${LABELS[$i]}/rep$rep
            rm -rf "$run_dir"
            mkdir -p "$run_dir"
            find "$case_dir" -maxdepth 1 -name '*.mol' -exec cp {} "$run_dir/" \;
            awk -v n="${SUITE_NITR:-$n_itr}" \
                '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next } { print }' \
                "$case_dir/$parm" > "$run_dir/$parm"

            start=$(now)
            ( cd "$run_dir" && "${BINARIES[$i]}" -f "$parm" -s "$SEED" > stdout.log 2> stderr.log )
            finish=$(now)

            printf '%s\t%s\t%s\t%s\t%s\n' "$id" "${SUITE_NITR:-$n_itr}" "${LABELS[$i]}" "$rep" \
                "$(awk -v a="$start" -v b="$finish" 'BEGIN { printf "%.3f", b - a }')" \
                >> "$WORK_DIR/timings.tsv"
        done
        printf '%-22s rep%-3s done\n' "$id" "$rep"
    done
done < "$SCRIPT_DIR/cases.tsv"

echo
python3 - "$WORK_DIR/timings.tsv" "${LABELS[@]}" <<'PY'
import statistics
import sys

path, *labels = sys.argv[1:]
rows = {}
n_itr = {}
with open(path) as handle:
    next(handle)
    for line in handle:
        case, iters, variant, _rep, seconds = line.rstrip("\n").split("\t")
        rows.setdefault(case, {}).setdefault(variant, []).append(float(seconds))
        n_itr[case] = iters

base = labels[0]
header = ["case", "nItr"] + [f"{l}_s" for l in labels] + [f"{l}/{base}" for l in labels[1:]]
print("\t".join(header))

totals = {l: 0.0 for l in labels}
for case in sorted(rows):
    med = {l: statistics.median(rows[case][l]) for l in labels if l in rows[case]}
    for l, v in med.items():
        totals[l] += v
    cells = [case, n_itr[case]] + [f"{med[l]:.3f}" for l in labels]
    cells += [f"{med[base] / med[l]:.3f}" for l in labels[1:]]
    print("\t".join(cells))

print()
for l in labels:
    speed = f"\tspeedup {totals[base] / totals[l]:.3f}x" if l != base else ""
    print(f"total_{l}_s\t{totals[l]:.3f}{speed}")
PY
