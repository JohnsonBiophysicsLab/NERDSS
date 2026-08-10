#!/usr/bin/env bash
# Compare two run_suite.sh result sets: per-case bitwise verdict and speedup.
#
# Usage: compare_suites.sh <baseline-label> <candidate-label> [result-root]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BASE=${1:?usage: compare_suites.sh <baseline-label> <candidate-label> [result-root]}
CAND=${2:?missing candidate label}
RESULT_ROOT=${3:-$SCRIPT_DIR/results}

BASE_DIR=$RESULT_ROOT/$BASE
CAND_DIR=$RESULT_ROOT/$CAND

# Median wall time per case, so one slow scheduling hiccup cannot move the
# reported speedup.
median_times() {
    awk -F'\t' 'NR > 1 { t[$1] = t[$1] " " $4 }
        END {
            for (c in t) {
                n = split(t[c], v, " ")
                # v[1] is empty because the accumulated string has a leading space.
                m = 0
                for (i = 1; i <= n; i++) if (v[i] != "") s[++m] = v[i] + 0
                for (i = 2; i <= m; i++) {
                    key = s[i]; j = i - 1
                    while (j > 0 && s[j] > key) { s[j+1] = s[j]; j-- }
                    s[j+1] = key
                }
                printf "%s\t%.3f\n", c, (m % 2) ? s[(m+1)/2] : (s[m/2] + s[m/2+1]) / 2
                delete s
            }
        }' "$1/timings.tsv" | LC_ALL=C sort
}

median_times "$BASE_DIR" > "$RESULT_ROOT/.median_$BASE"
median_times "$CAND_DIR" > "$RESULT_ROOT/.median_$CAND"

n_itr_of() { awk -F'\t' -v c="$1" 'NR > 1 && $1 == c { print $3; exit }' "$BASE_DIR/timings.tsv"; }

printf 'case\tnItr\tbitwise_identical\t%s_s\t%s_s\tspeedup\n' "$BASE" "$CAND"
while IFS=$'\t' read -r case_id base_s; do
    cand_s=$(awk -F'\t' -v c="$case_id" '$1 == c { print $2 }' "$RESULT_ROOT/.median_$CAND")
    [[ -z $cand_s ]] && cand_s=NA

    base_hashes=$(grep -E "^[0-9a-f]+  $case_id/" "$BASE_DIR/manifest.sha256" || true)
    cand_hashes=$(grep -E "^[0-9a-f]+  $case_id/" "$CAND_DIR/manifest.sha256" || true)
    if [[ -z $base_hashes || -z $cand_hashes ]]; then
        verdict=no_output
    elif [[ $base_hashes == "$cand_hashes" ]]; then
        verdict=yes
    else
        verdict=NO
    fi

    if [[ $cand_s == NA ]]; then
        speedup=NA
    else
        speedup=$(awk -v a="$base_s" -v b="$cand_s" \
            'BEGIN { if (b > 0) printf "%.3f", a / b; else printf "NA" }')
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$case_id" "$(n_itr_of "$case_id")" "$verdict" "$base_s" "$cand_s" "$speedup"
done < "$RESULT_ROOT/.median_$BASE"

# Aggregate over the per-case medians rather than over every repetition, so that
# one repetition disturbed by unrelated load on the machine cannot move the
# headline number.  Summing medians also weights cases by how much work they
# represent.
base_total=$(awk -F'\t' '{ s += $2 } END { printf "%.3f", s }' "$RESULT_ROOT/.median_$BASE")
cand_total=$(awk -F'\t' '{ s += $2 } END { printf "%.3f", s }' "$RESULT_ROOT/.median_$CAND")
printf '\ntotal_%s_s (sum of medians)\t%s\n' "$BASE" "$base_total"
printf 'total_%s_s (sum of medians)\t%s\n' "$CAND" "$cand_total"
awk -v a="$base_total" -v b="$cand_total" \
    'BEGIN { if (b > 0) printf "aggregate_speedup\t%.3f\n", a / b }'

rm -f "$RESULT_ROOT/.median_$BASE" "$RESULT_ROOT/.median_$CAND"
