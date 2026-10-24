#!/usr/bin/env bash
# Every candidate pair must be listed by exactly one rank.  At step 1 nothing
# has moved and every decomposition starts from rank 0's coordinates, so the
# number of candidate pairs the ranks can act on, summed over ranks, has to be
# the np=1 number.  A CHECK_STRUCTURE build prints that sum at step 1
# ("CANDIDATE PAIRS itr=1 total=N"; count_candidate_pairs() in
# src/mpi/region_split.cpp).
#
# Before measure_separations_to_identify_possible_reactions() paired ghosts
# with owned molecules, homoTrimer listed 1994/2016/1883 of its 2071 pairs at
# np=2/3/4 and hetTrimer 895/887/859 of 929: most pairs across each boundary
# were reached only from the ghost's side of the half-shell cell walk.
#
# Usage: pair_count_check.sh <nerdss_mpi built with -DCHECK_STRUCTURE=true> [work-dir]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
BINARY=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
WORK=${2:-$(mktemp -d)}
fail=0

# id <TAB> input directory under sample_inputs <TAB> parm file <TAB> largest np
while IFS=$'\t' read -r id dir parm most; do
    totals=()
    for np in $(seq 1 "$most"); do
        d=$WORK/${id}_np$np
        rm -rf "$d"
        mkdir -p "$d"
        cp "$ROOT_DIR/sample_inputs/$dir"/*.mol "$d/"
        awk '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 1"; next } { print }' \
            "$ROOT_DIR/sample_inputs/$dir/$parm" > "$d/$parm"
        (cd "$d" && mpirun --oversubscribe -np "$np" "$BINARY" -f "$parm" -s 1 \
            > stdout.log 2> stderr.log < /dev/null) || true
        totals+=("$(sed -n 's/^CANDIDATE PAIRS itr=1 total=//p' "$d/stderr.log")")
    done
    echo "$id: ${totals[*]} (np=1..$most)"
    for t in "${totals[@]}"; do
        if [[ -z $t || $t != "${totals[0]}" ]]; then
            echo "FAIL: $id lists a different set of candidate pairs at some np" >&2
            fail=1
            break
        fi
    done
done <<'EOF'
homoTrimer	VALIDATE_SUITE/homoTrimer	parmTri6.inp	4
hetTrimer	VALIDATE_SUITE/hetTrimer	parm_autodiff_hetTri.inp	4
EOF

[[ $fail == 0 ]] && echo "pair_count_check: ok"
exit $fail
