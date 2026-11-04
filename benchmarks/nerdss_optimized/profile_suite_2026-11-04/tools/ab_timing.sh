#!/usr/bin/env bash
# Interleaved A/B wall-time of two builds on selected cases (A, B, A, B, ...).
#   tools/ab_timing.sh <reps> <binA> <binB> <out.tsv> <case-table-lines-file>
# The case file has run_suite.sh's columns: id dir parm nItr [extra];
# INPUT_ROOT (default <repo>/sample_inputs) can point elsewhere per line with a
# sixth column (put - in the fifth when there is no extra line).  Wall time from date +%s.%N, not the 1 s poll of run_suite.sh.
set -euo pipefail
REPS=$1 A=$(readlink -f "$2") B=$(readlink -f "$3") OUT=$4 TABLE=$5
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)
WORK=$(mktemp -d)
printf 'case\trep\tbuild\twall_s\n' > "$OUT"
while IFS=$'\t' read -r id dir parm nitr extra root; do
  [[ -z "${id// }" || ${id:0:1} == "#" ]] && continue
  [[ ${extra:-} == - ]] && extra=
  src=${root:-$ROOT/sample_inputs}/$dir
  for rep in $(seq 1 "$REPS"); do
    for lab in A B; do
      bin=$A; [[ $lab == B ]] && bin=$B
      run=$WORK/$id.$lab; rm -rf "$run"; mkdir -p "$run"
      find "$src" -maxdepth 1 -name '*.mol' -exec cp {} "$run/" \;
      awk -v n="$nitr" -v extra="${extra:-}" '
        /^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next }
        /^[[:space:]]*end parameters/    { if (extra != "") print "    " extra; print; next }
        { print }' "$src/$parm" > "$run/$parm"
      t0=$(date +%s.%N)
      ( cd "$run" && "$bin" -f "$parm" -s 20260810 > /dev/null 2>&1 )
      t1=$(date +%s.%N)
      printf '%s\t%s\t%s\t%.3f\n' "$id" "$rep" "$lab" "$(echo "$t1 - $t0" | bc)" >> "$OUT"
    done
  done
done < "$TABLE"
rm -rf "$WORK"
