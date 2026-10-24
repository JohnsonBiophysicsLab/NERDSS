#!/usr/bin/env bash
# The reweighting history (Molecule::prevReweight) has to survive the
# exchanges.  A CHECK_STRUCTURE build counts, every step, how many of the pairs
# evaluated found last step's entry, for owned molecules in an edge column
# (0, 1, n-2, n-1: the ones the exchanges overwrite) and in the interior, and
# prints the running sums every 1000 steps ("REWEIGHT HISTORY ...";
# reweight_history_hits() in src/mpi/region_split.cpp).  The two hit rates
# must agree.
#
# Before the deserializers kept the local history, homoTrimer at np=2 over
# 3000 steps found it for 1.7% of edge pairs against 77.4% in the interior
# (np=1: 78.3%); after, 76.7% and 77.8%.
#
# Usage: reweight_history_check.sh <nerdss_mpi built with -DCHECK_STRUCTURE=true> [work-dir]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
BINARY=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
WORK=${2:-$(mktemp -d)}
CASE=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/homoTrimer

rm -rf "$WORK/run"
mkdir -p "$WORK/run"
cp "$CASE"/*.mol "$WORK/run/"
awk '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 3000"; next } { print }' \
    "$CASE/parmTri6.inp" > "$WORK/run/parmTri6.inp"
(cd "$WORK/run" && mpirun --oversubscribe -np 2 "$BINARY" -f parmTri6.inp -s 1 \
    > stdout.log 2> stderr.log < /dev/null)

line=$(grep 'REWEIGHT HISTORY' "$WORK/run/stderr.log" | tail -1)
echo "$line"
# edge=h/t interior=h/t -> pass when the edge hit rate is within 5 points
echo "$line" | awk '{
    split($4, e, "[=/]"); split($5, i, "[=/]")
    if (e[3] == 0 || i[3] == 0) { print "FAIL: no pairs counted"; exit 1 }
    re = e[2] / e[3]; ri = i[2] / i[3]
    printf "edge %.3f, interior %.3f\n", re, ri
    if (re < ri - 0.05) { print "FAIL: edge pairs lose their reweighting history"; exit 1 }
    print "reweight_history_check: ok"
}'
