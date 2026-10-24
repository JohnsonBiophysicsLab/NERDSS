#!/usr/bin/env bash
# Positive control for the exchange-conflict class: a model with fast enough
# association and dissociation that both neighbours of a rank boundary act on
# one complex in the same half-step often enough to count.  The two models in
# models/ are the trimer and hetTrimer validation inputs with onRate 1000 and
# offRatekb 1e5 (see README.md for what each found).
#
# Needs a nerdss_mpi built with -DCHECK_STRUCTURE=true (and, for the ownership
# column, -DCHECK_OWNERSHIP=true -DCHECK_OWNERSHIP_EVERY=1).  Runs one seed at
# a time, so it never puts more than <np> ranks on the machine.
#
# Usage: conflict_control.sh <nerdss_mpi> <model-dir> <np> <out-dir> <seed>...
#   ./conflict_control.sh ../../../bin/chk/nerdss_mpi models/het_fast 4 /tmp/cc 1 2 3 4
set -uo pipefail

BINARY=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
MODEL=$(cd "$2" && pwd)
NPROCS=$3
OUT=$4
shift 4
mkdir -p "$OUT"
parm=$(cd "$MODEL" && ls *.inp | head -1)

for seed in "$@"; do
    d=$OUT/$(basename "$MODEL")_np${NPROCS}_s$seed
    rm -rf "$d"
    cp -R "$MODEL" "$d"
    status=0
    (cd "$d" && mpirun --oversubscribe -np "$NPROCS" "$BINARY" -f "$parm" -s "$seed" \
        > stdout.log 2> stderr.log < /dev/null) || status=$?
    printf '%s np=%s seed=%s exit=%s conflicts=%s structure=%s ownership=%s\n' \
        "$(basename "$MODEL")" "$NPROCS" "$seed" "$status" \
        "$(grep -c 'EXCHANGE CONFLICT' "$d/stderr.log")" \
        "$(grep -c 'STRUCTURE VIOLATION' "$d/stderr.log")" \
        "$(grep -c 'OWNERSHIP VIOLATION' "$d/stderr.log")"
done
