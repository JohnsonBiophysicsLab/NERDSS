#!/usr/bin/env bash
# Runs the deferral model (models/deferral) for one configuration -- serial,
# or nerdss_mpi at a rank count -- at one off rate, one seed after another,
# with assocDissocWrite on, and records the wall seconds and exit status.
# The off rate replaces every offRatekb of the committed parm file, so the
# same model is run at several bond lifetimes; the run directory is
# <out>/kb<kb>/<config>_s<seed>.  deferral_lifetimes.py reads the results.
#
# NITR=<n> in the environment replaces nItr too (longer runs for slower
# dissociation).
#
# Usage: deferral_run.sh <binary> <np|serial> <kb> <out-dir> <seed>...
#   ./deferral_run.sh ../../../bin/nerdss serial 1e5 /tmp/def 301 302
#   ./deferral_run.sh ../../../bin/nerdss_mpi 4 1e5 /tmp/def 301 302
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BINARY=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
CONFIG=$2
KB=$3
OUT=$4
shift 4
MODEL=$SCRIPT_DIR/models/deferral
parm=parm_deferral.inp
label=$([[ $CONFIG == serial ]] && echo serial || echo "np$CONFIG")
mkdir -p "$OUT/kb$KB"
log=$OUT/kb$KB/timing_$label.tsv
[[ -f $log ]] || printf 'config\tkb\tseed\tseconds\texit\n' > "$log"

for seed in "$@"; do
    d=$OUT/kb$KB/${label}_s$seed
    rm -rf "$d"
    mkdir -p "$d"
    cp "$MODEL"/*.mol "$d/"
    sed -E "s/^([[:space:]]*offRatekb[[:space:]]*=[[:space:]]*)[^[:space:]#]+/\1$KB/" \
        "$MODEL/$parm" > "$d/$parm"
    if [[ -n ${NITR:-} ]]; then
        sed -E "s/^([[:space:]]*nItr[[:space:]]*=[[:space:]]*)[^[:space:]#]+/\1$NITR/" "$d/$parm" > "$d/$parm.tmp" \
            && mv "$d/$parm.tmp" "$d/$parm"
    fi
    start=$(date +%s)
    status=0
    if [[ $CONFIG == serial ]]; then
        (cd "$d" && "$BINARY" -f "$parm" -s "$seed" > stdout.log 2> stderr.log < /dev/null) || status=$?
    else
        (cd "$d" && mpirun --oversubscribe -np "$CONFIG" "$BINARY" -f "$parm" -s "$seed" \
            > stdout.log 2> stderr.log < /dev/null) || status=$?
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$KB" "$seed" "$(( $(date +%s) - start ))" "$status" >> "$log"
done
