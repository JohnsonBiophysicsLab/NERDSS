#!/usr/bin/env bash
# Profile every case of one or more case tables with the gperftools sampling
# profiler, one case at a time, and write a flat and a cumulative pprof table
# per case plus a wall-time table.
#
#   tools/profile_cases.sh <nerdss-binary> <out-dir> <case-table>...
#
# The binary should be built with symbols and frame pointers, e.g.
#   make serial -j4 CFLAGS='-O3 -g -fno-omit-frame-pointer'
# No relinking is needed: libprofiler is preloaded.
#
# Env: SEED (default 20260810), PROFILE_HZ (default 1000), NITR_SCALE (multiply
# each case's nItr, default 1), INPUT_ROOT (default <repo>/sample_inputs),
# LIBPROFILER (default /usr/lib/x86_64-linux-gnu/libprofiler.so), PPROF.
set -euo pipefail
BIN=$(readlink -f "${1:?binary}"); OUT=$(mkdir -p "${2:?out dir}" && cd "$2" && pwd); shift 2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)
INPUT_ROOT=${INPUT_ROOT:-$ROOT/sample_inputs}
SEED=${SEED:-20260810}
LIBPROFILER=${LIBPROFILER:-/usr/lib/x86_64-linux-gnu/libprofiler.so}
PPROF=${PPROF:-$(command -v google-pprof || command -v pprof)}
export CPUPROFILE_FREQUENCY=${PROFILE_HZ:-1000}
mkdir -p "$OUT/runs" "$OUT/profiles"
printf 'case\tnItr\twall_s\texit\tsamples\n' > "$OUT/timings.tsv"
for table in "$@"; do
  while IFS=$'\t' read -r id dir parm nitr extra; do
    [[ -z "${id// }" || ${id:0:1} == "#" ]] && continue
    [[ -f $INPUT_ROOT/$dir/$parm ]] || { echo "SKIP $id (no $dir/$parm)" >&2; continue; }
    n=$(awk -v a="$nitr" -v k="${NITR_SCALE:-1}" 'BEGIN{printf "%d", a*k}')
    run=$OUT/runs/$id; rm -rf "$run"; mkdir -p "$run"
    find "$INPUT_ROOT/$dir" -maxdepth 1 -name '*.mol' -exec cp {} "$run/" \;
    awk -v n="$n" -v extra="${extra:-}" '
      /^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next }
      /^[[:space:]]*end parameters/    { if (extra != "") print "    " extra; print; next }
      { print }' "$INPUT_ROOT/$dir/$parm" > "$run/$parm"
    t0=$(date +%s.%N); st=0
    ( cd "$run" && LD_PRELOAD=$LIBPROFILER CPUPROFILE=$OUT/profiles/$id.prof \
        "$BIN" -f "$parm" -s "$SEED" > stdout.log 2> stderr.log ) || st=$?
    t1=$(date +%s.%N)
    "$PPROF" --text "$BIN" "$OUT/profiles/$id.prof" 2>/dev/null > "$OUT/profiles/$id.flat.txt" || true
    "$PPROF" --text --cum "$BIN" "$OUT/profiles/$id.prof" 2>/dev/null > "$OUT/profiles/$id.cum.txt" || true
    samples=$(sed -n 's/^Total: \([0-9]*\) samples.*/\1/p' "$OUT/profiles/$id.flat.txt")
    printf '%s\t%s\t%.2f\t%s\t%s\n' "$id" "$n" "$(echo "$t1 - $t0" | bc)" "$st" "${samples:-0}" | tee -a "$OUT/timings.tsv"
    rm -rf "$run/DATA" "$run/PDB" "$run/RESTARTS" "$run"/*.xyz
  done < "$table"
done
