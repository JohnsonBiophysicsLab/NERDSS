#!/usr/bin/env bash
# Time one model under several binaries and rank counts, one run at a time,
# and record wall seconds, exit status, the computation / communication split
# each MPI build prints, and the final merged copy-number row.
#
# usage: bench_mpi.sh <result-root> <cases.tsv> <config>...
#   config = label=binary[:np]   (np absent: run the binary directly, serial)
#
# Every case in the table is staged like run_mpi_suite.sh does (the .mol
# files, nItr forced from the table) and the configs are run in the order
# given, so a machine that changes state between runs changes every config's
# timing the same way.  SUITE_SEED (default 20260810) is the seed, REPS
# (default 1) the number of passes over the configs.  Only the final rows and
# logs are kept.  A case dir that does not exist under sample_inputs is looked
# up next to this script (the wide-box models of SCALING_2026-10-30.md live in
# models/wide_box/).  Besides timing.tsv, ranks.tsv keeps every rank's own
# final computation / communication totals, which timing.tsv sums.  Every
# timing row ends with the 1-minute load average before and after the run,
# since on a shared machine a row taken at high load is worthless (the audit
# had to redo ten cells); REP_START (default 1) is the first pass number, so
# a batch can be resumed or interleaved with other work.
set -uo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=${ROOT_DIR:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
RESULT_ROOT=${1:?usage: bench_mpi.sh <result-root> <cases.tsv> <label=binary[:np]>...}
CASES_FILE=${2:?missing case table}
shift 2
CONFIGS=("$@")
SEED=${SUITE_SEED:-20260810}
REPS=${REPS:-1}
MPIRUN=${MPIRUN:-mpirun}
INPUT_ROOT=$ROOT_DIR/sample_inputs
CASE_TIMEOUT=${CASE_TIMEOUT:-1800}
mkdir -p "$RESULT_ROOT"
RES=$RESULT_ROOT/timing.tsv
[[ -f $RES ]] || printf 'case\tnItr\tconfig\trep\tseconds\texit\tcomputation_s\tcommunication_s\tfinal_row\tload_start\tload_end\n' > "$RES"
RANKS=$RESULT_ROOT/ranks.tsv
[[ -f $RANKS ]] || printf 'case\tnItr\tconfig\trep\trank\tcomputation_s\tcommunication_s\n' > "$RANKS"

run_with_timeout() {
    local limit=$1; shift
    "$@" & local child=$!
    ( local w=0; while [[ $w -lt $limit ]] && kill -0 "$child" 2>/dev/null; do sleep 1; w=$((w+1)); done
      if kill -0 "$child" 2>/dev/null; then kill -TERM "$child" 2>/dev/null; sleep 10; kill -9 "$child" 2>/dev/null; fi ) &
    local wd=$!; local st=0
    wait "$child" || st=$?
    wait "$wd" 2>/dev/null || true
    return $st
}

load1() {   # 1-minute load average, Linux or macOS
    if [[ -r /proc/loadavg ]]; then awk '{print $1}' /proc/loadavg
    else sysctl -n vm.loadavg 2>/dev/null | awk '{print $2}'; fi
}

while IFS=$'\t' read -r id dir parm n_itr extra; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    case_dir=$INPUT_ROOT/$dir
    [[ -d "$case_dir" ]] || case_dir=$SCRIPT_DIR/$dir
    [[ -f "$case_dir/$parm" ]] || { echo "SKIP $id: missing $case_dir/$parm" >&2; continue; }
    n=${SUITE_NITR:-$n_itr}
    for rep in $(seq "${REP_START:-1}" "$REPS"); do
        for cfg in "${CONFIGS[@]}"; do
            label=${cfg%%=*}; spec=${cfg#*=}
            bin=${spec%%:*}; np=""
            [[ $spec == *:* ]] && np=${spec##*:}
            run_dir=$RESULT_ROOT/runs/$id/$label/rep$rep
            rm -rf "$run_dir"; mkdir -p "$run_dir"
            find "$case_dir" -maxdepth 1 -name '*.mol' -exec cp {} "$run_dir/" \;
            awk -v n="$n" -v extra="${extra:-}" '
                /^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next }
                /^[[:space:]]*end parameters/    { if (extra != "") print "    " extra; print; next }
                { print }' "$case_dir/$parm" > "$run_dir/$parm"
            load0=$(load1); start=$(date +%s.%N); status=0
            if [[ -n $np ]]; then
                ( cd "$run_dir" && run_with_timeout "$CASE_TIMEOUT" $MPIRUN -np "$np" "$bin" -f "$parm" -s "$SEED" > stdout.log 2> stderr.log < /dev/null ) || status=$?
            else
                ( cd "$run_dir" && run_with_timeout "$CASE_TIMEOUT" "$bin" -f "$parm" -s "$SEED" > stdout.log 2> stderr.log < /dev/null ) || status=$?
            fi
            secs=$(echo "$(date +%s.%N) - $start" | bc); load9=$(load1)
            comp=""; comm=""; final=""
            if [[ -n $np ]]; then
                # 1.4's nerdss_mpi prints each rank's running totals every
                # few thousand steps and once at the end, to the stdout it
                # sends to DATA/output_<rank>.txt; the standalone branch
                # writes its into output_<rank>.txt.  The last line of each
                # file is the rank's total; the ranks are summed.
                comp=$(for f in "$run_dir/stdout.log" "$run_dir"/DATA/output_*.txt "$run_dir"/output_*.txt; do [[ -f $f ]] && grep -h 'total computation time' "$f" | tail -n 1; done | awk '{s+=$4} END{printf "%.1f", s}')
                comm=$(for f in "$run_dir/stdout.log" "$run_dir"/DATA/output_*.txt "$run_dir"/output_*.txt; do [[ -f $f ]] && grep -h 'total communication time' "$f" | tail -n 1; done | awk '{s+=$4} END{printf "%.1f", s}')
                final=$(tail -n 1 "$run_dir/mergeOUT/copy_numbers_time.dat" 2>/dev/null || true)
                for f in "$run_dir"/DATA/output_*.txt "$run_dir"/output_*.txt; do
                    [[ -f $f ]] || continue
                    r=${f##*output_}; r=${r%.txt}
                    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$n" "$label" "$rep" "$r" \
                        "$(grep -h 'total computation time' "$f" | tail -n 1 | awk '{print $4}')" \
                        "$(grep -h 'total communication time' "$f" | tail -n 1 | awk '{print $4}')" >> "$RANKS"
                done
            else
                final=$(tail -n 1 "$run_dir/DATA/copy_numbers_time.dat" 2>/dev/null || true)
            fi
            printf '%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "$n" "$label" "$rep" "$secs" "$status" "$comp" "$comm" "$final" "$load0" "$load9" >> "$RES"
            printf '%-20s %-12s rep%s %8.2fs exit=%s comp=%s comm=%s load=%s/%s %s\n' "$id" "$label" "$rep" "$secs" "$status" "$comp" "$comm" "$load0" "$load9" "${final:0:40}"
            rm -rf "$run_dir/DATA" "$run_dir/PDB" "$run_dir/RESTARTS" "$run_dir/mergePDB" "$run_dir"/*.dat "$run_dir"/*.xyz 2>/dev/null
        done
    done
done < "$CASES_FILE"
echo "results: $RES"
