#!/usr/bin/env bash
# Run every case in a case table with one nerdss_mpi build at one rank count
# and record, per case: the exit status, the number of OWNERSHIP and STRUCTURE
# violations the build reported on stderr (nonzero only when the build was
# made with CHECK_OWNERSHIP or CHECK_STRUCTURE on), and the final row of the
# merged copy_numbers_time.dat, which is what a mass-conservation check reads.
#
# Usage: run_mpi_suite.sh <path-to-nerdss_mpi> <label> <nprocs> [result-root]
#
# This is the np>1 counterpart of run_suite.sh and stages a case the same way:
# the .mol files and the parameter file with nItr forced from the table, and
# the optional fifth column injected into the parameters block.  Nothing is
# hashed: an MPI run is not expected to be bitwise reproducible, so the
# comparison across builds is on completion, on the violation counts and on
# the copy numbers.  Spherical models are skipped by reading isSphere/sphereR
# out of the parameter file, because the spherical system is known not to
# decompose; see benchmarks/mpi_validation/FINDINGS.md.
#
#   CASES_FILE=cases.tsv ./run_mpi_suite.sh ../../bin/nerdss_mpi fixed 2
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=${ROOT_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)}

BINARY=${1:?usage: run_mpi_suite.sh <path-to-nerdss_mpi> <label> <nprocs> [result-root]}
LABEL=${2:?missing label}
NPROCS=${3:?missing rank count}
RESULT_ROOT=${4:-$SCRIPT_DIR/results}

BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
SEED=${SUITE_SEED:-20260810}
INPUT_ROOT=$ROOT_DIR/sample_inputs
MPIRUN=${MPIRUN:-mpirun --oversubscribe}

CASES_FILE=${CASES_FILE:-cases.tsv}
[[ $CASES_FILE == /* ]] || CASES_FILE=$SCRIPT_DIR/$CASES_FILE
if [[ ! -f $CASES_FILE ]]; then
    echo "no such case table: $CASES_FILE" >&2
    exit 1
fi

OUT_DIR=$RESULT_ROOT/$LABEL
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

{
    printf 'label\t%s\n' "$LABEL"
    printf 'binary\t%s\n' "$BINARY"
    printf 'nprocs\t%s\n' "$NPROCS"
    printf 'seed\t%s\n' "$SEED"
    printf 'git_commit\t%s\n' "$(git -C "$ROOT_DIR" rev-parse --short HEAD)"
} > "$OUT_DIR/metadata.tsv"

printf 'case\tnItr\tseconds\texit_status\townership_violations\tstructure_violations\tfinal_row\n' \
    > "$OUT_DIR/results.tsv"

while IFS=$'\t' read -r id dir parm n_itr extra; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    extra=${extra:-}

    case_dir=$INPUT_ROOT/$dir
    if [[ ! -f "$case_dir/$parm" ]]; then
        echo "SKIP $id: missing $case_dir/$parm" >&2
        continue
    fi
    if grep -qiE '^[[:space:]]*(isSphere[[:space:]]*=[[:space:]]*true|sphereR[[:space:]]*=)' "$case_dir/$parm"; then
        printf '%-22s SKIP spherical\n' "$id"
        continue
    fi

    run_dir=$OUT_DIR/runs/$id
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    find "$case_dir" -maxdepth 1 -name '*.mol' -exec cp {} "$run_dir/" \;
    awk -v n="${SUITE_NITR:-$n_itr}" -v extra="$extra" '
        /^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next }
        /^[[:space:]]*end parameters/    { if (extra != "") print "    " extra; print; next }
        { print }
    ' "$case_dir/$parm" > "$run_dir/$parm"

    start=$(date +%s)
    status=0
    (
        cd "$run_dir"
        $MPIRUN -np "$NPROCS" "$BINARY" -f "$parm" -s "$SEED" > stdout.log 2> stderr.log < /dev/null
    ) || status=$?
    seconds=$(( $(date +%s) - start ))

    owned=$(grep -c 'OWNERSHIP VIOLATION' "$run_dir/stderr.log" || true)
    struct=$(grep -c 'STRUCTURE VIOLATION' "$run_dir/stderr.log" || true)
    final=$(tail -n 1 "$run_dir/mergeOUT/copy_numbers_time.dat" 2>/dev/null || true)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "${SUITE_NITR:-$n_itr}" "$seconds" "$status" \
        "$owned" "$struct" "$final" >> "$OUT_DIR/results.tsv"
    printf '%-22s %6ss exit=%-3s ownership=%-3s structure=%-3s %s\n' \
        "$id" "$seconds" "$status" "$owned" "$struct" "${final:0:60}"
    # The per-rank outputs are large and reproduced by any rerun; keep the
    # merged copy numbers, the logs and the parameter file.
    rm -rf "$run_dir/DATA" "$run_dir/PDB" "$run_dir/RESTARTS" "$run_dir/mergePDB"
done < "$CASES_FILE"

echo
echo "cases:   $CASES_FILE"
echo "results: $OUT_DIR/results.tsv"
