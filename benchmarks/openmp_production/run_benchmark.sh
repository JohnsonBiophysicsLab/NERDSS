#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
INPUT_DIR="$ROOT_DIR/benchmarks/openmp_production"
INPUT_FILE=${BENCH_INPUT:-parm.inp}
INPUT_PATH="$INPUT_DIR/$INPUT_FILE"
RESULT_ROOT=${1:-"$INPUT_DIR/results/$(date +%Y%m%d-%H%M%S)"}
REPETITIONS=${BENCH_REPS:-3}
THREAD_LIST=${BENCH_THREADS:-"1 2 4 8"}
SEED=${BENCH_SEED:-20260806}

mkdir -p "$RESULT_ROOT"

# Build in separate object directories so switching modes cannot reuse objects
# compiled with a different preprocessor or linker configuration.
make -C "$ROOT_DIR" serial
make -C "$ROOT_DIR" omp

printf 'variant\tthreads\trepetition\twall_seconds\n' > "$RESULT_ROOT/timings.tsv"
printf 'variant\tthreads\trepetition\traw_byte_identical\tsimulation_state_identical\traw_differing_entries\n' \
    > "$RESULT_ROOT/equivalence.tsv"

{
    printf 'git_commit=%s\n' "$(git -C "$ROOT_DIR" rev-parse HEAD)"
    printf 'benchmark_input=%s\n' "$INPUT_FILE"
    printf 'seed=%s\n' "$SEED"
    printf 'repetitions=%s\n' "$REPETITIONS"
    printf 'threads=%s\n' "$THREAD_LIST"
    printf 'OMP_WAIT_POLICY=active\nOMP_PROC_BIND=spread\nOMP_PLACES=cores\n'
    uname -a
    lscpu
    "${CXX:-g++}" --version
} > "$RESULT_ROOT/metadata.txt"

run_case() {
    local variant=$1
    local executable=$2
    local threads=$3
    local repetition=$4
    local run_dir="$RESULT_ROOT/${variant}_t${threads}_r${repetition}"

    mkdir -p "$run_dir"
    cp "$INPUT_PATH" "$run_dir/parm.inp"
    cp "$INPUT_DIR/A.mol" "$run_dir/"

    (
        cd "$run_dir"
        OMP_NUM_THREADS="$threads" \
        OMP_DYNAMIC=FALSE \
        OMP_WAIT_POLICY=active \
        OMP_PROC_BIND=spread \
        OMP_PLACES=cores \
        /usr/bin/time -f '%e' -o wall_seconds.txt \
            "$ROOT_DIR/bin/$executable" -f parm.inp -s "$SEED" \
            > stdout.log 2> stderr.log
    )

    local seconds
    seconds=$(cat "$run_dir/wall_seconds.txt")
    printf '%s\t%s\t%s\t%s\n' \
        "$variant" "$threads" "$repetition" "$seconds" \
        >> "$RESULT_ROOT/timings.tsv"

    (
        cd "$run_dir"
        find DATA PDB RESTARTS -type f -print0 2>/dev/null \
            | sort -z \
            | xargs -0 -r sha256sum
    ) > "$run_dir/output.sha256"

    # PDB line 1 contains the wall-clock creation time.  The normalized
    # manifest hashes its scientific body while all other outputs remain raw.
    (
        cd "$run_dir"
        find DATA RESTARTS -type f -print0 2>/dev/null \
            | sort -z \
            | xargs -0 -r sha256sum
        while IFS= read -r pdb_file; do
            printf '%s  %s\n' \
                "$(tail -n +2 "$pdb_file" | sha256sum | cut -d' ' -f1)" \
                "$pdb_file"
        done < <(find PDB -type f -print 2>/dev/null | sort)
    ) > "$run_dir/simulation_state.sha256"
}

for repetition in $(seq 1 "$REPETITIONS"); do
    run_case serial nerdss 1 "$repetition"
done

for threads in $THREAD_LIST; do
    for repetition in $(seq 1 "$REPETITIONS"); do
        run_case openmp nerdss_omp "$threads" "$repetition"
    done

    for repetition in $(seq 1 "$REPETITIONS"); do
        reference="$RESULT_ROOT/serial_t1_r1"
        candidate="$RESULT_ROOT/openmp_t${threads}_r${repetition}"
        diff_file="$RESULT_ROOT/openmp_t${threads}_r${repetition}_diff.txt"
        if diff -u "$reference/output.sha256" "$candidate/output.sha256" \
            > "$diff_file"; then
            raw_identical=yes
            differing=0
        else
            raw_identical=no
            differing=$(grep -c '^[+-][^+-]' "$diff_file" || true)
        fi
        if cmp -s "$reference/simulation_state.sha256" \
            "$candidate/simulation_state.sha256"; then
            state_identical=yes
        else
            state_identical=no
        fi
        printf 'openmp\t%s\t%s\t%s\t%s\t%s\n' \
            "$threads" "$repetition" "$raw_identical" \
            "$state_identical" "$differing" \
            >> "$RESULT_ROOT/equivalence.tsv"
    done
done

printf '%s\n' "$RESULT_ROOT"
