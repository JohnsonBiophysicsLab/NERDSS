#!/usr/bin/env bash
# Checks restarting nerdss_mpi on more than one rank from a checkpoint
# written by write_gathered_restart().
#
# 1. hetTrimer at np=2 for 1000 steps with a checkpoint at 500.  Restarting
#    that checkpoint at np=2, 3 and 4 has to start from the row the run wrote
#    at step 500 (the split of the file between the ranks is lossless), has
#    to run to the end with no checker report (build nerdss_mpi with
#    CHECK_STRUCTURE and CHECK_OWNERSHIP), and its final histogram block has
#    to hold every molecule of the file.
# 2. The same for implicit_lipid (free lipid sites split per state), with the
#    free sites plus bound sites summing to nSites in the final block.
# 3. Bound pairs at the end of the hetTrimer continuation, averaged over
#    SEEDS seeds (default 8) at np=2 and np=3, against the serial build's
#    continuation of the same checkpoints: printed as a Welch z.
# 4. The random streams.  A checkpoint written with rngwrite at np=2 saves
#    every rank's generator (RESTARTS/rng_state<itr>_<rank>_np2); restarted
#    at np=2 from that directory, each rank's generator at the top of the
#    first step after the checkpoint has to be the uninterrupted run's at the
#    same step (NERDSS_RNG_PROBE prints a digest of the state bytes to stderr
#    at the top of every step).  Without the state files, and at np=3 (whose
#    ranks the np=2 states cannot serve), the ranks re-seed with seed + rank
#    and the digests differ.
#
# Usage: restart_np_check.sh <nerdss_mpi> <nerdss> [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
MPI_BIN=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
SERIAL_BIN=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
WORK=${3:-$(mktemp -d)}
SEEDS=${SEEDS:-8}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
fail=0
mkdir -p "$WORK"

# stage <dir> <case-dir> <parm> <nItr> <checkPoint> <timeWrite>
stage() {
    mkdir -p "$1"
    cp "$2"/*.mol "$1/"
    awk -v n="$4" -v cp="$5" -v tw="$6" '
         /^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next }
         /^[[:space:]]*(restartWrite|checkPoint|rngwrite|timeWrite|pdbWrite|trajWrite)[[:space:]]*=/ { next }
         /^[[:space:]]*end parameters/ { print "    restartWrite = 0"; print "    checkPoint = " cp;
                                         print "    timeWrite = " tw; print "    pdbWrite = -1"; print "    trajWrite = -1"; print; next }
         { print }' "$2/$3" > "$1/$3"
}

reports() {  # <stderr.log>
    grep -c 'OWNERSHIP VIOLATION\|STRUCTURE VIOLATION\|EXCHANGE CONFLICT\|ERROR: AddressSanitizer\|dropping the member\|unresolved and\|owned by more than one rank' "$1" 2>/dev/null || true
}

# molecules per type summed over the last histogram block: "<type>=<n> ..."
last_block_totals() {  # <histogram file>
    python3 - "$1" <<'PY'
import re, sys
blocks = []
for line in open(sys.argv[1]):
    if line.startswith('Time'):
        blocks.append({}); continue
    if not blocks: continue
    parts = line.rstrip('\n').split('\t', 1)
    if len(parts) != 2: continue
    try: n = int(parts[0])
    except ValueError: continue
    for sp, k in re.findall(r'([A-Za-z_][A-Za-z0-9_]*)\s*:\s*(\d+)\.', parts[1]):
        blocks[-1][sp] = blocks[-1].get(sp, 0) + n * int(k)
print(' '.join(f'{k}={v}' for k, v in sorted(blocks[-1].items())) if blocks else '')
PY
}

# check_case <id> <case-dir> <parm> <nItr> <checkPoint> <timeWrite> <np list>
check_case() {
    local id=$1 cdir=$2 parm=$3 nItr=$4 cp=$5 tw=$6 nps=$7
    local d=$WORK/$id
    rm -rf "$d"; stage "$d/run" "$cdir" "$parm" "$nItr" "$cp" "$tw"
    (cd "$d/run" && $MPIRUN -np 2 "$MPI_BIN" -f "$parm" -s 1 > stdout.log 2> stderr.log < /dev/null)
    local ref
    ref=$(awk -F, -v s="$cp" 'NR > 1 && NR - 2 == s / '"$tw"' { print; exit }' "$d/run/mergeOUT/copy_numbers_time.dat")
    local refTotals
    refTotals=$(last_block_totals "$d/run/mergeOUT/histogram_complexes_time.dat")
    if [[ ! -f $d/run/RESTARTS/restart$cp.json ]]; then
        echo "$id FAIL: the np=2 run wrote no RESTARTS/restart$cp.json ($(reports "$d/run/stderr.log") reports)"; fail=1; return
    fi
    for np in $nps; do
        mkdir -p "$d/np$np"; cp "$d/run/RESTARTS/restart$cp.json" "$d/np$np/"
        (cd "$d/np$np" && $MPIRUN -np "$np" "$MPI_BIN" -r "restart$cp.json" -s 1 > stdout.log 2> stderr.log < /dev/null)
        local status=$? first rep totals
        first=$(sed -n 2p "$d/np$np/mergeOUT/copy_numbers_time.dat" 2>/dev/null)
        rep=$(reports "$d/np$np/stderr.log")
        totals=$(last_block_totals "$d/np$np/mergeOUT/histogram_complexes_time.dat" 2>/dev/null)
        if [[ $status == 0 && -n $ref && $first == "$ref" && $rep == 0 && $totals == "$refTotals" ]]; then
            echo "$id ok: np=$np restart starts from the checkpoint's row, runs clean, last block holds $totals"
        else
            echo "$id FAIL: np=$np restart exit=$status reports=$rep first='$first' expected='$ref' totals='$totals' expected='$refTotals'"
            fail=1
        fi
    done
}

HET=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/hetTrimer
IL=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/implicit_lipid
check_case hetTrimer "$HET" parm_autodiff_hetTri.inp 1000 500 10 "2 3 4"
check_case implicit_lipid "$IL" parms.inp 4000 2000 100 "2 3"

# --- 3: the continuation's end against the serial build's, over seeds
d=$WORK/hetTrimer_seeds
python_welch() {
    python3 - "$@" <<'PY'
import math, sys
def stats(xs):
    n = len(xs); m = sum(xs) / n
    v = sum((x - m) ** 2 for x in xs) / (n - 1) if n > 1 else 0.0
    return m, v, n
groups = {}
for arg in sys.argv[1:]:
    name, vals = arg.split(':', 1)
    groups[name] = [float(v) for v in vals.split(',') if v]
ms, vs, ns = stats(groups['serial'])
for name, xs in groups.items():
    if name == 'serial': continue
    m, v, n = stats(xs)
    se = math.sqrt(vs / ns + v / n) if ns > 1 and n > 1 else float('nan')
    z = (m - ms) / se if se else float('nan')
    print(f"  {name}: bound pairs {m:.2f} (n={n}) vs serial {ms:.2f} (n={ns}), z = {z:+.2f}")
PY
}
rm -rf "$d"; mkdir -p "$d"
serial=""; np2=""; np3=""
for s in $(seq 1 "$SEEDS"); do
    stage "$d/run$s" "$HET" parm_autodiff_hetTri.inp 1000 500 10
    (cd "$d/run$s" && $MPIRUN -np 2 "$MPI_BIN" -f parm_autodiff_hetTri.inp -s "$s" > stdout.log 2> stderr.log < /dev/null)
    for b in serial np2 np3; do
        mkdir -p "$d/run$s/$b"; cp "$d/run$s/RESTARTS/restart500.json" "$d/run$s/$b/"
        if [[ $b == serial ]]; then
            (cd "$d/run$s/$b" && "$SERIAL_BIN" -r restart500.json -s "$((s + 100))" > stdout.log 2> stderr.log < /dev/null)
            v=$(tail -n 1 "$d/run$s/$b/DATA/copy_numbers_time.dat" | awk -F, '{print $NF}')
        else
            (cd "$d/run$s/$b" && $MPIRUN -np "${b#np}" "$MPI_BIN" -r restart500.json -s "$((s + 100))" > stdout.log 2> stderr.log < /dev/null)
            v=$(tail -n 1 "$d/run$s/$b/mergeOUT/copy_numbers_time.dat" | awk -F, '{print $NF}')
        fi
        eval "$b=\"\${$b}\$v,\""
    done
done
echo "3: last copy-number column (bound species) at step 999, $SEEDS seeds each:"
python_welch "serial:$serial" "np2:$np2" "np3:$np3"

# --- 4: every rank's generator continues at the same rank count
d=$WORK/rng_streams
rm -rf "$d"; stage "$d/run" "$HET" parm_autodiff_hetTri.inp 600 500 100
awk '{ print } /^    restartWrite = 0$/ { print "    rngwrite = true" }' "$d/run/parm_autodiff_hetTri.inp" > "$d/run/parm.tmp" && mv "$d/run/parm.tmp" "$d/run/parm_autodiff_hetTri.inp"
(cd "$d/run" && $MPIRUN -x NERDSS_RNG_PROBE=1 -np 2 "$MPI_BIN" -f parm_autodiff_hetTri.inp -s 1 > stdout.log 2> stderr.log < /dev/null)
digests() {  # <stderr.log> <itr>: "rank0:<digest> rank1:<digest> ..."
    grep "^RNG STATE rank [0-9]* itr=$2 " "$1" | sort -t' ' -k4,4n | awk '{printf "rank%s:%s ", $4, $6}'
}
ref=$(digests "$d/run/stderr.log" 501)
if [[ ! -f $d/run/RESTARTS/rng_state500_0_np2 || ! -f $d/run/RESTARTS/rng_state500_1_np2 ]]; then
    echo "4 FAIL: the np=2 run did not write RESTARTS/rng_state500_<rank>_np2"; fail=1
fi
# (a) restarted at np=2 beside the checkpoint's own state files
mkdir -p "$d/same"; cp "$d/run/RESTARTS/restart500.json" "$d/run/RESTARTS/rng_state500_"*"_np2" "$d/same/"
(cd "$d/same" && $MPIRUN -x NERDSS_RNG_PROBE=1 -np 2 "$MPI_BIN" -r restart500.json -s 1 > stdout.log 2> stderr.log < /dev/null)
got=$(digests "$d/same/stderr.log" 501)
reads=$(grep -h -c "RNG state of rank [0-9]* read from " "$d/same/DATA"/output_*.txt | paste -sd+ - | bc)
if [[ -n $ref && $got == "$ref" && $reads == 2 ]]; then
    echo "4a ok: np=2 restart continues both ranks' streams (step 501 digests $got)"
else
    echo "4a FAIL: np=2 restart step-501 digests '$got', uninterrupted '$ref', $reads rank(s) read a state"; fail=1
fi
# (b) restarted at np=2 without the state files: seed + rank, different streams
mkdir -p "$d/none"; cp "$d/run/RESTARTS/restart500.json" "$d/none/"
(cd "$d/none" && $MPIRUN -x NERDSS_RNG_PROBE=1 -np 2 "$MPI_BIN" -r restart500.json -s 1 > stdout.log 2> stderr.log < /dev/null)
got=$(digests "$d/none/stderr.log" 501)
seeded=$(grep -h -c "no saved RNG state for rank [0-9]* of 2; seeding with " "$d/none/DATA"/output_*.txt | paste -sd+ - | bc)
r0=${got%% *}; r1=${got#* }
if [[ $got != "$ref" && $seeded == 2 && ${r0#rank0:} != "${r1#rank1:}" ]]; then
    echo "4b ok: without state files both ranks re-seed (seed + rank) and draw other streams"
else
    echo "4b FAIL: digests '$got' vs uninterrupted '$ref', $seeded rank(s) re-seeded"; fail=1
fi
# (c) restarted at np=3: the np=2 states do not apply, every rank re-seeds
mkdir -p "$d/np3"; cp "$d/run/RESTARTS/restart500.json" "$d/run/RESTARTS/rng_state500_"*"_np2" "$d/np3/"
(cd "$d/np3" && $MPIRUN -x NERDSS_RNG_PROBE=1 -np 3 "$MPI_BIN" -r restart500.json -s 1 > stdout.log 2> stderr.log < /dev/null)
seeded=$(grep -h -c "no saved RNG state for rank [0-9]* of 3; seeding with " "$d/np3/DATA"/output_*.txt | paste -sd+ - | bc)
if [[ $seeded == 3 ]]; then
    echo "4c ok: np=3 restart from the np=2 checkpoint re-seeds all three ranks"
else
    echo "4c FAIL: $seeded of 3 ranks re-seeded"; fail=1
fi

[[ $fail == 0 ]] && echo "restart_np_check: ok"
exit $fail
