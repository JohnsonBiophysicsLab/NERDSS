#!/usr/bin/env bash
# Checks the restart files nerdss_mpi writes (write_gathered_restart(): the
# whole system, gathered on rank 0, in the serial format).
#
# 1. Exact continuation at np=1.  hetTrimer for 1000 steps with a checkpoint
#    at 500 and rngwrite; a restart from the checkpoint at np=1 has to write
#    the same copy-number rows as the uninterrupted run, row for row.  This
#    is what caught bndlist being rebuilt in interface order: the step-top
#    dissociation pass draws a molecule's bonds in bndlist order.
# 2. A checkpoint written at np=<np> (default 2) at step 500 has to restart
#    with the serial build and with nerdss_mpi -np 1, and both continuations
#    have to start from the row the MPI run wrote at step 500.
# 3. Its bindPairList has to list every dissociable bond once, by the
#    molecule on the reaction's first product side, which is the rule
#    associate_box() keeps and check_for_unimolecular_reactions_population()
#    relies on.
#
# Usage: checkpoint_check.sh <nerdss_mpi> <nerdss> [np] [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
MPI_BIN=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
SERIAL_BIN=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
NPROCS=${3:-2}
WORK=${4:-$(mktemp -d)}
CASE=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/hetTrimer
PARM=parm_autodiff_hetTri.inp
fail=0

# stage <dir> <timeWrite>
stage() {
    mkdir -p "$1"
    cp "$CASE"/*.mol "$1/"
    awk -v tw="$2" '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 1000"; next }
         /^[[:space:]]*(restartWrite|checkPoint|rngwrite|timeWrite)[[:space:]]*=/ { next }
         /^[[:space:]]*end parameters/ { print "    restartWrite = 0"; print "    checkPoint = 500";
                                         print "    rngwrite = true"; print "    timeWrite = " tw; print; next }
         { print }' "$CASE/$PARM" > "$1/$PARM"
}

# --- 1
mkdir -p "$WORK"
d=$(mktemp -d "$WORK/exact.XXXX")
stage "$d/run" 10
(cd "$d/run" && mpirun -np 1 "$MPI_BIN" -f "$PARM" -s 1 > stdout.log 2> stderr.log < /dev/null)
mkdir -p "$d/cont"
cp "$d/run/RESTARTS/restart500.json" "$d/cont/"
cp "$d/run/RESTARTS/rng_state500" "$d/cont/rng_state_0"
(cd "$d/cont" && mpirun -np 1 "$MPI_BIN" -r restart500.json -s 1 > stdout.log 2> stderr.log < /dev/null)
awk -F, 'NR > 1 && $1 + 0 >= 0.00015' "$d/run/mergeOUT/copy_numbers_time.dat" > "$d/uninterrupted.txt"
tail -n +2 "$d/cont/mergeOUT/copy_numbers_time.dat" > "$d/continued.txt"
if [[ -s $d/continued.txt ]] && cmp -s "$d/uninterrupted.txt" "$d/continued.txt"; then
    echo "1 ok: np=1 continuation identical ($(wc -l < "$d/continued.txt" | tr -d ' ') rows)"
else
    echo "1 FAIL: np=1 continuation differs from the uninterrupted run"
    fail=1
fi

# --- 2
d=$(mktemp -d "$WORK/np$NPROCS.XXXX")
stage "$d/run" 500
(cd "$d/run" && mpirun --oversubscribe -np "$NPROCS" "$MPI_BIN" -f "$PARM" -s 1 > stdout.log 2> stderr.log < /dev/null)
ref=$(awk -F, 'NR > 1 && $1 + 0 == 0.00015' "$d/run/mergeOUT/copy_numbers_time.dat")
for b in serial np1; do
    mkdir -p "$d/$b"
    cp "$d/run/RESTARTS/restart500.json" "$d/$b/"
    if [[ $b == serial ]]; then
        (cd "$d/$b" && "$SERIAL_BIN" -r restart500.json -s 1 > stdout.log 2> stderr.log < /dev/null)
        status=$?; first=$(sed -n 2p "$d/$b/DATA/copy_numbers_time.dat" 2>/dev/null)
    else
        (cd "$d/$b" && mpirun -np 1 "$MPI_BIN" -r restart500.json -s 1 > stdout.log 2> stderr.log < /dev/null)
        status=$?; first=$(sed -n 2p "$d/$b/mergeOUT/copy_numbers_time.dat" 2>/dev/null)
    fi
    if [[ $status == 0 && -n $ref && $first == "$ref" ]]; then
        echo "2 ok: np=$NPROCS checkpoint continues with $b from the run's row at step 500"
    else
        echo "2 FAIL: $b continuation (exit $status) starts at '$first', the run wrote '$ref'"
        fail=1
    fi
done

# --- 3
python3 - "$d/run/RESTARTS/restart500.json" <<'EOF' || fail=1
import json, sys
j = json.load(open(sys.argv[1]))
mols = j['molecules']['list']
first = {}
for r in j['reactions']['forward']:
    if r['rxnType'] == 1 and r['productList']:      # bimolecular
        p = r['productList'][0]
        first[p['absIfaceIndex']] = (p['molTypeIndex'], p['relIfaceIndex'])
bad = 0
for state, listed in enumerate(j['counterArrays']['bindPairList']):
    if state not in first:
        continue
    t, rel = first[state]
    bonds = set()
    for m in listed:
        f = mols[m]['interfaceList'][rel]
        bonds.add(frozenset([m, f['interaction']['partnerIndex']]) if f['isBound'] else ('unbound', m))
    expected = {frozenset([m['index'], f['interaction']['partnerIndex']])
                for m in mols if not m['isImplicitLipid'] and m['molTypeIndex'] == t
                for f in m['interfaceList'] if f['isBound'] and f['index'] == state and f['relIndex'] == rel}
    if bonds != expected or len(listed) != len(bonds):
        bad += 1
print('3 ok: bindPairList lists every dissociable bond once' if not bad
      else f'3 FAIL: bindPairList wrong for {bad} product(s)')
sys.exit(1 if bad else 0)
EOF

[[ $fail == 0 ]] && echo "checkpoint_check: ok"
exit $fail
