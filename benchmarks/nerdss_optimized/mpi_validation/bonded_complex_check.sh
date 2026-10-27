#!/usr/bin/env bash
# nerdss_mpi has to write the serial build's DATA/COMPLEXES/<step>.json
# (bondedComplexWrite) for the whole system.  write_gathered_bonded_complexes()
# gathers the system the way the checkpoint writer does, so the file it writes
# at a step can be compared with what the serial build's own writer makes of
# the checkpoint written at the same step: the serial build restarted from it
# with nothing left to run (nItr = step + 1) writes DATA/COMPLEXES/<step>.json
# at its end, and the two have to be the same file.
#
# Usage: bonded_complex_check.sh <nerdss_mpi> <nerdss> [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
MPI_BIN=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
SERIAL_BIN=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
WORK=${3:-$(mktemp -d)}
STEP=500
fail=0

# input directory under sample_inputs, parm file, rank counts
while read -r dir parm nps; do
    for np in $nps; do
        d=$WORK/$(basename "$dir")_np$np
        mkdir -p "$d/run" "$d/serial"
        cp "$ROOT_DIR/sample_inputs/$dir"/*.mol "$d/run/"
        awk -v k="$STEP" '/^[[:space:]]*(nItr|checkPoint|restartWrite|bondedComplexWrite)[[:space:]]*=/ { next }
             /^[[:space:]]*end parameters/ { print "    nItr = " 2 * k + 1; print "    checkPoint = " k;
                                             print "    restartWrite = 0"; print "    bondedComplexWrite = " k;
                                             print; next }
             { print }' "$ROOT_DIR/sample_inputs/$dir/$parm" > "$d/run/$parm"
        (cd "$d/run" && mpirun --oversubscribe -np "$np" "$MPI_BIN" -f "$parm" -s 1 \
            > stdout.log 2> stderr.log < /dev/null)
        status=$?
        # the checkpoint, set to end where it starts
        python3 - "$d/run/RESTARTS/restart$STEP.json" "$d/serial/restart.json" "$STEP" <<'EOF' || status=1
import json, sys
j = json.load(open(sys.argv[1]))
j['parameters']['nItr'] = int(sys.argv[3]) + 1
json.dump(j, open(sys.argv[2], 'w'))
EOF
        (cd "$d/serial" && "$SERIAL_BIN" -r restart.json -s 1 > stdout.log 2> stderr.log < /dev/null) || status=1
        mpiFile=$d/run/DATA/COMPLEXES/$STEP.json
        serialFile=$d/serial/DATA/COMPLEXES/$STEP.json
        last=$d/run/DATA/COMPLEXES/$((2 * STEP)).json
        n=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" "$mpiFile" 2>/dev/null)
        if [[ $status == 0 && -n $n && $n -gt 0 && -s $last ]] && cmp -s "$mpiFile" "$serialFile"; then
            echo "$(basename "$dir") np=$np: step $STEP, $n bonded complexes, identical to serial's; step $((2 * STEP)) written"
        else
            echo "FAIL: $(basename "$dir") np=$np (exit $status, ${n:-no} complexes; serial file $( [[ -s $serialFile ]] && echo present || echo missing ), last $( [[ -s $last ]] && echo present || echo missing ))"
            fail=1
        fi
    done
done <<'CASES'
VALIDATE_SUITE/hetTrimer parm_autodiff_hetTri.inp 1 2 3
VALIDATE_SUITE/homoTrimer parmTri6.inp 2 4
CASES

[[ $fail == 0 ]] && echo "bonded_complex_check: ok"
exit $fail
