#!/usr/bin/env bash
# With assocDissocWrite on, every rank logs its own associations and
# dissociations to DATA/assoc_dissoc_time_<rank>.dat, molecules named by id.
# Replayed together in step order (BOND adds a bond, BREAK removes it), the
# logs have to make no bond twice, break none that was not made, and end with
# as many bonds as the restart file the run writes at its end.
#
# Before: the reaction passes each opened a shared DATA/assoc_dissoc_time.dat
# themselves, truncating it every pass on every rank, and the per-rank files
# stayed empty (hetTrimer, 2000 steps, np=1: 0 bytes; serial 1525).
#
# Usage: assoc_dissoc_check.sh <nerdss_mpi> [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
BINARY=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
WORK=${2:-$(mktemp -d)}
fail=0

# input directory under sample_inputs, parm file, steps, rank counts
while read -r dir parm steps nps; do
    for np in $nps; do
        d=$WORK/$(basename "$dir")_np$np
        mkdir -p "$d"
        cp "$ROOT_DIR/sample_inputs/$dir"/*.mol "$d/"
        awk -v n="$steps" '/^[[:space:]]*(nItr|assocDissocWrite)[[:space:]]*=/ { next }
             /^[[:space:]]*end parameters/ { print "    nItr = " n; print "    assocDissocWrite = true"; print; next }
             { print }' "$ROOT_DIR/sample_inputs/$dir/$parm" > "$d/$parm"
        (cd "$d" && mpirun --oversubscribe -np "$np" "$BINARY" -f "$parm" -s 1 \
            > stdout.log 2> stderr.log < /dev/null)
        status=$?
        python3 - "$d" "$np" "$status" "$(basename "$dir")" <<'EOF' || fail=1
import collections, json, sys
d, np, status, case = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
events = []
for r in range(np):
    for line in open(f'{d}/DATA/assoc_dissoc_time_{r}.dat'):
        itr, kind, n1, i1, f1, n2, i2, f2 = line.strip().split(',')
        key = frozenset([(n1, int(i1), int(f1)), (n2, int(i2), int(f2))])
        events.append((int(itr.split(':')[1]), kind == 'BOND', key))
events.sort(key=lambda e: (e[0], not e[1]))     # within a step, BOND first
restart = json.load(open(f'{d}/DATA/restart.json'))
mols = restart['molecules']['list']
templates = restart['molTemplates']['templates']
lipid = {t['molName'] for t in templates if t.get('isImplicitLipid')}
live = collections.Counter()
twice = early = 0
for _, bond, key in events:
    # the implicit lipid is one molecule bound many times over
    shared = any(name in lipid for name, _, _ in key)
    if bond:
        if live[key] and not shared: twice += 1
        live[key] += 1
    elif live[key] == 0:
        early += 1
    else:
        live[key] -= 1
final = sum(1 for m in mols if not m.get('isImplicitLipid')
            for f in m['interfaceList'] if f['isBound'])
il = sum(1 for m in mols if not m.get('isImplicitLipid') for f in m['interfaceList']
         if f['isBound'] and mols[f['interaction']['partnerIndex']].get('isImplicitLipid'))
expected = (final - il) // 2 + il
got = sum(live.values())
bad = status != 0 or not events or twice or early or got != expected
print(f"{case} np={np}: {len(events)} events, {got} bonds left, restart has {expected}, "
      f"made twice {twice}, broken before made {early}{' FAIL' if bad else ''}")
sys.exit(1 if bad else 0)
EOF
    done
done <<'CASES'
VALIDATE_SUITE/hetTrimer parm_autodiff_hetTri.inp 2000 1 2 3
VALIDATE_SUITE/implicit_lipid parms.inp 5000 1 2
CASES

[[ $fail == 0 ]] && echo "assoc_dissoc_check: ok"
exit $fail
