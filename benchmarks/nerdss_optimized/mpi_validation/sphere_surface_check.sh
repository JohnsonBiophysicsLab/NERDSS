#!/usr/bin/env bash
# Membrane-bound complexes on a sphere have to keep finite coordinates under
# nerdss_mpi.  When an overlap sends a surface complex back for a new move,
# create_complex_propagation_vectors_on_sphere() projects it onto a sphere of
# radius Parameters::sphereR, which nerdss_mpi left at 0: arc = dl / 0, and
# the complex got NaN coordinates.  cluster_gagsphere (gag on a 70 nm sphere,
# clusterOverlapCheck on) at np=1 had 12 of its 200 gag at NaN by step 100.
#
# Runs that model for 200 steps at np=1 and np=2 and reads the restart file
# the run writes at its end (NaN is written to JSON as null).
#
# Usage: sphere_surface_check.sh <nerdss_mpi> [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
BINARY=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
WORK=${2:-$(mktemp -d)}
CASE=$ROOT_DIR/sample_inputs/gagsphere
fail=0

for np in 1 2; do
    d=$WORK/np$np
    mkdir -p "$d"
    cp "$CASE"/*.mol "$d/"
    awk '/^[[:space:]]*(nItr|restartWrite|checkPoint)[[:space:]]*=/ { next }
         /^[[:space:]]*end parameters/ { print "    nItr = 200"; print "    restartWrite = 0";
                                         print "    checkPoint = 0";
                                         print "    clusterOverlapCheck = true"; print; next }
         /^[[:space:]]*clusterOverlapCheck[[:space:]]*=/ { next }
         { print }' "$CASE/parms.inp" > "$d/parms.inp"
    (cd "$d" && mpirun --oversubscribe -np "$np" "$BINARY" -f parms.inp -s 1 \
        > stdout.log 2> stderr.log < /dev/null)
    status=$?
    python3 - "$d/DATA/restart.json" "$np" "$status" <<'EOF' || fail=1
import json, math, sys
path, np, status = sys.argv[1], sys.argv[2], int(sys.argv[3])
try:
    mols = json.load(open(path))['molecules']['list']
except (OSError, ValueError) as e:
    print(f'np={np}: FAIL, exit {status}, no restart file ({e})')
    sys.exit(1)
def crds(m):
    c = m['comCoord']
    return c.values() if isinstance(c, dict) else c
bad = sum(1 for m in mols if not m.get('isImplicitLipid') and
          any(v is None or not math.isfinite(v) for v in crds(m)))
n = sum(1 for m in mols if not m.get('isImplicitLipid'))
print(f"np={np}: {bad} of {n} molecules with non-finite coordinates at step 200"
      f"{'' if bad == 0 and status == 0 else ' FAIL'}")
sys.exit(1 if bad or status else 0)
EOF
done

[[ $fail == 0 ]] && echo "sphere_surface_check: ok"
exit $fail
