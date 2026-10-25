#!/usr/bin/env bash
# A 1D (fiber) model has to bind at the serial rate at every rank count.
# fiber_b8_3/model_kin: 50 A and 200 B on a 2000 nm fiber, A + B -> A.B
# irreversible.  The script runs it for 30 steps, over the same seeds, with the
# serial build and with nerdss_mpi at np = 1, 2, 4.  Each run's number of bound
# pairs is averaged over steps 1-30; the mean of that over seeds at each np has
# to agree with the serial mean within 4 combined standard errors, and A and B
# have to be conserved in every row.
#
# Three defects stopped this model under MPI:
#   - The MPI update_memberMolLists() took D.z == 0 to mean "on the membrane".
#     A fiber molecule (D = [D.x, 0, 0]) has that too, so every fiber model
#     exited at set-up ("Molecule 0 of type B is off the membrane"), at any np.
#   - Molecule::isPromoter, Complex::onFiber and MolTemplate::isPromoter were
#     not serialized.  The start-up distribution from rank 0 goes through
#     serialize(), so even np=1 lost them, and a fiber molecule moved and
#     reacted as a 3D one.
#   - RxnBase::area3Dto1D was not serialized.  The 1D rate is
#     rate / area3Dto1D, and every rank had the default -1.
# With only the first one fixed, 15-24 of the 50 A were bound at step 200;
# serial binds all 50 by step ~30.
#
# Usage: fiber_check.sh <nerdss_mpi> <nerdss> [seeds] [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MODEL=$(cd "$SCRIPT_DIR/../fiber_b8_3" && pwd)
MPI_BIN=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
SERIAL_BIN=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
SEEDS=${3:-12}
WORK=${4:-$(mktemp -d)}
fail=0

for run in serial 1 2 4; do
    for s in $(seq 1 "$SEEDS"); do
        d=$WORK/$run/s$s
        mkdir -p "$d"
        cp "$MODEL"/model/*.mol "$d/"
        awk '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 30"; next } { print }' \
            "$MODEL/model_kin/parms.inp" > "$d/parms.inp"
        if [[ $run == serial ]]; then
            (cd "$d" && "$SERIAL_BIN" -f parms.inp -s "$s" > stdout.log 2> stderr.log < /dev/null)
        else
            (cd "$d" && mpirun --oversubscribe -np "$run" "$MPI_BIN" -f parms.inp -s "$s" \
                > stdout.log 2> stderr.log < /dev/null)
        fi
        status=$?
        if [[ $status != 0 ]]; then
            echo "FAIL: $run seed $s exited $status" >&2
            fail=1
        fi
    done
done

python3 - "$WORK" "$SEEDS" <<'EOF' || fail=1
import math, statistics as st, sys
work, seeds = sys.argv[1], int(sys.argv[2])
def per_seed(run):
    out = []
    for s in range(1, seeds + 1):
        f = f'{work}/{run}/s{s}/' + ('DATA' if run == 'serial' else 'mergeOUT') + '/copy_numbers_time.dat'
        try:
            rows = [[float(x) for x in l.split(',')] for l in open(f).read().split('\n')[1:] if l]
        except OSError:
            return None
        if len(rows) < 31 or any(r[1] + r[3] != 200 or r[2] + r[3] != 50 for r in rows):
            print(f'FAIL: {run} seed {s}: {len(rows)} rows, or A/B not conserved')
            return None
        out.append(sum(r[3] for r in rows[1:31]) / 30)
    return out
ref = per_seed('serial')
if ref is None:
    sys.exit(1)
bad = 0
for run in ('1', '2', '4'):
    v = per_seed(run)
    if v is None:
        bad += 1
        continue
    se = math.sqrt(st.variance(v) / len(v) + st.variance(ref) / len(ref))
    ok = abs(st.mean(v) - st.mean(ref)) < 4 * se
    bad += not ok
    print(f"np={run}: {st.mean(v):.2f} bound pairs over steps 1-30, serial {st.mean(ref):.2f}"
          f" (4 SE = {4 * se:.2f}) {'ok' if ok else 'FAIL'}")
sys.exit(1 if bad else 0)
EOF

[[ $fail == 0 ]] && echo "fiber_check: ok"
exit $fail
