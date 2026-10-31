#!/usr/bin/env bash
# Checks merge_outputs(): the per-rank tallies of an np=N run, merged into
# mergeOUT/, against a serial run of the same model (hetTrimer with the A+B
# reaction observed, 1000 steps, rows every 100, transition blocks every 200,
# the association log on).  The values differ -- the runs draw different
# streams -- but
#
# 1. every merged file has the serial file's header, row (or block) count,
#    fields per row and time column, and the row at time 0 is the serial
#    row (the ranks start from one placement);
# 2. every tally column is the sum of the ranks' columns: bound pairs,
#    monomers and dimers, copy numbers, the event counters per label, the
#    transition matrices element by element, with the lifetime lists laid
#    end to end;
# 3. the merged observable AB equals the merged bound-pair column 'A,B' and
#    the merged copy-number column A(a1!1).B(b2!1) on every row: three
#    tallies of the same bonds, kept by three different writers;
# 4. the association log holds every rank's lines, in step order, in the
#    serial format;
# 5. a restart of the np=N checkpoint at np=N keeps 3 (each rank starts
#    from the whole system's observable, which the merge corrects for) and
#    writes the bound-pair header once.
#
# Usage: merge_check.sh <nerdss_mpi> <nerdss> [np] [work-dir]
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
MPI_BIN=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
SERIAL_BIN=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
NPROCS=${3:-2}
WORK=${4:-$(mktemp -d)}
MPIRUN=${MPIRUN:-mpirun --oversubscribe}
CASE=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/hetTrimer
PARM=parm_autodiff_hetTri.inp
fail=0
mkdir -p "$WORK"

stage() {  # <dir>
    mkdir -p "$1"
    cp "$CASE"/*.mol "$1/"
    awk 'BEGIN { labelled = 0 }
         /^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 1000"; next }
         /^[[:space:]]*(restartWrite|checkPoint|rngwrite|timeWrite|transitionWrite|pdbWrite|trajWrite|assocDissocWrite)[[:space:]]*=/ { next }
         /^[[:space:]]*end parameters/ { print "    restartWrite = 0"; print "    checkPoint = 500"; print "    timeWrite = 100";
                                         print "    transitionWrite = 200"; print "    pdbWrite = -1"; print "    trajWrite = -1";
                                         print "    assocDissocWrite = true"; print; next }
         { print }
         /^[[:space:]]*onRate3Dka/ && !labelled { print "    observeLabel = AB"; labelled = 1 }' "$CASE/$PARM" > "$1/$PARM"
}

d=$WORK/np$NPROCS
rm -rf "$d"; stage "$d/serial"; stage "$d/mpi"
(cd "$d/serial" && "$SERIAL_BIN" -f "$PARM" -s 1 > stdout.log 2> stderr.log < /dev/null) || { echo "serial run failed"; exit 1; }
(cd "$d/mpi" && $MPIRUN -np "$NPROCS" "$MPI_BIN" -f "$PARM" -s 1 > stdout.log 2> stderr.log < /dev/null) || { echo "np=$NPROCS run failed"; exit 1; }
mkdir -p "$d/restart"; cp "$d/mpi/RESTARTS/restart500.json" "$d/restart/"
(cd "$d/restart" && $MPIRUN -np "$NPROCS" "$MPI_BIN" -r restart500.json -s 7 > stdout.log 2> stderr.log < /dev/null) || { echo "np=$NPROCS restart failed"; exit 1; }

python3 - "$d" "$NPROCS" <<'PY' || fail=1
import re, sys
d, n = sys.argv[1], int(sys.argv[2])
bad = 0
def fail(msg):
    global bad
    bad += 1
    print("FAIL:", msg)
def lines(path):
    return open(path).read().split('\n')[:-1]
def serial(name): return lines(f"{d}/serial/DATA/{name}.dat")
def merged(name, run="mpi"): return lines(f"{d}/{run}/mergeOUT/{name}.dat")
def ranks(name, run="mpi"): return [lines(f"{d}/{run}/DATA/{name}_{r}.dat") for r in range(n)]
def num(s):
    try: return int(s)
    except ValueError: return float(s)

# --- 1 + 2: row files
cols = {}
for name, sep in (("copy_numbers_time", ","), ("bound_pair_time", "\t"), ("mono_dimer_time", "\t"), ("observables_time", ",")):
    s, m, rk = serial(name), merged(name), ranks(name)
    if not s or not m: fail(f"{name}: empty file (serial {len(s)} lines, merged {len(m)})"); continue
    if m[0] != s[0]: fail(f"{name}: header '{m[0]}' differs from serial '{s[0]}'")
    if len(m) != len(s): fail(f"{name}: {len(m)} lines, serial {len(s)}")
    if any(len(r) != len(m) for r in rk): fail(f"{name}: a rank's file has another length than the merged one")
    table = []
    for i in range(1, min(len(m), len(s))):
        mf, sf = m[i].split(sep), s[i].split(sep)
        if len(mf) != len(sf): fail(f"{name} row {i}: {len(mf)} fields, serial {len(sf)}"); continue
        if mf[0] != sf[0]: fail(f"{name} row {i}: time {mf[0]}, serial {sf[0]}")
        if i == 1 and name != "observables_time" and mf != sf:
            fail(f"{name}: the time-0 row differs from serial's")
        rf = [r[i].split(sep) for r in rk]
        for j in range(1, len(mf)):
            if mf[j] == "" or name == "observables_time": continue
            total = sum(num(f[j]) for f in rf)
            if abs(num(mf[j]) - total) > 1e-9:
                fail(f"{name} row {i} field {j}: merged {mf[j]}, ranks sum to {total}")
        table.append([num(x) if x not in ("",) else None for x in mf])
    cols[name] = (m[0].split(sep), table)
    print(f"{name}: {len(m) - 1} rows, {len(m[1].split(sep))} fields, header and time column as serial, every column the ranks' sum")

# --- 3: three tallies of the A-B bonds
def column(name, header):
    hdr, table = cols[name]
    j = hdr.index(header)
    return [row[j] for row in table]
ab_pairs = column("bound_pair_time", "'A,B'")
ab_copy = column("copy_numbers_time", "A(a1!1).B(b2!1)")
ab_obs = column("observables_time", "AB")
if ab_pairs == ab_copy == ab_obs and ab_obs[-1] > 0:
    print(f"observable AB == bound pairs 'A,B' == copy numbers A(a1!1).B(b2!1) on all {len(ab_obs)} rows (last {ab_obs[-1]})")
else:
    fail(f"the three A-B tallies differ: pairs {ab_pairs}, copy numbers {ab_copy}, observable {ab_obs}")

# --- 2: event counters per label
def blocks(ls, head):
    out = []
    for l in ls:
        if l.startswith(head): out.append((l, {}))
        else:
            k, v = l.rsplit(": ", 1); out[-1][1][k] = int(v)
    return out
s, m, rk = blocks(serial("event_counters_time"), "time (s):"), blocks(merged("event_counters_time"), "time (s):"), [blocks(r, "time (s):") for r in ranks("event_counters_time")]
if [b[0] for b in m] != [b[0] for b in s]: fail("event_counters_time: block headers differ from serial's")
order = lambda k: ({"3D": 0, "3Dto2D": 1, "2D": 2}[k.split()[0]], 0 if "dimers" in k else int(re.search(r"n[=>](\d+)", k).group(1)) * (2 if ">" in k else 1))
for i, (head, vals) in enumerate(m):
    total = {}
    for r in rk:
        for k, v in r[i][1].items(): total[k] = total.get(k, 0) + v
    if vals != total: fail(f"event_counters_time block {head}: merged {vals}, ranks sum to {total}")
    keys = list(vals)
    if keys != sorted(keys, key=order): fail(f"event_counters_time block {head}: labels out of the writer's order: {keys}")
print(f"event_counters_time: {len(m)} blocks as serial, every label the ranks' sum, in the writer's order")

# --- 2: transition matrix
s, m, rk = serial("transition_matrix_time"), merged("transition_matrix_time"), ranks("transition_matrix_time")
text = lambda ls: [l for l in ls if not (l == "" or l.startswith(" "))]
if text(m) != text(s): fail(f"transition_matrix_time: the text lines differ from serial's ({len(text(m))} vs {len(text(s))})")
if len(m) != len(s): fail(f"transition_matrix_time: {len(m)} lines, serial {len(s)}")
lifetimes = False
summed = concatenated = 0
for i, l in enumerate(m):
    if l.startswith("time:"): lifetimes = False; continue
    if l.startswith("lifetime"): lifetimes = True; continue
    if not (l == "" or l.startswith(" ")): continue
    if lifetimes:
        if sorted(l.split()) != sorted(sum((r[i].split() for r in rk), [])): fail(f"transition_matrix_time line {i + 1}: lifetimes are not the ranks' laid end to end")
        concatenated += 1
    else:
        total = [sum(x) for x in zip(*[[int(v) for v in r[i].split()] for r in rk])]
        if [int(v) for v in l.split()] != total: fail(f"transition_matrix_time line {i + 1}: merged row is not the ranks' sum")
        summed += 1
print(f"transition_matrix_time: {m.count('time: 0') and sum(1 for l in m if l.startswith('time:'))} blocks as serial, {summed} matrix rows summed, {concatenated} lifetime rows concatenated")

# --- 4: the association log
s, m, rk = serial("assoc_dissoc_time"), merged("assoc_dissoc_time"), ranks("assoc_dissoc_time")
pat = re.compile(r"^ITR:\d+,(BOND|BREAK),\w+,\d+,\d+,\w+,\d+,\d+$")
if not all(pat.match(l) for l in s): fail("assoc_dissoc_time: the serial lines do not match the expected format (probe is stale)")
if not all(pat.match(l) for l in m): fail("assoc_dissoc_time: a merged line is not in the serial format")
if len(m) != sum(len(r) for r in rk): fail(f"assoc_dissoc_time: {len(m)} merged lines, the ranks wrote {sum(len(r) for r in rk)}")
steps = [int(l.split(",")[0][4:]) for l in m]
if steps != sorted(steps): fail("assoc_dissoc_time: merged lines are not in step order")
if sorted(m) != sorted(sum(rk, [])): fail("assoc_dissoc_time: merged lines are not the ranks' lines")
print(f"assoc_dissoc_time: {len(m)} lines, every rank's, in step order, serial format (serial wrote {len(s)})")

# --- 5: the restart continuation
m = merged("bound_pair_time", "restart")
if sum(1 for l in m if l.startswith("TIME")) != 1 or m[0] != serial("bound_pair_time")[0]: fail(f"restart: bound_pair_time header written {sum(1 for l in m if l.startswith('TIME'))} times or differs from serial's")
# by time: the observables file of a restart has no header and no row at the
# restart step (as the serial build's, which appends to the run's file)
def by_time(name, sep, header):
    ls = merged(name, "restart")
    hdr = ls[0].split(sep) if ls[0].startswith(("Time", "TIME")) else None
    j = hdr.index(header) if hdr else 1
    return {l.split(sep)[0]: num(l.split(sep)[j]) for l in ls[1 if hdr else 0:]}
ab_pairs, ab_copy, ab_obs = by_time("bound_pair_time", "\t", "'A,B'"), by_time("copy_numbers_time", ",", "A(a1!1).B(b2!1)"), by_time("observables_time", ",", "AB")
times = [t for t in ab_obs if t in ab_pairs and t in ab_copy]
if len(times) >= 3 and all(ab_pairs[t] == ab_copy[t] == ab_obs[t] for t in times):
    print(f"restart at np={n}: observable AB == bound pairs == copy numbers on all {len(times)} rows (first {ab_obs[times[0]]}, last {ab_obs[times[-1]]}), one bound-pair header")
else:
    fail(f"restart: the three A-B tallies differ: pairs {ab_pairs}, copy numbers {ab_copy}, observable {ab_obs}")
sys.exit(1 if bad else 0)
PY

[[ $fail == 0 ]] && echo "merge_check: ok"
exit $fail
