#!/usr/bin/env bash
# Timing pass that interleaves the builds instead of running each build's whole
# suite in turn.
#
# run_suite.sh measures one build at a time, which is fine for producing the
# output hashes but leaves the timings exposed to anything else happening on the
# machine: a desktop that becomes busy halfway through penalises whichever build
# happened to be running then.  This script runs case A build 1, case A build 2,
# case A build 3, then repeats, so every build meets the same conditions, and
# reports the median over repetitions.
#
# Usage: interleaved_timing.sh <reps> <label>=<binary> [<label>=<binary> ...]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# ROOT_DIR is overridable so a case table can point at generated models that do
# not live under the repository's own sample_inputs -- scale_cases.tsv does.
ROOT_DIR=${ROOT_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)}

REPS=${1:?usage: interleaved_timing.sh <reps> <label>=<binary> ...}
shift
[[ $# -ge 2 ]] || { echo "need at least two builds to compare" >&2; exit 1; }

LABELS=()
BINARIES=()
for spec in "$@"; do
    LABELS+=("${spec%%=*}")
    path=${spec#*=}
    BINARIES+=("$(cd "$(dirname "$path")" && pwd)/$(basename "$path")")
done

SEED=${SUITE_SEED:-20260810}

# Relative names are resolved against this directory, so CASES_FILE can be given
# as a bare table name, matching run_suite.sh.  This script hardcoded cases.tsv
# until now, which is what cost section 11.6 a measurement.
CASES_FILE=${CASES_FILE:-cases.tsv}
[[ $CASES_FILE == /* ]] || CASES_FILE=$SCRIPT_DIR/$CASES_FILE
if [[ ! -f $CASES_FILE ]]; then
    echo "no such case table: $CASES_FILE" >&2
    exit 1
fi

WORK_DIR=${WORK_DIR:-$SCRIPT_DIR/results/interleaved}
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"

# `time -l` is a BSD flag, and on Apple silicon it is what reports retired
# instructions and elapsed cycles.  GNU time rejects it and some Linux images
# ship no /usr/bin/time at all, either of which would abort this script under
# `set -e`.  Probe once and degrade: -l gives everything, POSIX -p gives CPU
# time, and bare gives wall time only.  The summary below already skips a
# metric that any build is missing.
TIME_PREFIX=()
if /usr/bin/time -l true >/dev/null 2>&1; then
    TIME_PREFIX=(/usr/bin/time -l)
elif /usr/bin/time -p true >/dev/null 2>&1; then
    TIME_PREFIX=(/usr/bin/time -p)
    echo "note: /usr/bin/time -l unavailable; no instruction or cycle counts" >&2
else
    echo "note: no usable /usr/bin/time; wall clock only" >&2
fi

if command -v perl >/dev/null 2>&1; then
    now() { perl -MTime::HiRes=time -e 'printf "%.6f\n", time'; }
else
    now() { python3 -c 'import time; print("%.6f" % time.time())'; }
fi

printf 'case\tnItr\tvariant\trepetition\tcpu_seconds\tinstructions\tcycles\trss_bytes\twall_seconds\n' \
    > "$WORK_DIR/timings.tsv"

while IFS=$'\t' read -r id dir parm n_itr; do
    [[ -z "${id// }" || "${id:0:1}" == "#" ]] && continue
    case_dir=$ROOT_DIR/sample_inputs/$dir
    [[ -f "$case_dir/$parm" ]] || continue

    for rep in $(seq 1 "$REPS"); do
        for i in "${!LABELS[@]}"; do
            run_dir=$WORK_DIR/$id/${LABELS[$i]}/rep$rep
            rm -rf "$run_dir"
            mkdir -p "$run_dir"
            find "$case_dir" -maxdepth 1 -name '*.mol' -exec cp {} "$run_dir/" \;
            awk -v n="${SUITE_NITR:-$n_itr}" \
                '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = " n; next } { print }' \
                "$case_dir/$parm" > "$run_dir/$parm"

            # Five observables, because the first three are only as good as the
            # machine is idle and this one rarely is.  Section 11.6 lost a whole
            # pass to wall clock under load; CPU time is better but still counts
            # a performance core and an efficiency core differently.
            #
            # `time -l` on Apple silicon also reports retired instructions and
            # elapsed cycles, and those are nearly load-proof: measured against
            # three unrelated processes at 99% CPU, instructions varied by 0.1%
            # and cycles by 1.5%, while wall time was inflated threefold.  They
            # also separate the two things a change can do -- fewer instructions
            # means less work, fewer cycles at the same instruction count means
            # fewer stalls -- which is exactly the distinction a layout change
            # turns on.  Both are empty on platforms whose `time -l` omits them.
            start=$(now)
            ( cd "$run_dir" && "${TIME_PREFIX[@]}" "${BINARIES[$i]}" -f "$parm" -s "$SEED" \
                > stdout.log 2> stderr.log )
            finish=$(now)

            # Two layouts: `time -l` puts real/user/sys on one line, POSIX
            # `time -p` puts them on three.  Either may be absent entirely.
            counters=$(awk '
                / real .* user .* sys/ { cpu = $3 + $5; haveCpu = 1 }
                /^user /               { u = $2; havePosix = 1 }
                /^sys /                { s = $2; havePosix = 1 }
                /instructions retired/ { ins = $1 }
                /cycles elapsed/       { cyc = $1 }
                /maximum resident set size/ { rss = $1 }
                END {
                    if (!haveCpu && havePosix) { cpu = u + s; haveCpu = 1 }
                    printf "%s\t%s\t%s\t%s", haveCpu ? sprintf("%.3f", cpu) : "", ins, cyc, rss
                }
            ' "$run_dir/stderr.log")

            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$id" "${SUITE_NITR:-$n_itr}" "${LABELS[$i]}" "$rep" \
                "$counters" \
                "$(awk -v a="$start" -v b="$finish" 'BEGIN { printf "%.3f", b - a }')" \
                >> "$WORK_DIR/timings.tsv"
        done
        printf '%-22s rep%-3s done\n' "$id" "$rep"
    done
done < "$CASES_FILE"

echo
python3 - "$WORK_DIR/timings.tsv" "${LABELS[@]}" <<'PYEOF'
import statistics
import sys

# Four tables, because they answer different questions.  CPU time is what a user
# feels but is only trustworthy on an idle host.  Cycles say the same thing far
# more robustly.  Instructions say whether the change did less *work*, which is
# how a stall reduction is told apart from a work reduction: a layout change
# should move cycles a long way and instructions barely at all.  RSS says
# whether the memory actually went where it was supposed to.
path, *labels = sys.argv[1:]
METRICS = [("cpu_seconds", "s", 3), ("cycles", "cyc", 0),
           ("instructions", "ins", 0), ("rss_bytes", "rss", 0)]

rows = {}
n_itr = {}
with open(path) as handle:
    header = next(handle).rstrip("\n").split("\t")
    col = {name: i for i, name in enumerate(header)}
    for line in handle:
        f = line.rstrip("\n").split("\t")
        case, variant = f[col["case"]], f[col["variant"]]
        n_itr[case] = f[col["nItr"]]
        bucket = rows.setdefault(case, {}).setdefault(variant, {})
        for name, _short, _dp in METRICS:
            raw = f[col[name]] if col[name] < len(f) else ""
            if raw:
                bucket.setdefault(name, []).append(float(raw))

base = labels[0]


def emit(metric, short, dp):
    present = [c for c in sorted(rows)
               if all(rows[c].get(l, {}).get(metric) for l in labels)]
    if not present:
        return
    print("### " + metric)
    # Standard deviations are printed because a ratio without them is not a
    # result: one pass of the size sweep showed a 1.07x that a second pass with
    # more repetitions turned into 0.90x, and only the spread said which to
    # believe.
    head = ["case", "nItr"]
    for l in labels:
        head += ["%s_%s" % (l, short), "%s_sd" % l]
    head += ["%s/%s" % (l, base) for l in labels[1:]]
    print("\t".join(head))

    totals = dict((l, 0.0) for l in labels)
    for case in present:
        med = dict((l, statistics.median(rows[case][l][metric])) for l in labels)
        for l, v in med.items():
            totals[l] += v
        cells = [case, n_itr[case]]
        for l in labels:
            samples = rows[case][l][metric]
            sd = statistics.stdev(samples) if len(samples) > 1 else 0.0
            cells += ["%.*f" % (dp, med[l]), "%.*f" % (dp, sd)]
        cells += ["%.3f" % (med[base] / med[l]) if med[l] else "-"
                  for l in labels[1:]]
        print("\t".join(cells))
    print()
    for l in labels:
        ratio = "\t%.3fx" % (totals[base] / totals[l]) if l != base and totals[l] else ""
        print("total_%s\t%.*f%s" % (l, dp, totals[l], ratio))
    print()


for metric, short, dp in METRICS:
    emit(metric, short, dp)

# Instructions per cycle, which is where a layout change shows up directly.
ipc_cases = [c for c in sorted(rows)
             if all(rows[c].get(l, {}).get("cycles") and rows[c][l].get("instructions")
                    for l in labels)]
if ipc_cases:
    print("### instructions per cycle")
    print("\t".join(["case"] + labels))
    for case in ipc_cases:
        cells = [case]
        for l in labels:
            ins = statistics.median(rows[case][l]["instructions"])
            cyc = statistics.median(rows[case][l]["cycles"])
            cells.append("%.3f" % (ins / cyc) if cyc else "-")
        print("\t".join(cells))
PYEOF
