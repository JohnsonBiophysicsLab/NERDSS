#!/usr/bin/env bash
# Check histogramMolStates, the switch that makes DATA/histogram_complexes_time.dat
# tell molecules of one type apart by their interface states.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# Three state-change models are each run three times with one seed: as is
# (plain), with `histogramMolStates = false` (off) and with it true (on).  Then:
#
#   * off == plain on every output file: the explicit default changes nothing.
#   * on == plain on every output file but the histogram, and its restart files
#     differ from plain's only by the parameter: the switch changes the
#     histogram and nothing else.
#   * check_histogram.py: the state-aware histogram collapses onto the plain one
#     frame by frame, every state-bearing type carries its states, and its
#     counts agree with copy_numbers_time.dat.
#   * a restart from the on run's checkpoint keeps the switch: its histogram
#     carries states and equals the on run's at every shared time.
#
# Exit status 0 when every case passes.
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../../.." && pwd)

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
SEED=${SEED:-20260810}
mkdir -p "$WORK"

echo "binary: $BINARY"
echo "work:   $WORK"
echo

# id, model directory, parm file, nItr, timeWrite, checkpoint step, binding
# sites among the state interfaces (space separated Type.iface, or -).
#
# auto_phos: A(ser~U/~P) binds through ser, so the ser columns count monomers.
# unimol_state_change_reversible: A carries two state interfaces (tag, ser),
# which exercises the A(tag~U,ser~P) naming; neither binds.
# enzyme: the motivating model, a lipid pip2(head~U/~P) phosphorylated by syn
# and bound through head by ap2, on a membrane.  The first P appears around
# step 3000 with this seed, so 10000 steps see both states (about 15 s a leg).
CASES='
auto_phos|sample_inputs/auto_phos|autophos_D10.inp|20000|1000|10000|A.ser
unimol_state_change|sample_inputs/VALIDATE_SUITE/unimol_state_change_reversible|uni_state_rev.inp|20000|1000|10000|-
enzyme|sample_inputs/enzyme|parms_clat_enzyme.inp|10000|500|5000|pip2.head
'

# Writes the parameter file for one leg: nItr, timeWrite, checkPoint and RNGwrite
# set for the case, and the switch line added when asked.
write_parms() { # <source> <dest> <nItr> <timeWrite> <checkpoint> <switch line or empty>
    awk -v n="$3" -v tw="$4" -v cp="$5" -v extra="$6" '
        /^[[:space:]]*nItr[[:space:]]*=/        { print "    nItr = " n; next }
        /^[[:space:]]*timeWrite[[:space:]]*=/   { print "    timeWrite = " tw; next }
        /^[[:space:]]*checkPoint[[:space:]]*=/  { print "    checkPoint = " cp; next }
        /^[[:space:]]*end parameters/ {
            print "    RNGwrite = true"
            if (extra != "") print "    " extra
            print; next
        }
        { print }' "$1" > "$2"
}

run_leg() { # <case dir> <leg dir> <parm> <nItr> <timeWrite> <checkpoint> <switch line or empty>
    rm -rf "$2"; mkdir -p "$2"
    cp "$1"/*.mol "$2/"
    write_parms "$1/$3" "$2/$3" "$4" "$5" "$6" "$7"
    (cd "$2" && "$BINARY" -f "$3" -s "$SEED" > stdout.log 2> stderr.log)
}

# Lists the output files of a run directory, relative, sorted; PDB frames carry
# a wall-clock stamp on line 1 and are left out, as run_suite.sh hashes them.
output_files() { (cd "$1" && find DATA RESTARTS -type f 2>/dev/null | grep -v '^PDB' | LC_ALL=C sort); }

# True when every file in both directories is byte-identical, excluding the
# paths matched by the regex in $3 ('' excludes nothing).
dirs_identical() { # <a> <b> <exclude regex>
    local a=$1 b=$2 exclude=$3 ok=0
    diff <(output_files "$a") <(output_files "$b") > /dev/null || { echo "    file lists differ"; return 1; }
    while read -r f; do
        [[ -n $exclude && $f =~ $exclude ]] && continue
        cmp -s "$a/$f" "$b/$f" || { echo "    differs: $f"; ok=1; }
    done < <(output_files "$a")
    return $ok
}

failures=0
pass_case() { printf '%-22s PASS\n' "$1"; }
fail_case() { printf '%-22s FAIL: %s\n' "$1" "$2"; failures=$((failures + 1)); }

while IFS='|' read -r id dir parm n_itr time_write checkpoint sites; do
    [[ -z "${id// }" ]] && continue
    case_dir=$ROOT_DIR/$dir
    case_work=$WORK/$id
    echo "== $id ($dir/$parm, nItr $n_itr)"

    run_leg "$case_dir" "$case_work/plain" "$parm" "$n_itr" "$time_write" "$checkpoint" "" \
        || { fail_case "$id" "plain run failed (see $case_work/plain/stdout.log)"; continue; }
    run_leg "$case_dir" "$case_work/off" "$parm" "$n_itr" "$time_write" "$checkpoint" "histogramMolStates = false" \
        || { fail_case "$id" "off run failed"; continue; }
    run_leg "$case_dir" "$case_work/on" "$parm" "$n_itr" "$time_write" "$checkpoint" "histogramMolStates = true" \
        || { fail_case "$id" "on run failed"; continue; }

    bad=""
    echo "  off vs plain, every output file:"
    dirs_identical "$case_work/plain" "$case_work/off" "" && echo "    identical" || bad="off differs from plain"

    echo "  on vs plain, every output file but the histogram and the restart files:"
    dirs_identical "$case_work/plain" "$case_work/on" 'histogram_complexes_time\.dat$|restart[0-9]*\.json$' \
        && echo "    identical" || bad="${bad:+$bad; }on changes more than the histogram"

    echo "  on vs plain, restart files:"
    for f in $(output_files "$case_work/on" | grep 'restart[0-9]*\.json$'); do
        if python3 - "$case_work/plain/$f" "$case_work/on/$f" <<'PYEOF'
import json, sys
plain = json.load(open(sys.argv[1])); on = json.load(open(sys.argv[2]))
if 'histogramMolStates' in plain['parameters']:
    sys.exit('plain run wrote histogramMolStates')
if on['parameters'].pop('histogramMolStates', None) is not True:
    sys.exit('on run did not write histogramMolStates = true')
if plain != on:
    sys.exit('differs beyond histogramMolStates')
PYEOF
        then echo "    $f: differs only by histogramMolStates"
        else bad="${bad:+$bad; }$f differs beyond the switch"
        fi
    done

    echo "  states histogram vs plain histogram and copy numbers:"
    site_args=()
    [[ $sites != - ]] && for s in $sites; do site_args+=(--binding-site "$s"); done
    python3 "$SCRIPT_DIR/check_histogram.py" \
        --plain "$case_work/plain/DATA/histogram_complexes_time.dat" \
        --states "$case_work/on/DATA/histogram_complexes_time.dat" \
        --copy-numbers "$case_work/on/DATA/copy_numbers_time.dat" ${site_args[@]+"${site_args[@]}"} 2>&1 | sed 's/^/    /'
    [[ ${PIPESTATUS[0]} -eq 0 ]] || bad="${bad:+$bad; }check_histogram.py"

    echo "  restart from the on run's checkpoint:"
    cont=$case_work/on_restart
    rm -rf "$cont"; mkdir -p "$cont"
    if [[ -s "$case_work/on/RESTARTS/restart$checkpoint.json" ]]; then
        cp "$case_work/on/RESTARTS/restart$checkpoint.json" "$cont/restart.json"
        cp "$case_work/on/RESTARTS/rng_state$checkpoint" "$cont/rng_state"
        if (cd "$cont" && "$BINARY" -r restart.json -s 1 > stdout.log 2> stderr.log); then
            if grep -q '~' "$cont/DATA/histogram_complexes_time.dat" \
                && python3 "$SCRIPT_DIR/compare_frames.py" "$case_work/on/DATA/histogram_complexes_time.dat" \
                    "$cont/DATA/histogram_complexes_time.dat" 2>&1 | sed 's/^/    /' \
                && [[ ${PIPESTATUS[0]} -eq 0 ]]; then
                echo "    the restart keeps the switch and reproduces the on run's frames"
            else
                bad="${bad:+$bad; }restart lost the switch or its frames differ"
            fi
        else
            bad="${bad:+$bad; }restart failed (see $cont/stdout.log)"
        fi
    else
        bad="${bad:+$bad; }no RESTARTS/restart$checkpoint.json"
    fi

    if [[ -z $bad ]]; then pass_case "$id"; else fail_case "$id" "$bad"; fi
    echo
done <<< "$CASES"

if [[ $failures -eq 0 ]]; then
    echo "PASS: every case"
else
    echo "FAIL: $failures case(s)"
    exit 1
fi
