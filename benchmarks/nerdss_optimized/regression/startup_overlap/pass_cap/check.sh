#!/usr/bin/env bash
# Check that a startup overlap pass which runs out of iterations says so on
# stderr, and that the run still exits 0.
#
# Usage: check.sh <path-to-nerdss> [work-dir]
#
# P and Q each carry one interface 0.5 nm from the center, and P(p) + Q(q) has
# sigma = 5 nm.  Three runs of that model:
#
#   full   1.5 nm box, no -c.  Two centers are never more than 2.6 nm apart, so
#          two interfaces never more than 3.6 nm, and the overlap pass over all
#          molecules finds the pair on every one of its 50 iterations.
#   coord  3 nm box, P put at the center by -c.  Q's center is never more than
#          2.6 nm from P's, so its interface is never more than 3.6 nm from P's,
#          and the -c pass finds the pair on every iteration.  The full pass
#          before it is not checked here.  (The harness this comes from used a
#          4 nm box, but a Q center in a corner is 3.46 nm from P's, past the
#          3 nm that are_molecules_in_vicinity() looks within, so some seeds
#          escape the overlap there instead of exhausting the pass.)
#   clear  200 nm box, no -c.  The pair is resolved in a pass or two, so no
#          warning may appear.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

BINARY=${1:?usage: check.sh <path-to-nerdss> [work-dir]}
BINARY=$(cd "$(dirname "$BINARY")" && pwd)/$(basename "$BINARY")
WORK=${2:-$(mktemp -d)}
SEED=${SEED:-7}

FULL_WARNING='WARNING: the startup overlap pass over all molecules gave up after 50 iterations'
COORD_WARNING='WARNING: the overlap pass over the molecules placed from the coordinate file (-c) gave up after 50 iterations'

echo "model:  $SCRIPT_DIR"
echo "binary: $BINARY"
echo "seed:   $SEED"
echo "work:   $WORK"
echo

failures=0
fail() {
    echo "  FAIL  $1"
    failures=$((failures + 1))
}
pass() {
    echo "  PASS  $1"
}

# run <name> <box edge in nm> [coordinate file contents]
run() {
    local dir="$WORK/$1"
    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$SCRIPT_DIR"/P.mol "$SCRIPT_DIR"/Q.mol "$dir"
    sed "s/@BOX@/$2/g" "$SCRIPT_DIR"/parms.inp > "$dir"/parms.inp
    local args=(-f parms.inp -s "$SEED")
    if [[ $# -ge 3 ]]; then
        printf '%s\n' "$3" > "$dir"/coords.pdb
        args+=(-c coords.pdb)
    fi
    local status=0
    (cd "$dir" && "$BINARY" "${args[@]}" > stdout.log 2> stderr.log) || status=$?
    echo "$status" > "$dir"/exit_status
}

# passes <name> <prefix>: how many iterations a pass reported ('' full, '|-' -c)
passes() {
    # '|' is a literal in a basic regular expression
    grep -c "^$2Iteration: " "$WORK/$1/stdout.log" || true
}

echo "full: 1.5 nm box, no -c"
run full 1.5
[[ $(cat "$WORK/full/exit_status") == 0 ]] && pass "exits 0" || fail "exits $(cat "$WORK/full/exit_status")"
n=$(passes full '')
[[ $n == 50 ]] && pass "the full pass ran 50 iterations" || fail "the full pass ran $n iterations, expected 50"
grep -qF "$FULL_WARNING" "$WORK/full/stderr.log" && pass "stderr warns about the full pass" \
    || fail "no full-pass warning on stderr"
grep -qF "gave up" "$WORK/full/stdout.log" && fail "the warning went to stdout" || pass "stdout has no warning"

echo "coord: 3 nm box, P at the center from -c"
run coord 3 'ATOM      1 COM  P               0.000   0.000   0.000'
[[ $(cat "$WORK/coord/exit_status") == 0 ]] && pass "exits 0" || fail "exits $(cat "$WORK/coord/exit_status")"
n=$(passes coord '|-')
[[ $n == 50 ]] && pass "the -c pass ran 50 iterations" || fail "the -c pass ran $n iterations, expected 50"
grep -qF "$COORD_WARNING" "$WORK/coord/stderr.log" && pass "stderr warns about the -c pass" \
    || fail "no -c pass warning on stderr"

echo "clear: 200 nm box, no -c"
run clear 200
[[ $(cat "$WORK/clear/exit_status") == 0 ]] && pass "exits 0" || fail "exits $(cat "$WORK/clear/exit_status")"
n=$(passes clear '')
[[ $n -lt 50 ]] && pass "the full pass finished after $n iteration(s)" || fail "the full pass ran $n iterations"
grep -qF "gave up" "$WORK/clear/stderr.log" && fail "a warning appeared although the overlap was resolved" \
    || pass "no warning"

echo
if [[ $failures -ne 0 ]]; then
    echo "FAIL: $failures check(s)" >&2
    exit 1
fi
echo "PASS"
