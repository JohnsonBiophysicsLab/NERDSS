#!/usr/bin/env bash
# Regression check for nerdss_mpi -c (initial coordinates).  The MPI parse path
# used to drop the file name, so a run given -c started from generated
# coordinates and said nothing.
#
# Takes rev_3D, writes a coordinate file that puts every A molecule at
# x = -400 nm (the left edge of the 940 nm box, in the centred frame the -c
# reader expects; NERDSS's own PDB output is in a 0..L frame), and runs it at
# np=2: rank 1 owns the right half of the box, so its time-0 copy-number row
# must hold no free A.  Without the fix both ranks hold about half.
#
# The coordinate file is made from a one-step serial run's PDB frame
# (nerdss_mpi writes CIF), so the serial build is needed too.
#
# Usage: coordinate_file_check.sh <nerdss_mpi> [work-dir]
#   SERIAL_BIN=<nerdss> (default: bin/nerdss of this tree)
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
BINARY=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
SERIAL_BIN=${SERIAL_BIN:-$ROOT_DIR/bin/nerdss}
WORK=${2:-$(mktemp -d)}
CASE=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D

rm -rf "$WORK/gen" "$WORK/run"
mkdir -p "$WORK/gen" "$WORK/run"
cp "$CASE"/*.mol "$WORK/gen/"
awk '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 1"; next }
     /^[[:space:]]*end parameters/    { print "    pdbWrite = 1"; print; next }
     { print }' "$CASE/parms3d.inp" > "$WORK/gen/parms3d.inp"
(cd "$WORK/gen" && "$SERIAL_BIN" -f parms3d.inp -s 1 > stdout.log 2>&1 < /dev/null)

# COM lines only; shift into the centred frame and move every A to x = -400.
awk '{
    if (substr($0, 1, 4) == "ATOM" && substr($0, 13, 4) ~ /COM/) {
        x = substr($0, 31, 8) - 470; y = substr($0, 39, 8) - 470; z = substr($0, 47, 8) - 470
        if (substr($0, 18, 4) ~ /A/) x = -400
        printf "%s%8.1f%8.1f%8.1f%s\n", substr($0, 1, 30), x, y, z, substr($0, 55)
    } else print
}' "$WORK/gen/PDB/0.pdb" > "$WORK/run/left.pdb"

cp "$CASE"/*.mol "$WORK/run/"
awk '/^[[:space:]]*nItr[[:space:]]*=/ { print "    nItr = 2"; next } { print }' \
    "$CASE/parms3d.inp" > "$WORK/run/parms3d.inp"
(cd "$WORK/run" && mpirun --oversubscribe -np 2 "$BINARY" -f parms3d.inp -c left.pdb -s 1 \
    > stdout.log 2> stderr.log < /dev/null)

# Column 2 of copy_numbers_time is free A; row 2 is time 0.
a0=$(sed -n 2p "$WORK/run/DATA/copy_numbers_time_0.dat" | cut -d, -f2)
a1=$(sed -n 2p "$WORK/run/DATA/copy_numbers_time_1.dat" | cut -d, -f2)
echo "free A at time 0: rank 0 $a0, rank 1 $a1"
if [[ $a1 != 0 || $a0 -le 0 ]]; then
    echo "FAIL: the coordinate file was not applied" >&2
    exit 1
fi
echo "coordinate_file_check: ok"
