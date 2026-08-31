#!/usr/bin/env bash
# Generate constant-density size variants of rev_3D, for scale_cases.tsv.
#
# cases.tsv tops out at 3,955 molecules, and section 22 of RESULTS.md measured
# the layout effect appearing between 6,410 and 10,000.  A suite that cannot
# reach those sizes reports "no effect" for a change worth 1.21x at 40,000, so
# any layout work needs models above the suite's range.
#
# Copy numbers scale by k and the box side by the cube root of k, so the number
# density, the reactions, the diffusion constants and the timestep are all
# unchanged and the only variable is how much memory the molecule array
# occupies.  That is what makes these a size sweep rather than four unrelated
# models.
#
# The generated inputs are not checked in: they are two numbers away from
# rev_3D and regenerating them is cheaper than reviewing them.
#
# Usage: make_scale_inputs.sh [output-root]
#   default output root: results/scale_inputs, which .gitignore already covers
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "$SCRIPT_DIR/../.." && pwd)
SRC=$ROOT_DIR/sample_inputs/VALIDATE_SUITE/bimolecular_reversible/rev_3D
OUT=${1:-$SCRIPT_DIR/results/scale_inputs}

[[ -f $SRC/parms3d.inp ]] || { echo "missing $SRC/parms3d.inp" >&2; exit 1; }

# rev_3D is 1000 A + 1000 R in a 939.993 nm cube.
BASE_SIDE=939.993
BASE_COPIES=1000

mkdir -p "$OUT/sample_inputs"

emit() {
    local name=$1 factor=$2
    local dir=$OUT/sample_inputs/$name
    mkdir -p "$dir"
    find "$SRC" -maxdepth 1 -name '*.mol' -exec cp {} "$dir/" \;

    local copies side
    copies=$((BASE_COPIES * factor))
    side=$(awk -v s="$BASE_SIDE" -v k="$factor" 'BEGIN { printf "%.3f", s * (k ^ (1.0 / 3.0)) }')

    sed -e "s/WaterBox = \[$BASE_SIDE,$BASE_SIDE,$BASE_SIDE\]/WaterBox = [$side,$side,$side]/" \
        -e "s/^\( *\)A : $BASE_COPIES$/\1A : $copies/" \
        -e "s/^\( *\)R : $BASE_COPIES$/\1R : $copies/" \
        "$SRC/parms3d.inp" > "$dir/parms3d.inp"

    # A silent sed miss would produce a model of the wrong size that still runs,
    # which is the worst failure mode for a size sweep.
    grep -q "A : $copies" "$dir/parms3d.inp" || { echo "$name: copy substitution failed" >&2; exit 1; }
    grep -q "WaterBox = \[$side," "$dir/parms3d.inp" || { echo "$name: box substitution failed" >&2; exit 1; }

    printf '%-12s %7d molecules  side %9s nm  Molecule bytes %6.1f MB\n' \
        "$name" "$((copies * 2))" "$side" \
        "$(awk -v n="$((copies * 2))" 'BEGIN { printf "%.1f", n * 656 / 1e6 }')"
}

emit scale_2k 1
emit scale_10k 5
emit scale_40k 20

# enzyme is not generated -- it is a real model, and the largest one in
# sample_inputs that runs.  clathrin_coat/flat_clat-ap2-pip2.dir declares
# 500,400 explicit molecules but never reaches the timestep loop: the call after
# coordinate generation is fixOverlappingMolecules(), whose fallback branch at
# generate_coordinates.cpp:354 is an all-pairs O(N^2) double loop repeated up to
# 50 times.  clathrin_pioneer declares 10,000 pip2 but ships no pip2.mol.
mkdir -p "$OUT/sample_inputs/enzyme"
find "$ROOT_DIR/sample_inputs/enzyme" -maxdepth 1 -name '*.mol' -exec cp {} "$OUT/sample_inputs/enzyme/" \;
cp "$ROOT_DIR/sample_inputs/enzyme/parms_clat_enzyme.inp" "$OUT/sample_inputs/enzyme/"
printf '%-12s %7d molecules  (copied from sample_inputs)\n' enzyme 6410

echo
echo "Run the sweep with:"
echo "  ROOT_DIR=$OUT CASES_FILE=scale_cases.tsv \\"
echo "    $SCRIPT_DIR/interleaved_timing.sh 7 base=<bin> cand=<bin>"
