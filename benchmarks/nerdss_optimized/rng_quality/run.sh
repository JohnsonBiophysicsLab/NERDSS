#!/usr/bin/env bash
# Build and run the RNG quality comparison for issues #10 and #12.
# Usage: run.sh [orientation-samples] [gauss-samples]
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BIN=$SCRIPT_DIR/rng_quality

"${CXX:-g++}" -O3 -std=c++0x $(gsl-config --cflags) \
    -o "$BIN" "$SCRIPT_DIR/rng_quality.cpp" $(gsl-config --libs)

"$BIN" "${1:-4000000}" "${2:-20000000}"
