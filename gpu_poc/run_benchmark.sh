#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
make
mkdir -p results
stamp=$(date -u +%Y%m%dT%H%M%SZ)
{
  echo "# generated_utc=$stamp"
  nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv,noheader
  nvcc --version | tail -n 1
  g++ --version | head -n 1
  ./reaction_propagation_bench "$@"
} | tee "results/benchmark_${stamp}.csv"
