#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
make clean all
mkdir -p results
stamp=$(date -u +%Y%m%dT%H%M%SZ)
{
  echo "# generated_utc=$stamp"
  lscpu | grep -E '^(Model name|CPU\(s\)|Core\(s\) per socket|Thread\(s\) per core):'
  g++ --version | head -n 1
  echo "# OMP_BENCH_THREADS=1,2,4,8,OMP_WAIT_POLICY=active,OMP_PROC_BIND=spread,OMP_PLACES=cores"
  OMP_BENCH_THREADS=1,2,4,8 OMP_WAIT_POLICY=active \
    OMP_PROC_BIND=spread OMP_PLACES=cores \
    ./reaction_propagation_omp_bench "$@"
} | tee "results/benchmark_${stamp}.csv"
