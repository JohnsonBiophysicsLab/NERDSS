#!/usr/bin/env bash
# run_time.sh <build> <rep>: 100 000-step run, nothing else of ours running, wall time only
S=${FIBER_SCRATCH:?set FIBER_SCRATCH to the directory holding build/ and fiber/}
b=$1; rep=$2; d=$S/fiber/timing/$b/rep$rep
rm -rf "$d"; mkdir -p "$d"; cp $S/fiber/model_time/*.mol $S/fiber/model_time/parms.inp "$d"/
cd "$d" || exit 1
start=$(date +%s.%N); $S/build/$b/bin/nerdss -f parms.inp -s 20260811 > stdout.log 2> stderr.log; echo $? > exit_status; end=$(date +%s.%N)
echo "$end - $start" | bc > wall_seconds; echo "done time $b rep$rep $(cat wall_seconds)s exit=$(cat exit_status)"
