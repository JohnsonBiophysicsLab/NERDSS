#!/usr/bin/env bash
# run_one.sh <build> <seed>: run the fiber model, record wall time, md5 of the trajectory, overlaps
S=${FIBER_SCRATCH:?set FIBER_SCRATCH to the directory holding build/ and fiber/}
b=$1; seed=$2; d=$S/fiber/runs/$b/seed$seed
rm -rf "$d"; mkdir -p "$d"; cp $S/fiber/model/*.mol $S/fiber/model/parms.inp "$d"/
cd "$d" || exit 1
start=$(date +%s.%N)
$S/build/$b/bin/nerdss -f parms.inp -s $seed > stdout.log 2> stderr.log; echo $? > exit_status
end=$(date +%s.%N)
echo "$end - $start" | bc > wall_seconds
md5 -q DATA/trajectory.xyz > trajectory.md5 2>/dev/null
python3 $S/fiber/scripts/analyze_ab.py "$d" BA 2.0 > overlaps.txt
echo "done $b seed$seed $(cat wall_seconds)s exit=$(cat exit_status)"
