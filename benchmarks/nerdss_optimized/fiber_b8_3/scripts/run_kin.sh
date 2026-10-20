#!/usr/bin/env bash
# run_kin.sh <build> <seed>: short binding-kinetics run, copy numbers every 5 steps
S=${FIBER_SCRATCH:?set FIBER_SCRATCH to the directory holding build/ and fiber/}
b=$1; seed=$2; d=$S/fiber/kin/$b/seed$seed
rm -rf "$d"; mkdir -p "$d"; cp $S/fiber/model_kin/*.mol $S/fiber/model_kin/parms.inp "$d"/
cd "$d" || exit 1
$S/build/$b/bin/nerdss -f parms.inp -s $seed > stdout.log 2> stderr.log; echo $? > exit_status
echo "done kin $b seed$seed exit=$(cat exit_status)"
