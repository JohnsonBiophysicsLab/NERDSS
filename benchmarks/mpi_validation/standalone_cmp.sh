#!/bin/bash
# Is the multi-rank failure a defect of the MPI implementation itself, or something
# the merge introduced?
#
# Controlled comparison: ONE input (the standalone branch's own bimolecular_3D
# benchmark, A + B <-> A.B, so conservation and "no homo-dimer" are both checkable),
# the SAME seeds, run under:
#   merged serial      -- reference; defines the correct answer
#   merged  mpi np=1/2 -- unified tree, branch mpi-serial-build-modes @015789c
#   standalone mpi np=1/2 -- origin/mpi, built from that branch alone
#
# If the standalone build shows the same conservation violations and crashes, the
# defect predates the merge and lives in the MPI implementation.  If it does not,
# the merge introduced it.
export PATH="/opt/homebrew/bin:$PATH"
SP=/private/tmp/claude-501/-Users-yue-Workspace-NERDSS/2decf821-1783-4fbb-b3f6-5b053862dc82/scratchpad
OUT=$SP/standalone_cmp.tsv
R=$SP/scmp_runs; rm -rf $R; mkdir -p $R
NITR=${NITR:-5000}
SEEDS="5001 5002 5003 5004 5005 5006 5007 5008 5009 5010"
printf "build\tnp\tseed\tstatus\tfinal_row\n" > $OUT

run () {  # run <label> <binary> <np>
  local lbl=$1 bin=$2 np=$3
  for s in $SEEDS; do
    local d=$SP/scmp/${lbl}_$s; rm -rf $d; mkdir -p $d; cp $SP/xin/bi3D/* $d/
    sed -i '' -E "s/nItr = [0-9]+/nItr = $NITR/" $d/parms.inp
    cd $d
    if [ "$np" = 0 ]; then "$bin" -f parms.inp -s $s > o.txt 2>&1 </dev/null
    else mpirun --oversubscribe -np $np "$bin" -f parms.inp -s $s > o.txt 2>&1 </dev/null; fi
    local f=DATA/copy_numbers_time.dat; [ -f "$f" ] || f=mergeOUT/copy_numbers_time.dat
    local h=DATA/histogram_complexes_time.dat; [ -f "$h" ] || h=mergeOUT/histogram_complexes_time.dat
    local row=$(tail -1 "$f" 2>/dev/null)
    local st
    if grep -qiE 'received signal|exited on signal' o.txt; then st=CRASH
    elif [ -z "$row" ]; then st=NOOUT; else st=ok; fi
    printf "%s\t%s\t%s\t%s\t%s\n" "$lbl" "$np" "$s" "$st" "$row" >> $OUT
    [ "$st" = ok ] && [ -s "$h" ] && { od=$R/${lbl}__np${np}__${s}; mkdir -p $od; cp "$h" $od/hist.dat; }
    rm -rf $d/DATA $d/PDB $d/RESTARTS $d/mergeOUT $d/mergePDB 2>/dev/null
  done
  echo "  $lbl done"
}

run merged_serial   $SP/bench_bins/nerdss_serial 0
run merged_mpi      $SP/bench_bins/nerdss_mpi    1
run merged_mpi      $SP/bench_bins/nerdss_mpi    2
run standalone_mpi  $SP/bins/standalone_mpi      1
run standalone_mpi  $SP/bins/standalone_mpi      2
echo "=== wrote $OUT"
