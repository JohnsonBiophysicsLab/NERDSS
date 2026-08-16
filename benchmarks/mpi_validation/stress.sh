#!/bin/bash
# Higher-stress replicate of the standalone-vs-merged comparison.
# usage: stress.sh <label> <binary> <nitr>
# Longer runs raise the merged build's multi-rank failure rate (~37% at nItr=10000
# in earlier probing), so the contrast with standalone is measurable with fewer runs.
export PATH="/opt/homebrew/bin:$PATH"
SP=/private/tmp/claude-501/-Users-yue-Workspace-NERDSS/2decf821-1783-4fbb-b3f6-5b053862dc82/scratchpad
lbl=$1; bin=$2; nitr=$3
OUT=$SP/stress_$lbl.tsv
R=$SP/stress_runs; mkdir -p $R
printf "build\tnp\tseed\tstatus\tfinal_row\n" > $OUT
for s in 6001 6002 6003 6004 6005 6006 6007 6008 6009 6010 6011 6012 6013 6014 6015; do
  d=$SP/stress/${lbl}_$s; rm -rf $d; mkdir -p $d; cp $SP/xin/bi3D/* $d/
  sed -i '' -E "s/nItr = [0-9]+/nItr = $nitr/" $d/parms.inp
  cd $d
  mpirun --oversubscribe -np 2 "$bin" -f parms.inp -s $s > o.txt 2>&1 </dev/null
  f=DATA/copy_numbers_time.dat; [ -f "$f" ] || f=mergeOUT/copy_numbers_time.dat
  h=DATA/histogram_complexes_time.dat; [ -f "$h" ] || h=mergeOUT/histogram_complexes_time.dat
  row=$(tail -1 "$f" 2>/dev/null)
  if grep -qiE 'received signal|exited on signal' o.txt; then st=CRASH
  elif [ -z "$row" ]; then st=NOOUT; else st=ok; fi
  printf "%s\t2\t%s\t%s\t%s\n" "$lbl" "$s" "$st" "$row" >> $OUT
  [ "$st" = ok ] && [ -s "$h" ] && { od=$R/${lbl}__np2__${s}; mkdir -p $od; cp "$h" $od/hist.dat; }
  rm -rf $d/DATA $d/PDB $d/RESTARTS $d/mergeOUT $d/mergePDB 2>/dev/null
done
echo "$lbl done"
