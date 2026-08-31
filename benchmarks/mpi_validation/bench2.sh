#!/bin/bash
# Serial vs MPI benchmark over non-spherical sample inputs.
#
# Reads $SP/cases2.txt, one line per case:  name|/abs/path/to/parms.inp|nItr
# Spherical systems are refused here as well as in the screen, so a sphere case
# cannot reach the benchmark by accident.
#
# Records, per run: the final copy_numbers row (species) into a TSV, and the
# whole histogram_complexes_time.dat (assembly pathway) into its own directory.
# The pathway file is the reason the DATA tree is not simply deleted afterwards.
export PATH="/opt/homebrew/bin:$PATH"
SP=/private/tmp/claude-501/-Users-yue-Workspace-NERDSS/2decf821-1783-4fbb-b3f6-5b053862dc82/scratchpad
B=$SP/bench2; rm -rf $B; mkdir -p $B
R=$SP/bench2_runs; rm -rf $R; mkdir -p $R
rm -rf $SP/bench2_headers
OUT=$SP/bench2_species.tsv
SEEDS="3001 3002 3003 3004 3005 3006 3007 3008 3009 3010 3011 3012 3013 3014 3015 3016 3017 3018 3019 3020"

is_sphere () { grep -qiE '^[[:space:]]*(isSphere[[:space:]]*=[[:space:]]*true|sphereR[[:space:]]*=)' "$1"; }

printf "case\tconfig\tseed\tstatus\tfinal_row\n" > $OUT

while IFS='|' read -r name inp nitr; do
  [ -z "$name" ] && continue
  case "$name" in '#'*) continue;; esac
  if is_sphere "$inp"; then echo "=== $name SKIPPED (spherical)"; continue; fi
  src=$(dirname "$inp")
  echo "=== $name (nItr=$nitr)"
  for cfg in serial np1 np2; do
    for s in $SEEDS; do
      d=$B/${name}_${cfg}_$s; mkdir -p $d
      find "$src" -maxdepth 1 -type f ! -name '*.inp' ! -name '*~' -exec cp {} $d/ \; 2>/dev/null
      sed -E "s/[nN]Itr[[:space:]]*=[[:space:]]*[0-9]+/nItr = $nitr/" "$inp" > $d/parms.inp
      cd $d
      case $cfg in
        serial) $SP/bench_bins/nerdss_serial -f parms.inp -s $s > out.txt 2>&1 </dev/null ;;
        np1) mpirun --oversubscribe -np 1 $SP/bench_bins/nerdss_mpi -f parms.inp -s $s > out.txt 2>&1 </dev/null ;;
        np2) mpirun --oversubscribe -np 2 $SP/bench_bins/nerdss_mpi -f parms.inp -s $s > out.txt 2>&1 </dev/null ;;
      esac
      rc=$?
      f=DATA/copy_numbers_time.dat;        [ -f "$f" ] || f=mergeOUT/copy_numbers_time.dat
      h=DATA/histogram_complexes_time.dat; [ -f "$h" ] || h=mergeOUT/histogram_complexes_time.dat
      row=$(tail -1 "$f" 2>/dev/null)
      if grep -qiE 'received signal|exited on signal' out.txt; then st=CRASH
      elif [ -z "$row" ]; then st=NOOUT
      elif [ $rc -ne 0 ]; then st="rc$rc"
      else st=ok; fi
      printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$cfg" "$s" "$st" "$row" >> $OUT
      if [ "$st" = ok ] && [ -s "$h" ]; then
        od=$R/${name}__${cfg}__${s}; mkdir -p $od; cp "$h" $od/hist.dat
      fi
      # species column names, once per case, so the analysis can label columns
      if [ "$st" = ok ] && [ ! -f "$SP/bench2_headers/$name.txt" ]; then
        mkdir -p $SP/bench2_headers; head -1 "$f" > $SP/bench2_headers/$name.txt
      fi
      rm -rf $d/DATA $d/PDB $d/RESTARTS $d/mergeOUT $d/mergePDB 2>/dev/null
    done
    echo "    $cfg done"
  done
done < $SP/cases2.txt
echo "=== species -> $OUT ($(wc -l < $OUT) lines);  pathways -> $R ($(ls $R | wc -l) runs)"
