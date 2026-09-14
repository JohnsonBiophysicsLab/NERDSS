#!/bin/bash
# Feasibility screen over non-spherical sample inputs.
#
# Spherical systems are excluded by reading the parameter file, not by name:
# RefinedGagSphere is not a sphere system and testAdd/sphere is, so the name is
# not a reliable signal.  A system is spherical iff its .inp sets isSphere or
# sphereR.
source /private/tmp/claude-501/-Users-yue-Workspace-NERDSS/2decf821-1783-4fbb-b3f6-5b053862dc82/scratchpad/timeout.sh
export PATH="/opt/homebrew/bin:$PATH"
SP=/private/tmp/claude-501/-Users-yue-Workspace-NERDSS/2decf821-1783-4fbb-b3f6-5b053862dc82/scratchpad
SI=/Users/yue/Workspace/NERDSS/sample_inputs
SC=$SP/screen2; rm -rf $SC; mkdir -p $SC
NITR=2000

is_sphere () { grep -qiE '^[[:space:]]*(isSphere[[:space:]]*=[[:space:]]*true|sphereR[[:space:]]*=)' "$1"; }

CANDIDATES="
$SI/implicit_lipid/parms.inp
$SI/clathrin_coat/flat_clat.dir/clath_kon1uM.inp
$SI/clathrin_coat/flat_clat-ap2-pip2.dir/parms_clath_ap_pip.inp
$SI/clathrin_coat/puckered_clat.dir/parms_pucker.inp
$SI/auto_phos/autophos_D10.inp
$SI/genetic_oscillator/clock_model.inp
$SI/enzyme/parms_clat_enzyme.inp
$SI/pucadyil/parms_pucadyil_ka3d.inp
$SI/gag_coat/solution/parms.inp
$SI/VALIDATE_SUITE/bimolecular_reversible/rev_3D/parms3d.inp
$SI/VALIDATE_SUITE/bimolecular_reversible/rev_2D/parms2D.inp
$SI/VALIDATE_SUITE/bimolecular_reversible/rev_3Dto2D/parms3dto2d.inp
$SI/VALIDATE_SUITE/create_destroy/create.inp
$SI/VALIDATE_SUITE/unimol_state_change_reversible/uni_state_rev.inp
$SI/VALIDATE_SUITE/michaelis_menten/michaelis.inp
$SI/VALIDATE_SUITE/implicit_lipid/parms.inp
$SI/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/IL/parms.inp
$SI/VALIDATE_SUITE/mem_localization/SmallBox/FastDsol/EL/parms.inp
$SI/membrane_rev_localization/SmallBox/FastDsol/IL/parms.inp
$SI/VALIDATE_SUITE/clathrin/parms_clath_kon1uM.inp
$SI/clathrin_invitro_invivo/invitro/parms.inp
$SI/gagsphere/gagOriginalModelSolution.inp
"

printf "%-52s %-8s %-8s %-8s %-7s %s\n" CASE SERIAL MPI_np1 MPI_np2 SEC HIST
i=0
for inp in $CANDIDATES; do
  [ -f "$inp" ] || continue
  if is_sphere "$inp"; then
    printf "%-52s SKIPPED (spherical system per .inp)\n" "$(echo $inp | sed "s|$SI/||")"; continue
  fi
  i=$((i+1)); name=$(printf "c%02d" $i)
  src=$(dirname "$inp"); res=""; sec=""; hist=""
  for mode in serial np1 np2; do
    d=$SC/${name}_$mode; mkdir -p $d
    find "$src" -maxdepth 1 -type f ! -name '*.inp' ! -name '*~' -exec cp {} $d/ \; 2>/dev/null
    sed -E "s/[nN]Itr[[:space:]]*=[[:space:]]*[0-9]+/nItr = $NITR/" "$inp" > $d/parms.inp
    cd $d; s0=$(date +%s)
    case $mode in
      serial) tmo 240 $SP/bench_bins/nerdss_serial -f parms.inp -s 7 > out.txt 2>&1 </dev/null ;;
      np1) tmo 240 mpirun --oversubscribe -np 1 $SP/bench_bins/nerdss_mpi -f parms.inp -s 7 > out.txt 2>&1 </dev/null ;;
      np2) tmo 240 mpirun --oversubscribe -np 2 $SP/bench_bins/nerdss_mpi -f parms.inp -s 7 > out.txt 2>&1 </dev/null ;;
    esac
    rc=$?; s1=$(date +%s)
    [ $mode = serial ] && sec=$((s1-s0))
    f=DATA/copy_numbers_time.dat; [ -f "$f" ] || f=mergeOUT/copy_numbers_time.dat
    h=DATA/histogram_complexes_time.dat; [ -f "$h" ] || h=mergeOUT/histogram_complexes_time.dat
    [ $mode = serial ] && { [ -s "$h" ] && hist=yes || hist=NO; }
    if [ $rc -eq 0 ] && [ -s "$f" ]; then res="$res ok"
    elif grep -qiE 'received signal|exited on signal' out.txt; then res="$res CRASH"
    elif [ $rc -eq 124 ]; then res="$res TMOUT"
    else res="$res fail"; fi
  done
  printf "%-52s %-8s %-8s %-8s %-7s %s\n" "$(echo $inp | sed "s|$SI/||")" $res "$sec" "$hist"
done
