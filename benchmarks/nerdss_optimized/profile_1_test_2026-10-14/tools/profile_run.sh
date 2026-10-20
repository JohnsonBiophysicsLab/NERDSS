#!/bin/zsh
# profile_run.sh <rundir> <binary> <parms> [window_s]
# Runs nerdss under tsrun.py (timestamped stdout) and takes back-to-back
# sample(1) windows until the process exits.
set -u
RUN=$1; BIN=$2; PARMS=$3; WIN=${4:-30}
TOOLS=$(cd "$(dirname "$0")" && pwd)
cd "$RUN" || exit 1
rm -rf DATA PDB RESTARTS samples; mkdir -p samples
/usr/bin/time -l python3 "$TOOLS/tsrun.py" "$RUN/stdout.ts" "$BIN" -f "$PARMS" 2> "$RUN/time.txt" &
WRAP=$!
for i in {1..50}; do [ -f "$RUN/stdout.ts.pid" ] && break; sleep 0.1; done
PID=$(cat "$RUN/stdout.ts.pid")
echo "nerdss pid $PID"
n=0
while kill -0 "$PID" 2>/dev/null; do
  n=$((n+1))
  t=$(python3 -c 'import time; print(f"{time.monotonic():.3f}")')
  sample "$PID" "$WIN" 1 -mayDie -file "$RUN/samples/w$(printf %03d $n).txt" >/dev/null 2>&1
  echo "$n $t" >> "$RUN/samples/windows.txt"
done
wait $WRAP
echo "done: $n windows"
