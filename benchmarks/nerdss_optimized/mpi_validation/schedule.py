#!/usr/bin/env python3
"""Run run_mpi_suite.sh jobs under a core budget, strictly in file order.

usage: schedule.py <jobs-file> <result-root> <cores>

jobs file: one job per line, tab separated:
  label  binary  np  seed  casefile  [serial_bin  [result_root]]
casefile is relative to benchmarks/nerdss_optimized unless absolute.  A job
with serial_bin set runs that serial binary through serial_as_mpirun.sh (and
counts as one core); a binary path containing 'asan' gets
ASAN_OPTIONS=detect_leaks=0:halt_on_error=0.  Results go to
<result-root>/<label>; name labels <build>_<table>_np<N>_s<seed> or
serial_<table>_s<seed> for analyze_audit.py.
"""
import os, subprocess, sys, time
jobs_file, root, budget = sys.argv[1], sys.argv[2], int(sys.argv[3])
tools = os.path.dirname(os.path.abspath(__file__))
here = os.path.dirname(tools)                      # benchmarks/nerdss_optimized
# Copy run_mpi_suite.sh somewhere else and point SUITE at the copy if the tree
# may be edited while this runs: bash reads a script as it goes.
suite = os.environ.get('SUITE', os.path.join(here, 'run_mpi_suite.sh'))
jobs = []
for line in open(jobs_file):
    line = line.rstrip('\n')
    if not line or line.startswith('#'): continue
    f = line.split('\t')
    jobs.append(f)
running = []
os.makedirs(root, exist_ok=True)
log = open(os.path.join(root, 'schedule.log'), 'a')
def cost(f): return 1 if len(f) > 5 and f[5] else int(f[2])
while jobs or running:
    still = []
    for (p, f, t) in running:
        if p.poll() is None:
            still.append((p, f, t))
        else:
            log.write(f"{time.strftime('%H:%M:%S')} done {f[0]} rc={p.returncode} {time.time()-t:.0f}s\n"); log.flush()
    running = still
    used = sum(cost(f) for _, f, _ in running)
    started = False
    # strictly in order: a wide job at the head waits for room rather than
    # being overtaken by narrow ones forever
    for f in list(jobs)[:1]:
        if used + cost(f) <= budget or not running:
            env = dict(os.environ)
            env['SUITE_SEED'] = f[3]
            env['CASES_FILE'] = f[4] if f[4].startswith('/') else os.path.join(os.environ.get('CASE_DIR', here), f[4])
            env['ROOT_DIR'] = os.environ.get('ROOT_DIR', os.path.dirname(os.path.dirname(here)))
            if len(f) > 5 and f[5]:
                env['SERIAL_BIN'] = f[5]
                env['MPIRUN'] = os.path.join(tools, 'serial_as_mpirun.sh')
            if 'asan' in f[1]:
                env['ASAN_OPTIONS'] = 'detect_leaks=0:halt_on_error=0'
            jroot = f[6] if len(f) > 6 and f[6] else root
            os.makedirs(jroot, exist_ok=True)
            out = open(os.path.join(jroot, f[0] + '.out'), 'w')
            p = subprocess.Popen([suite, f[1], f[0], f[2], jroot], env=env, stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
            running.append((p, f, time.time()))
            log.write(f"{time.strftime('%H:%M:%S')} start {f[0]}\n"); log.flush()
            jobs.remove(f)
            used += cost(f)
            started = True
            break
    if not started:
        time.sleep(2)
