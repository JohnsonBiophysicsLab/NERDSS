#!/usr/bin/env python3
"""Markdown tables comparing two audits made with schedule.py.

usage: audit_tables.py <tip-root> <branch-root>

Both roots hold run_mpi_suite.sh label directories as analyze_audit.py reads
them, including the serial reference (serial_* labels; symlink them into the
second root).  Prints two tables:

1. per case and np: runs, failed runs, and the summed ownership / structure /
   ASan / exchange-conflict reports and mass-conservation failures, for both
   audits (checker, -O3 and ASan runs together);
2. per case and np: the largest |t| over copy-number columns, second half of
   the run (equilibrium) and first half (relaxation), MPI against serial,
   full-length runs only, for both audits.
"""
import os, subprocess, sys
from collections import defaultdict

here = os.path.dirname(os.path.abspath(__file__))

def rows(root):
    out = subprocess.run([sys.executable, os.path.join(here, 'analyze_audit.py'), '--tsv', root],
                         capture_output=True, text=True).stdout.splitlines()
    res = {}
    for line in out:
        f = line.split('\t')
        if len(f) < 14 or f[1] == 'serial':
            continue
        res[(f[0], f[1])] = dict(runs=int(f[2]), fail=int(f[3]), own=int(f[4]), struct=int(f[5]),
                                 asan=int(f[6]), confl=int(f[7]), mass=f[9], t1=f[10], t2=f[11], p=f[12])
    return res

tip, br = rows(sys.argv[1]), rows(sys.argv[2])
cases = sorted({c for c, _ in list(tip) + list(br)})

def counts(d, case, np):
    tot = defaultdict(int); mass = []
    for cfg in (f'np{np}', f'np{np}-asan'):
        r = d.get((case, cfg))
        if not r: continue
        for k in ('runs', 'fail', 'own', 'struct', 'asan', 'confl'): tot[k] += r[k]
        mass.append(r['mass'])
    if not tot: return '-'
    bad = [m for m in mass if m.endswith('bad')]
    s = f"{tot['runs']}/{tot['fail']}"
    rep = [f"{tot[k]} {k}" for k in ('own', 'struct', 'asan', 'confl') if tot[k]]
    if bad: rep.append('mass ' + ','.join(bad))
    return s + (' ' + ', '.join(rep) if rep else '')

print('| case | np | tip: runs/failed, reports | branch: runs/failed, reports |')
print('|---|---|---|---|')
for case in cases:
    for np in (2, 3, 4):
        a, b = counts(tip, case, np), counts(br, case, np)
        if a == '-' and b == '-': continue
        print(f'| {case} | {np} | {a} | {b} |')

print()
print('| case | np | tip: max\\|t\\| 2nd half (p) / 1st half | branch: max\\|t\\| 2nd half (p) / 1st half |')
print('|---|---|---|---|')
def stat(d, case, np):
    r = d.get((case, f'np{np}'))
    if not r or not r['t2']: return '-'
    return f"{r['t2']} ({r['p']}) / {r['t1'] or '-'}"
for case in cases:
    for np in (2, 3, 4):
        a, b = stat(tip, case, np), stat(br, case, np)
        if a == '-' and b == '-': continue
        print(f'| {case} | {np} | {a} | {b} |')
