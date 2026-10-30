#!/usr/bin/env python3
"""Print bench_mpi.sh's timing.tsv as a Markdown table: one row per case,
one column per config, wall seconds (minimum over reps) and, in brackets,
the ratio to the first config's time.  A config that failed shows its exit
status.

usage: bench_table.py <timing.tsv> [config ...]
"""
import sys
from collections import defaultdict
rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1]) if l.strip()]
hdr, rows = rows[0], rows[1:]
col = {k: i for i, k in enumerate(hdr)}
best = defaultdict(dict)   # case -> config -> (secs, exit)
order = []
cases = []
for r in rows:
    case, cfg, secs, ex = r[col['case']], r[col['config']], float(r[col['seconds']]), r[col['exit']]
    if case not in cases: cases.append(case)
    if cfg not in order: order.append(cfg)
    cur = best[case].get(cfg)
    if ex != '0':
        if cur is None: best[case][cfg] = (secs, ex)
    elif cur is None or cur[1] != '0' or secs < cur[0]:
        best[case][cfg] = (secs, ex)
configs = sys.argv[2:] or order
print('| case | ' + ' | '.join(configs) + ' |')
print('|---|' + '---|' * len(configs))
for case in cases:
    ref = best[case].get(configs[0])
    cells = []
    for cfg in configs:
        v = best[case].get(cfg)
        if v is None: cells.append('-')
        elif v[1] != '0': cells.append(f'exit {v[1]}')
        elif ref and ref[1] == '0' and cfg != configs[0]: cells.append(f'{v[0]:.1f} ({v[0]/ref[0]:.2f}x)')
        else: cells.append(f'{v[0]:.1f}')
    print(f'| {case} | ' + ' | '.join(cells) + ' |')
