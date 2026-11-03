#!/usr/bin/env python3
"""Print analyze_audit.py --tsv output as the two Markdown tables of the
audit write-up: failures and reports per case and np, and the copy-number
comparison against serial.

usage: analyze_audit.py --tsv <root>... | audit_tables.py [--label NAME]
"""
import sys
label = 'build'
if '--label' in sys.argv: label = sys.argv[sys.argv.index('--label') + 1]
rows = [l.rstrip('\n').split('\t') for l in sys.stdin if l.strip()]
rows = [r for r in rows if r[1] != 'serial' and not r[1].endswith('-asan')]
print(f'| case | np | runs/failed | own | struct | ASan | conflicts | mass | max\\|t\\| 1st | max\\|t\\| 2nd (p) | column (2nd half) |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    case, cfg, runs, fail, own, st, asan, cf, dr, mass, t1, t2, p, col = r
    t2p = f'{t2} ({p})' if t2 else '-'
    print(f'| {case} | {cfg[2:]} | {runs}/{fail} | {own} | {st} | {asan} | {cf} | {mass} | {t1 or "-"} | {t2p} | {col} |')
