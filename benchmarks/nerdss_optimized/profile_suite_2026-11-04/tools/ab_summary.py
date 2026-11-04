#!/usr/bin/env python3
"""Median wall time per build and speed-up (A/B) from ab_timing.sh output.
usage: ab_summary.py <ab.tsv>"""
import csv, sys, statistics as st, collections
t = collections.defaultdict(lambda: collections.defaultdict(list))
order = []
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    if r['case'] not in order: order.append(r['case'])
    t[r['case']][r['build']].append(float(r['wall_s']))
print('| case | reps | baseline median (s) | prefilter median (s) | speed-up | per-rep ratios |')
print('| --- | ---: | ---: | ---: | ---: | --- |')
ga = gb = 0
for c in order:
    a, b = t[c]['A'], t[c]['B']; n = min(len(a), len(b))
    if n == 0: continue
    ma, mb = st.median(a[:n]), st.median(b[:n]); ga += ma; gb += mb
    ratios = ', '.join(f'{x / y:.3f}' for x, y in zip(a, b))
    print(f'| {c} | {n} | {ma:.2f} | {mb:.2f} | {ma / mb:.3f}x | {ratios} |')
print(f'| **sum of medians** | | {ga:.2f} | {gb:.2f} | **{ga / gb:.3f}x** | |')
