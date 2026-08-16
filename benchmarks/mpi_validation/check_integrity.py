#!/usr/bin/env python3
"""Two integrity checks on the pathway data that need no distributional assumption
and work even when only a few MPI runs survive.

1. Mass conservation: summing every complex composition must return the molecule
   count the input declares, in every run, exactly.  Serial defines the truth.
2. Novel compositions: any complex composition that MPI produces and serial never
   produces in any of its runs.  This needs no knowledge of the reaction network --
   serial's own output defines what the network can make.
"""
import os, re, sys
from collections import defaultdict

root = sys.argv[1]

def final_block(p):
    blocks, cur = [], None
    for line in open(p):
        if 'Time (s)' in line:
            cur = []; blocks.append(cur); continue
        if cur is not None and line.strip():
            m = re.match(r'\s*(\d+)\s+(.*\S)\s*$', line)
            if m: cur.append((int(m.group(1)), m.group(2)))
    return blocks[-1] if blocks else []

def canon(c):
    p = re.findall(r'([A-Za-z_][A-Za-z0-9_~!\-\.]*)\s*:\s*(\d+)', c)
    return ','.join(f'{n}:{k}' for n, k in sorted(p)) if p else c.strip()

totals = defaultdict(lambda: defaultdict(list))   # (case,cfg) -> species -> [totals]
comps  = defaultdict(set)                          # (case,cfg) -> {composition}
nruns  = defaultdict(int)
for d in sorted(os.listdir(root)):
    p = os.path.join(root, d, 'hist.dat')
    if not os.path.exists(p): continue
    case, cfg, seed = d.split('__')
    nruns[(case, cfg)] += 1
    cnt = defaultdict(int)
    for n, comp in final_block(p):
        comps[(case, cfg)].add(canon(comp))
        for sp, k in re.findall(r'([A-Za-z_][A-Za-z0-9_~!\-\.]*)\s*:\s*(\d+)', comp):
            cnt[sp] += n * int(k)
    for sp, v in cnt.items():
        totals[(case, cfg)][sp].append(v)

cases = sorted({c for c, _ in totals})
print('=' * 92)
print('1. MASS CONSERVATION  (serial defines the expected count; exact match required)')
print('=' * 92)
print(f"{'case':<17}{'cfg':<8}{'species':<12}{'expected':>9}{'mean':>10}{'min':>7}{'max':>7}{'runs':>6}  status")
for case in cases:
    exp = {}
    for sp, v in totals[(case, 'serial')].items():
        exp[sp] = v[0] if len(set(v)) == 1 else None
    for cfg in ('serial', 'np1', 'np2'):
        if (case, cfg) not in totals: continue
        for sp in sorted(totals[(case, cfg)]):
            v = totals[(case, cfg)][sp]
            e = exp.get(sp)
            ok = (e is not None and min(v) == max(v) == e)
            if e is None: st = 'serial itself varies'
            elif ok: st = 'conserved'
            else: st = '*** VIOLATED ***'
            print(f"{case:<17}{cfg:<8}{sp:<12}{str(e):>9}{sum(v)/len(v):>10.2f}{min(v):>7}{max(v):>7}"
                  f"{len(v):>6}  {st}")

print()
print('=' * 92)
print('2. COMPOSITIONS MPI PRODUCES THAT SERIAL NEVER DOES')
print('=' * 92)
any_novel = False
for case in cases:
    base = comps.get((case, 'serial'), set())
    for cfg in ('np1', 'np2'):
        if (case, cfg) not in comps: continue
        novel = comps[(case, cfg)] - base
        if novel:
            any_novel = True
            print(f"  {case:<17}{cfg:<6}{len(novel)} novel over {nruns[(case,cfg)]} runs: "
                  + ', '.join(sorted(novel)[:8]))
        else:
            print(f"  {case:<17}{cfg:<6}none")
if not any_novel:
    print("  (none anywhere)")
