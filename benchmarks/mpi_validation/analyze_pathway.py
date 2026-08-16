#!/usr/bin/env python3
"""Compare assembly pathways (histogram_complexes_time.dat) between serial and MPI.

copy_numbers_time.dat says how much of each species exists; it does not say what
the material is assembled INTO.  histogram_complexes_time.dat does: at each
recorded time it lists, for every distinct complex composition, how many copies
of that composition exist.  Two runs can agree on every species count and still
reach it by different routes -- many small complexes versus few large ones -- so
the pathway is tested separately here.

File format, one block per recorded time:

    Time (s): 1e-05
    565     IL: 1.
    35      IL: 1. A: 1.
    465     A: 1.

The MPI build emits compositions in a different order than serial, so each
composition is canonicalised to a sorted "A:1,IL:1" key before comparison.

Two tests per case:

  * per-composition Welch t-test, Benjamini-Hochberg corrected across the whole
    family, saying WHICH complex types differ.
  * a global label-permutation test on the L1 distance between the serial and
    MPI mean histograms.  This asks whether the pathway as a whole differs, makes
    no distributional assumption, and needs no correction because it is one test
    per case.  A permutation test is the right tool here: the per-composition
    counts are neither independent nor normal (they are constrained by mass
    conservation), which is exactly the situation where a t-test per bin and a
    naive chi-square both mislead.
"""
import math, os, re, sys, json
from collections import defaultdict

# ---------- parsing ----------
def canon(comp):
    """'IL: 1. A: 1.' -> 'A:1,IL:1' so serial and MPI orderings agree."""
    parts = re.findall(r'([A-Za-z_][A-Za-z0-9_~!\-\.]*)\s*:\s*(\d+)', comp)
    if not parts:
        return comp.strip()
    return ','.join(f'{n}:{c}' for n, c in sorted(parts))

def parse_hist(path):
    """Return {time: {canonical_composition: count}} for the whole file."""
    if not os.path.exists(path):
        return {}
    out, t = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            m = re.match(r'\s*Time \(s\):\s*(\S+)', line)
            if m:
                t = float(m.group(1)); out[t] = {}
                continue
            if t is None:
                continue
            m = re.match(r'\s*(\d+)\s+(.*\S)\s*$', line)
            if m:
                out[t][canon(m.group(2))] = out[t].get(canon(m.group(2)), 0) + int(m.group(1))
    return out

def final_hist(path):
    h = parse_hist(path)
    if not h:
        return None
    return h[max(h)]

# ---------- statistics (no scipy available) ----------
def _betacf(a, b, x):
    MAXIT, EPS, FPMIN = 300, 3e-16, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0 / d; h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d; de = d * c; h *= de
        if abs(de - 1.0) < EPS: break
    return h

def betai(a, b, x):
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    lb = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lb + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b

def t_p(t, df):
    if df <= 0 or t != t: return float('nan')
    return betai(df / 2.0, 0.5, df / (df + t * t))

def mean(v): return sum(v) / len(v) if v else float('nan')

def sd(v):
    if len(v) < 2: return float('nan')
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))

def welch(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2: return None
    va, vb = sd(a) ** 2, sd(b) ** 2
    se = math.sqrt(va / na + vb / nb)
    diff = mean(b) - mean(a)
    if se == 0:
        return dict(diff=diff, se=0.0, p=(0.0 if diff != 0 else 1.0), degenerate=True)
    t = diff / se
    num = (va / na + vb / nb) ** 2
    den = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    df = num / den if den > 0 else na + nb - 2
    return dict(diff=diff, se=se, p=t_p(t, df), degenerate=False)

def bh(pairs):
    live = sorted([(k, p) for k, p in pairs if p == p], key=lambda kp: kp[1])
    n, out, prev = len(live), {}, 1.0
    for i in range(n - 1, -1, -1):
        k, p = live[i]
        prev = min(prev, min(1.0, p * n / (i + 1))); out[k] = prev
    for k, p in pairs:
        if p != p: out[k] = float('nan')
    return out

# deterministic PRNG: the harness forbids seeding from the clock, and a fixed
# stream also makes the permutation p-values reproducible run to run.
class LCG:
    def __init__(self, seed): self.s = seed & 0xFFFFFFFF
    def next(self, n):
        self.s = (1103515245 * self.s + 12345) & 0x7FFFFFFF
        return self.s % n

def perm_test(A, B, keys, nperm=20000, seed=12345):
    """Label-permutation test on the L1 distance between group mean histograms."""
    def meanvec(rows):
        return [sum(r[k] for r in rows) / len(rows) for k in keys]
    def l1(u, v):
        return sum(abs(x - y) for x, y in zip(u, v))
    obs = l1(meanvec(A), meanvec(B))
    pool = A + B
    na = len(A)
    rng = LCG(seed)
    ge = 0
    for _ in range(nperm):
        idx = list(range(len(pool)))
        for i in range(len(idx) - 1, 0, -1):        # Fisher-Yates
            j = rng.next(i + 1)
            idx[i], idx[j] = idx[j], idx[i]
        pa = [pool[i] for i in idx[:na]]
        pb = [pool[i] for i in idx[na:]]
        if l1(meanvec(pa), meanvec(pb)) >= obs - 1e-12:
            ge += 1
    return obs, (ge + 1) / (nperm + 1)          # add-one: never reports p = 0

# ---------- load runs ----------
# expects <root>/<case>__<cfg>__<seed>/hist.dat
root = sys.argv[1] if len(sys.argv) > 1 else 'bench2_runs'
runs = defaultdict(list)     # (case,cfg) -> [(seed, {comp:count})]
for d in sorted(os.listdir(root)):
    p = os.path.join(root, d, 'hist.dat')
    if not os.path.exists(p): continue
    try: case, cfg, seed = d.split('__')
    except ValueError: continue
    h = final_hist(p)
    if h: runs[(case, cfg)].append((seed, h))

cases = []
for (c, _cfg) in runs:
    if c not in cases: cases.append(c)
cases.sort()

print('=' * 100)
print('PATHWAY COMPARISON  (histogram_complexes_time.dat, final recorded time)')
print('=' * 100)

all_tests, summary = [], []
for case in cases:
    ser = runs.get((case, 'serial'), [])
    if len(ser) < 2: continue
    for cfg in ('np1', 'np2'):
        mpi = runs.get((case, cfg), [])
        if len(mpi) < 2: continue
        keys = sorted({k for _s, h in ser for k in h} | {k for _s, h in mpi for k in h})
        A = [{k: h.get(k, 0) for k in keys} for _s, h in ser]
        B = [{k: h.get(k, 0) for k in keys} for _s, h in mpi]
        obs, pperm = perm_test(A, B, keys)
        tot = sum(mean([r[k] for r in A]) for k in keys)
        summary.append(dict(case=case, cfg=cfg, nser=len(A), nmpi=len(B),
                            ncomp=len(keys), l1=obs,
                            l1rel=(100.0 * obs / tot if tot else float('nan')),
                            pperm=pperm))
        for k in keys:
            w = welch([r[k] for r in A], [r[k] for r in B])
            if w: all_tests.append(dict(case=case, cfg=cfg, comp=k,
                                        a=[r[k] for r in A], b=[r[k] for r in B], **w))

padj = bh([((t['case'], t['cfg'], t['comp']), t['p']) for t in all_tests])

print(f"\n{'case':<22}{'cfg':<6}{'n_ser':>6}{'n_mpi':>6}{'#comps':>8}"
      f"{'L1(mean hist)':>15}{'L1 as % total':>15}{'p_perm':>10}  verdict")
for s in summary:
    v = 'PATHWAY DIFFERS' if s['pperm'] < 0.05 else 'indistinguishable'
    print(f"{s['case']:<22}{s['cfg']:<6}{s['nser']:>6}{s['nmpi']:>6}{s['ncomp']:>8}"
          f"{s['l1']:>15.2f}{s['l1rel']:>15.2f}{s['pperm']:>10.5f}  {v}")

print()
print('=' * 100)
print('PER-COMPOSITION DETAIL  (only compositions differing at FDR 5%)')
print('=' * 100)
shown = defaultdict(int)
for t in sorted(all_tests, key=lambda x: (x['case'], x['cfg'], x['p'])):
    k = (t['case'], t['cfg'], t['comp'])
    if padj[k] >= 0.05: continue
    if shown[(t['case'], t['cfg'])] == 0:
        print(f"\n--- {t['case']}  {t['cfg']}")
        print(f"  {'composition':<40}{'serial':>16}{'mpi':>16}{'diff':>10}{'p_BH':>10}")
    shown[(t['case'], t['cfg'])] += 1
    if shown[(t['case'], t['cfg'])] > 12:
        continue
    sa, sb = sd(t['a']) / math.sqrt(len(t['a'])), sd(t['b']) / math.sqrt(len(t['b']))
    print(f"  {t['comp']:<40}{mean(t['a']):9.2f}+/-{sa:<5.2f}{mean(t['b']):9.2f}+/-{sb:<5.2f}"
          f"{t['diff']:>10.2f}{padj[k]:>10.4f}")
for key, n in shown.items():
    if n > 12:
        print(f"  ... {n - 12} further differing compositions in {key[0]} {key[1]} not shown")
for s in summary:
    if shown[(s['case'], s['cfg'])] == 0:
        print(f"\n--- {s['case']}  {s['cfg']}: no individual composition differs at FDR 5%")

json.dump({'summary': summary}, open(os.path.join(root, '..', 'pathway_summary.json'), 'w'), indent=1)
print('\nwrote pathway_summary.json')
