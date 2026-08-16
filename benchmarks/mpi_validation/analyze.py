#!/usr/bin/env python3
"""Serial vs MPI mean-shift analysis over every species in every sample input.

MPI draws random numbers in a different order than serial, so no seed reproduces
a serial trajectory step for step and bitwise comparison is meaningless.  What is
testable is whether the DISTRIBUTION of the end state is shifted.  So: many seeds
per configuration, then Welch's t-test on each species count at the final time.

Welch rather than Student because the MPI groups are visibly noisier than serial
and there is no reason to assume equal variances.  Every species of every case is
tested, so p-values are corrected two ways: Holm (controls the chance of ANY
false positive across the whole family -- strict) and Benjamini-Hochberg (controls
the expected FRACTION of false positives among the flagged -- the usual choice
when screening many endpoints).
"""
import csv, math, sys, os
from collections import defaultdict

SP = os.path.dirname(os.path.abspath(__file__))

# ---------- Student t distribution via the regularized incomplete beta ----------
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

def t_two_tailed_p(t, df):
    if df <= 0 or t != t: return float('nan')
    return betai(df / 2.0, 0.5, df / (df + t * t))

def t_crit_975(df):
    lo, hi = 0.0, 200.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if t_two_tailed_p(mid, df) > 0.05: lo = mid
        else: hi = mid
    return (lo + hi) / 2.0

def mean(v): return sum(v) / len(v)

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
        return dict(diff=diff, se=0.0, t=float('nan'), df=float('nan'),
                    p=(0.0 if diff != 0 else 1.0), degenerate=True)
    t = diff / se
    num = (va / na + vb / nb) ** 2
    den = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    df = num / den if den > 0 else na + nb - 2
    return dict(diff=diff, se=se, t=t, df=df, p=t_two_tailed_p(t, df), degenerate=False)

def holm(pairs):
    live = sorted([(k, p) for k, p in pairs if p == p], key=lambda kp: kp[1])
    n, out, run = len(live), {}, 0.0
    for i, (k, p) in enumerate(live):
        run = max(run, min(1.0, (n - i) * p)); out[k] = run
    for k, p in pairs:
        if p != p: out[k] = float('nan')
    return out

def bh(pairs):
    live = sorted([(k, p) for k, p in pairs if p == p], key=lambda kp: kp[1])
    n, out, prev = len(live), {}, 1.0
    for i in range(n - 1, -1, -1):
        k, p = live[i]
        prev = min(prev, min(1.0, p * n / (i + 1))); out[k] = prev
    for k, p in pairs:
        if p != p: out[k] = float('nan')
    return out

# ---------- headers ----------
def load_headers(case):
    for f in (os.path.join(SP, 'bench2_headers', case + '.txt'),
              os.path.join(SP, 'calib', case, 'DATA', 'copy_numbers_time.dat')):
        if os.path.exists(f):
            with open(f) as fh:
                return [c.strip() for c in fh.readline().strip().split(',')]
    return None

# ---------- load results ----------
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SP, 'bench_results.tsv')
status = defaultdict(lambda: defaultdict(int))
rows = defaultdict(list)          # (case,cfg) -> [ [floats...] ]
case_order = []
with open(path) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        case, cfg = r['case'], r['config']
        if case not in case_order: case_order.append(case)
        status[(case, cfg)][r['status']] += 1
        if r['status'] != 'ok' or not r['final_row']: continue
        try: rows[(case, cfg)].append([float(x) for x in r['final_row'].split(',')])
        except ValueError: status[(case, cfg)]['BADROW'] += 1

# ---------- completion ----------
print('=' * 104)
print('RUN COMPLETION')
print('=' * 104)
print(f"{'case':<16}{'serial':<26}{'mpi np=1':<26}{'mpi np=2':<26}")
for case in case_order:
    cells = []
    for cfg in ('serial', 'np1', 'np2'):
        st = status[(case, cfg)]; tot = sum(st.values())
        bad = ' '.join(f'{k}:{v}' for k, v in sorted(st.items()) if k != 'ok')
        cells.append(f"ok {st.get('ok',0)}/{tot}" + (f'  [{bad}]' if bad else ''))
    print(f'{case:<16}' + ''.join(f'{c:<26}' for c in cells))

# ---------- build the test family ----------
tests = []
for case in case_order:
    hdr = load_headers(case)
    ser = rows[(case, 'serial')]
    if not ser: continue
    ncol = min(len(r) for r in ser)
    for col in range(1, ncol):                       # column 0 is time
        name = hdr[col] if hdr and col < len(hdr) else f'col{col}'
        a = [r[col] for r in ser if len(r) > col]
        for cfg in ('np1', 'np2'):
            b = [r[col] for r in rows[(case, cfg)] if len(r) > col]
            if len(a) < 2 or len(b) < 2: continue
            if sd(a) == 0 and sd(b) == 0 and mean(a) == mean(b):
                continue                              # constant and identical: nothing to test
            w = welch(a, b)
            if w is None: continue
            tests.append(dict(case=case, col=col, name=name, cfg=cfg, a=a, b=b, **w))

p_holm = holm([((t['case'], t['col'], t['cfg']), t['p']) for t in tests])
p_bh   = bh([((t['case'], t['col'], t['cfg']), t['p']) for t in tests])

print()
print('=' * 104)
print(f'PER-SPECIES TEST  ({len(tests)} tests across {len(case_order)} cases; diff = MPI mean - serial mean)')
print('=' * 104)
for case in case_order:
    ts = [t for t in tests if t['case'] == case]
    if not ts: continue
    print(f'\n--- {case}')
    print(f"  {'species':<26}{'cfg':<6}{'serial mean+/-SEM':<22}{'mpi mean+/-SEM':<22}"
          f"{'diff':>9}{'rel%':>8}{'p':>9}{'p_BH':>9}{'p_Holm':>9}  verdict")
    for t in sorted(ts, key=lambda x: (x['col'], x['cfg'])):
        a, b = t['a'], t['b']
        sa, sb = sd(a) / math.sqrt(len(a)), sd(b) / math.sqrt(len(b))
        base = mean(a)
        rel = 100.0 * t['diff'] / base if base else float('nan')
        k = (t['case'], t['col'], t['cfg'])
        ph, pb = p_holm[k], p_bh[k]
        if t['degenerate'] and t['diff'] != 0:
            verdict = 'SHIFTED (deterministic)'
        elif pb < 0.05:
            verdict = 'SHIFTED' + ('' if ph < 0.05 else ' (BH only)')
        else:
            verdict = 'no shift'
        print(f"  {t['name']:<26}{t['cfg']:<6}"
              f"{mean(a):9.2f} +/-{sa:<8.2f}{mean(b):9.2f} +/-{sb:<8.2f}"
              f"{t['diff']:>9.2f}{rel:>8.2f}{t['p']:>9.4f}{pb:>9.4f}{ph:>9.4f}  {verdict}")

print()
print('=' * 104)
print('SENSITIVITY: smallest shift each comparison could have resolved')
print('(half-width of the 95% CI on the difference, as % of the serial mean)')
print('=' * 104)
for case in case_order:
    ts = [t for t in tests if t['case'] == case and not t['degenerate']]
    if not ts: continue
    for cfg in ('np1', 'np2'):
        sub = [t for t in ts if t['cfg'] == cfg]
        if not sub: continue
        vals = []
        for t in sub:
            base = mean(t['a'])
            if base: vals.append(100.0 * t_crit_975(t['df']) * t['se'] / base)
        if vals:
            print(f'  {case:<16}{cfg:<6}median +/-{sorted(vals)[len(vals)//2]:5.2f}%   '
                  f'best +/-{min(vals):5.2f}%   worst +/-{max(vals):5.2f}%')

print()
print('=' * 104)
print('SUMMARY')
print('=' * 104)
for case in case_order:
    for cfg in ('np1', 'np2'):
        ts = [t for t in tests if t['case'] == case and t['cfg'] == cfg]
        if not ts: continue
        nsh = sum(1 for t in ts if p_bh[(t['case'], t['col'], t['cfg'])] < 0.05
                  or (t['degenerate'] and t['diff'] != 0))
        print(f'  {case:<16}mpi {cfg:<5}{nsh}/{len(ts)} species shifted at FDR 5%')
