#!/usr/bin/env python3
"""Summarise run_mpi_suite.sh result directories: failures, checker reports,
mass conservation, and copy numbers against a serial seed ensemble.

usage: analyze_audit.py [--tsv] <root> [<root> ...]

Every <root>/<label>/results.tsv is read.  The label must end in
_np<N>_s<seed> (MPI) or be serial_..._s<seed>; the configuration is the np
(or 'serial').  Per case and configuration:

* runs, nonzero exits, and the summed ownership / structure / ASan /
  exchange-conflict counts, plus 'dropping the member' lines in stderr;
* mass conservation: in every block of the merged complex histogram, the
  number of each molecule type summed over compositions must equal the
  serial runs' constant total (types whose serial total changes over time,
  i.e. created or destroyed ones, are skipped);
* copy numbers: each column of copy_numbers_time.dat averaged over the
  second half of the rows (equilibrium) and over the first half after
  time 0 (relaxation, where a rate error shows), per run; Welch t of the
  MPI ensemble against the serial one, reported as the largest |t| over
  columns (uncorrected for their number) and, for the second half, the
  column and its p.
"""
import math, os, re, sys
from collections import defaultdict

def tcdf_two_sided(t, df):
    # p = I_{df/(df+t^2)}(df/2, 1/2), by continued fraction (Numerical Recipes)
    x = df / (df + t * t)
    a, b = df / 2.0, 0.5
    def betacf(a, b, x):
        qab, qap, qam = a + b, a + 1, a - 1
        c, d = 1.0, 1 - qab * x / qap
        d = 1 / (d if abs(d) > 1e-30 else 1e-30); h = d
        for m in range(1, 200):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-30 else 1e-30)
            c = 1 + aa / c; c = c if abs(c) > 1e-30 else 1e-30
            h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-30 else 1e-30)
            c = 1 + aa / c; c = c if abs(c) > 1e-30 else 1e-30
            de = d * c; h *= de
            if abs(de - 1) < 3e-12: break
        return h
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return front * betacf(a, b, x) / a
    return 1 - front * betacf(b, a, 1 - x) / b

def welch(x, y):
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2: return None
    mx, my = sum(x) / nx, sum(y) / ny
    vx = sum((v - mx) ** 2 for v in x) / (nx - 1)
    vy = sum((v - my) ** 2 for v in y) / (ny - 1)
    se2 = vx / nx + vy / ny
    if se2 == 0:
        return (0.0, 1.0) if mx == my else (float('inf'), 0.0)
    t = (mx - my) / math.sqrt(se2)
    df = se2 ** 2 / ((vx / nx) ** 2 / (nx - 1 if nx > 1 else 1) + (vy / ny) ** 2 / (ny - 1 if ny > 1 else 1) + 1e-300)
    return t, tcdf_two_sided(t, df)

def cn_means(path, half=2):
    try:
        lines = [l.strip() for l in open(path) if l.strip()]
    except OSError:
        return None, None
    if len(lines) < 3: return None, None
    head = lines[0].split(',')[1:]
    timed = []
    for l in lines[1:]:
        try: timed.append([float(v) for v in l.split(',')])
        except ValueError: pass
    # Keep the timeWrite rows only, so that serial and MPI average the same
    # times.  The serial build runs steps 1 to nItr - 1 and writes a row at its
    # last step whatever timeWrite is; nerdss_mpi runs to nItr and writes the
    # timeWrite rows only.  When nItr is not a multiple of timeWrite that row
    # is extra: hexamer's "second half" was {10, 12.5 ms} for serial and
    # {10 ms} for MPI, 11% more bonds for serial.  When it is, serial's last
    # row stands one step short of MPI's (rev_2D: step 199, 1.99e-05, against
    # 200), which is the same row for this purpose, so a row is regular when
    # it lands within 5% of a step of where one is due.  The duplicate filter
    # is for nItr - 1 itself a multiple of timeWrite, where serial writes the
    # row twice.
    timed = [r for k, r in enumerate(timed) if k == 0 or r[0] != timed[k - 1][0]]
    if len(timed) >= 3:
        step = timed[1][0] - timed[0][0]
        if abs(timed[-1][0] - timed[-2][0] - step) > 0.05 * step:
            timed.pop()
    rows = [r[1:] for r in timed]
    # half 2: the second half of the rows (equilibrium); half 1: the first
    # half after the time-0 row (relaxation, where a rate error shows)
    rows = rows[len(rows) // 2:] if half == 2 else rows[1:len(rows) // 2]
    if not rows: return head, None
    return head, [sum(r[i] for r in rows) / len(rows) for i in range(len(head))]

SPEC = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\s*:\s*(\d+)\.')
def hist_totals(path):
    """list of {type: molecules} per block"""
    out = []
    try:
        f = open(path)
    except OSError:
        return None
    cur = None
    for line in f:
        if line.startswith('Time'):
            cur = defaultdict(int); out.append(cur); continue
        if cur is None or not line.strip(): continue
        parts = line.split('\t', 1)
        if len(parts) != 2: continue
        try: n = int(parts[0])
        except ValueError: continue
        for sp, k in SPEC.findall(parts[1]):
            cur[sp] += n * int(k)
    return out

tsv = '--tsv' in sys.argv[1:]
runs = []   # dicts
for root in [a for a in sys.argv[1:] if a != '--tsv']:
    for label in sorted(os.listdir(root)):
        res = os.path.join(root, label, 'results.tsv')
        if not os.path.exists(res): continue
        m = re.search(r'_np(\d+)_s(\d+)$', label)
        if label.startswith('serial'):
            cfg = 'serial'; seed = re.search(r'_s(\d+)$', label).group(1)
        elif m:
            cfg = 'np' + m.group(1); seed = m.group(2)
        else:
            continue
        build = label.split('_')[0]
        # runs of a different length are not comparable with the serial ones
        if build == 'asan': cfg += '-asan'
        lines = open(res).read().splitlines()
        hdr = lines[0].split('\t')
        for l in lines[1:]:
            f = dict(zip(hdr, l.split('\t')))
            rd = os.path.join(root, label, 'runs', f['case'])
            drops = 0
            try:
                drops = sum(1 for x in open(os.path.join(rd, 'stderr.log'), errors='replace') if 'dropping the member' in x)
            except OSError: pass
            runs.append(dict(case=f['case'], cfg=cfg, seed=seed, build=build, label=label,
                             exit=f['exit_status'], own=int(f.get('ownership_violations') or 0),
                             struct=int(f.get('structure_violations') or 0),
                             asan=int(f.get('asan_reports') or 0),
                             confl=int(f.get('exchange_conflicts') or 0), drops=drops,
                             secs=float(f['seconds']), dir=rd))

by = defaultdict(list)
for r in runs: by[(r['case'], r['cfg'])].append(r)
cases = sorted({r['case'] for r in runs})
cfgs = sorted({r['cfg'] for r in runs}, key=lambda c: (c != 'serial', c))

# serial references
ser_cn, ser_mass = {}, {}
for case in cases:
    vals = []
    totals = defaultdict(set)
    for r in by.get((case, 'serial'), []):
        h, mns = cn_means(os.path.join(r['dir'], 'mergeOUT', 'copy_numbers_time.dat'))
        if mns: vals.append((h, mns))
        ht = hist_totals(os.path.join(r['dir'], 'mergeOUT', 'histogram_complexes_time.dat')) or []
        for blk in ht:
            for sp, v in blk.items(): totals[sp].add(v)
    ser_cn[case] = vals
    ser_mass[case] = {sp: next(iter(v)) for sp, v in totals.items() if len(v) == 1}

if not tsv: print(f"{'case':<20}{'cfg':<7}{'runs':>5}{'fail':>5}{'own':>6}{'struct':>7}{'asan':>6}{'confl':>7}{'drops':>6}{'mass':>7}{'|t|1st':>8}{'|t|2nd':>8}{'p 2nd':>9}  column (2nd half)")
for case in cases:
    for cfg in cfgs:
        rs = by.get((case, cfg))
        if not rs: continue
        fail = sum(1 for r in rs if r['exit'] != '0')
        own = sum(r['own'] for r in rs); st = sum(r['struct'] for r in rs)
        asan = sum(r['asan'] for r in rs); cf = sum(r['confl'] for r in rs)
        dr = sum(r['drops'] for r in rs)
        # mass
        massbad = 0
        if cfg != 'serial' and ser_mass.get(case):
            for r in rs:
                ht = hist_totals(os.path.join(r['dir'], 'mergeOUT', 'histogram_complexes_time.dat')) or []
                for blk in ht:
                    if any(blk.get(sp, 0) != v for sp, v in ser_mass[case].items()):
                        massbad += 1; break
        mass = '-' if cfg == 'serial' or not ser_mass.get(case) else (f'{massbad}bad' if massbad else 'ok')
        # stats
        tmax, pmin, col = '', '', ''
        if cfg != 'serial' and len(ser_cn.get(case, [])) >= 2:
            mv = []
            for r in rs:
                if r['exit'] != '0': continue
                h, mns = cn_means(os.path.join(r['dir'], 'mergeOUT', 'copy_numbers_time.dat'))
                if mns: mv.append((h, mns))
            if len(mv) >= 2:
                h0 = ser_cn[case][0][0]
                best = None
                for i, name in enumerate(h0):
                    x = [m[i] for hh, m in mv if hh == h0]
                    y = [m[i] for hh, m in ser_cn[case] if hh == h0]
                    w = welch(x, y)
                    if w and (best is None or abs(w[0]) > abs(best[0])):
                        best = (w[0], w[1], name, sum(x)/len(x), sum(y)/len(y))
                if best:
                    tmax = f'{abs(best[0]):.2f}'; pmin = f'{best[1]:.2g}'
                    col = f'{best[2]} mpi={best[3]:.1f} ser={best[4]:.1f}'
        # the same for the relaxation half
        t1 = ''
        if cfg != 'serial' and len(ser_cn.get(case, [])) >= 2:
            mv1 = [cn_means(os.path.join(r['dir'], 'mergeOUT', 'copy_numbers_time.dat'), 1) for r in rs if r['exit'] == '0']
            sv1 = [cn_means(os.path.join(r['dir'], 'mergeOUT', 'copy_numbers_time.dat'), 1) for r in by.get((case, 'serial'), [])]
            mv1 = [m for m in mv1 if m[1]]; sv1 = [m for m in sv1 if m[1]]
            if len(mv1) >= 2 and len(sv1) >= 2:
                h0 = sv1[0][0]; best = 0.0
                for i in range(len(h0)):
                    w = welch([m[i] for hh, m in mv1 if hh == h0], [m[i] for hh, m in sv1 if hh == h0])
                    if w and abs(w[0]) > best: best = abs(w[0])
                t1 = f'{best:.2f}'
        if tsv:
            print('\t'.join(map(str, [case, cfg, len(rs), fail, own, st, asan, cf, dr, mass, t1, tmax, pmin, col])))
            continue
        print(f"{case:<20}{cfg:<7}{len(rs):>5}{fail:>5}{own:>6}{st:>7}{asan:>6}{cf:>7}{dr:>6}{mass:>7}{t1:>8}{tmax:>8}{pmin:>9}  {col}")
