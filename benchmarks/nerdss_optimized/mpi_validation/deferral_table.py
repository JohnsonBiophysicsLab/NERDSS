#!/usr/bin/env python3
"""Tables for DEFERRAL_2026-10-30.md from deferral_run.sh output.

For each off rate under <out-dir> (kb*/), reads every <config>_s<seed> run
with deferral_lifetimes.py's analysis and prints, per configuration, the
mean over seeds of the per-run mean bond lifetime (all bonds, gainers,
others), the bound fraction, with standard errors over seeds, and a Welch t
of each MPI configuration against serial on the per-seed means.  The pooled
column treats every counted bond as one sample (the per-run standard errors
combined), which is what the per-seed t should agree with when the seeds are
independent.

Usage: deferral_table.py [--burn 500] [--tail 400] <out-dir>
"""
import argparse
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deferral_lifetimes as dl   # noqa: E402


def betacf(a, b, x):
    # continued fraction for the incomplete beta function (Numerical Recipes)
    MAXIT, EPS, FPMIN = 300, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < FPMIN: d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN: d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS: break
    return h


def betai(a, b, x):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * betacf(a, b, x) / a
    return 1.0 - bt * betacf(b, a, 1 - x) / b


def t_two_sided_p(t, df):
    if math.isnan(t) or df <= 0: return float('nan')
    return betai(df / 2, 0.5, df / (df + t * t))


def welch(xs, ys):
    """t, df, p of mean(xs) - mean(ys)."""
    mx, sex, nx = dl.mean_se(xs)
    my, sey, ny = dl.mean_se(ys)
    if nx < 2 or ny < 2: return float('nan'), float('nan'), float('nan')
    vx, vy = sex ** 2, sey ** 2
    t = (mx - my) / math.sqrt(vx + vy)
    df = (vx + vy) ** 2 / (vx ** 2 / (nx - 1) + vy ** 2 / (ny - 1))
    return t, df, t_two_sided_p(t, df)


def pooled_z(runs_a, runs_b, key_m, key_se, key_n):
    """All counted bonds as one sample per arm: weighted mean, combined SE."""
    def pool(runs):
        n = sum(r[key_n] for r in runs)
        m = sum(r[key_m] * r[key_n] for r in runs) / n
        var = sum((r[key_se] ** 2) * r[key_n] ** 2 for r in runs) / n ** 2
        return m, math.sqrt(var), n
    ma, sa, na = pool(runs_a)
    mb, sb, nb = pool(runs_b)
    return (ma - mb) / math.sqrt(sa ** 2 + sb ** 2), ma, sa, na


def fmt(m, se, digits=3):
    return f"{m:.{digits}f} ± {se:.{digits}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--burn', type=int, default=500)
    ap.add_argument('--tail', type=int, default=None,
                    help='default: 20 mean bond lifetimes (20 / (kb dt)) of the off rate')
    ap.add_argument('--hub', default='A')
    ap.add_argument('--sites', type=int, default=4)
    ap.add_argument('out')
    args = ap.parse_args()
    for kbdir in sorted(glob.glob(os.path.join(args.out, 'kb*'))):
        kb = os.path.basename(kbdir)[2:]
        runs = {}
        if args.tail is None:
            parm = glob.glob(os.path.join(kbdir, '*_s*', '*.inp'))
            if not parm: continue
            dt = float(re.search(r'timestep\s*=\s*([^\s#]+)', open(parm[0]).read(), re.I).group(1))
            tail = int(round(20.0 / (float(kb) * dt * 1e-6)))
        else:
            tail = args.tail
        args_kb = argparse.Namespace(**vars(args)); args_kb.tail = tail
        # only runs deferral_run.sh recorded as finished with exit 0
        done = set()
        for tsv in glob.glob(os.path.join(kbdir, 'timing_*.tsv')):
            for line in open(tsv).read().strip().split('\n')[1:]:
                f = line.split('\t')
                if len(f) == 5 and f[4].strip() == '0': done.add((f[0], f[2]))
        for d in sorted(glob.glob(os.path.join(kbdir, '*_s*'))):
            m = re.match(r'(.+)_s(\d+)$', os.path.basename(d))
            if not m or not os.path.isdir(d) or (m.group(1), m.group(2)) not in done: continue
            try:
                r = dl.analyse(d, args_kb)
            except Exception as e:   # a run that did not finish
                print(f"skipping {d}: {e}", file=sys.stderr)
                continue
            if r['made_twice'] or r['broken_early'] or \
               r['alive_at_end'] != int(r['last_copy_number_bound']):
                print(f"replay anomaly in {d}: {r}", file=sys.stderr)
            runs.setdefault(m.group(1), []).append(r)
        if 'serial' not in runs: continue
        configs = ['serial'] + sorted(c for c in runs if c != 'serial')
        print(f"\n### offRatekb = {kb} (bonds formed in steps [{args.burn}, nItr - {tail}])\n")
        print("| config | seeds | bonds/run | lifetime (steps) | Welch t (p) | pooled z | gainers' lifetime | t (p) | others' lifetime | t (p) | gainer share | gains per bond-step | bound fraction | t (p) |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        ser = runs['serial']
        for c in configs:
            rs = runs[c]
            cells = [c, str(len(rs)), f"{sum(r['n'] for r in rs) / len(rs):.0f}"]
            for key, nk in (('mean', 'n'), ('mean_gain', 'n_gain'), ('mean_nogain', 'n_nogain')):
                m, se, _ = dl.mean_se([r[key] for r in rs])
                cells.append(fmt(m, se))
                if c == 'serial':
                    cells.append('—')
                    if key == 'mean': cells.append('—')
                else:
                    t, df, p = welch([r[key] for r in rs], [r[key] for r in ser])
                    cells.append(f"{t:+.2f} ({p:.3g})")
                    if key == 'mean':
                        z, _, _, _ = pooled_z(rs, ser, key, 'se' if key == 'mean' else key.replace('mean', 'se'), nk)
                        cells.append(f"{z:+.2f}")
            m, se, _ = dl.mean_se([r['gain_fraction'] for r in rs])
            cells.append(fmt(m, se, 4))
            m, se, _ = dl.mean_se([r['gains_per_bond_step'] for r in rs])
            cells.append(fmt(m, se, 4))
            m, se, _ = dl.mean_se([r['bound_fraction'] for r in rs])
            cells.append(fmt(m, se, 4))
            if c == 'serial':
                cells.append('—')
            else:
                t, df, p = welch([r['bound_fraction'] for r in rs], [r['bound_fraction'] for r in ser])
                cells.append(f"{t:+.2f} ({p:.3g})")
            print('| ' + ' | '.join(cells) + ' |')
        # relative shifts against serial
        ms, _, _ = dl.mean_se([r['mean'] for r in ser])
        mgs, _, _ = dl.mean_se([r['mean_gain'] for r in ser])
        mns, _, _ = dl.mean_se([r['mean_nogain'] for r in ser])
        for c in configs[1:]:
            rs = runs[c]
            m, se, _ = dl.mean_se([r['mean'] for r in rs])
            mg, seg, _ = dl.mean_se([r['mean_gain'] for r in rs])
            mn, sen, _ = dl.mean_se([r['mean_nogain'] for r in rs])
            print(f"\n{c} vs serial: lifetime {100 * (m / ms - 1):+.2f}% ({m - ms:+.3f} steps), "
                  f"gainers {100 * (mg / mgs - 1):+.2f}% ({mg - mgs:+.3f} steps), "
                  f"others {100 * (mn / mns - 1):+.2f}% ({mn - mns:+.3f} steps)")
        tsvs = glob.glob(os.path.join(kbdir, 'timing_*.tsv'))
        if tsvs:
            print()
            for tsv in sorted(tsvs):
                rows = [l.split('\t') for l in open(tsv).read().strip().split('\n')[1:]]
                if not rows: continue
                secs = [int(r[3]) for r in rows]
                bad = sum(1 for r in rows if r[4].strip() != '0')
                print(f"{rows[0][0]}: {len(rows)} runs, {sum(secs) / len(secs):.0f} s each (wall, loaded machine), {bad} non-zero exits")


if __name__ == '__main__':
    main()
