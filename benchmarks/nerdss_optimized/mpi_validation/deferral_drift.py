#!/usr/bin/env python3
"""Left-right imbalance of hubs and free partners in the deferral runs, from
the RESTARTS/restart<itr>.json snapshots (see DEFERRAL_2026-10-30.md).

    python3 deferral_drift.py <kb dir> <config> [<config> ...] [--from 1200]

For each configuration (serial, np1, np2, ...): per seed, the mean over the
snapshots at itr >= --from of (hubs in the left half of the box) - (hubs in
the right half), the same for unbound partner molecules, then the mean and
SE over seeds with t = mean/SE; and the mean count of hubs and free partners
per MPI column (8 columns over the box) at each snapshot time, to show
whether an imbalance is there from the start or builds up.
"""
import argparse
import collections
import glob
import json
import math
import os


def snapshot_counts(path, hub, lx, ncol):
    r = json.load(open(path))
    hub_type = next(t['molTypeIndex'] for t in r['molTemplates']['templates']
                    if t['molName'] == hub)
    hubs = [0] * ncol
    free = [0] * ncol
    for m in r['molecules']['list']:
        if m.get('isEmpty'):
            continue
        k = min(ncol - 1, max(0, int((m['comCoord'][0] + lx / 2) / lx * ncol)))
        if m['molTypeIndex'] == hub_type:
            hubs[k] += 1
        elif not any(f['isBound'] for f in m['interfaceList']):
            free[k] += 1
    return hubs, free


def mean_se(xs):
    n = len(xs)
    m = sum(xs) / n
    se = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1) / n) if n > 1 else float('nan')
    return m, se


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('kbdir')
    ap.add_argument('configs', nargs='+')
    ap.add_argument('--from', dest='start', type=int, default=1200)
    ap.add_argument('--hub', default='A')
    ap.add_argument('--lx', type=float, default=270.0)
    ap.add_argument('--columns', type=int, default=8)
    args = ap.parse_args()
    half = args.columns // 2
    for cfg in args.configs:
        rows = []
        prof_h = collections.defaultdict(lambda: [0.0] * args.columns)
        prof_f = collections.defaultdict(lambda: [0.0] * args.columns)
        nprof = collections.Counter()
        for run in sorted(glob.glob(os.path.join(args.kbdir, f'{cfg}_s*'))):
            dh = db = 0.0
            snaps = []
            for p in glob.glob(os.path.join(run, 'RESTARTS', 'restart*.json')):
                itr = int(os.path.basename(p)[7:-5])
                hubs, free = snapshot_counts(p, args.hub, args.lx, args.columns)
                for k in range(args.columns):
                    prof_h[itr][k] += hubs[k]
                    prof_f[itr][k] += free[k]
                nprof[itr] += 1
                if itr >= args.start:
                    snaps.append((hubs, free))
            if not snaps:
                continue
            for hubs, free in snaps:
                dh += (sum(hubs[:half]) - sum(hubs[half:])) / len(snaps)
                db += (sum(free[:half]) - sum(free[half:])) / len(snaps)
            rows.append((os.path.basename(run), dh, db))
        if not rows:
            print(f'{cfg}: no finished runs')
            continue
        mh, seh = mean_se([r[1] for r in rows])
        mb, seb = mean_se([r[2] for r in rows])
        print(f'{cfg}: {len(rows)} seeds, snapshots at itr >= {args.start}')
        print(f'  hubs left - right   = {mh:+.1f} +- {seh:.1f} (t = {mh / seh:+.1f})'
              '  per seed: ' + ' '.join(f'{r[1]:+.0f}' for r in rows))
        print(f'  free B left - right = {mb:+.1f} +- {seb:.1f} (t = {mb / seb:+.1f})'
              '  per seed: ' + ' '.join(f'{r[2]:+.0f}' for r in rows))
        print(f'  per column (left..right), mean over seeds: hubs | free B')
        for itr in sorted(prof_h):
            n = nprof[itr]
            print(f'    itr {itr:5d}: ' + ' '.join(f'{v / n:5.1f}' for v in prof_h[itr]) +
                  ' | ' + ' '.join(f'{v / n:5.1f}' for v in prof_f[itr]))


if __name__ == '__main__':
    main()
