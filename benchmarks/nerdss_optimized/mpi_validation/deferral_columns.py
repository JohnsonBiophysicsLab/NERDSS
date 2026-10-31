#!/usr/bin/env python3
"""Bound fraction by x column from the restart snapshots of the deferral runs.

Every run writes RESTARTS/restart<itr>.json ten times over the run (nerdss_mpi
gathers one serial-format file), with every molecule's position and bonds.
The x axis is cut into the same columns nerdss_mpi decomposes it into:
floor(Lx / rMaxLimit) of them, rMaxLimit read from the snapshot.  For each
column this prints the hub sites' bound fraction (mean over the seeds of the
per-seed mean over snapshots after the burn-in, standard error over seeds)
per configuration, a Welch t against serial's same column, and the mean hub
count per column so that an ownership artefact cannot be mistaken for a
spatial one.

At np=2 on 8 columns the ranks own 0-3 | 4-7, so column 3 is the one the
right neighbour binds into; at np=4 (1,3,3,1) those are columns 0, 3 and 6.

Usage: deferral_columns.py [--burn 500] [--bins N] <out-dir>
"""
import argparse
import collections
import glob
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deferral_lifetimes as dl   # noqa: E402
import deferral_table as dt       # noqa: E402


def snapshot_columns(path, hub_type, ncol, lx, bins):
    r = json.load(open(path))
    mols = r['molecules']['list']
    width = lx / bins
    hubs = [0] * bins
    bound = [0] * bins
    free_b = [0] * bins
    for m in mols:
        if m.get('isEmpty'): continue
        x = m['comCoord'][0]
        k = min(bins - 1, max(0, int((x + lx / 2) // width)))
        if m['molTypeIndex'] == hub_type:
            hubs[k] += 1
            bound[k] += sum(1 for f in m['interfaceList'] if f['isBound'])
        elif not any(f['isBound'] for f in m['interfaceList']):
            free_b[k] += 1
    return hubs, bound, free_b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--burn', type=int, default=500)
    ap.add_argument('--bins', type=int, default=None, help='x bins (default: the MPI columns)')
    ap.add_argument('--hub', default='A')
    ap.add_argument('--sites', type=int, default=4)
    ap.add_argument('out')
    args = ap.parse_args()
    for kbdir in sorted(glob.glob(os.path.join(args.out, 'kb*'))):
        kb = os.path.basename(kbdir)[2:]
        done = set()
        for tsv in glob.glob(os.path.join(kbdir, 'timing_*.tsv')):
            for line in open(tsv).read().strip().split('\n')[1:]:
                f = line.split('\t')
                if len(f) == 5 and f[4].strip() == '0': done.add((f[0], f[2]))
        per = collections.defaultdict(list)   # config -> per seed [frac per bin], [hubs per bin], [free B per bin]
        ncol = bins = None
        for d in sorted(glob.glob(os.path.join(kbdir, '*_s*'))):
            m = re.match(r'(.+)_s(\d+)$', os.path.basename(d))
            if not m or (m.group(1), m.group(2)) not in done: continue
            parm = glob.glob(os.path.join(d, '*.inp'))[0]
            text = open(parm).read()
            lx = float(re.search(r'WaterBox\s*=\s*\[\s*([^,\]]+)', text).group(1))
            snaps = []
            for p in glob.glob(os.path.join(d, 'RESTARTS', 'restart*.json')):
                itr = int(re.search(r'restart(\d+)\.json', p).group(1))
                if itr >= args.burn: snaps.append(p)
            if not snaps: continue
            if ncol is None:
                r = json.load(open(snaps[0]))
                ncol = int(lx // r['parameters']['rMaxLimit'])
                bins = args.bins or ncol
                hub_type = next(t['molTypeIndex'] for t in r['molTemplates']['templates'] if t['molName'] == args.hub)
            hubs = [0] * bins; bound = [0] * bins; freeb = [0] * bins
            for p in snaps:
                h, b, fb = snapshot_columns(p, hub_type, ncol, lx, bins)
                for k in range(bins):
                    hubs[k] += h[k]; bound[k] += b[k]; freeb[k] += fb[k]
            n = len(snaps)
            per[m.group(1)].append(([bound[k] / (args.sites * hubs[k]) if hubs[k] else float('nan') for k in range(bins)],
                                    [hubs[k] / n for k in range(bins)], [freeb[k] / n for k in range(bins)]))
        if 'serial' not in per: continue
        print(f"\n### offRatekb = {kb}: bound fraction by x bin ({bins} bins over {lx:g} nm, {ncol} MPI columns; snapshots at itr >= {args.burn})\n")
        configs = ['serial'] + sorted(c for c in per if c != 'serial')
        print('| bin | ' + ' | '.join(f"{c} ({len(per[c])} seeds)" + ('' if c == 'serial' else ' / t vs serial') for c in configs) + ' |')
        print('|---|' + '---|' * len(configs))
        for k in range(bins):
            cells = [str(k)]
            ser = [s[0][k] for s in per['serial']]
            for c in configs:
                xs = [s[0][k] for s in per[c]]
                mm, se, _ = dl.mean_se(xs)
                if c == 'serial':
                    cells.append(dt.fmt(mm, se, 4))
                else:
                    t, df, p = dt.welch(xs, ser)
                    cells.append(f"{dt.fmt(mm, se, 4)} / {t:+.1f}")
            print('| ' + ' | '.join(cells) + ' |')
        print(f"\nmean {args.hub} per bin per snapshot, and free B per bin:\n")
        print('| bin | ' + ' | '.join(configs) + ' |')
        print('|---|' + '---|' * len(configs))
        for k in range(bins):
            cells = [str(k)]
            for c in configs:
                h, _, _ = dl.mean_se([s[1][k] for s in per[c]])
                fb, _, _ = dl.mean_se([s[2][k] for s in per[c]])
                cells.append(f"{h:.1f} / {fb:.1f}")
            print('| ' + ' | '.join(cells) + ' |')


if __name__ == '__main__':
    main()
