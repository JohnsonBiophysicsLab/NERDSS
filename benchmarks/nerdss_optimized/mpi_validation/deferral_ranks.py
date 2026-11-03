#!/usr/bin/env python3
"""Per-rank view of the deferral runs: where in the box a difference sits.

Each rank's copy_numbers_time_<rank>.dat counts the molecules it owns, so at
np=4 on the 8-column deferral box (columns split 1,3,3,1) rank 0's file is
exactly the edge column that rank 1 binds into during its first half, and
rank 3's file is a column no neighbour binds into.  Each rank's
assoc_dissoc_time_<rank>.dat holds the reactions it performed, so a bond's
BREAK rank says roughly where it lived (a complex moves about a column width
in a 20-step lifetime).

Prints, for every rank count found under <out-dir>/kb*/, per rank over the
seeds: the bound fraction of the rank's hub sites after the burn-in, the mean
lifetime of the bonds whose BREAK the rank logged, and the share of gains it
logged; plus the serial values for reference.

Usage: deferral_ranks.py [--burn 500] [--tail 400] <out-dir>
"""
import argparse
import collections
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deferral_lifetimes as dl   # noqa: E402
import deferral_table as dt       # noqa: E402


def rank_bound_fraction(run, rank, hub, sites_per_hub, burn_time):
    path = os.path.join(run, 'DATA', f'copy_numbers_time_{rank}.dat') if rank >= 0 \
        else os.path.join(run, 'DATA', 'copy_numbers_time.dat')
    with open(path) as fh:
        header = fh.readline().strip().split(',')
        bound_cols = [i for i, h in enumerate(header) if '!' in h and h.startswith(hub + '(')]
        # the hub's first site, free or bound, counts every owned hub once
        first = header[1:][0] if False else None
        free_first = next(i for i, h in enumerate(header) if h.startswith(hub + '(') and '!' not in h)
        bound_first = bound_cols[0]
        fr = []
        for line in fh:
            row = line.strip().split(',')
            if len(row) != len(header): continue
            if float(row[0]) < burn_time: continue
            hubs = float(row[free_first]) + float(row[bound_first])
            if hubs == 0: continue
            fr.append(sum(float(row[i]) for i in bound_cols) / (sites_per_hub * hubs))
    return sum(fr) / len(fr) if fr else float('nan')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--burn', type=int, default=500)
    ap.add_argument('--tail', type=int, default=None)
    ap.add_argument('--hub', default='A')
    ap.add_argument('--sites', type=int, default=4)
    ap.add_argument('out')
    args = ap.parse_args()
    for kbdir in sorted(glob.glob(os.path.join(args.out, 'kb*'))):
        kb = os.path.basename(kbdir)[2:]
        parm = glob.glob(os.path.join(kbdir, '*_s*', '*.inp'))
        if not parm: continue
        text = open(parm[0]).read()
        dt_us = float(re.search(r'timestep\s*=\s*([^\s#]+)', text, re.I).group(1))
        nitr = int(float(re.search(r'nItr\s*=\s*([^\s#]+)', text).group(1)))
        tail = args.tail if args.tail is not None else int(round(20.0 / (float(kb) * dt_us * 1e-6)))
        lo, hi = args.burn, nitr - tail
        done = set()
        for tsv in glob.glob(os.path.join(kbdir, 'timing_*.tsv')):
            for line in open(tsv).read().strip().split('\n')[1:]:
                f = line.split('\t')
                if len(f) == 5 and f[4].strip() == '0': done.add((f[0], f[2]))
        per = collections.defaultdict(lambda: collections.defaultdict(list))  # config -> rank -> per-seed values
        for d in sorted(glob.glob(os.path.join(kbdir, '*_s*'))):
            m = re.match(r'(.+)_s(\d+)$', os.path.basename(d))
            if not m or (m.group(1), m.group(2)) not in done: continue
            config = m.group(1)
            events, nranks = dl.read_events(d)
            bonds, _, _, _ = dl.replay(events, args.hub)
            window = [b for b in bonds if lo <= b[0] <= hi]
            gains_total = sum(b[2] for b in window) or 1
            for r in range(nranks):
                rr = r if nranks > 1 else -1
                f = rank_bound_fraction(d, rr, args.hub, args.sites, lo * dt_us * 1e-6)
                mine = [b for b in window if b[4] == r]
                life = [b[1] - b[0] for b in mine]
                gain = [b[1] - b[0] for b in mine if b[2] > 0]
                other = [b[1] - b[0] for b in mine if b[2] == 0]
                # gains this rank logged: BOND events on a hub that already held a bond
                per[config][r].append((f, dl.mean_se(life)[0], dl.mean_se(gain)[0], dl.mean_se(other)[0],
                                       len(mine), sum(1 for b in mine if b[2] > 0) / max(1, len(mine))))
        if 'serial' not in per: continue
        print(f"\n### offRatekb = {kb}, per rank (bonds formed in [{lo}, {hi}])\n")
        print("| config | rank | seeds | bound fraction | t vs serial | lifetime of bonds broken here | t | gainers | t | others | t | bonds/run broken here |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        ser = per['serial'][0]
        for config in ['serial'] + sorted(c for c in per if c != 'serial'):
            for r in sorted(per[config]):
                vals = per[config][r]
                cells = [config, 'all' if config == 'serial' else str(r), str(len(vals))]
                for k, digits in ((0, 4), (1, 3), (2, 3), (3, 3)):
                    m, se, _ = dl.mean_se([v[k] for v in vals])
                    cells.append(dt.fmt(m, se, digits))
                    if config == 'serial':
                        cells.append('—')
                    else:
                        t, df, p = dt.welch([v[k] for v in vals], [v[k] for v in ser])
                        cells.append(f"{t:+.1f}")
                cells.append(f"{sum(v[4] for v in vals) / len(vals):.0f}")
                print('| ' + ' | '.join(cells) + ' |')


if __name__ == '__main__':
    main()
