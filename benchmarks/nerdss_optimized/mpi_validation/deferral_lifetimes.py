#!/usr/bin/env python3
"""Bond lifetimes from the association/dissociation log, for the deferral
measurement (proposal item 10, DEFERRAL_2026-10-30.md).

A run directory holds DATA/assoc_dissoc_time.dat (serial) or
DATA/assoc_dissoc_time_<rank>.dat for every rank (nerdss_mpi; the ranks'
records are concatenated and replayed in step order, BREAK before BOND within
a step, as a serial step orders them).  A bond is the record pair BOND ..
BREAK on the same (molecule, interface) pair; its lifetime is the step
difference.  While a bond is alive, every further BOND on the same hub
molecule is a member its complex gained; the deferred first-order pass skips
one dissociation draw for the older bonds of a complex the right neighbour
bound into, so the bonds that gained a member while alive ("gainers") are
where any bias has to show.

Only bonds formed in [burn, nItr - tail] are counted, so the ensemble is at
equilibrium and no counted bond is cut off by the end of the run.  The bound
fraction is the mean over the copy-number rows after the burn-in of the bound
hub sites over all hub sites.

Usage: deferral_lifetimes.py --burn 500 --tail 400 --hub A --sites 4
                             [--json] RUN_DIR...
"""
import argparse
import collections
import glob
import json
import math
import os
import re
import sys


def read_events(run):
    files = sorted(glob.glob(os.path.join(run, 'DATA', 'assoc_dissoc_time_*.dat')))
    if not files:
        files = [os.path.join(run, 'DATA', 'assoc_dissoc_time.dat')]
    events = []
    for rank, path in enumerate(files):
        with open(path) as fh:
            for line in fh:
                parts = line.strip().split(',')
                if len(parts) != 8:
                    continue
                itr = int(parts[0].split(':')[1])
                bond = parts[1] == 'BOND'
                m1 = (parts[2], int(parts[3]), int(parts[4]))
                m2 = (parts[5], int(parts[6]), int(parts[7]))
                events.append((itr, 0 if not bond else 1, m1, m2, rank))
    # BREAK (0) before BOND (1) within a step: a serial step draws its
    # first-order reactions before it binds, and a site freed at the top of
    # a step can bind again in the same step.
    events.sort(key=lambda e: (e[0], e[1]))
    return events, len(files)


def replay(events, hub):
    """Returns the list of finished bonds as (start, end, gains, rank that
    logged the BOND, rank that logged the BREAK) and the anomalies (a BOND on
    a live pair, a BREAK on a dead one)."""
    live = {}                      # key -> [start, gains, bond rank]
    hub_bonds = collections.defaultdict(set)   # hub (name, id) -> live keys
    bonds = []
    twice = early = 0
    for itr, kind, m1, m2, rank in events:
        key = (m1, m2) if m1 <= m2 else (m2, m1)
        hubs = [m[:2] for m in (m1, m2) if m[0] == hub]
        if kind == 1:
            if key in live:
                twice += 1
                continue
            for h in hubs:
                for other in hub_bonds[h]:
                    live[other][1] += 1
                hub_bonds[h].add(key)
            live[key] = [itr, 0, rank]
        else:
            rec = live.pop(key, None)
            if rec is None:
                early += 1
                continue
            for h in hubs:
                hub_bonds[h].discard(key)
            bonds.append((rec[0], itr, rec[1], rec[2], rank))
    return bonds, live, twice, early


def copy_numbers(run, nranks, hub, sites, burn_time):
    path = os.path.join(run, 'mergeOUT', 'copy_numbers_time.dat') if nranks > 1 \
        else os.path.join(run, 'DATA', 'copy_numbers_time.dat')
    with open(path) as fh:
        header = fh.readline().strip().split(',')
        bound_cols = [i for i, h in enumerate(header)
                      if '!' in h and h.startswith(hub + '(')]
        fracs = []
        last = None
        for line in fh:
            row = line.strip().split(',')
            if len(row) != len(header):
                continue
            t = float(row[0])
            bound = sum(float(row[i]) for i in bound_cols)
            last = bound
            if t >= burn_time:
                fracs.append(bound / sites)
    return (sum(fracs) / len(fracs) if fracs else float('nan')), len(fracs), last


def mean_se(xs):
    n = len(xs)
    if n == 0:
        return float('nan'), float('nan'), 0
    m = sum(xs) / n
    v = sum((x - m) ** 2 for x in xs) / (n - 1) if n > 1 else float('nan')
    return m, math.sqrt(v / n) if n > 1 else float('nan'), n


def analyse(run, args):
    events, nranks = read_events(run)
    bonds, live, twice, early = replay(events, args.hub)
    with open(os.path.join(run, glob.glob(os.path.join(run, '*.inp'))[0])) as fh:
        text = fh.read()
    nitr = int(float(re.search(r'nItr\s*=\s*([^\s#]+)', text).group(1)))
    dt = float(re.search(r'timestep\s*=\s*([^\s#]+)', text, re.I).group(1))
    hub_count = int(re.search(r'^\s*%s\s*:\s*(\d+)' % re.escape(args.hub), text, re.M).group(1))
    sites = hub_count * args.sites
    lo, hi = args.burn, nitr - args.tail
    window = [b for b in bonds if lo <= b[0] <= hi]
    life = [b[1] - b[0] for b in window]
    gain = [b[1] - b[0] for b in window if b[2] > 0]
    nogain = [b[1] - b[0] for b in window if b[2] == 0]
    gains = sum(b[2] for b in window)
    # gains per bond-step of exposure: the rate at which a live bond's complex
    # takes in a member
    exposure = sum(life)
    frac, nrows, last_bound = copy_numbers(run, nranks, args.hub, sites, lo * dt * 1e-6)
    m, se, n = mean_se(life)
    mg, seg, ng = mean_se(gain)
    mn, sen, nn = mean_se(nogain)
    return {
        'run': run, 'ranks': nranks, 'events': len(events), 'nItr': nitr,
        'bonds_total': len(bonds), 'alive_at_end': len(live),
        'last_copy_number_bound': last_bound, 'made_twice': twice, 'broken_early': early,
        'n': n, 'mean': m, 'se': se,
        'n_gain': ng, 'mean_gain': mg, 'se_gain': seg,
        'n_nogain': nn, 'mean_nogain': mn, 'se_nogain': sen,
        'gain_fraction': ng / n if n else float('nan'),
        'gains_per_bond_step': gains / exposure if exposure else float('nan'),
        'bound_fraction': frac, 'copy_rows': nrows,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--burn', type=int, default=500, help='first step a counted bond may form at')
    ap.add_argument('--tail', type=int, default=400, help='no counted bond forms within this many steps of the end')
    ap.add_argument('--hub', default='A')
    ap.add_argument('--sites', type=int, default=4, help='binding sites per hub molecule')
    ap.add_argument('--json', action='store_true')
    ap.add_argument('runs', nargs='+')
    args = ap.parse_args()
    out = [analyse(r, args) for r in args.runs]
    if args.json:
        json.dump(out, sys.stdout, indent=1)
        print()
        return
    for r in out:
        print(f"{r['run']}: ranks {r['ranks']}, {r['events']} events, {r['bonds_total']} bonds, "
              f"{r['alive_at_end']} alive at end (copy numbers {r['last_copy_number_bound']:.0f}), "
              f"twice {r['made_twice']} early {r['broken_early']}")
        print(f"  window n {r['n']}: lifetime {r['mean']:.3f} +- {r['se']:.3f} steps; "
              f"gainers {r['n_gain']} ({100*r['gain_fraction']:.1f}%) {r['mean_gain']:.3f} +- {r['se_gain']:.3f}; "
              f"others {r['mean_nogain']:.3f} +- {r['se_nogain']:.3f}; "
              f"gains per bond-step {r['gains_per_bond_step']:.4f}; bound fraction {r['bound_fraction']:.4f} "
              f"over {r['copy_rows']} rows")


if __name__ == '__main__':
    main()
