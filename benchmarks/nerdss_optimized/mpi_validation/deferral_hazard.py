#!/usr/bin/env python3
"""Dissociation hazard by bond age, and lifetimes by the rank that logged the
BREAK, for the deferral runs (see DEFERRAL_2026-10-30.md).

    python3 deferral_hazard.py <runs root> [--burn 500] [--ages 1,2,3,5,8,12]

For every kb*/ directory and every configuration with exit 0 in its timing
file, pools the finished bonds of all seeds and prints
  * hazard(age) = bonds broken at that age / bonds alive at that age, for all
    bonds, for bonds whose hub gained another partner while they were alive
    (gain), for the rest (nogain), and at np > 1 for bonds whose BOND and
    BREAK were logged by different ranks (moved);
  * per rank that logged the BREAK: bonds per run, mean lifetime +- SE over
    seeds, and the hazard at the chosen ages;
  * the number of hub-steps with two BREAKs or two BONDs on the same hub,
    which the serial rule (one dissociation per molecule per step, one
    association per molecule per step) makes zero.
Bonds count if formed in [burn, nItr - 400].  Hazards are pooled over seeds
(no error bars); the per-rank lifetimes carry the seed-to-seed SE.
"""
import argparse
import collections
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deferral_lifetimes as dl  # noqa: E402


def runs_ok(kbdir):
    ok = collections.defaultdict(list)
    for t in glob.glob(os.path.join(kbdir, 'timing_*.tsv')):
        for line in open(t).read().splitlines()[1:]:
            parts = line.split('\t')
            if len(parts) == 5 and parts[4] == '0':
                ok[parts[0]].append(parts[2])
    return ok


def mean_se(xs):
    n = len(xs)
    m = sum(xs) / n
    se = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1) / n) if n > 1 else float('nan')
    return m, se


def hazards(counter, maxage):
    alive = sum(counter.values())
    out = []
    for a in range(maxage + 1):
        out.append(counter[a] / alive if alive else float('nan'))
        alive -= counter[a]
    return out


def same_step_doubles(run, hub):
    """Hub-steps with >= 2 BREAKs and with >= 2 BONDs on one hub."""
    breaks = collections.Counter()
    bonds = collections.Counter()
    for f in sorted(glob.glob(os.path.join(run, 'DATA', 'assoc_dissoc_time*.dat'))):
        for line in open(f):
            p = line.strip().split(',')
            if len(p) < 8 or not p[0].startswith('ITR:'):
                continue
            key = (int(p[0][4:]), p[3] if p[2] == hub else p[6])
            (breaks if p[1] == 'BREAK' else bonds)[key] += 1
    return (sum(1 for v in breaks.values() if v > 1),
            sum(1 for v in bonds.values() if v > 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root')
    ap.add_argument('--burn', type=int, default=500)
    ap.add_argument('--hub', default='A')
    ap.add_argument('--ages', default='1,2,3,5,8,12')
    args = ap.parse_args()
    ages = [int(a) for a in args.ages.split(',')]
    maxage = max(ages)
    for kbdir in sorted(glob.glob(os.path.join(args.root, 'kb*'))):
        print(f'## {os.path.basename(kbdir)}')
        for cfg, seeds in sorted(runs_ok(kbdir).items()):
            cnt = {k: collections.Counter() for k in ('all', 'gain', 'nogain', 'moved')}
            by_rank = collections.defaultdict(collections.Counter)
            rank_means = collections.defaultdict(list)
            rank_n = collections.defaultdict(list)
            dbl_break = dbl_bond = 0
            for seed in seeds:
                run = os.path.join(kbdir, f'{cfg}_s{seed}')
                events, nranks = dl.read_events(run)
                nitr = max(e[0] for e in events)
                bonds, live, twice, early = dl.replay(events, args.hub)
                acc = collections.defaultdict(list)
                for start, end, gains, bond_rank, break_rank in bonds:
                    if start < args.burn or start > nitr - 400:
                        continue
                    life = end - start
                    cnt['all'][life] += 1
                    cnt['gain' if gains else 'nogain'][life] += 1
                    if bond_rank != break_rank:
                        cnt['moved'][life] += 1
                    by_rank[break_rank][life] += 1
                    acc[break_rank].append(life)
                for r, lifes in acc.items():
                    rank_means[r].append(sum(lifes) / len(lifes))
                    rank_n[r].append(len(lifes))
                db, dn = same_step_doubles(run, args.hub)
                dbl_break += db
                dbl_bond += dn
            print(f'{cfg}: {len(seeds)} seeds, {sum(cnt["all"].values())} bonds; '
                  f'hub-steps with two BREAKs {dbl_break}, with two BONDs {dbl_bond}')
            for k, c in cnt.items():
                if not c:
                    continue
                h = hazards(c, maxage)
                print(f'  {k:7s} n={sum(c.values()):9d} hazard ' +
                      ' '.join(f'a{a}={h[a]:.4f}' for a in ages))
            for r in sorted(by_rank):
                m, se = mean_se(rank_means[r])
                h = hazards(by_rank[r], maxage)
                print(f'  break-rank {r}: bonds/run {sum(rank_n[r]) / len(rank_n[r]):.0f}, '
                      f'lifetime {m:.3f} +- {se:.3f}, hazard ' +
                      ' '.join(f'a{a}={h[a]:.4f}' for a in ages))


if __name__ == '__main__':
    main()
