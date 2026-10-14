#!/usr/bin/env python3
"""Check a histogram_complexes_time.dat written with histogramMolStates = true.

Usage:
  check_histogram.py --plain HIST_OFF --states HIST_ON --copy-numbers COPY_NUMBERS_ON
                     [--binding-site Type.iface ...]

Three things must hold, frame by frame:

1. Collapsing the state-aware histogram onto molecule types -- dropping the
   "(iface~S,...)" part of every name and adding the counts of members that
   then share a name -- gives exactly the plain histogram.  The parameter only
   changes how the histogram is written, not the trajectory, so the two runs
   (same seed) see the same complexes.

2. Every molecule type with a multi-state interface (read off the species
   columns of copy_numbers_time.dat) is always written with its states.

3. The state-aware histogram agrees with copy_numbers_time.dat.  A column
   Type(iface~S) counts the molecules of that type whose interface is in state
   S and unbound.  Summing count x members over the histogram's lines gives the
   same number, provided members whose interface is bound are left out.  The
   histogram does not record bonds, so --binding-site Type.iface declares that
   the interface binds, and members of such a type are then only counted on
   monomer lines; that is exact when the interface is the type's only binding
   interface, which holds for the models check.sh runs.
"""
import argparse
import re
import sys
from collections import Counter

ENTRY = re.compile(r'(\S+): (\d+)\.')
STATE_COLUMN = re.compile(r'^(\w+)\((\w+)~(\w+)\)$')


def read_histogram(path):
    """-> list of (time, Counter{composition: count}), composition = tuple of (name, n)."""
    frames = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('Time (s): '):
                frames.append((float(line[len('Time (s): '):]), Counter()))
                continue
            if not line.strip():
                continue
            count, composition = line.split('\t', 1)
            entries = tuple((name, int(n)) for name, n in ENTRY.findall(composition))
            if not entries:
                sys.exit(f'{path}: cannot parse line {line!r}')
            frames[-1][1][entries] += int(count)
    return frames


def read_copy_numbers(path):
    """-> (columns, list of (time, {column: value}))."""
    with open(path) as f:
        columns = f.readline().rstrip('\n').split(',')[1:]
        rows = []
        for line in f:
            if not line.strip():
                continue
            values = line.rstrip('\n').split(',')
            rows.append((float(values[0]), dict(zip(columns, map(int, values[1:])))))
    return columns, rows


def base_type(name):
    return name.split('(', 1)[0]


def states_of(name):
    """'A(tag~U,ser~P)' -> {'tag': 'U', 'ser': 'P'}; 'A' -> {}."""
    if '(' not in name:
        return {}
    inner = name[name.index('(') + 1:-1]
    return dict(part.split('~', 1) for part in inner.split(',') if part)


def collapse(composition):
    tally = Counter()
    for name, n in composition:
        tally[base_type(name)] += n
    return tuple(sorted(tally.items()))


def same_time(a, b):
    return abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--plain', required=True)
    ap.add_argument('--states', required=True)
    ap.add_argument('--copy-numbers', required=True)
    ap.add_argument('--binding-site', action='append', default=[], metavar='Type.iface')
    args = ap.parse_args()

    plain = read_histogram(args.plain)
    states = read_histogram(args.states)
    columns, rows = read_copy_numbers(args.copy_numbers)
    binding_sites = {tuple(s.split('.', 1)) for s in args.binding_site}

    failures = []

    # 1. collapse(states) == plain
    if len(plain) != len(states):
        failures.append(f'frame count: plain {len(plain)}, states {len(states)}')
    for (t_plain, c_plain), (t_states, c_states) in zip(plain, states):
        if not same_time(t_plain, t_states):
            failures.append(f'frame time: plain {t_plain}, states {t_states}')
            continue
        collapsed = Counter()
        for composition, count in c_states.items():
            collapsed[collapse(composition)] += count
        expected = Counter({tuple(sorted(composition)): count for composition, count in c_plain.items()})
        if collapsed != expected:
            missing = expected - collapsed
            extra = collapsed - expected
            failures.append(f't={t_plain}: collapsed states histogram differs from plain: '
                            f'missing {dict(missing)}, extra {dict(extra)}')

    # 2. every state-bearing type is written with its states
    state_ifaces = {}
    for column in columns:
        m = STATE_COLUMN.match(column)
        if m:
            state_ifaces.setdefault(m.group(1), set()).add(m.group(2))
    for t, frame in states:
        for composition in frame:
            for name, _ in composition:
                mol_type = base_type(name)
                if mol_type in state_ifaces:
                    have = set(states_of(name))
                    if have != state_ifaces[mol_type]:
                        failures.append(f't={t}: {name!r} should carry states for {sorted(state_ifaces[mol_type])}')

    # 3. copy-number columns
    checked = 0
    for t, frame in states:
        row = next((r for rt, r in rows if same_time(rt, t)), None)
        if row is None:
            failures.append(f't={t}: no copy_numbers_time.dat row at this time')
            continue
        for column in columns:
            m = STATE_COLUMN.match(column)
            if not m:
                continue
            mol_type, iface, state = m.groups()
            monomer_only = (mol_type, iface) in binding_sites
            total = 0
            for composition, count in frame.items():
                is_monomer = len(composition) == 1 and composition[0][1] == 1
                if monomer_only and not is_monomer:
                    continue
                for name, n in composition:
                    if base_type(name) == mol_type and states_of(name).get(iface) == state:
                        total += count * n
            checked += 1
            if total != row[column]:
                failures.append(f't={t}: histogram gives {total} for {column}, copy_numbers_time.dat has {row[column]}')

    print(f'frames: {len(states)}; state columns checked: {checked}; '
          f'state-bearing types: {", ".join(sorted(state_ifaces)) or "none"}')
    if not state_ifaces:
        failures.append('copy_numbers_time.dat has no Type(iface~S) column: this model has no states to check')
    if failures:
        for f in failures[:20]:
            print('FAIL:', f, file=sys.stderr)
        if len(failures) > 20:
            print(f'... and {len(failures) - 20} more', file=sys.stderr)
        sys.exit(1)
    print('PASS: states histogram collapses to the plain one and matches copy_numbers_time.dat')


if __name__ == '__main__':
    main()
