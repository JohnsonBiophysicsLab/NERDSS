#!/usr/bin/env python3
"""Assert that every bond NERDSS reports is a bond between interfaces that are
actually touching.

A bond is written at the binding radius sigma, and a complex moves as a rigid
body, so the two interface sites of a bonded pair should stay sigma apart for as
long as the bond lives.  The one legitimate exception is a loop closure: a
reaction between two molecules that are already in the same complex bonds them
in place, without moving anything, and is allowed at any separation below
`bindRadSameCom * sigma` (bindRadSameCom defaults to 1.1).  So the invariant this
checks is: no bonded pair is further apart than `--tolerance` sigma, default 1.15,
which is just above the largest separation a legitimate loop closure can produce.

This reads `DATA/COMPLEXES/*.json`, written when `bondedComplexWrite` is set.
Each frame gives, per bonded complex, the member molecules' centres of mass, the
quaternion (w, x, y, z) that rotates each one's .mol template into its current
orientation, and the list of bonds as (member index, site name) pairs.  A site's
position is therefore

    site = com + R(q) . template_site

where `template_site` is the site's coordinate in the .mol file, taken relative
to that file's COM line.

Usage:
    check_bond_site_separation.py [run-dir] [--tolerance 1.15] [--quiet]

`run-dir` holds `DATA/COMPLEXES/` and the `.mol` files, and defaults to the
working directory.  Exits 0 when every bond is within tolerance, 1 otherwise, and
2 when there is nothing to check (no frames, which would otherwise pass
vacuously).
"""

import argparse
import json
import math
import os
import re
import sys


def parse_mol_template(path):
    """Return (molecule name, {site name: (x, y, z) relative to COM}).

    The coordinate block of a .mol file runs from its `COM` line to the `bonds`
    line, one `name x y z` per line.
    """
    name = None
    coords = {}
    in_coords = False
    for line in open(path):
        stripped = line.split('#')[0].strip()
        if not stripped:
            continue
        lowered = stripped.lower()
        if lowered.startswith('name'):
            name = stripped.split('=', 1)[1].strip()
            continue
        if lowered.startswith('com'):
            in_coords = True
        elif lowered.startswith('bonds'):
            in_coords = False
            continue
        if not in_coords:
            continue
        fields = stripped.split()
        if len(fields) != 4:
            continue
        try:
            coords[fields[0]] = tuple(float(v) for v in fields[1:])
        except ValueError:
            continue
    if name is None:
        raise ValueError(f'{path}: no "Name =" line')
    if 'COM' not in coords:
        raise ValueError(f'{path}: no COM coordinate')
    com = coords.pop('COM')
    return name, {s: tuple(c[i] - com[i] for i in range(3)) for s, c in coords.items()}


def rotate(quat, vec):
    """Rotate `vec` by the unit quaternion `quat`, given as (w, x, y, z)."""
    w, x, y, z = quat
    norm = math.sqrt(w * w + x * x + y * y + z * z)
    if norm == 0.0:
        raise ValueError('zero-length quaternion')
    w, x, y, z = (c / norm for c in (w, x, y, z))
    return (
        (1 - 2 * (y * y + z * z)) * vec[0] + 2 * (x * y - w * z) * vec[1] + 2 * (x * z + w * y) * vec[2],
        2 * (x * y + w * z) * vec[0] + (1 - 2 * (x * x + z * z)) * vec[1] + 2 * (y * z - w * x) * vec[2],
        2 * (x * z - w * y) * vec[0] + 2 * (y * z + w * x) * vec[1] + (1 - 2 * (x * x + y * y)) * vec[2],
    )


def frame_files(complexes_dir):
    """Frame files, ordered by the iteration in the file name."""
    frames = []
    for entry in os.listdir(complexes_dir):
        match = re.fullmatch(r'(\d+)\.json', entry)
        if match:
            frames.append((int(match.group(1)), os.path.join(complexes_dir, entry)))
    return sorted(frames)


def check_frame(iteration, path, templates):
    """Yield (iteration, description, gap in sigma) for every bond in one frame."""
    for complex_index, one_complex in enumerate(json.load(open(path))):
        names = one_complex['names']
        coords = one_complex['coords']
        rotations = one_complex['rotations']
        for bond in one_complex['bonds']:
            members = (bond['molindex1'], bond['molindex2'])
            sites = (bond['site1'], bond['site2'])
            sigma = one_complex['bond_types'][bond['type']]['sigma']
            positions = []
            for member, site in zip(members, sites):
                template = templates[names[member]]
                if site not in template:
                    raise KeyError(f'{path}: {names[member]} has no site "{site}"')
                offset = rotate(rotations[member], template[site])
                positions.append([coords[member][d] + offset[d] for d in range(3)])
            description = (
                f'complex {complex_index}: '
                f'{names[members[0]]}[{members[0]}].{sites[0]} -- '
                f'{names[members[1]]}[{members[1]}].{sites[1]}'
            )
            yield iteration, description, math.dist(*positions) / sigma


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run_dir', nargs='?', default='.',
                        help='directory holding DATA/COMPLEXES and the .mol files')
    parser.add_argument('--tolerance', type=float, default=1.15,
                        help='largest allowed bond site separation, in sigma (default: 1.15)')
    parser.add_argument('--quiet', action='store_true',
                        help='print only failures')
    args = parser.parse_args()

    complexes_dir = os.path.join(args.run_dir, 'DATA', 'COMPLEXES')
    if not os.path.isdir(complexes_dir):
        print(f'{complexes_dir}: not a directory -- was bondedComplexWrite set?', file=sys.stderr)
        return 2

    templates = {}
    for entry in sorted(os.listdir(args.run_dir)):
        if entry.endswith('.mol'):
            name, coords = parse_mol_template(os.path.join(args.run_dir, entry))
            templates[name] = coords
    if not templates:
        print(f'{args.run_dir}: no .mol files', file=sys.stderr)
        return 2

    frames = frame_files(complexes_dir)
    if not frames:
        print(f'{complexes_dir}: no frames to check', file=sys.stderr)
        return 2

    worst = 0.0
    n_bonds = 0
    offenders = []
    for iteration, path in frames:
        for iteration, description, gap in check_frame(iteration, path, templates):
            n_bonds += 1
            worst = max(worst, gap)
            if gap > args.tolerance:
                offenders.append((iteration, description, gap))

    if offenders:
        iteration, description, gap = offenders[0]
        print(f'FAIL {args.run_dir}: {len(offenders)} of {n_bonds} bonds exceed '
              f'{args.tolerance} sigma over {len(frames)} frames; worst {worst:.3f} sigma')
        print(f'  first at iteration {iteration}, {description}, {gap:.3f} sigma')
        return 1

    if not args.quiet:
        print(f'ok {args.run_dir}: {n_bonds} bonds over {len(frames)} frames, '
              f'worst {worst:.3f} sigma (limit {args.tolerance})')
    return 0


if __name__ == '__main__':
    sys.exit(main())
