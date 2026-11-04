#!/usr/bin/env python3
"""Aggregate the per-case pprof tables written by profile_cases.sh.

usage: aggregate.py <out.md> <profiles-dir>...

Writes, for every function: self% and inclusive% per case, the unweighted mean
over cases (every case counts once, whatever its length), and in how many cases
it is among the top five by self time.  libgsl is stripped, so its static
functions are reported under the nearest exported symbol; that is why
gsl_ran_wishart appears (it is mt_get / the Mersenne Twister, see REPORT.md).
"""
import re, sys, pathlib, collections
out = pathlib.Path(sys.argv[1]); dirs = [pathlib.Path(d) for d in sys.argv[2:]]
row = re.compile(r'^\s*(\d+)\s+([\d.]+)%\s+[\d.]+%\s+(\d+)\s+([\d.]+)%\s+(.*)$')
self_, cum, cases = collections.defaultdict(dict), collections.defaultdict(dict), []
for d in dirs:
    for f in sorted(d.glob('*.flat.txt')):
        case = f.name[:-len('.flat.txt')]
        cases.append(case)
        for line in f.read_text().splitlines():
            m = row.match(line)
            if m:
                fn = m.group(5).strip()
                self_[fn][case] = self_[fn].get(case, 0) + float(m.group(2))
                cum[fn][case] = max(cum[fn].get(case, 0), float(m.group(4)))
n = len(cases)
top5 = collections.Counter()
for c in cases:
    for fn, _ in sorted(((fn, v.get(c, 0)) for fn, v in self_.items()), key=lambda t: -t[1])[:5]:
        top5[fn] += 1
mean = lambda d: sum(d.values()) / n
fns = sorted(self_, key=lambda fn: -mean(self_[fn]))
with out.open('w') as o:
    o.write(f'# Aggregated profile over {n} cases\n\nUnweighted mean over cases. "top5" = number of cases in which the function is among the five largest by self time.\n\n')
    o.write('| function | mean self % | mean inclusive % | max self % (case) | top5 |\n| --- | --- | --- | --- | --- |\n')
    for fn in fns[:45]:
        mx = max(self_[fn].items(), key=lambda t: t[1])
        o.write(f'| `{fn}` | {mean(self_[fn]):.1f} | {mean(cum[fn]):.1f} | {mx[1]:.1f} ({mx[0]}) | {top5[fn]} |\n')
    key = [fn for fn in fns[:14]]
    o.write('\n## Self % per case, largest functions\n\n| case | ' + ' | '.join(f'`{k[:28]}`' for k in key) + ' |\n|' + ' --- |' * (len(key) + 1) + '\n')
    for c in cases:
        o.write(f'| {c} | ' + ' | '.join(f'{self_[k].get(c, 0):.1f}' for k in key) + ' |\n')
    incl = ['measure_separations_to_identify_possible_reactions', 'check_bimolecular_reactions',
            'Complex::propagate', 'create_complex_propagation_vectors', 'SimulVolume::update_memberMolLists',
            'Complex::update_properties', 'gsl_ran_gaussian_ziggurat', 'QuatRotation::rotate (inline)',
            'determine_3D_bimolecular_reaction_probability', 'reflect_traj_complex_rad_rot',
            'SimulVolume::add_member (inline)', 'ShellIndex::rebin', 'remove_empty_slots']
    o.write('\n## Inclusive % per case, selected frames\n\n| case | ' + ' | '.join(f'`{k[:28]}`' for k in incl) + ' |\n|' + ' --- |' * (len(incl) + 1) + '\n')
    for c in cases:
        o.write(f'| {c} | ' + ' | '.join(f'{cum[k].get(c, 0):.1f}' for k in incl) + ' |\n')
print(f'{n} cases -> {out}')
