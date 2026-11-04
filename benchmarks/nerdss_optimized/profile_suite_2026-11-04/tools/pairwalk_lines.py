#!/usr/bin/env python3
"""Share of all samples spent on chosen source lines, per case.

usage: pairwalk_lines.py <binary> <prof>...

Line numbers refer to the source at 9200537:
  measure_separations_to_identify_possible_reactions.cpp
    175      refresh_occupied_cells()
    199      moleculeList[target].isImplicitLipid (first load of the target Molecule)
    303-307  the neighbour scan: 13 SubBox typeMask loads per target molecule
  class_SimulVolume.cpp / .hpp: update_memberMolLists + add_member (inclusive)
"""
import re, subprocess, sys, pathlib
binary, profs = sys.argv[1], sys.argv[2:]
row = re.compile(r'^\s+(\S+)\s+(\S+)\s+(\d+):')
print('| case | samples | neighbour scan 303-307 | target load 199 | refresh 175 | update_memberMolLists incl. |')
print('| --- | ---: | ---: | ---: | ---: | ---: |')
for p in profs:
    case = pathlib.Path(p).stem
    def pp(*a):
        return subprocess.run(['google-pprof', *a, binary, p], capture_output=True, text=True).stdout
    total = int(re.search(r'Total: (\d+) samples', pp('--text')).group(1))
    scan = tgt = ref = 0
    for line in pp('--list=^measure_separations_to_identify_possible_reactions').splitlines():
        m = row.match(line)
        if not m: continue
        f = 0 if m.group(1) == '.' else int(m.group(1)); ln = int(m.group(3))
        if 303 <= ln <= 307: scan += f
        elif ln == 199: tgt += f
        elif ln == 175: ref += f
    cum = 0
    for line in pp('--text', '--cum').splitlines():
        if line.rstrip().endswith('SimulVolume::update_memberMolLists'):
            cum = float(line.split()[4].rstrip('%')); break
    pct = lambda x: f'{100.0 * x / total:.1f}%'
    print(f'| {case} | {total} | {pct(scan)} | {pct(tgt)} | {pct(ref)} | {cum:.1f}% |')
