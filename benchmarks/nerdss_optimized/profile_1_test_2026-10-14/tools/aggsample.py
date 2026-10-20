#!/usr/bin/env python3
"""Aggregate sample(1) reports: inclusive and self counts per function.
Usage: aggsample.py <sample files...>  [--top N] [--only nerdss]"""
import sys, re, glob, collections
args = [a for a in sys.argv[1:] if not a.startswith('--')]
top = 25
only = None
for i, a in enumerate(sys.argv):
    if a == '--top': top = int(sys.argv[i+1]); args.remove(sys.argv[i+1])
    if a == '--only': only = sys.argv[i+1]; args.remove(sys.argv[i+1])
node = re.compile(r'^\s*([\+\!\:\| ]*?)\s*(\d+) (.*?)  \(in ([^)]+)\)')
incl = collections.Counter(); self_ = collections.Counter(); total = 0
def fname(s):
    s = s.split('(')[0].strip()
    return s[:70]
for f in args:
    txt = open(f).read()
    cg = txt.split('Call graph:')[1].split('Total number in stack')[0]
    lines = cg.splitlines()
    stack = []  # (depth, name)
    for ln in lines:
        mt = re.match(r'^\s+(\d+) Thread_', ln)
        if mt: total += int(mt.group(1)); stack = []; continue
        m = node.match(ln)
        if not m: continue
        depth = len(m.group(1)); cnt = int(m.group(2)); nm = fname(m.group(3)); lib = m.group(4)
        if nm.startswith('Thread_'): total += cnt; stack=[]; continue
        while stack and stack[-1][0] >= depth: stack.pop()
        if nm not in [s[1] for s in stack]:
            incl[(nm, lib)] += cnt
        stack.append((depth, nm))
        # self = cnt - children: computed later via 'Sort by top of stack' instead
    tos = txt.split('Sort by top of stack, same collapsed (when >= 5):')[1]
    for ln in tos.splitlines():
        m = re.match(r'^\s+(.*?)  \(in ([^)]+)\)\s+(\d+)$', ln)
        if m: self_[(fname(m.group(1)), m.group(2))] += int(m.group(3))
print(f"total samples {total}")
print("\n-- inclusive (nerdss functions) --")
rows = [(k, v) for k, v in incl.items() if (only is None or only in k[1])]
for (nm, lib), v in sorted(rows, key=lambda kv: -kv[1])[:top]:
    print(f"{100*v/total:6.2f}%  {v:7d}  {nm}  [{lib}]")
print("\n-- self (top of stack, all libraries) --")
for (nm, lib), v in sorted(self_.items(), key=lambda kv: -kv[1])[:top]:
    print(f"{100*v/total:6.2f}%  {v:7d}  {nm}  [{lib}]")
