#!/usr/bin/env python3
"""Find verbatim duplicated runs of >= W significant lines between any two places
(across files, or across functions in the same file).  Lines are comment-stripped and
whitespace-normalized; lines that are only braces / 'else' are ignored."""
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.environ.get('NERDSS_ROOT', _os.path.abspath(_os.path.join(HERE, '..', '..', '..')))
OUTDIR = _os.environ.get('DUPSCAN_OUT', _os.path.join(HERE, 'out'))
_os.makedirs(OUTDIR, exist_ok=True)
import os, re, sys, collections
sys.path.insert(0, HERE)
from extract import strip_comments, ROOT, DIRS

W = int(sys.argv[1]) if len(sys.argv) > 1 else 8
TRIV = re.compile(r'^[{}();\s]*$|^(else|break;|continue;|return;|\}\s*else\s*\{|default:)$')

files = {}
for d in DIRS:
    for dp, dn, fn in os.walk(os.path.join(ROOT, d)):
        for f in fn:
            if f.endswith(('.cpp', '.hpp')) and f != 'json.hpp' and 'Faddeeva' not in f:
                p = os.path.join(dp, f)
                s = strip_comments(open(p, errors='replace').read())
                sig = []
                for i, l in enumerate(s.split('\n'), 1):
                    n = ' '.join(l.split())
                    if n and not TRIV.match(n) and not n.startswith('#include'):
                        sig.append((i, n))
                files[os.path.relpath(p, ROOT)] = sig

idx = collections.defaultdict(list)
for f, sig in files.items():
    for k in range(len(sig) - W + 1):
        key = tuple(x[1] for x in sig[k:k + W])
        idx[key].append((f, k))

# pair (fA,kA)-(fB,kB) -> extend to maximal runs
pairs = collections.defaultdict(set)
for key, locs in idx.items():
    if len(locs) < 2 or len(locs) > 8:
        continue
    for a in range(len(locs)):
        for b in range(a + 1, len(locs)):
            (fa, ka), (fb, kb) = locs[a], locs[b]
            if fa == fb and abs(ka - kb) < W:
                continue
            pairs[(fa, fb, ka - kb if fa != fb else ka - kb)].add((ka, kb))
runs = []
for (fa, fb, off), starts in pairs.items():
    st = sorted(starts)
    cur = None
    for ka, kb in st:
        if cur and ka <= cur[1] + 1:
            cur[1] = ka
        else:
            if cur: runs.append((fa, fb, cur[0], cur[1], cur[2]))
            cur = [ka, ka, kb]
    if cur: runs.append((fa, fb, cur[0], cur[1], cur[2]))
out = []
for fa, fb, s0, s1, kb in runs:
    n = s1 - s0 + W
    la = files[fa][s0][0]; la2 = files[fa][s1 + W - 1][0]
    lb = files[fb][kb][0]; lb2 = files[fb][kb + (s1 - s0) + W - 1][0]
    out.append((n, f"{fa}:{la}-{la2}", f"{fb}:{lb}-{lb2}"))
out.sort(key=lambda x: -x[0])
tot = 0
for n, a, b in out:
    print(f"{n:4d}  {a}  ==  {b}")
