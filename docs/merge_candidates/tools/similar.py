#!/usr/bin/env python3
"""Pairwise similarity of function bodies.

  norm  = difflib ratio on identifier/number-normalized token streams (type-2 clone)
  raw   = difflib ratio on raw token streams (type-1 clone, names matter)
  lines = difflib ratio on whitespace-normalized source lines
"""
import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
ROOT = _os.environ.get('NERDSS_ROOT', _os.path.abspath(_os.path.join(HERE, '..', '..', '..')))
OUTDIR = _os.environ.get('DUPSCAN_OUT', _os.path.join(HERE, 'out'))
_os.makedirs(OUTDIR, exist_ok=True)
import json, difflib, sys
F = json.load(open(_os.path.join(OUTDIR, 'funcs.json')))
MINL = int(sys.argv[1]) if len(sys.argv) > 1 else 12
THR = float(sys.argv[2]) if len(sys.argv) > 2 else 0.55
F = [f for f in F if f['nlines'] >= MINL]
K = 6
sh = []
for f in F:
    t = f['ntoks']
    sh.append(set(tuple(t[i:i + K]) for i in range(max(1, len(t) - K + 1))))
res = []
for i in range(len(F)):
    for j in range(i + 1, len(F)):
        a, b = sh[i], sh[j]
        la, lb = len(a), len(b)
        if min(la, lb) / max(la, lb) < 0.25:
            continue
        inter = len(a & b)
        jac = inter / (la + lb - inter)
        cont = inter / min(la, lb)
        if jac < 0.25 and cont < 0.5:
            continue
        sm = difflib.SequenceMatcher(None, F[i]['ntoks'], F[j]['ntoks'], autojunk=False)
        norm = sm.ratio()
        if norm < THR and cont < 0.6:
            continue
        raw = difflib.SequenceMatcher(None, F[i]['rtoks'], F[j]['rtoks'], autojunk=False).ratio()
        lin = difflib.SequenceMatcher(None, F[i]['lines'], F[j]['lines'], autojunk=False).ratio()
        res.append((norm, raw, lin, cont, F[i], F[j]))
res.sort(key=lambda r: -(r[0] * min(r[4]['nlines'], r[5]['nlines'])))
print(f"{'norm':>5} {'raw':>5} {'lines':>5} {'cont':>5} {'nA':>4} {'nB':>4}  A  |  B")
for norm, raw, lin, cont, a, b in res:
    print(f"{norm:5.2f} {raw:5.2f} {lin:5.2f} {cont:5.2f} {a['nlines']:4d} {b['nlines']:4d}  "
          f"{a['file']}:{a['start']} {a['name']}  |  {b['file']}:{b['start']} {b['name']}")
