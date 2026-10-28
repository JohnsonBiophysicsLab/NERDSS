#!/usr/bin/env python3
"""Check that the -c overlap pass gives the same answer whichever molecule came from the file.

After a coordinate file is applied, fixOverlappingMolecules(checkOnlyChanged)
compares each file molecule, as the first molecule, with every other one, and
re-places the other one on an overlap.  checkAndResolveOverlap() used to test
only the first molecule's interfaces that react with the second's type, each
with the sigma of the first such reaction, against every interface of the
second molecule.  So the same two molecules in the same place overlapped or not
depending on which of them the file held.

Model: one P (interface p) and one Q (q binds p at sigma 1, listed first; s
binds p at sigma 2.5; r reacts with nothing).  A run without -c tells us where
the random molecules start; a run with -c puts one molecule from the file at a
chosen distance from an interface of the other, which is random and is the one
that moves on an overlap.  Each placement is run both ways round:

  p 0.3 nm from r     no overlap  (before: file P moved Q, file Q left P)
  p 1.8 nm from s     overlap     (before: file P left Q, as p held s to q's
                                   sigma 1; file Q moved P)
  p 0.3 nm from q     overlap     (control: both moved)

Every other interface pair stays outside its sigma in every case.

usage: check.py <path-to-nerdss> [work-dir]     (SEEDS="1 2 ..." overrides 1-3)
"""
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

NERDSS = Path(sys.argv[1]).resolve()
WORK = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.mkdtemp())
MODEL = Path(__file__).resolve().parent
SEEDS = [int(s) for s in os.environ.get("SEEDS", "1 2 3").split()]
IFACES = {"P": ["p"], "Q": ["q", "s", "r"]}
CASES = [  # (Q interface, distance from p, overlap?)
    ("r", 0.3, False),
    ("s", 1.8, True),
    ("q", 0.3, True),
]

failures = []


def check(cond, msg):
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        failures.append(msg)


def read_xyz(path):
    """{name: (com, {iface: xyz})} from DATA/initial_crds.xyz (one molecule per type here)."""
    rows = path.read_text().splitlines()[2:]
    mols, i = {}, 0
    while i < len(rows):
        name, *com = rows[i].split()
        if name == "EMTY":
            break
        ifaces = {f: tuple(map(float, rows[i + 1 + k].split()[1:])) for k, f in enumerate(IFACES[name])}
        mols[name] = (tuple(map(float, com)), ifaces)
        i += 1 + len(IFACES[name])
    return mols


def pdb_line(serial, resname, xyz):
    x, y, z = xyz
    return f"ATOM  {serial:5d} {'COM':<4} {resname:<4}{'':9}{x:8.3f}{y:8.3f}{z:8.3f}\n"


def run(tag, seed, coords=None):
    d = WORK / tag
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    for f in ("P.mol", "Q.mol", "parms.inp"):
        shutil.copy(MODEL / f, d)
    cmd = [str(NERDSS), "-f", "parms.inp", "-s", str(seed)]
    if coords is not None:
        (d / "coords.pdb").write_text("".join(pdb_line(i + 1, n, c) for i, (n, c) in enumerate(coords)))
        cmd += ["-c", "coords.pdb"]
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=600)
    (d / "log.txt").write_text(p.stdout + p.stderr)
    if p.returncode != 0:
        print(f"  FAIL  {tag}: nerdss exited {p.returncode} (see {d / 'log.txt'})")
        sys.exit(1)
    return read_xyz(d / "DATA" / "initial_crds.xyz")


def add(a, b, k=1.0):
    return tuple(x + k * y for x, y in zip(a, b))


def unit(v):
    n = math.sqrt(sum(x * x for x in v))
    return tuple(x / n for x in v)


print(f"binary: {NERDSS}\nwork:   {WORK}\n")
for seed in SEEDS:
    ref = run(f"seed_{seed}/reference", seed)
    pCom, pIf = ref["P"]
    qCom, qIf = ref["Q"]
    for iface, dist, overlap in CASES:
        # File P: its p goes `dist` beyond Q's interface, on the far side from Q's center.
        target = add(qIf[iface], unit(add(qIf[iface], qCom, -1)), dist)
        fileP = add(pCom, add(target, pIf["p"], -1))
        # File Q: its interface goes `dist` beyond P's p, on the far side from P's center.
        target = add(pIf["p"], unit(add(pIf["p"], pCom, -1)), dist)
        fileQ = add(qCom, add(target, qIf[iface], -1))

        outP = run(f"seed_{seed}/file_P_{iface}", seed, [("P", fileP)])
        outQ = run(f"seed_{seed}/file_Q_{iface}", seed, [("Q", fileQ)])
        movedQ = math.dist(outP["Q"][0], qCom) > 1e-6
        movedP = math.dist(outQ["P"][0], pCom) > 1e-6
        want = "moved" if overlap else "left in place"
        print(f"seed {seed}: p {dist} nm from {iface}: file P -> Q {'moved' if movedQ else 'stayed'}, "
              f"file Q -> P {'moved' if movedP else 'stayed'}")
        check(movedQ == overlap, f"seed {seed}, p {dist} nm from {iface}, file P: random Q {want}")
        check(movedP == overlap, f"seed {seed}, p {dist} nm from {iface}, file Q: random P {want}")

print()
if failures:
    print(f"FAIL ({len(failures)})")
    sys.exit(1)
print("PASS")
