#!/usr/bin/env python3
"""Check that the startup overlap pass excludes interface pairs, not molecule types.

checkAndResolveOverlap() used to take, for each interface of the first
molecule, the first reaction between the two molecule *types*, and test that
interface against every interface of the second molecule with that reaction's
sigma.  An interface that reacts with nothing then pushed partners away as if
it were the binding site.

Model: 300 P and 300 Q in a 60 nm box.  P has one interface p; Q has q, which
binds p at sigma 2 nm, and r, which reacts with nothing, 3 nm from Q's center
on the side away from q.  For every seed the start must
  * leave no p-q pair closer than sigma (the overlaps that do count), and,
over all seeds together,
  * leave p-r pairs closer than sigma about as often as chance places them
    there.  Before the fix r excluded p like q does and there were none.

usage: check.py <path-to-nerdss> [work-dir]     (SEEDS="1 2 ..." overrides 1-5)
"""
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

NERDSS = Path(sys.argv[1]).resolve()
WORK = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.mkdtemp())
MODEL = Path(__file__).resolve().parent
SEEDS = [int(s) for s in os.environ.get("SEEDS", "1 2 3 4 5").split()]
SIGMA = 2.0
IFACES = {"P": ["p"], "Q": ["q", "r"]}
# Chance alone puts about 12 p within sigma of an r per seed (125 over 10
# seeds in the review's measurement); the old criterion left none.
MIN_PR_PER_SEED = 4

failures = []


def check(cond, msg):
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        failures.append(msg)


def read_xyz(path):
    """[(name, com, {iface: xyz})] in molecule order from DATA/initial_crds.xyz."""
    rows = path.read_text().splitlines()[2:]
    mols, i = [], 0
    while i < len(rows):
        name, *com = rows[i].split()
        if name == "EMTY":
            break
        ifaces = {}
        for k, iface in enumerate(IFACES[name]):
            ifaces[iface] = tuple(map(float, rows[i + 1 + k].split()[1:]))
        mols.append((name, tuple(map(float, com)), ifaces))
        i += 1 + len(IFACES[name])
    return mols


def run(seed):
    d = WORK / f"seed_{seed}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    for f in ("P.mol", "Q.mol", "parms.inp"):
        shutil.copy(MODEL / f, d)
    p = subprocess.run([str(NERDSS), "-f", "parms.inp", "-s", str(seed)], cwd=d,
                       capture_output=True, text=True, timeout=600)
    log = p.stdout + p.stderr
    for f in sorted((d / "DATA").glob("output_*.txt")):
        log += f.read_text()
    (d / "log.txt").write_text(log)
    if p.returncode != 0:
        print(f"  FAIL  seed {seed}: nerdss exited {p.returncode} (see {d / 'log.txt'})")
        sys.exit(1)
    return read_xyz(d / "DATA" / "initial_crds.xyz"), log


print(f"binary: {NERDSS}\nwork:   {WORK}\n")
total_pr = 0
for seed in SEEDS:
    mols, log = run(seed)
    P = [m for m in mols if m[0] == "P"]
    Q = [m for m in mols if m[0] == "Q"]
    pq = sum(math.dist(a[2]["p"], b[2]["q"]) < SIGMA for a in P for b in Q)
    pr = sum(math.dist(a[2]["p"], b[2]["r"]) < SIGMA for a in P for b in Q)
    moved = sum(int(n) for n in re.findall(r"Overlapping Proteins \(true\): (\d+)", log))
    print(f"seed {seed}: {len(P)} P, {len(Q)} Q, re-placements {moved}, "
          f"p-q pairs < sigma {pq}, p-r pairs < sigma {pr}")
    check(len(P) == 300 and len(Q) == 300, f"seed {seed}: all 600 molecules read")
    check(pq == 0, f"seed {seed}: no p-q pair within sigma")
    total_pr += pr

check(total_pr >= MIN_PR_PER_SEED * len(SEEDS),
      f"p-r pairs within sigma over {len(SEEDS)} seeds: {total_pr} "
      f"(>= {MIN_PR_PER_SEED * len(SEEDS)}; r reacts with nothing and must not exclude p)")

print()
if failures:
    print(f"FAIL ({len(failures)})")
    sys.exit(1)
print("PASS")
