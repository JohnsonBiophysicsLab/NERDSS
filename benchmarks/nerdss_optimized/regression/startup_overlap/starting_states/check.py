#!/usr/bin/env python3
"""Check that the startup overlap pass sees the states molecules start in.

Per-state starting copy numbers ("Q : 150 (q~U), 150 (q~P)") are applied by
initialize_states(), which used to run only after generate_coordinates(), so
every molecule was placed as if its interfaces were in their first state.  A
pair that can bind only in a later state was never tested, and started inside
sigma.  A re-placed molecule also comes back from create_random_coords() in the
first states, so it has to get its own back before the next pass tests it.

Model: 300 P, all p~P, and 300 Q, the first 150 q~U and the last 150 q~P, in a
40 nm box.  Only P(p~P) + Q(q~P) binds, at sigma 2.  For every seed the start
must
  * leave no P closer than sigma to a phosphorylated Q, and, over all seeds,
  * leave P closer than sigma to an unphosphorylated Q about as often as
    chance places them there: those cannot bind and must not be excluded.

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
NQ_U = 150  # initialize_states() gives the first 150 Q, by index, q~U
# Chance alone puts about 18 P within sigma of an unphosphorylated Q per seed
# (36-38 P/Q pairs in all at this density, half of them with a q~U).
MIN_PU_PER_SEED = 6

failures = []


def check(cond, msg):
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        failures.append(msg)


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
    # DATA/initial_crds.xyz: two header lines, then per molecule, in index
    # order, a COM line and its one interface
    rows = (d / "DATA" / "initial_crds.xyz").read_text().splitlines()[2:]
    mols = [(c.split()[0], tuple(map(float, i.split()[1:]))) for c, i in zip(rows[0::2], rows[1::2])]
    return [m for m in mols if m[0] != "EMTY"], log


print(f"binary: {NERDSS}\nwork:   {WORK}\n")
total_pu = 0
for seed in SEEDS:
    mols, log = run(seed)
    P = [m[1] for m in mols if m[0] == "P"]
    Q = [m[1] for m in mols if m[0] == "Q"]
    pu = sum(math.dist(a, b) < SIGMA for a in P for b in Q[:NQ_U])
    pp = sum(math.dist(a, b) < SIGMA for a in P for b in Q[NQ_U:])
    moved = sum(int(n) for n in re.findall(r"Overlapping Proteins \(true\): (\d+)", log))
    print(f"seed {seed}: {len(P)} P, {len(Q)} Q, re-placements {moved}, "
          f"P-Q(q~P) pairs < sigma {pp}, P-Q(q~U) pairs < sigma {pu}")
    check(len(P) == 300 and len(Q) == 300, f"seed {seed}: all 600 molecules read")
    check(pp == 0, f"seed {seed}: no P within sigma of a phosphorylated Q")
    total_pu += pu

check(total_pu >= MIN_PU_PER_SEED * len(SEEDS),
      f"P-Q(q~U) pairs within sigma over {len(SEEDS)} seeds: {total_pu} "
      f"(>= {MIN_PU_PER_SEED * len(SEEDS)}; q~U binds nothing and must not exclude P)")

print()
if failures:
    print(f"FAIL ({len(failures)})")
    sys.exit(1)
print("PASS")
