#!/usr/bin/env python3
"""Check that `nerdss -c` re-checks the random molecules it re-places.

After a coordinate file is applied, fixOverlappingMolecules(checkOnlyChanged)
compares every file molecule with every other molecule and re-places at random
any random molecule that collides with one.  The re-placed molecule then has to
be checked against everything else as well, or it can land on a random molecule
and the run starts with an overlap the log calls resolved.

Model: a crowded 30 nm box with 120 P and 20 Q, one interface each 0.5 nm from
the center, P(p) + Q(q) at sigma 2.  A run without -c and a run with -c use the
same random stream up to the moment the coordinate file is read, so a reference
run tells us where every random molecule will be.  The file puts the first 20 P
exactly on the 20 Q of the reference, forcing every Q to be re-placed.

For every seed the -c run must
  * leave no P/Q pair closer than sigma,
  * keep the 20 file P exactly where the file put them,
  * keep the 100 random P exactly where the reference put them (only the
    re-placed molecule moves when a re-placed molecule collides again),
  * have moved all 20 Q, and end its -c pass on "Overlapping Proteins (false)".

A last case puts two file molecules on top of each other: that stays a warning,
and neither molecule moves.

usage: check.py <path-to-nerdss> [work-dir]     (SEEDS=n scans seeds 1..n, default 30)

For nerdss_mpi, pass a wrapper script that runs `mpirun -np N nerdss_mpi "$@"`.
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
SEEDS = range(1, int(os.environ.get("SEEDS", 30)) + 1)
SIGMA = 2
BOX = 30
NPF, NPR, NQ = 20, 100, 20   # file P, random P, Q

failures = []


def check(cond, msg):
    print(("  PASS  " if cond else "  FAIL  ") + msg)
    if not cond:
        failures.append(msg)


def pdb_line(serial, resname, xyz):
    # cols 1-6 ATOM, 13-16 atom name, 18-21 residue (= molecule) name, 31-54 xyz
    x, y, z = xyz
    return f"ATOM  {serial:5d} {'COM':<4} {resname:<4}{'':9}{x:8.3f}{y:8.3f}{z:8.3f}\n"


def run(tag, np_, nq, seed, coords=None):
    d = WORK / tag
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    for f in ("P.mol", "Q.mol"):
        shutil.copy(MODEL / f, d)
    parms = (MODEL / "parms.template").read_text()
    for k, v in {"@BOX@": BOX, "@NP@": np_, "@NQ@": nq, "@SIGMA@": SIGMA}.items():
        parms = parms.replace(k, str(v))
    (d / "parms.inp").write_text(parms)
    cmd = [str(NERDSS), "-f", "parms.inp", "-s", str(seed)]
    if coords is not None:
        (d / "coords.pdb").write_text(
            "".join(pdb_line(i + 1, n, c) for i, (n, c) in enumerate(coords)))
        cmd += ["-c", "coords.pdb"]
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=600)
    log = p.stdout + p.stderr
    # nerdss_mpi sends each rank's stdout to DATA/output_<rank>.txt
    for f in sorted((d / "DATA").glob("output_*.txt")):
        log += f.read_text()
    (d / "log.txt").write_text(log)
    if p.returncode != 0:
        print(f"  FAIL  {tag}: nerdss exited {p.returncode} (see {d / 'log.txt'})")
        sys.exit(1)
    # DATA/initial_crds.xyz: two header lines, then per molecule a COM line and
    # one line per interface (one interface per molecule here)
    xyz = (d / "DATA" / "initial_crds.xyz").read_text().splitlines()[2:]
    mols = []  # (name, com, iface)
    for com_line, if_line in zip(xyz[0::2], xyz[1::2]):
        n, *c = com_line.split()
        _, *i = if_line.split()
        mols.append((n, tuple(map(float, c)), tuple(map(float, i))))
    return mols, log


def overlaps(mols):
    """(i, j) molecule-list indices of P/Q pairs whose interfaces are < sigma apart."""
    P = [i for i, m in enumerate(mols) if m[0] == "P"]
    Q = [i for i, m in enumerate(mols) if m[0] == "Q"]
    return [(i, j) for i in P for j in Q if math.dist(mols[i][2], mols[j][2]) < SIGMA]


def same(a, b, tol=1e-5):
    return all(abs(x - y) <= tol for x, y in zip(a, b))


def last_c_pass(log):
    pat = re.compile(r"^\|-Overlapping Proteins \((true|false)\): (\d+)$", re.M)
    found = pat.findall(log)
    return (found[-1][0] == "true", int(found[-1][1])) if found else None


print(f"binary: {NERDSS}\nwork:   {WORK}\n")
print(f"{BOX} nm box, {NPF} file P on top of the {NQ} random Q, {NPR} random P, sigma {SIGMA}")
left_total = 0
for seed in SEEDS:
    ref, _ = run(f"ref_{seed}", NPF + NPR, NQ, seed)
    qs = [k for k, m in enumerate(ref) if m[0] == "Q"]
    coords = [("P", tuple(round(v, 3) for v in ref[k][1])) for k in qs]
    cur, log = run(f"coord_{seed}", NPF + NPR, NQ, seed, coords=coords)
    ov = overlaps(cur)
    left_total += len(ov)
    moved_q = [k for k in qs if not same(cur[k][1], ref[k][1])]
    file_ok = all(same(cur[k][1], coords[k][1]) for k in range(NPF))
    rand_ok = all(same(cur[k][1], ref[k][1]) for k in range(NPF, NPF + NPR))
    last = last_c_pass(log)
    print(f"  seed {seed:2d}: overlaps left {len(ov)}  Q moved {len(moved_q):2d}/{NQ}  "
          f"file P in place {file_ok}  random P untouched {rand_ok}  last -c pass {last}")
    if ov:
        failures.append(f"seed {seed}: P/Q overlaps left {ov}")
    if not file_ok:
        failures.append(f"seed {seed}: a file P is not at its file position")
    if not rand_ok:
        failures.append(f"seed {seed}: a random P moved")
    if len(moved_q) != NQ:
        failures.append(f"seed {seed}: {NQ - len(moved_q)} Q under a file P were not moved")
    if last is None or last[0]:
        failures.append(f"seed {seed}: last -c pass is {last}, expected no overlaps")
check(left_total == 0, f"no P/Q overlap left over {len(SEEDS)} seeds (found {left_total})")

print("\nTwo file molecules on top of each other")
cur, log = run("file_vs_file", 1, 1, 7, coords=[("P", (10, 10, 10)), ("Q", (10, 10, 10))])
check("WARNING!!!" in log and "overlaps. Please check your input!" in log,
      "a WARNING line is printed")
check(same(cur[0][1], (10, 10, 10)) and same(cur[1][1], (10, 10, 10)),
      "both file molecules stay at (10, 10, 10)")
check(last_c_pass(log) == (False, 0), "the -c pass does not count it as an overlap")

print()
for f in failures:
    print("FAIL: " + f)
print("PASS" if not failures else f"{len(failures)} FAILURE(S)")
sys.exit(1 if failures else 0)
