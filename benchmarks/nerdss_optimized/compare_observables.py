#!/usr/bin/env python3
"""Compare seed-averaged species copy numbers between two NERDSS builds.

Reads the run tree produced by statistical_check.sh:

    <root>/<case>/{baseline,candidate}/seed<N>/DATA/copy_numbers_time.dat

For each run the second half of the trajectory is averaged, which drops the
approach to steady state and leaves the plateau.  Those per-seed averages are
then combined into a mean and standard error per species, and the two builds are
compared with a Welch z-score.

|z| below about 2 means the two builds are statistically indistinguishable for
that species at this number of seeds.  A change in the physics shows up as a
|z| that grows with the number of seeds; a pure change of random stream does not.

    compare_observables.py <root> [--window LO:HI]

--window averages rows LO*n to HI*n of each trajectory instead (fractions,
default 0.5:1).  An early window such as 0:0.1 compares the approach to steady
state, which reaction probabilities move more directly than the plateau does.
"""

import math
import statistics
import sys
from pathlib import Path

VARIANTS = ("baseline", "candidate")


def read_plateau_means(path, window=(0.5, 1.0)):
    """Average each numeric column over a window of the trajectory (the second half by default)."""
    with path.open() as handle:
        header = handle.readline().strip()
        rows = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.replace("\t", ",").split(",")
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue

    if not rows:
        return None, None

    names = [n.strip() for n in header.replace("\t", ",").split(",")]
    lo = int(len(rows) * window[0])
    hi = max(lo + 1, int(round(len(rows) * window[1])))
    half = rows[lo:hi]
    width = min(len(r) for r in half)
    # Column 0 is time.
    means = [statistics.fmean(r[c] for r in half) for c in range(1, width)]
    labels = names[1:width] if len(names) >= width else [f"col{c}" for c in range(1, width)]
    return labels, means


def welch_z(mean_a, sem_a, mean_b, sem_b):
    spread = math.sqrt(sem_a ** 2 + sem_b ** 2)
    if spread == 0.0:
        return 0.0 if mean_a == mean_b else float("inf")
    return (mean_b - mean_a) / spread


def summarize(values):
    mean = statistics.fmean(values)
    if len(values) < 2:
        return mean, 0.0
    return mean, statistics.stdev(values) / math.sqrt(len(values))


def main():
    root = Path(sys.argv[1])
    window = (0.5, 1.0)
    if len(sys.argv) > 3 and sys.argv[2] == "--window":
        lo, hi = sys.argv[3].split(":")
        window = (float(lo), float(hi))
        print(f"averaging rows {window[0]:g} to {window[1]:g} of each trajectory")
    worst = 0.0
    any_case = False

    for case_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        per_variant = {}
        labels = None
        for variant in VARIANTS:
            seed_dirs = sorted((case_dir / variant).glob("seed*"))
            series = []
            for seed_dir in seed_dirs:
                data = seed_dir / "DATA" / "copy_numbers_time.dat"
                if not data.is_file():
                    continue
                cols, means = read_plateau_means(data, window)
                if means is None:
                    continue
                labels = labels or cols
                series.append(means)
            per_variant[variant] = series

        if not all(per_variant.get(v) for v in VARIANTS):
            print(f"\n### {case_dir.name}: no usable output, skipped")
            continue

        any_case = True
        width = min(min(len(row) for row in per_variant[v]) for v in VARIANTS)

        print(f"\n### {case_dir.name}   "
              f"seeds: baseline={len(per_variant['baseline'])}, candidate={len(per_variant['candidate'])}")
        print(f"{'species':<28}{'baseline mean+-SEM':>26}{'candidate mean+-SEM':>26}{'z':>9}")

        for col in range(width):
            base_mean, base_sem = summarize([row[col] for row in per_variant["baseline"]])
            cand_mean, cand_sem = summarize([row[col] for row in per_variant["candidate"]])
            z = welch_z(base_mean, base_sem, cand_mean, cand_sem)
            if math.isfinite(z):
                worst = max(worst, abs(z))
            label = labels[col] if labels and col < len(labels) else f"col{col + 1}"
            print(f"{label:<28}{base_mean:>14.3f} +-{base_sem:>8.3f}"
                  f"{cand_mean:>14.3f} +-{cand_sem:>8.3f}{z:>9.2f}")

    if any_case:
        print(f"\nlargest |z| across all species and cases: {worst:.2f}")
        print("|z| < 2 means the two builds agree within seed-to-seed scatter.")


if __name__ == "__main__":
    main()
