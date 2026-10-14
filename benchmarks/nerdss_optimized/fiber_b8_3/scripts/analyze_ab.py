"""Realized overlaps in the two-species fiber model (fast version: sorted + bisect).
usage: analyze_ab.py <run dir> <AB|BA> [sigma]
Molecules: 2 atoms each (COM, site).  Order BA = 200 B then 50 A; AB = 50 A then 200 B.
An A-site/B-site pair closer than 0.9 sigma is an overlap (bound pairs sit at exactly sigma).
"""
import sys, pathlib, collections, bisect
run = pathlib.Path(sys.argv[1]); order = sys.argv[2]; sigma = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
nA, nB = 50, 200; thr = 0.9 * sigma; hmax = 3 * sigma
frames = 0; overlaps = 0; frames_with = 0; hist = collections.Counter()
with (run / "DATA" / "trajectory.xyz").open() as f:
    while True:
        line = f.readline()
        if not line: break
        n = int(line.strip()); f.readline()
        xs = [float(f.readline().split()[1]) for _ in range(n)]
        sites = xs[1::2]
        ax, bx = (sites[nB:], sites[:nB]) if order == "BA" else (sites[:nA], sites[nA:])
        bx = sorted(bx); frames += 1; fo = 0
        for a in ax:
            lo = bisect.bisect_left(bx, a - hmax); hi = bisect.bisect_right(bx, a + hmax)
            for b in bx[lo:hi]:
                d = abs(a - b)
                if d < thr: fo += 1
                hist[int(d / (0.1 * sigma))] += 1
        overlaps += fo; frames_with += (fo > 0)
status = (run / "exit_status").read_text().strip() if (run / "exit_status").exists() else "?"
print(f"{run.parent.name}/{run.name}\texit={status}\tframes={frames}\toverlap_pair_frames={overlaps}\tframes_with_overlap={frames_with}\t"
      + "hist(0.1sigma bins, 0..3sigma)=" + ",".join(str(hist[k]) for k in range(30)))
