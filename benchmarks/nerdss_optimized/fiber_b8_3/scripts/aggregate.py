"""Aggregate the B8.3 fiber benchmark into results.json.
usage: aggregate.py <scratch fiber dir> <out results.json>
Reads runs/<build>/seed*/ (300 000-step runs: overlaps.txt, trajectory.md5, wall_seconds, stderr.log counters)
and kin/<build>/seed*/DATA/copy_numbers_time.dat (3000-step binding-kinetics runs)."""
import sys, pathlib, re, json, statistics as st
root = pathlib.Path(sys.argv[1]); out = {"runs": {}, "kin": {}}
for bdir in sorted((root / "runs").iterdir()):
    seeds = {}
    for sd in sorted(bdir.glob("seed*")):
        r = {}
        ov = (sd / "overlaps.txt").read_text() if (sd / "overlaps.txt").exists() else ""
        m = re.search(r"frames=(\d+)\toverlap_pair_frames=(\d+)\tframes_with_overlap=(\d+)\thist\(0.1sigma bins, 0..3sigma\)=([\d,]+)", ov)
        if not m: continue
        r["frames"], r["overlap_pair_frames"], r["frames_with_overlap"] = map(int, m.groups()[:3])
        r["hist"] = list(map(int, m.group(4).split(",")))
        r["md5"] = (sd / "trajectory.md5").read_text().strip() if (sd / "trajectory.md5").exists() else None
        r["wall_s"] = float((sd / "wall_seconds").read_text()) if (sd / "wall_seconds").exists() else None
        err = (sd / "stderr.log").read_text() if (sd / "stderr.log").exists() else ""
        c = re.search(r"B83_COUNT rMaxLimit=([\d.]+) fiber_seen=(\d+) fiber_beyond=(\d+) react_admit=(\d+) react_admit_beyond=(\d+) excl_admit=(\d+) excl_admit_beyond=(\d+)", err)
        if c:
            r["counts"] = dict(zip(["rMaxLimit", "fiber_seen", "fiber_beyond", "react_admit", "react_admit_beyond", "excl_admit", "excl_admit_beyond"],
                                   [float(c.group(1))] + list(map(int, c.groups()[1:]))))
            for key in ("SEEN", "REACT", "EXCL"):
                h = re.search(r"B83_HIST_" + key + r"((?: \d+)+)", err)
                r["hist_" + key.lower()] = list(map(int, h.group(1).split()))
        seeds[sd.name] = r
    if not seeds: continue
    per1e5 = [v["overlap_pair_frames"] / v["frames"] * 5000 for v in seeds.values()]  # a frame every 20 steps
    summ = {"seeds": seeds, "overlaps_per_1e5_mean": st.fmean(per1e5),
            "overlaps_per_1e5_sem": st.stdev(per1e5) / len(per1e5) ** 0.5 if len(per1e5) > 1 else 0.0}
    walls = [v["wall_s"] for v in seeds.values() if v["wall_s"]]
    if walls: summ["wall_mean"] = st.fmean(walls); summ["wall_sem"] = st.stdev(walls) / len(walls) ** 0.5 if len(walls) > 1 else 0.0
    out["runs"][bdir.name] = summ
# byte-exactness: which builds match master seed for seed
ms = out["runs"].get("master", {}).get("seeds", {})
for b, summ in out["runs"].items():
    summ["matches_master"] = {s: (ms.get(s, {}).get("md5") == v["md5"]) for s, v in summ["seeds"].items()} if ms else {}
# kinetics
if (root / "kin").exists():
    for bdir in sorted((root / "kin").iterdir()):
        curves = []
        for sd in sorted(bdir.glob("seed*")):
            f = sd / "DATA" / "copy_numbers_time.dat"
            if not f.exists(): continue
            rows = [l.strip().split(",") for l in f.read_text().splitlines()[1:] if l.strip()]
            curves.append([(float(t), int(b)) for t, *_, b in rows])
        if not curves: continue
        n = min(len(c) for c in curves); times = [curves[0][i][0] for i in range(n)]
        mean = [st.fmean(c[i][1] for c in curves) for i in range(n)]
        sem = [st.stdev([c[i][1] for c in curves]) / len(curves) ** 0.5 if len(curves) > 1 else 0 for i in range(n)]
        out["kin"][bdir.name] = {"n_seeds": len(curves), "t_us": [t * 1e6 for t in times], "bound_mean": mean, "bound_sem": sem}
# sequential timing runs (one at a time, 100 000 steps, no trajectory)
out["timing"] = {}
if (root / "timing").exists():
    for bdir in sorted((root / "timing").iterdir()):
        w = [float((r / "wall_seconds").read_text()) for r in sorted(bdir.glob("rep*")) if (r / "wall_seconds").exists()]
        if w: out["timing"][bdir.name] = {"wall": w, "mean": st.fmean(w), "sem": st.stdev(w) / len(w) ** 0.5 if len(w) > 1 else 0.0}
pathlib.Path(sys.argv[2]).write_text(json.dumps(out, indent=1))
for b, s in out["runs"].items():
    print(f"{b:24s} seeds={len(s['seeds'])} overlaps/1e5={s['overlaps_per_1e5_mean']:.1f}+-{s['overlaps_per_1e5_sem']:.1f} "
          f"raw={[v['overlap_pair_frames'] for v in s['seeds'].values()]} wall={s.get('wall_mean', 0):.0f}s "
          f"=master:{sum(s['matches_master'].values())}/{len(s['matches_master'])}")
for b, t in out["timing"].items(): print(f"timing {b}: {t['wall']}")
for b, k in out["kin"].items(): print(f"kin {b}: {k['n_seeds']} seeds, bound at t[10]={k['bound_mean'][10]:.1f}")
