"""Figures for B8 (B8.2 + B8.3) on the two-species fiber model: master against the 1.3 branch.
usage: make_figures.py <results.json> <out dir>      (needs matplotlib)"""
import json, sys, pathlib, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = json.load(open(sys.argv[1])); out = pathlib.Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#9a9892"
TXT, TXT2, GRID = "#0b0b0b", "#52514e", "#d9d8d3"
plt.rcParams.update({"font.family": "Helvetica", "font.size": 10, "axes.edgecolor": TXT2, "axes.labelcolor": TXT,
                     "xtick.color": TXT2, "ytick.color": TXT2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "axes.titleweight": "semibold"})
runs = R["runs"]; kin = R["kin"]
BUILDS = [  # key, label, colour
    ("master", "master (c507ae30)", BLUE),
    ("b13_before", "1.3 before: B8.2 only, pre-cut drops fiber pairs", ORANGE),
    ("b13_fix", "1.3 fixed: B8.2 + B8.3", AQUA),
    ("b13_grid", "1.3 fixed + grid sized by the 1D reach", AQUA),
    ("b13_fix_legacy_nob82", "legacy stream, B8.2 reverted, B8.3 fixed (= master byte for byte)", GRAY),
    ("b13_before_legacy", "legacy stream, B8.2 only", GRAY),
    ("b13_fix_legacy", "legacy stream, B8.2 + B8.3", GRAY),
]
BUILDS = [b for b in BUILDS if b[0] in runs]
STEPS = 300000
def style(ax):
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)

# ---------------- Figure 1: master vs the fixed branch --------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5)); fig.subplots_adjust(hspace=0.45, wspace=0.35)
# A. realized overlaps per 300 000-step run
ax = axs[0, 0]; style(ax)
ys = range(len(BUILDS))[::-1]
for y, (k, lab, col) in zip(ys, BUILDS):
    s = runs[k]; vals = [v["overlap_pair_frames"] for v in s["seeds"].values()]
    m = sum(vals) / len(vals); sem = s["overlaps_per_1e5_sem"] * STEPS / 1e5
    ax.barh(y, m, color=col, height=0.6); ax.errorbar(m, y, xerr=sem, fmt="none", ecolor=TXT, capsize=3, lw=1.2)
    ax.scatter(vals, [y] * len(vals), s=18, color="white", edgecolor=TXT, zorder=3, lw=0.8)
    ax.text(max(max(vals), m + sem) + 5, y, f"{m:.1f} ± {sem:.1f}", va="center", fontsize=9, color=TXT)
ax.set_yticks(list(ys)); ax.set_yticklabels([b[1] for b in BUILDS], fontsize=8.5)
ax.set_xlabel("A-site / B-site pairs closer than 0.9 σ, per 300 000-step run (4 seeds, mean ± SEM, dots = seeds)")
ax.set_title("A. Realized overlaps of a free B site with a bound A site", loc="left", fontsize=10)
ax.set_xlim(0, max(v["overlap_pair_frames"] for b in BUILDS for v in runs[b[0]]["seeds"].values()) * 1.35)
# B. separation histogram near sigma
ax = axs[0, 1]
sig = 2.0; edges = [i * 0.1 for i in range(31)]
for k, lab, col in [b for b in BUILDS if b[0] in ("master", "b13_fix")]:
    hs = [v["hist"] for v in runs[k]["seeds"].values()]; frames = sum(v["frames"] for v in runs[k]["seeds"].values())
    mean = [sum(h[i] for h in hs) / frames for i in range(30)]
    ax.stairs(mean, edges, color=col, lw=2, label=lab.split(":")[0])
ax.set_yscale("log"); ax.set_xlabel("A-site to B-site separation along the fiber (units of σ = 2 nm)")
ax.set_ylabel("pairs per frame (every 20th step)"); ax.axvspan(0, 0.9, color=ORANGE, alpha=0.12)
ax.text(0.45, ax.get_ylim()[1] * 0.4, "overlap\n(< 0.9 σ)", ha="center", fontsize=8.5, color=TXT2)
ax.annotate("bound pairs sit at exactly σ", xy=(1.05, 50), xytext=(1.5, 20), fontsize=8.5, color=TXT2,
            arrowprops=dict(arrowstyle="-", color=TXT2, lw=0.8))
ax.legend(frameon=False, fontsize=9, loc="lower right"); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.set_title("B. Separation distribution next to the exclusion zone", loc="left", fontsize=10)
# C. binding kinetics
ax = axs[1, 0]
for k, lab, col in [b for b in BUILDS if b[0] in kin]:
    d = kin[k]; t = d["t_us"]; m = d["bound_mean"]; e = d["bound_sem"]
    ax.plot(t, m, color=col, lw=2, label=f"{lab.split(':')[0]} ({d['n_seeds']} seeds)")
    ax.fill_between(t, [a - b for a, b in zip(m, e)], [a + b for a, b in zip(m, e)], color=col, alpha=0.18, lw=0)
ax.set_xlabel("time (µs; dt = 1 µs, copy numbers every step)"); ax.set_ylabel("A·B bonds (of 50 A)"); ax.set_xlim(0, 60)
ax.legend(frameon=False, fontsize=9, loc="lower right"); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
ax.set_title("C. Association kinetics over the first 60 steps, mean ± SEM over seeds", loc="left", fontsize=10)
# D. wall time, sequential 100 000-step runs
ax = axs[1, 1]; style(ax)
T = R.get("timing", {}); tb = [b for b in BUILDS if b[0] in T]
for y, (k, lab, col) in zip(range(len(tb))[::-1], tb):
    m = T[k]["mean"]; e = T[k]["sem"]
    ax.barh(y, m, color=col, height=0.6); ax.errorbar(m, y, xerr=e, fmt="none", ecolor=TXT, capsize=3)
    ax.scatter(T[k]["wall"], [y] * len(T[k]["wall"]), s=18, color="white", edgecolor=TXT, zorder=3, lw=0.8)
    ax.text(max(T[k]["wall"]) + 1.5, y, f"{m:.1f} s", va="center", fontsize=9)
ax.set_yticks(range(len(tb))[::-1]); ax.set_yticklabels([b[1].split(":")[0].replace(" + grid", "\n+ grid") for b in tb], fontsize=9)
ax.set_xlabel("wall time of a 100 000-step run (s), one run at a time, 2 repeats (dots)")
ax.set_title("D. Run time", loc="left", fontsize=10)
if tb: ax.set_xlim(0, max(max(T[b[0]]["wall"]) for b in tb) * 1.25)
fig.suptitle("B8 on the two-species fiber model (50 A + 200 B sliding on a 2000 nm fiber, excludeVolumeBound): master vs nerdss-optimized-1.3",
             fontsize=11, fontweight="semibold", x=0.02, ha="left")
fig.savefig(out / "fig1_b8_master_vs_1.3.png", dpi=150, bbox_inches="tight"); fig.savefig(out / "fig1_b8_master_vs_1.3.svg", bbox_inches="tight")

# ---------------- Figure 2: what the pre-cut dropped (B8.3) ----------------------------------
if "b13_fix_count" in runs:
    seeds = runs["b13_fix_count"]["seeds"]; c0 = next(iter(seeds.values()))
    rmax = c0["counts"]["rMaxLimit"]; nb = len(c0["hist_seen"])
    tot = lambda key: [sum(v[key][i] for v in seeds.values()) for i in range(nb)]
    seen, react, excl = tot("hist_seen"), tot("hist_react"), tot("hist_excl")
    C = {k: sum(v["counts"][k] for v in seeds.values()) for k in c0["counts"] if k != "rMaxLimit"}
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6)); fig.subplots_adjust(wspace=0.3)
    x = [i + 0.5 for i in range(nb)]
    ax = axs[0]
    ax.stairs([v / 4 for v in seen], range(nb + 1), color=GRAY, lw=1.5, label="fiber pairs offered by the cell list (per molecule pair, every step)")
    ax.stairs([v / 4 for v in excl], range(nb + 1), color=ORANGE, lw=2, label="admitted by the 1D exclusion rule, bound site vs free site (per interface pair, both blocks; the spike at σ is the bound pairs themselves)")
    ax.stairs([v / 4 for v in react], range(nb + 1), color=AQUA, lw=2, label="admitted by the 1D reaction rule, free A vs free B (only while free A remain, the first ~30 steps)")
    ax.set_yscale("log"); ax.set_ylim(1, 1e8)
    ax.axvline(rmax, color=TXT, lw=1.2, ls="--"); ax.text(rmax - 1, 5e7, f"rMaxLimit = {rmax:.1f} nm\n(3D estimate)", ha="right", va="top", fontsize=8.5)
    hi = max(i for i in range(nb) if react[i] > 0) + 1
    ax.axvspan(rmax, hi, color=AQUA, alpha=0.12)
    ax.text((rmax + hi) / 2, 5e7, "dropped by the\npre-cut before\nthe fix", ha="center", va="top", fontsize=8.5, color=TXT2)
    ax.set_xlim(0, 100); ax.set_xlabel("pair separation along the fiber |Δx| (nm), 1 nm bins"); ax.set_ylabel("pair evaluations per 300 000-step run")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.16)); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    ax.set_title("A. Where the fiber pairs the 1D rules admit sit (fixed build, instrumented, 4 seeds)", loc="left", fontsize=10)
    ax = axs[1]; style(ax)
    rows = [("fiber pairs beyond rMaxLimit\n(of all fiber pairs offered)", C["fiber_beyond"], C["fiber_seen"]),
            ("admitted by the 1D reaction rule\nwhile beyond rMaxLimit", C["react_admit_beyond"], C["react_admit"]),
            ("admitted by the 1D exclusion rule\nwhile beyond rMaxLimit", C["excl_admit_beyond"], C["excl_admit"])]
    for y, (lab, a, b) in zip(range(len(rows))[::-1], rows):
        ax.barh(y, 100 * a / b, color=AQUA, height=0.55); ax.text(100 * a / b + 0.5, y, f"{100*a/b:.1f}%  ({a:,} of {b:,})", va="center", fontsize=9)
    ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[0] for r in rows], fontsize=9); ax.set_xlim(0, 60); ax.set_xlabel("% of pair evaluations, 4 runs")
    ax.set_title("B. Share of admitted pairs the pre-cut dropped before the fix", loc="left", fontsize=10)
    fig.suptitle("B8.3: the rMaxLimit pre-cut against the 1D reach, two-species fiber model", fontsize=11, fontweight="semibold", x=0.02, y=1.04, ha="left")
    fig.savefig(out / "fig2_b8_3_precut.png", dpi=150, bbox_inches="tight"); fig.savefig(out / "fig2_b8_3_precut.svg", bbox_inches="tight")
print("wrote", sorted(p.name for p in out.iterdir()))

# ---------------- Figure 3: the cell grid sized by the 1D reach ------------------------------
if "b13_fix_count" in runs and "b13_grid_count" in runs:
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6)); fig.subplots_adjust(wspace=0.3)
    panels = [("b13_fix_count", "A. rMaxLimit from the 3D estimate (62 nm, 32 cells of 62.5 nm)"),
              ("b13_grid_count", "B. rMaxLimit from the 1D reach (82 nm, 24 cells of 83.3 nm)")]
    for ax, (k, title) in zip(axs, panels):
        seeds = runs[k]["seeds"]; c0 = next(iter(seeds.values())); rmax = c0["counts"]["rMaxLimit"]; nb = len(c0["hist_seen"])
        tot = lambda key: [sum(v[key][i] for v in seeds.values()) / len(seeds) for i in range(nb)]
        seen, react, excl = tot("hist_seen"), tot("hist_react"), tot("hist_excl")
        ax.stairs(seen, range(nb + 1), color=GRAY, lw=1.5, label="fiber pairs offered by the cell list (per molecule pair)")
        ax.stairs(excl, range(nb + 1), color=ORANGE, lw=2, label="admitted by the 1D exclusion rule")
        ax.stairs(react, range(nb + 1), color=AQUA, lw=2, label="admitted by the 1D reaction rule")
        ax.set_yscale("log"); ax.set_ylim(1, 1e8); ax.set_xlim(0, 130)
        ax.axvline(rmax, color=TXT, lw=1.2, ls="--"); ax.text(rmax + 1, 5e7, f"rMaxLimit = {rmax:.0f} nm", ha="left", va="top", fontsize=8.5)
        ax.axvline(82, color=AQUA, lw=1, ls=":"); ax.text(83, 3e3, "1D reach\n82 nm", fontsize=8, color=TXT2)
        ax.set_xlabel("pair separation along the fiber |Δx| (nm), 1 nm bins"); ax.set_ylabel("pair evaluations per 300 000-step run")
        ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True); ax.set_title(title, loc="left", fontsize=10)
        C = {kk: sum(v["counts"][kk] for v in seeds.values()) for kk in c0["counts"] if kk != "rMaxLimit"}
        ax.text(0.02, 0.05, f"admitted beyond rMaxLimit: reaction {C['react_admit_beyond']:,} of {C['react_admit']:,}, exclusion {C['excl_admit_beyond']:,} of {C['excl_admit']:,}",
                transform=ax.transAxes, fontsize=8, color=TXT2)
    axs[0].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(1.1, -0.16), ncol=3)
    fig.suptitle("Cell grid sized by the 1D reach: what the cell list offers (instrumented builds, 4 seeds, per run)", fontsize=11, fontweight="semibold", x=0.02, y=1.04, ha="left")
    fig.savefig(out / "fig3_grid_1d_reach.png", dpi=150, bbox_inches="tight"); fig.savefig(out / "fig3_grid_1d_reach.svg", bbox_inches="tight")
    print("wrote fig3")
