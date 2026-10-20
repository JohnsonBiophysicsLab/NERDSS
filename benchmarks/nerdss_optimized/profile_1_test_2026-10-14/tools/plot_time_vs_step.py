#!/usr/bin/env python3
"""Plot wall-clock time vs iteration from tsrun.py-timestamped nerdss stdout.
Usage: plot_time_vs_step.py out.png label1=stdout.ts [label2=stdout.ts ...]
Also writes out.csv next to out.png."""
import sys, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse(path):
    t_begin = None; t_exit = None; pts = []; avg = []
    for ln in open(path):
        if ln.startswith('#'): continue
        t, rest = ln.split(None, 1) if ' ' in ln.strip() else (ln.strip(), '')
        t = float(t); rest = rest.strip()
        if 'BEGIN SIMULATION' in rest: t_begin = t
        m = re.match(r'End iteration: (\d+),', rest)
        if m: pts.append((int(m.group(1)), t))
        m = re.match(r'Avg timestep duration: ([0-9.e+-]+),', rest)
        if m: avg.append(float(m.group(1)))
        if rest.startswith('[exit status'): t_exit = t
    return t_begin, pts, avg, t_exit

out = sys.argv[1]
runs = [(a.rsplit('=',1)[0], a.rsplit('=',1)[1]) for a in sys.argv[2:]]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
csv = open(out.rsplit('.',1)[0] + '.csv', 'w')
csv.write('run,iteration,wall_s_since_launch,avg_step_ms_last_1000\n')
colors = ['#c0392b', '#2471a3', '#1e8449', '#7d3c98']
for (label, path), c in zip(runs, colors):
    t_begin, pts, avg, t_exit = parse(path)
    if not pts: continue
    it = [0] + [p[0] for p in pts]; tt = [t_begin] + [p[1] for p in pts]
    ax1.plot(it, tt, '-', color=c, lw=1.8, label=f'{label}')
    ax1.plot([0], [t_begin], 'o', color=c, ms=5)
    ax1.axhline(t_begin, color=c, lw=0.6, ls=':')
    for (i, t), a in zip(pts, avg):
        csv.write(f'{label},{i},{t:.3f},{a*1e3:.4f}\n')
    ax2.plot([p[0] for p in pts], [a*1e3 for a in avg], '-', color=c, lw=1.6, label=label)
    end = f', total {t_exit:.0f} s' if t_exit else ''
    ax1.annotate(f'startup {t_begin:.2f} s{end}', (0, t_begin), xytext=(6, 8 if c==colors[0] else -14),
                 textcoords='offset points', fontsize=8, color=c)
ax1.set_xlabel('iteration'); ax1.set_ylabel('wall-clock time since launch (s)')
ax1.set_title('1_test (parms_0.inp), nerdss-optimized v1.2: real time vs step')
ax1.grid(alpha=0.3); ax1.legend(fontsize=9, loc='upper left')
ax2.set_xlabel('iteration'); ax2.set_ylabel('mean step duration, last 1000 steps (ms)')
ax2.set_yscale('log'); ax2.grid(alpha=0.3, which='both'); ax2.legend(fontsize=9)
ax2.set_title('per-step cost (log scale)')
fig.tight_layout(); fig.savefig(out, dpi=150)
print('wrote', out)
