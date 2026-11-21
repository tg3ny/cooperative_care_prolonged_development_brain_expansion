"""Figure S7: Simulation 1 with baseline mortality of about 10% per generation, for three forms of the social return."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os; HERE = os.path.dirname(os.path.abspath(__file__))
rows = [json.loads(l) for l in open(os.path.join(HERE, 'out_main.jsonl'))]
INK, INK2, GRID = '#111111', '#555555', '#e3e3e3'
COL = {'care': '#0072B2', 'control': '#E69F00'}      # the colors used for treatment and control in the main figures
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'font.size': 10, 'axes.labelsize': 11, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
                     'axes.edgecolor': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK})
def series(sel, idx):
    runs = [r for r in rows if sel(r)]
    gens = np.array([t[0] for t in runs[0]['traj']] + [99999]) / 1000.0
    key = 'final_B' if idx == 1 else 'final_L'
    arr = np.array([[t[idx] for t in r['traj']] + [r[key]] for r in runs])
    return gens, arr.mean(0), arr.std(0, ddof=1) / np.sqrt(len(runs)), len(runs)
forms = [('superlinear', 'Superlinear return (default)'), ('linear', 'Linear return'), ('saturating', 'Saturating return')]
fig, axes = plt.subplots(2, 3, figsize=(9.9, 5.6), sharex=True, sharey='row')
ctrl = lambda r: r['arm'] == 'A' and r['group'] == 'control'
for j, (form, title) in enumerate(forms):
    care = lambda r, f=form: r['arm'] == 'A' and r['group'] == 'treatment' and r['form'] == f
    for i, (idx, ylab, ylim) in enumerate([(2, 'Childhood Duration (L)', (0, 21.5)), (1, 'Brain Size (B)', (36, 46))]):
        ax = axes[i, j]
        ax.axvline(10, color=INK2, lw=0.9, ls=(0, (3, 3)), zorder=1)
        if i == 0: ax.axhline(20, color=INK2, lw=0.9, ls=(0, (1, 2)), zorder=1)
        for name, sel in (('control', ctrl), ('care', care)):
            g, m, se, n = series(sel, idx)
            ax.fill_between(g, m - se, m + se, color=COL[name], alpha=0.22, lw=0, zorder=2)
            ax.plot(g, m, color=COL[name], lw=2, solid_capstyle='round', zorder=3, label={'care': 'Treatment (cooperative care)', 'control': 'Control'}[name])
        ax.set_ylim(*ylim); ax.set_xlim(0, 100)
        ax.grid(axis='y', color=GRID, lw=0.8); ax.set_axisbelow(True)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        if j == 0: ax.set_ylabel(ylab)
        if i == 0: ax.set_title(title, fontsize=11, color=INK, pad=8)
        if i == 1: ax.set_xlabel(r'Generation ($\times10^{3}$)')
axes[0, 0].text(-0.2, 1.08, 'A', transform=axes[0, 0].transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
axes[1, 0].text(-0.2, 1.04, 'B', transform=axes[1, 0].transAxes, fontsize=14, fontweight='bold', va='bottom', ha='left')
axes[0, 2].annotate(r'upper limit of $L$', (97, 20), xytext=(0, -11), textcoords='offset points', color=INK2, fontsize=8, ha='right')
axes[0, 0].annotate('care and added\nmortality begin', (10, 13.2), xytext=(5, 0), textcoords='offset points', color=INK2, fontsize=8, va='center')
h, l = axes[0, 0].get_legend_handles_labels()
axes[1, 2].legend(h[::-1], l[::-1], loc='lower right', frameon=False, fontsize=9, handlelength=1.8)
fig.subplots_adjust(left=0.075, right=0.985, top=0.93, bottom=0.095, wspace=0.10, hspace=0.20)
fig.savefig(os.path.join(HERE, 'Figure_S7.png'), dpi=300, facecolor='white'); print('saved S7')
