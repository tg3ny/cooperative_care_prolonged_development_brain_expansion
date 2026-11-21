"""
Figure-generation script for Figures S1, S2, S3, and S5 of the Supporting Information.

Produces, from the deposited data only:
  Figure S1  (full trajectories of L, B, S, and K across Simulation 1)   <- Data_S3_Trajectories.csv.gz
  Figure S2  (distributions of replicate means at the end of Simulation 1) <- Data_S1_Summary.csv
  Figure S3  (individual replicate trajectories, replicates 1 to 20)      <- Data_S3_Trajectories.csv.gz
  Figure S5  (parameter sensitivity of decoupling)                        <- Data_S2_Sensitivity.csv

In the naming of earlier manuscript versions used elsewhere in this repository, these are Figure 1-figure
supplements 1, 2, 3, and 5. Figure S4 is produced by visualize_alternative_models.py, Figure S6 by
generate_figures.py, and Figure S7 by selection_strength_check/make_figure_s7.py.

Usage, from the repository root:
    python code/generate_supplementary_figures.py --data-dir . --output-dir figures/generated

Color scheme: orange (#E69F00) = control, blue (#0072B2) = treatment.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

COLORS = {'control': '#E69F00', 'treatment': '#0072B2', 'neutral': '#666666', 'zone': '#009E73', 'outside': '#E6E6E6'}
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'savefig.facecolor': 'white', 'axes.linewidth': 1.0,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.facecolor': 'white', 'figure.facecolor': 'white',
})
LW = 2
AB = 0.25
PNL = dict(loc='left', fontweight='bold', fontsize=14)
T_PHASE1 = 10_000
GROUPS = [('control', 'Control'), ('treatment', 'Treatment')]
TRAITS = [('L_mean', 'Childhood Duration (L)'), ('B_mean', 'Brain Size (B)'),
          ('S_mean', 'Social Learning (S)'), ('K_mean', 'Knowledge (K)')]


def aggregate(traj, group, col):
    g = traj[traj['group'] == group].groupby('generation')[col]
    m = g.mean()
    sem = g.std(ddof=1) / np.sqrt(g.count())
    return m.index.values, m.values, sem.values


def save(out, name):
    Path(out).mkdir(parents=True, exist_ok=True)
    plt.savefig(f'{out}/{name}.pdf', format='pdf')
    plt.savefig(f'{out}/{name}.png', format='png')
    plt.close()
    print(f'  {name} saved')


def figure_s1(out, traj):
    """Mean trajectories (+/- SEM across 100 replicates) of the four quantities across Phases 1 and 2."""
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.2), gridspec_kw=dict(wspace=0.22, hspace=0.25))
    for ax, (col, label), panel in zip(axes.ravel(), TRAITS, 'ABCD'):
        for grp, lab in GROUPS:
            x, m, s = aggregate(traj, grp, col)
            ax.fill_between(x, m - s, m + s, alpha=AB, color=COLORS[grp], lw=0)
            ax.plot(x, m, color=COLORS[grp], lw=1.5, label=lab)
        ax.axvline(T_PHASE1, color=COLORS['neutral'], ls='--', alpha=0.6, lw=1.5)
        ax.set_xlabel('Generation')
        ax.set_ylabel(label)
        ax.set_title(panel, **PNL)
    axes[0, 0].legend(loc='lower right', frameon=False)
    save(out, 'FigureS1_simulation1_trajectories')


def figure_s2(out, summ):
    """Distributions of the replicate means at the end of Simulation 1."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), gridspec_kw=dict(wspace=0.18))
    for ax, (col, label), panel in zip(axes, [('phase2_L', 'Childhood Duration (L)'), ('phase2_B', 'Brain Size (B)')], 'AB'):
        for grp, lab in GROUPS:
            ax.hist(summ.loc[summ['group'] == grp, col], bins=15, color=COLORS[grp], alpha=0.7,
                    edgecolor='white', lw=1, label=lab)
        ax.set_xlabel(label)
        ax.set_ylabel('Frequency')
        ax.legend(loc='upper right', frameon=False)
        ax.set_title(panel, **PNL)
    save(out, 'FigureS2_simulation1_distributions')


def figure_s3(out, traj, n_shown=20):
    """Trajectories of individual replicates (replicates 1 to n_shown of each condition)."""
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), gridspec_kw=dict(wspace=0.18))
    for ax, (col, label), panel in zip(axes, TRAITS[:2], 'AB'):
        for grp, _ in GROUPS:
            sub = traj[(traj['group'] == grp) & (traj['replicate'] <= n_shown)]
            for _, run in sub.groupby('replicate'):
                ax.plot(run['generation'], run[col], color=COLORS[grp], lw=0.6, alpha=0.15)
        ax.axvline(T_PHASE1, color=COLORS['neutral'], ls='--', alpha=0.6, lw=1.5)
        ax.set_xlabel('Generation')
        ax.set_ylabel(label)
        ax.set_title(panel, **PNL)
    axes[0].legend(handles=[Line2D([0], [0], color=COLORS[g], lw=LW, label=lab) for g, lab in GROUPS],
                   loc='lower right', frameon=False)
    save(out, 'FigureS3_simulation1_replicates')


def figure_s5(out, sweep, default=(0.30, 1.8)):
    """Sensitivity sweep: childhood extension (A) and cells meeting the decoupling criterion (B)."""
    rates = np.sort(sweep['r'].unique())
    exps = np.sort(sweep['e'].unique())
    dL = sweep.pivot(index='e', columns='r', values='delta_L_percent').loc[exps, rates].values
    met = ((sweep['delta_L_percent'] >= 10) & (sweep['delta_B_percent'].abs() < 2))
    zone = sweep.assign(met=met.astype(int)).pivot(index='e', columns='r', values='met').loc[exps, rates].values
    n_met, n_all = int(met.sum()), len(sweep)
    ci = int(np.argmin(np.abs(rates - default[0])))
    ri = int(np.argmin(np.abs(exps - default[1])))

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), gridspec_kw=dict(wspace=0.3))
    ax = axes[0]
    im = ax.imshow(dL, origin='lower', aspect='auto', cmap='Blues', vmin=0)
    fig.colorbar(im, ax=ax, shrink=0.8, pad=0.04).set_label('$\\Delta$L (%)')
    ax.set_title('A  Childhood Extension ($\\Delta$L%)', loc='left', fontweight='bold', fontsize=12)
    ax = axes[1]
    ax.imshow(zone, origin='lower', aspect='auto', cmap=ListedColormap([COLORS['outside'], COLORS['zone']]), vmin=0, vmax=1)
    ax.set_title(f'B  Decoupling Zone ({n_met}/{n_all} = {100 * n_met / n_all:.1f}%)', loc='left', fontweight='bold', fontsize=12)
    ax.legend(handles=[Patch(facecolor=COLORS['outside'], edgecolor='gray', label='No decoupling'),
                       Patch(facecolor=COLORS['zone'], edgecolor='gray', label='Decoupling'),
                       Line2D([0], [0], marker='s', color='none', markerfacecolor='none', markeredgecolor='black',
                              markeredgewidth=1.5, markersize=8, label=f'Default (r = {default[0]:.2f}, e = {default[1]:.1f})')],
              loc='upper right', fontsize=7.5, framealpha=0.8)
    for ax in axes:
        ax.set_xticks(range(len(rates)))
        ax.set_xticklabels([f'{r:.2f}' for r in rates], rotation=45, ha='right')
        ax.set_yticks(range(len(exps)))
        ax.set_yticklabels([f'{e:.1f}' for e in exps])
        ax.set_xlabel('Social benefit rate (r)')
        ax.set_ylabel('Convexity exponent (e)')
        ax.plot(ci, ri, marker='s', markersize=14, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.5)
        ax.plot(ci, ri, marker='+', markersize=26, color='black', markeredgewidth=1.2)
    save(out, 'FigureS5_sensitivity')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Regenerate Figures S1, S2, S3, and S5 from the deposited data.')
    parser.add_argument('--data-dir', default='.',
                        help='Directory containing Data_S1_Summary.csv, Data_S2_Sensitivity.csv, and Data_S3_Trajectories.csv.gz')
    parser.add_argument('--output-dir', default='figures/generated', help='Directory for generated figures')
    args = parser.parse_args()
    data_dir = Path(args.data_dir)
    print('Generating supplementary figures...')
    trajectories = pd.read_csv(data_dir / 'Data_S3_Trajectories.csv.gz')
    summary = pd.read_csv(data_dir / 'Data_S1_Summary.csv')
    figure_s1(args.output_dir, trajectories)
    figure_s2(args.output_dir, summary)
    figure_s3(args.output_dir, trajectories)
    figure_s5(args.output_dir, pd.read_csv(data_dir / 'Data_S2_Sensitivity.csv'))
    print('Done.')
