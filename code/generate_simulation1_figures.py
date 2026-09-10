"""
Figure-generation script for Figures 1 and 2 (Simulation 1).

Produces:
  Figure 1  (developmental decoupling under cooperative care)
    <- Data_S3_Trajectories.csv.gz (panel A) and Data_S1_Summary.csv (panels B, C)
  Figure 2  (social-learning propensity and the model-defined learning-storage gap)
    <- Data_S3_Trajectories.csv.gz

Figures 3, 4, and Figure 3-figure supplement 1 are produced by generate_figures.py.
Figure 1-figure supplement 4 is produced by visualize_alternative_models.py.

Terminology: Phase 1 (burn-in) / Phase 2 (care divergence).
Color scheme: orange (#E69F00) = control, blue (#0072B2) = treatment.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

COLORS = {'control': '#E69F00', 'treatment': '#0072B2',
          'highlight': '#D55E00', 'neutral': '#666666'}
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


def cohens_d(a, b):
    """Cohen's d with pooled sample SD (ddof = 1) and its 95% CI (large-sample SE)."""
    na, nb = len(a), len(b)
    sp = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    d = (a.mean() - b.mean()) / sp
    se = np.sqrt((na + nb) / (na * nb) + d ** 2 / (2 * (na + nb)))
    return d, 1.96 * se


def aggregate(traj, group, col):
    g = traj[traj['group'] == group].groupby('generation')[col]
    m = g.mean()
    sem = g.std(ddof=1) / np.sqrt(g.count())
    return m.index.values, m.values, sem.values


def figure1(out, trajectory_file, summary_file):
    traj = pd.read_csv(trajectory_file)
    summ = pd.read_csv(summary_file)
    ctrl = summ[summ['group'] == 'control']
    treat = summ[summ['group'] == 'treatment']

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2),
                             gridspec_kw=dict(width_ratios=[1.3, 1.0, 0.8], wspace=0.35))

    # Panel A: childhood-duration trajectories across Phases 1 and 2
    ax = axes[0]
    for grp, lab in [('control', 'Control'), ('treatment', 'Treatment')]:
        x, m, s = aggregate(traj, grp, 'L_mean')
        ax.fill_between(x, m - s, m + s, alpha=AB, color=COLORS[grp], lw=0)
        ax.plot(x, m, color=COLORS[grp], lw=LW, label=lab)
    ax.axvline(T_PHASE1, color=COLORS['neutral'], ls='--', alpha=0.6, lw=1)
    ax.text(T_PHASE1 * 0.15, ax.get_ylim()[1] * 0.97, 'Phase 1', color=COLORS['neutral'],
            fontsize=9, va='top')
    ax.text(T_PHASE1 * 6.0, ax.get_ylim()[1] * 0.97, 'Phase 2', color=COLORS['neutral'],
            fontsize=9, va='top')
    ax.set_xlabel('Generation')
    ax.set_ylabel('Childhood Duration (L)')
    ax.legend(loc='lower right', frameon=False)
    ax.set_title('A', **PNL)

    # Panel B: end-of-Simulation-1 trait space
    ax = axes[1]
    ax.scatter(ctrl['phase2_B'], ctrl['phase2_L'], c=COLORS['control'], alpha=0.5, s=30,
               edgecolors='white', lw=0.5)
    ax.scatter(treat['phase2_B'], treat['phase2_L'], c=COLORS['treatment'], alpha=0.5, s=30,
               edgecolors='white', lw=0.5)
    cB, cL = ctrl['phase2_B'].mean(), ctrl['phase2_L'].mean()
    tB, tL = treat['phase2_B'].mean(), treat['phase2_L'].mean()
    ax.plot([cB, tB], [cL, tL], color=COLORS['highlight'], lw=2, zorder=4)
    ax.scatter([cB], [cL], c=COLORS['control'], s=180, edgecolors='black', lw=1.5, zorder=5,
               label='Control')
    ax.scatter([tB], [tL], c=COLORS['treatment'], s=180, edgecolors='black', lw=1.5, zorder=5,
               label='Treatment')
    dL = (tL - cL) / cL * 100
    dB = (tB - cB) / cB * 100
    ax.text(0.04, 0.96, f'$\\Delta L$ = {dL:+.1f}%\n$\\Delta B$ = {dB:+.1f}%',
            transform=ax.transAxes, ha='left', va='top', fontsize=9, color=COLORS['highlight'],
            fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                          edgecolor=COLORS['highlight'], lw=0.8))
    ax.set_xlabel('Brain Size (B)')
    ax.set_ylabel('Childhood Duration (L)')
    ax.legend(loc='upper right', frameon=True, fontsize=8)
    ax.set_title('B', **PNL)

    # Panel C: effect sizes with 95% CI
    ax = axes[2]
    dB_d, dB_ci = cohens_d(treat['phase2_B'], ctrl['phase2_B'])
    dL_d, dL_ci = cohens_d(treat['phase2_L'], ctrl['phase2_L'])
    ax.bar([0, 1], [dB_d, dL_d], yerr=[dB_ci, dL_ci], capsize=5,
           color=[COLORS['neutral'], COLORS['treatment']], edgecolor='black', lw=1, width=0.6)
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(0.2, color=COLORS['neutral'], ls=':', lw=1, alpha=0.7)
    ax.axhline(-0.2, color=COLORS['neutral'], ls=':', lw=1, alpha=0.7)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Brain Size\n(equivalent)', 'Childhood\n***'], fontsize=9)
    ax.set_ylabel("Cohen's d")
    ax.set_title('C', **PNL)

    Path(out).mkdir(parents=True, exist_ok=True)
    plt.savefig(f'{out}/Figure1_simulation1.pdf', format='pdf')
    plt.savefig(f'{out}/Figure1_simulation1.png', format='png')
    plt.close()
    print('  Fig 1 saved')


def figure2(out, trajectory_file):
    traj = pd.read_csv(trajectory_file)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw=dict(wspace=0.3))

    # Panel A: social-learning propensity trajectories
    ax = axes[0]
    for grp, lab in [('control', 'Control'), ('treatment', 'Treatment')]:
        x, m, s = aggregate(traj, grp, 'S_mean')
        ax.fill_between(x / 1000, m - s, m + s, alpha=AB, color=COLORS[grp], lw=0)
        ax.plot(x / 1000, m, color=COLORS[grp], lw=LW, label=lab)
    ax.axvline(T_PHASE1 / 1000, color=COLORS['neutral'], ls='--', alpha=0.6, lw=1)
    ax.text(T_PHASE1 / 1000 + 1.0, 0.47, 'Care\nintroduced', color=COLORS['neutral'],
            fontsize=8, ha='left', va='top')
    ax.set_xlabel(r'Generation ($\times 10^3$)')
    ax.set_ylabel('Social learning propensity (S)')
    ax.legend(loc='upper right', frameon=False)
    ax.set_title('A', **PNL)

    # Panel B: percentage change from the Phase 1 baseline, treatment population only
    ax = axes[1]
    tr = traj[traj['group'] == 'treatment']
    piv_S = tr.pivot(index='generation', columns='replicate', values='S_mean')
    piv_K = tr.pivot(index='generation', columns='replicate', values='K_mean')
    piv_S = piv_S[piv_S.index >= T_PHASE1]
    piv_K = piv_K[piv_K.index >= T_PHASE1]
    base_S = piv_S.iloc[0].mean()
    base_K = piv_K.iloc[0].mean()
    pct_S = (piv_S - base_S) / base_S * 100
    pct_K = (piv_K - base_K) / base_K * 100
    x = piv_S.index.values / 1000
    mS, sS = pct_S.mean(axis=1).values, pct_S.std(axis=1, ddof=1).values / np.sqrt(pct_S.shape[1])
    mK = pct_K.mean(axis=1).values
    ax.fill_between(x, mS - sS, mS + sS, alpha=AB, color=COLORS['treatment'], lw=0)
    ax.fill_between(x, mK, mS, where=mS > mK, alpha=0.12, color=COLORS['treatment'], lw=0)
    ax.plot(x, mS, color=COLORS['treatment'], lw=LW, label='Social learning (S)')
    ax.plot(x, mK, color=COLORS['highlight'], lw=LW, ls='--', label='Knowledge (K)')
    ax.axhline(0, color=COLORS['neutral'], lw=0.8, alpha=0.6)
    ax.annotate('Learning-storage\ngap', xy=(x[-120], (mS[-120] + mK[-120]) / 2 + 12),
                xytext=(x[-120], mK[-120] + 22), ha='center', va='bottom', fontsize=9,
                color=COLORS['neutral'], fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=COLORS['neutral'], lw=1))
    ax.set_xlabel(r'Generation ($\times 10^3$)')
    ax.set_ylabel('Change from Phase 1 baseline (%)')
    ax.legend(loc='upper left', frameon=False)
    ax.set_title('B', **PNL)

    Path(out).mkdir(parents=True, exist_ok=True)
    plt.savefig(f'{out}/Figure2_simulation1.pdf', format='pdf')
    plt.savefig(f'{out}/Figure2_simulation1.png', format='png')
    plt.close()
    print('  Fig 2 saved')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Regenerate Figures 1 and 2 from the deposited Simulation 1 data.')
    parser.add_argument('--data-dir', default='.',
                        help='Directory containing Data_S1_Summary.csv and Data_S3_Trajectories.csv.gz')
    parser.add_argument('--output-dir', default='figures/generated',
                        help='Directory for generated figures')
    args = parser.parse_args()
    data_dir = Path(args.data_dir)
    print('Generating Simulation 1 figures...')
    figure1(args.output_dir, data_dir / 'Data_S3_Trajectories.csv.gz', data_dir / 'Data_S1_Summary.csv')
    figure2(args.output_dir, data_dir / 'Data_S3_Trajectories.csv.gz')
    print('Done.')
