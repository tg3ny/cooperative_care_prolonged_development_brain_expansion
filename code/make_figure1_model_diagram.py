"""Figure 1 of the revised manuscript: structure of the agent-based model and the manipulations across phases.
Usage: python code/make_figure1_model_diagram.py [output.png]   (added in release v1.1.0)"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
INK, INK2, EDGE = '#111111', '#4d4d4d', '#555555'
BLUE, BLUE_F = '#0072B2', '#e3f0f8'
VERM, VERM_F = '#D55E00', '#fbe8dc'
GRAY_F = '#f3f3f1'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'DejaVu Sans'], 'font.size': 9,
                     'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic', 'mathtext.bf': 'Liberation Sans:bold'})
fig = plt.figure(figsize=(7.2, 5.3)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 73.6); ax.axis('off')
def rbox(x, y, w, h, fc=GRAY_F, ec=EDGE, lw=1.0, r=1.2, z=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=%s' % r, fc=fc, ec=ec, lw=lw, zorder=z))
def title(x, y, s, size=10): ax.text(x, y, s, ha='left', va='top', fontsize=size, fontweight='bold', color=INK, zorder=5)
def txt(x, y, s, size=8.6, color=INK, **kw): ax.text(x, y, s, ha=kw.pop('ha', 'left'), va=kw.pop('va', 'top'), fontsize=size, color=color, zorder=5, linespacing=1.22, **kw)
def arrow(p0, p1, color=INK2, lw=1.4, style='-|>', ls='-', ms=11):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, lw=lw, color=color, linestyle=ls, zorder=4, shrinkA=0, shrinkB=0))

# ---------------- one generation ----------------
title(2, 72.4, 'One generation', 11); txt(24.5, 72.1, r'$N$ = 1,000 agents', 9, INK2)
Y1, H1 = 53.0, 14.6          # row 1
rbox(2, Y1, 26, H1); title(3.8, Y1 + H1 - 1.7, 'Heritable traits')
txt(3.8, Y1 + H1 - 6.0, r'$B$    brain size' + '\n' + r'$L$    childhood duration' + '\n' + r'$S$    social-learning propensity', 8.6)
rbox(33, Y1, 31, H1); title(34.8, Y1 + H1 - 1.7, 'Learning')
txt(34.8, Y1 + H1 - 6.0, r'Asocial:  $a\,B\,L^{q}$' + '\nSocial: knowledge of 3 demonstrators\nfrom the previous generation,\nweighted by ' + r'$S$ and $L$', 8.6)
rbox(69, Y1, 29, H1); title(70.8, Y1 + H1 - 1.7, r'Knowledge ($K$)')
txt(70.8, Y1 + H1 - 6.0, 'Sum of asocial and social\nlearning, capped at ' + r'0.6$\,B$' + '\n(storage limit)', 8.6)
arrow((28, Y1 + H1 / 2), (33, Y1 + H1 / 2)); arrow((64, Y1 + H1 / 2), (69, Y1 + H1 / 2))
Y2, H2 = 24.5, 23.5          # row 2
rbox(2, Y2, 26, H2); title(3.8, Y2 + H2 - 1.7, 'Reproduction')
txt(3.8, Y2 + H2 - 6.0, 'Survivors are parents of\n1,000 offspring; each trait\nis copied from 1 of 2\nparents, with mutation', 8.6)
rbox(33, Y2, 16, H2); title(34.8, Y2 + H2 - 1.7, 'Survival')
txt(34.8, Y2 + H2 - 6.0, 'Probability rises\nwith net fitness\n(bounded at\n0.001 and 0.999)', 8.6)
NX, NW = 54, 44
rbox(NX, Y2, NW, H2); title(NX + 1.8, Y2 + H2 - 1.7, 'Net fitness')
cy = Y2 + H2 - 6.2
txt(NX + 1.9, cy, r'+  knowledge benefit  $\lambda K$', 8.6)
txt(NX + 22.7, cy, r'$-$  brain cost  $\alpha B + \gamma B^{2}$', 8.6)
txt(NX + 22.7, cy - 4.1, r'$-$  time cost of childhood $\beta L$', 8.4)
rbox(NX + 1.0, Y2 + 1.2, 19.8, 7.6, fc=BLUE_F, ec=BLUE, lw=1.0, r=0.8, z=3)
txt(NX + 2.0, Y2 + 7.7, r'+  social return on $L$' + '\n    (only under care)', 8.6)
rbox(NX + 21.8, Y2 + 1.2, 21.2, 7.6, fc=VERM_F, ec=VERM, lw=1.0, r=0.8, z=3)
txt(NX + 22.7, Y2 + 7.7, r'$-$  penalty as $B$ nears' + '\n    the ceiling ' + r'$B_{max}$', 8.6)
arrow((83.5, Y1), (83.5, Y2 + H2)); arrow((NX, Y2 + H2 / 2), (49, Y2 + H2 / 2)); arrow((33, Y2 + H2 / 2), (28, Y2 + H2 / 2))
arrow((15, Y2 + H2), (15, Y1)); txt(16.3, (Y2 + H2 + Y1) / 2, 'next generation', 8.2, INK2, va='center')

# ---------------- manipulations across phases ----------------
title(2, 20.0, 'Manipulations across phases', 11)
y0, h0 = 10.0, 6.0
for x, w, name, gens, fc, ec in [(2, 22, 'Phase 1', 'generations 0 to 10,000', GRAY_F, EDGE), (25.5, 46, 'Phase 2', 'generations 10,001 to 100,000', BLUE_F, BLUE), (73, 25, 'Phase 3', '100,001 to 200,000', VERM_F, VERM)]:
    rbox(x, y0, w, h0, fc=fc, ec=ec, r=0.8)
    title(x + 1.3, y0 + h0 - 1.1, name, 9.2); txt(x + 1.3, y0 + 0.9, gens, 7.6, INK2, va='bottom')
txt(3.3, y0 - 1.0, 'Burn-in without care;\n' + r'$B_{max}$ = 50', 8.0)
txt(26.8, y0 - 1.0, 'Cooperative care switched on in treatment\npopulations (Simulation 1); ' + r'$B_{max}$ = 50', 8.0)
txt(74.3, y0 - 1.0, 'Ceiling raised under care\n(Simulation 2) or also\nwithout care (Simulation 3)', 8.0)
arrow((64.9, y0 + h0), (64.9, Y2 + 1.2), color=BLUE, ls=(0, (4, 2.5)), lw=1.2, ms=10)
arrow((86.4, y0 + h0), (86.4, Y2 + 1.2), color=VERM, ls=(0, (4, 2.5)), lw=1.2, ms=10)
import os, sys
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'figures', 'generated', 'Figure_1_model_diagram.png')
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
fig.savefig(out, dpi=300, facecolor='white'); print('saved', out)
