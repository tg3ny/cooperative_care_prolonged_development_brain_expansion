"""Summarize the selection-strength analysis (Supporting Information, text S1(f), Table S5).
Run from anywhere: python selection_strength_check/analyze.py   (reads the result files in this folder and ../Data_S1_Summary.csv)."""
import json, os, numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
def load(name): return [json.loads(l) for l in open(os.path.join(HERE, name))]
rows = load('out_main.jsonl') + load('out_default.jsonl')
d = pd.DataFrame([{k: v for k, v in r.items() if k not in ('traj', 'phase1_end')} for r in rows])
pub = pd.read_csv(os.path.join(ROOT, 'Data_S1_Summary.csv')); pub = pub[pub.replicate <= 20]

# 1. The patched code with default switches reproduces the published replicates 1-20 (identical up to rounding in the last digit).
worst = 0.0
for g in ('control', 'treatment'):
    x = d[(d.arm == 'B') & (d.form == 'superlinear') & (d.group == g)].sort_values('rep'); y = pub[pub.group == g].sort_values('replicate')
    assert len(x) == 20 and list(x.seed) == list(y.seed)
    for a, b in (('final_B', 'phase2_B'), ('final_L', 'phase2_L'), ('final_K', 'phase2_K'), ('final_S', 'phase2_S')):
        rel = np.abs(x[a].values - y[b].values) / np.abs(y[b].values); worst = max(worst, rel.max())
        assert rel.max() < 1e-12, (g, a, rel.max())
print('Default-form runs reproduce published replicates 1-20 to floating-point rounding (largest relative difference %.1e).' % worst)

def tost(t, c, frac=0.02):
    n1, n2 = len(t), len(c); diff = t.mean() - c.mean()
    sp = np.sqrt(((n1 - 1) * t.var(ddof=1) + (n2 - 1) * c.var(ddof=1)) / (n1 + n2 - 2)); se = sp * np.sqrt(1 / n1 + 1 / n2); df = n1 + n2 - 2
    b = frac * c.mean()
    return max(1 - stats.t.cdf((diff + b) / se, df), stats.t.cdf((diff - b) / se, df))
def sem(x): return x.std(ddof=1) / np.sqrt(len(x))

out = []
for off, arm, regime in ((0.0, 'B', 'default (offset 0)'), (3.5, 'A', 'baseline mortality ~10% (offset 3.5)')):
    c = d[(d.arm == arm) & (d.group == 'control')]; assert len(c) == 20
    print(f"\n{regime}: control L {c.final_L.mean():.2f} ± {sem(c.final_L):.2f}, B {c.final_B.mean():.2f} ± {sem(c.final_B):.2f}, survival {c.final_surv_mean.mean():.3f}, "
          f"agents at bound {100 * c.final_frac_at_cap.mean():.0f}%, S {c.final_S.mean():.3f}")
    for form in ('superlinear', 'linear', 'saturating'):
        t = d[(d.arm == arm) & (d.group == 'treatment') & (d.form == form)]; assert len(t) == 20
        dL = 100 * (t.final_L.mean() / c.final_L.mean() - 1); dB = 100 * (t.final_B.mean() / c.final_B.mean() - 1)
        tt = stats.ttest_ind(t.final_L, c.final_L); tb = stats.ttest_ind(t.final_B, c.final_B); pt = tost(t.final_B, c.final_B)
        rec = dict(selection=regime, form=form, n_control=len(c), n_treatment=len(t), L_control=c.final_L.mean(), L_control_sem=sem(c.final_L), L_care=t.final_L.mean(), L_care_sem=sem(t.final_L),
                   delta_L_pct=dL, t_L=tt.statistic, p_L=tt.pvalue, B_control=c.final_B.mean(), B_control_sem=sem(c.final_B), B_care=t.final_B.mean(), B_care_sem=sem(t.final_B),
                   delta_B_pct=dB, p_diff_B=tb.pvalue, p_tost_B=pt, criterion_met=bool(dL >= 10 and abs(dB) < 2), survival_care=t.final_surv_mean.mean(),
                   agents_at_bound_pct=100 * t.final_frac_at_cap.mean(), S_care=t.final_S.mean(), survival_control=c.final_surv_mean.mean())
        out.append(rec)
        print(f"  {form:<12} care L {rec['L_care']:5.2f} ± {rec['L_care_sem']:.2f} (dL {dL:+6.1f}%)  care B {rec['B_care']:5.2f} ± {rec['B_care_sem']:.2f} (dB {dB:+5.2f}%, TOST p {pt:.2g})  "
              f"criterion met {rec['criterion_met']}  survival {rec['survival_care']:.4f}  at bound {rec['agents_at_bound_pct']:.0f}%  S {rec['S_care']:.3f}")
pd.DataFrame(out).to_csv(os.path.join(HERE, 'table_S5.csv'), index=False)
d.sort_values(['arm', 'form', 'group', 'rep']).to_csv(os.path.join(HERE, 'replicate_results.csv'), index=False)
print('\nsaved table_S5.csv and replicate_results.csv;', len(d), 'runs')
