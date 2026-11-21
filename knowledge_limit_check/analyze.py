"""Summarize the knowledge-limit analysis (Supporting Information, text S1(g), Table S6).
Run from anywhere: python knowledge_limit_check/analyze.py
Reads out.jsonl in this folder; rows for the limit 0.6B come from ../selection_strength_check/out_main.jsonl (same seeds)."""
import json, os, sys, numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
def load(path): return [json.loads(l) for l in open(path)]
rows = load(os.path.join(HERE, 'out.jsonl'))
ref = [dict(r, cap=0.6) for r in load(os.path.join(ROOT, 'selection_strength_check', 'out_main.jsonl')) if r['arm'] == 'A' and r['form'] in ('superlinear', 'saturating')]
d = pd.DataFrame([{k: v for k, v in r.items() if k not in ('traj', 'phase1_end')} for r in rows + ref])
d['L_burnin'] = [r['phase1_end']['L_mean'] for r in rows + ref]
def tost(t, c, frac=0.02):
    n1, n2 = len(t), len(c); diff = t.mean() - c.mean()
    sp = np.sqrt(((n1 - 1) * t.var(ddof=1) + (n2 - 1) * c.var(ddof=1)) / (n1 + n2 - 2)); se = sp * np.sqrt(1 / n1 + 1 / n2); df = n1 + n2 - 2
    b = frac * c.mean()
    return max(1 - stats.t.cdf((diff + b) / se, df), stats.t.cdf((diff - b) / se, df))
def sem(x): return x.std(ddof=1) / np.sqrt(len(x))
out = []; complete = True
for cap in (0.4, 0.6, 0.8):
    c = d[(d.cap == cap) & (d.group == 'control')]
    if len(c) < 3: complete = False; print(f'limit {cap}B: only {len(c)} control runs so far'); continue
    print(f"\nlimit {cap}B (constant {c.offset.iloc[0]}): control n={len(c)}  L {c.final_L.mean():.2f} ± {sem(c.final_L):.2f}  B {c.final_B.mean():.2f} ± {sem(c.final_B):.2f}  "
          f"survival {c.final_surv_mean.mean():.3f}  at bound {100 * c.final_frac_at_cap.mean():.0f}%  S {c.final_S.mean():.3f}  K/limit {(c.final_K / (cap * c.final_B)).mean():.4f}  L at end of burn-in {c.L_burnin.mean():.2f}")
    for form in ('superlinear', 'saturating'):
        t = d[(d.cap == cap) & (d.group == 'treatment') & (d.form == form)]
        if len(t) < 3: complete = False; print(f'  {form}: only {len(t)} care runs so far'); continue
        if len(t) != 20 or len(c) != 20: complete = False
        dL = 100 * (t.final_L.mean() / c.final_L.mean() - 1); dB = 100 * (t.final_B.mean() / c.final_B.mean() - 1)
        tt = stats.ttest_ind(t.final_L, c.final_L); tb = stats.ttest_ind(t.final_B, c.final_B); pt = tost(t.final_B, c.final_B)
        rec = dict(limit=cap, constant=float(c.offset.iloc[0]), form=form, n_control=len(c), n_care=len(t), L_control=c.final_L.mean(), L_control_sem=sem(c.final_L), L_care=t.final_L.mean(), L_care_sem=sem(t.final_L),
                   delta_L_pct=dL, t_L=tt.statistic, p_L=tt.pvalue, B_control=c.final_B.mean(), B_control_sem=sem(c.final_B), B_care=t.final_B.mean(), B_care_sem=sem(t.final_B), delta_B_pct=dB, p_diff_B=tb.pvalue, p_tost_B=pt,
                   criterion_met=bool(dL >= 10 and abs(dB) < 2), survival_control=c.final_surv_mean.mean(), survival_care=t.final_surv_mean.mean(), agents_at_bound_pct=100 * t.final_frac_at_cap.mean(),
                   S_control=c.final_S.mean(), S_care=t.final_S.mean(), K_over_limit_care=(t.final_K / (cap * t.final_B)).mean(), L_care_min=t.final_L.min(), L_care_max=t.final_L.max())
        out.append(rec)
        print(f"  {form:<12} n={len(t)}  care L {rec['L_care']:5.2f} ± {rec['L_care_sem']:.2f} (dL {dL:+6.1f}%, t={tt.statistic:.2f}, p={tt.pvalue:.2g})  care B {rec['B_care']:5.2f} ± {rec['B_care_sem']:.2f} (dB {dB:+5.2f}%, TOST p {pt:.2g})  "
              f"criterion met {rec['criterion_met']}  survival {rec['survival_care']:.3f}  at bound {rec['agents_at_bound_pct']:.0f}%  S {rec['S_care']:.3f}  final L range {rec['L_care_min']:.1f}-{rec['L_care_max']:.1f}")
print('\nruns read from out.jsonl:', len(rows), 'of 120', '(complete)' if complete and len(rows) == 120 else '(INCOMPLETE)')
if complete and len(rows) == 120:
    pd.DataFrame(out).to_csv(os.path.join(HERE, 'table_S6.csv'), index=False)
    d[d.cap != 0.6].sort_values(['cap', 'form', 'group', 'rep']).drop(columns=['arm']).to_csv(os.path.join(HERE, 'replicate_results.csv'), index=False)
    print('saved table_S6.csv and replicate_results.csv')
