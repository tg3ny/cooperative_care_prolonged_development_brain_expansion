#!/usr/bin/env python3
"""Recompute the statistics reported in the article and its Supporting Information from the deposited data.

Run from the repository root:

    python code/reported_statistics.py

The script reads only files deposited in this repository and reruns no simulation. It prints the values in the
order in which they appear in the article, with the place where each is reported. Conventions, as in the article:
SEM and SD use the sample standard deviation (n - 1); Cohen's d uses the pooled sample SD; t tests are two-sided
independent-samples tests with equal variances; equivalence of brain size is tested with two one-sided tests
(TOST) with bounds of +/-2% of the control mean. Net fitness and survival "at the replicate means" evaluate the
model equations at the population means deposited for each replicate.

Statistics already printed by other scripts are not repeated here: Tables S5 and S6 and the text that summarizes
them come from selection_strength_check/analyze.py and knowledge_limit_check/analyze.py; replicate-level details
of those two analyses that the text reports in addition are recomputed below.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BOUND_NF = float(np.log(999.0))          # net fitness at which survival reaches its upper bound of 0.999


def path(*parts):
    return os.path.join(ROOT, *parts)


def sem(x):
    x = np.asarray(x, float)
    return x.std(ddof=1) / np.sqrt(len(x))


def cohens_d(t, c):
    t, c = np.asarray(t, float), np.asarray(c, float)
    sp = np.sqrt(((len(t) - 1) * t.var(ddof=1) + (len(c) - 1) * c.var(ddof=1)) / (len(t) + len(c) - 2))
    return (t.mean() - c.mean()) / sp


def diff_se_df(t, c):
    t, c = np.asarray(t, float), np.asarray(c, float)
    n1, n2 = len(t), len(c)
    sp = np.sqrt(((n1 - 1) * t.var(ddof=1) + (n2 - 1) * c.var(ddof=1)) / (n1 + n2 - 2))
    return t.mean() - c.mean(), sp * np.sqrt(1 / n1 + 1 / n2), n1 + n2 - 2


def tost(t, c, frac=0.02):
    diff, se, df = diff_se_df(t, c)
    b = frac * np.mean(c)
    return max(1 - stats.t.cdf((diff + b) / se, df), stats.t.cdf((diff - b) / se, df))


def ci(t, c, level):
    diff, se, df = diff_se_df(t, c)
    h = stats.t.ppf(0.5 + level / 2, df) * se
    return diff - h, diff + h


def penalty(B, Bmax):
    return 10.0 / (1.0 + np.exp(-0.5 * (B - Bmax)))


def social_benefit(L, care, form="superlinear", r=0.30, e=1.8, Lmax=20.0, B=None):
    if form == "superlinear":
        return care * r * (L / Lmax) ** e * Lmax
    if form == "linear":
        return care * r * L
    if form == "saturating":
        return care * r * Lmax * L / (L + 5.0)
    if form == "weak":
        return care * r * (L / Lmax) ** e * Lmax * (B / 50.0) ** 0.2
    if form == "null":                    # linear form with L replaced by a uniform draw on [0, 20]: expected value
        return care * r * (Lmax / 2.0)
    raise ValueError(form)


def net_fitness(B, L, K, care, Bmax=50.0, form="superlinear"):
    return 0.5 * K + social_benefit(L, care, form, B=B) - 0.1 * B - 0.001 * B ** 2 - penalty(B, Bmax) - 0.12 * L


def survival(nf):
    return np.clip(1.0 / (1.0 + np.exp(-nf)), 0.001, 0.999)


def compare(label, t, c, tost_too=False):
    tt = stats.ttest_ind(t, c)
    line = (f"  {label:<34} control {np.mean(c):.3f} +/- {sem(c):.3f} (SD {np.std(c, ddof=1):.3f})   "
            f"treatment {np.mean(t):.3f} +/- {sem(t):.3f} (SD {np.std(t, ddof=1):.3f})   "
            f"difference {100 * (np.mean(t) / np.mean(c) - 1):+.1f}%   t({len(t) + len(c) - 2}) = {tt.statistic:.2f}, "
            f"p = {tt.pvalue:.2g}, d = {cohens_d(t, c):.2f}")
    if tost_too:
        line += f"   TOST p = {tost(t, c):.3f}"
    print(line)


def heading(text):
    print("\n" + text + "\n" + "-" * len(text))


# ----------------------------------------------------------------------------------------------------------------------
# Simulation 1
# ----------------------------------------------------------------------------------------------------------------------
heading("Simulation 1, end of Phase 2 (Results; Table 2; Table S2; Figure 2)")
s1 = pd.read_csv(path("Data_S1_Summary.csv"))
c1, t1 = s1[s1.group == "control"], s1[s1.group == "treatment"]
for trait, name in (("L", "Childhood duration (L)"), ("B", "Brain size (B)"), ("S", "Social-learning propensity (S)"), ("K", "Knowledge (K)")):
    compare(name, t1["phase2_" + trait].values, c1["phase2_" + trait].values, tost_too=(trait == "B"))
b_t, b_c = t1.phase2_B.values, c1.phase2_B.values
print(f"  Brain size: equivalence bounds +/-{0.02 * b_c.mean():.2f}; 90% CI of the difference {ci(b_t, b_c, 0.90)[0]:+.2f} to {ci(b_t, b_c, 0.90)[1]:+.2f}; "
      f"99% CI {ci(b_t, b_c, 0.99)[0]:+.2f} to {ci(b_t, b_c, 0.99)[1]:+.2f}  (Methods, Statistical analysis; Results)")

heading("Net fitness and survival at the replicate means, end of Simulation 1 (Methods; text S1(f))")
print(f"  Net fitness at which survival reaches the bound of 0.999: {BOUND_NF:.2f}")
for g, care, x in (("control", 0.0, c1), ("care", 0.6, t1)):
    nf = net_fitness(x.phase2_B.values, x.phase2_L.values, x.phase2_K.values, care)
    print(f"  {g:<8} mean net fitness {nf.mean():.2f}; survival at that value {1 / (1 + np.exp(-nf.mean())):.3f}; "
          f"mean social benefit {social_benefit(x.phase2_L.values, care).mean():.2f}")

heading("Simulation 1 trajectories (Figure 3 legend; text S1(c))")
tr = pd.read_csv(path("Data_S3_Trajectories.csv.gz"))
piv = {v: tr.pivot_table(index="generation", columns=["group", "replicate"], values=v) for v in ("L_mean", "S_mean", "B_mean")}
gens = piv["L_mean"].index.values
S_c, S_t = piv["S_mean"]["control"].mean(axis=1), piv["S_mean"]["treatment"].mean(axis=1)
dS = 100 * (S_t / S_c - 1)
ph2 = gens >= 10000
print(f"  Difference in S between conditions: peak {dS[ph2].max():+.1f}% (generation {int(dS[ph2].idxmax()):,}); at the last sampled generation {dS.iloc[-1]:+.1f}%")
L_c, L_t = piv["L_mean"]["control"], piv["L_mean"]["treatment"]
pvals = np.array([stats.ttest_ind(L_t.loc[g].values, L_c.loc[g].values).pvalue for g in gens])
dL = 100 * (L_t.mean(axis=1) / L_c.mean(axis=1) - 1).values
sig = pvals < 0.005
after = gens >= 10000
first_sig = gens[after & sig][0]
last_nonsig = gens[after & ~sig][-1]
last_below10 = gens[after & (dL < 10)][-1]
print(f"  Difference in L first significant (p < 0.005) at generation {first_sig:,} ({first_sig - 10000:,} generations after care is introduced); "
      f"significant at every sampled generation after {last_nonsig:,} ({last_nonsig - 10000:,} after introduction)")
print(f"  Difference in L at least 10% at every sampled generation after {last_below10:,} ({last_below10 - 10000:,} after introduction)")
print(f"  Difference in L at generation 90,000: {dL[list(gens).index(90000)]:+.1f}%; smallest difference at sampled generations 90,000 to 99,900: {dL[gens >= 90000].min():+.1f}%; "
      f"at the final generation (Data_S1_Summary.csv): {100 * (t1.phase2_L.mean() / c1.phase2_L.mean() - 1):+.1f}%")

# ----------------------------------------------------------------------------------------------------------------------
# Alternative specifications and sensitivity sweep
# ----------------------------------------------------------------------------------------------------------------------
heading("Alternative-model analysis (Results; text S1(d) and S1(f); Figure S4)")
am = pd.read_csv(path("Data_S4_Alternative_Models.csv"))
def form_of(model):
    name = model.upper()
    for key, form in (("WEAK", "weak"), ("NULL", "null"), ("SATURATING", "saturating"), ("SUPERLINEAR", "superlinear"), ("LINEAR", "linear")):
        if key in name:
            return form
    raise ValueError(model)
for _, row in am.iterrows():
    form = form_of(row.model)
    L, B = row.L_treatment_mean, row.B_treatment_mean
    ben = social_benefit(L, 0.6, form, B=B)
    nf = 0.5 * 0.6 * B + ben - 0.1 * B - 0.001 * B ** 2 - penalty(B, 50.0) - 0.12 * L       # knowledge at its limit, 0.6B
    print(f"  {row.model:<34} dL {row.delta_L_pct:+.1f}%   dB {row.delta_B_pct:+.2f}%   TOST p = {row.p_tost_B:.4f}   "
          f"benefit at the mean traits of care populations {ben:.2f}; net fitness {nf:.2f}")

heading("Sensitivity sweep (Results; text S1(b) and S1(d); Figure S5; Table S3)")
sw = pd.read_csv(path("Data_S2_Sensitivity.csv"))
okL, okB = sw.delta_L_percent >= 10, sw.delta_B_percent.abs() < 2
print(f"  Cells meeting the criterion (dL >= 10% and |dB| < 2%): {int((okL & okB).sum())} of {len(sw)} ({100 * (okL & okB).mean():.1f}%)")
print(f"  Of the {int((~(okL & okB)).sum())} others: {int((~okL & okB).sum())} fail on childhood extension only, {int((okL & ~okB).sum())} on brain size only, {int((~okL & ~okB).sum())} on both")
print(f"  Brain size lower under care in {int((sw.delta_B_percent[~okB] < 0).sum())} of the {int((~okB).sum())} cells that fail on brain size, and in {int((sw.delta_B_percent < 0).sum())} of {len(sw)} cells overall")
print(f"  Cells with e > 1 meeting the criterion: {int((okL & okB & (sw.e > 1.0)).sum())} of {int((sw.e > 1.0).sum())}; cells with e = 1.0 meeting it: {int((okL & okB & np.isclose(sw.e, 1.0)).sum())}")
for bound in (1.5, 3.0):
    print(f"  With |dB| < {bound}%: {int((okL & (sw.delta_B_percent.abs() < bound)).sum())} cells meet the criterion")

# ----------------------------------------------------------------------------------------------------------------------
# Simulations 2 and 3
# ----------------------------------------------------------------------------------------------------------------------
heading("Simulation 2, end of Phase 3 (Results; Table 3; Figure 4)")
s2 = pd.read_csv(path("Data_S6_Pure_Summary.csv"))
c2, t2 = s2[s2.group == "control"], s2[s2.group == "treatment"]
a, b = s2.sort_values(["group", "replicate"]), s1.sort_values(["group", "replicate"])
worst = max(np.max(np.abs(a["phase2_" + v].values - b["phase2_" + v].values) / np.abs(b["phase2_" + v].values)) for v in "BLSK")
print(f"  Phase 2 endpoints equal to those of Simulation 1 in every replicate (same seeds: {list(a.seed) == list(b.seed)}; largest relative difference {worst:.1e})")
print(f"  Treatment B: {t2.phase2_B.mean():.2f} -> {t2.phase3_B.mean():.2f} (change {t2.delta_B.mean():+.2f}, {100 * t2.delta_B.mean() / t2.phase2_B.mean():+.1f}%)")
tt = stats.ttest_ind(t2.delta_B, c2.delta_B)
print(f"  Change in B, treatment vs control: t({len(s2) - 2}) = {tt.statistic:.2f}, p = {tt.pvalue:.2g}, d = {cohens_d(t2.delta_B, c2.delta_B):.2f}")
one = stats.ttest_1samp(c2.delta_B, 0.0)
print(f"  Control change in B: {c2.delta_B.mean():+.2f} (one-sample t test against zero, p = {one.pvalue:.2f}); control change in L: {c2.delta_L.mean():+.2f} +/- {sem(c2.delta_L):.2f} "
      f"({c2.phase2_L.mean():.2f} to {c2.phase3_L.mean():.2f})")
for trait, name in (("B", "Brain size (B)"), ("L", "Childhood duration (L)"), ("S", "Social-learning propensity (S)"), ("K", "Knowledge (K)")):
    compare(name, t2["phase3_" + trait].values, c2["phase3_" + trait].values)
print(f"  K / 0.6B at the endpoint: control {(c2.phase3_K / (0.6 * c2.phase3_B)).mean():.4f}, treatment {(t2.phase3_K / (0.6 * t2.phase3_B)).mean():.4f}")
tr2 = pd.concat([pd.read_csv(path(f"experiment2_{g}_trajectories.csv.gz"), usecols=["group", "replicate", "generation", "B_mean"]) for g in ("control", "treatment")])
m2 = tr2.groupby(["group", "generation"]).B_mean.mean().unstack(0)
last = m2[m2.index >= 190000]
print(f"  Treatment B exceeds control B at every sampled generation from 190,000 on: {bool((last.treatment > last.control).all())}; "
      f"mean treatment B at generations 150,000, 175,000, and 199,900: {m2.treatment.loc[150000]:.2f}, {m2.treatment.loc[175000]:.2f}, {m2.treatment.loc[199900]:.2f}")

heading("Gradual ceiling lift (Results; Figure S6)")
ramp = json.load(open(path("Data_S10_Ramp_Results.json")))["analysis"]
n = ramp["n_control"]
print(f"  Treatment change in B {ramp['treat_delta_B']['mean']:+.2f}; treatment vs control: t = {ramp['dB_t_statistic']:.2f}, p = {ramp['dB_p_value']:.2g}, "
      f"d = {ramp['dB_cohen_d'] * np.sqrt((n - 1) / n):.2f} (pooled sample SD; the deposited JSON value, {ramp['dB_cohen_d']:.2f}, uses the population SD)")

heading("Simulation 3, end of Phase 3 (Results; Table 4; Figure 5; text S1(e))")
s3 = pd.read_csv(path("Data_S7_Mixed_Summary.csv"))
c3, t3 = s3[s3.group == "control"], s3[s3.group == "treatment"]      # control = care-independent route; treatment = cooperative-care route
identical = all(np.array_equal(t3.sort_values("replicate")["phase3_" + v].values, t2.sort_values("replicate")["phase3_" + v].values) for v in "BLSK")
print(f"  Cooperative-care rows identical to the treatment rows of Simulation 2: {identical}")
print(f"  Change in B: care-independent {c3.delta_B.mean():+.2f} +/- {sem(c3.delta_B):.2f}, cooperative-care {t3.delta_B.mean():+.2f} +/- {sem(t3.delta_B):.2f}; p = {stats.ttest_ind(t3.delta_B, c3.delta_B).pvalue:.2f}")
print(f"  Change in L during Phase 3: care-independent {c3.delta_L.mean():+.2f} +/- {sem(c3.delta_L):.2f}, cooperative-care {t3.delta_L.mean():+.2f} +/- {sem(t3.delta_L):.2f}; p = {stats.ttest_ind(t3.delta_L, c3.delta_L).pvalue:.2f}")
for trait, name in (("B", "Brain size (B)"), ("L", "Childhood duration (L)"), ("S", "Social-learning propensity (S)"), ("K", "Knowledge (K)")):
    compare(name, t3["phase3_" + trait].values, c3["phase3_" + trait].values)
print(f"  K / 0.6B at the endpoint: care-independent {(c3.phase3_K / (0.6 * c3.phase3_B)).mean():.4f}, cooperative-care {(t3.phase3_K / (0.6 * t3.phase3_B)).mean():.4f}")
for label, x, care in (("cooperative-care (Simulations 2 and 3)", t3, 0.6), ("care-independent (Simulation 3)", c3, 0.0), ("control (Simulation 2)", c2, 0.0)):
    nf = net_fitness(x.phase3_B.values, x.phase3_L.values, x.phase3_K.values, care, Bmax=x.phase3_B_max.values)
    r, p = stats.pearsonr(x.delta_L, x.delta_B)
    print(f"  {label:<40} populations with net fitness at the replicate means at or above {BOUND_NF:.2f}: {int((nf >= BOUND_NF).sum())} of {len(x)}; "
          f"r(change in L, change in B) = {r:.2f}, p = {p:.2f}")

# ----------------------------------------------------------------------------------------------------------------------
# Details of the two analyses added in revision that the text reports in addition to Tables S5 and S6
# ----------------------------------------------------------------------------------------------------------------------
heading("Selection-strength analysis, replicate-level details (text S1(c) and S1(f); Figure S7)")
rows = [json.loads(line) for line in open(path("selection_strength_check", "out_main.jsonl"))]
def group(arm, form, grp):
    return [r for r in rows if r["arm"] == arm and r["form"] == form and r["group"] == grp]
def mean_L_at(rs, generation):
    return float(np.mean([dict((int(a[0]), a[2]) for a in r["traj"])[generation] for r in rs]))
ctrl = group("A", "superlinear", "control")
print(f"  Mean survival in the final generation: control {np.mean([r['final_surv_mean'] for r in ctrl]):.2f}; care " +
      ", ".join(f"{np.mean([r['final_surv_mean'] for r in group('A', f, 'treatment')]):.2f} ({f})" for f in ("superlinear", "linear", "saturating")))
print("  Control L at generations 10,000, 20,000, 30,000, and 99,000: " + ", ".join(f"{mean_L_at(ctrl, g):.2f}" for g in (10000, 20000, 30000, 99000)))
for form in ("superlinear", "linear", "saturating"):
    rs = group("A", form, "treatment")
    final = np.array([r["final_L"] for r in rs])
    print(f"  {form:<12} care L at generations 30,000, 50,000, 75,000, and 99,000: " + ", ".join(f"{mean_L_at(rs, g):.2f}" for g in (30000, 50000, 75000, 99000)) +
          f"; final L {final.mean():.2f}; replicates with final L above 19: {int((final > 19).sum())} of {len(final)}; mean S {np.mean([r['final_S'] for r in rs]):.2f}")
print(f"  Control S {np.mean([r['final_S'] for r in ctrl]):.2f}")
start = [r["phase1_end"]["L_mean"] for r in group("A", "superlinear", "treatment")]
print(f"  L of care populations at the last sample of the burn-in (generation 9,900): {min(start):.1f} to {max(start):.1f}")
print(f"  Duration above which the marginal superlinear return exceeds the temporal cost: L = {20 * (0.12 / (0.6 * 0.30 * 1.8)) ** (1 / 0.8):.2f}; "
      f"optimum under the saturating return: L = {np.sqrt(0.6 * 0.30 * 20 * 5 / 0.12) - 5:.2f}; marginal linear return {0.6 * 0.30:.2f}")

heading("Knowledge-limit analysis, replicate-level details (text S1(g); Table S6 note)")
kl = pd.read_csv(path("knowledge_limit_check", "replicate_results.csv"))
c04 = kl[(kl.cap == 0.4) & (kl.group == "control")].final_L.values
t04 = np.sort(kl[(kl.cap == 0.4) & (kl.group == "treatment") & (kl.form == "superlinear")].final_L.values)
tt = stats.ttest_ind(t04, c04)
print(f"  0.4B, superlinear: {int((t04 <= 5.0).sum())} of {len(t04)} replicates with final L of {t04[t04 <= 5.0].min():.1f} to {t04[t04 <= 5.0].max():.1f}; the others {', '.join(f'{v:.1f}' for v in t04[t04 > 5.0])}; "
      f"median {np.median(t04):.2f} (controls {np.median(c04):.2f}); t({len(t04) + len(c04) - 2}) = {tt.statistic:.2f}, p = {tt.pvalue:.3f}; "
      f"{int((t04 > c04.max()).sum())} of {len(t04)} above the largest control value")
pb = []
for cap, constant in ((0.4, -0.3), (0.8, 7.9)):
    c = kl[(kl.cap == cap) & (kl.group == "control")]
    for form in ("superlinear", "saturating"):
        t = kl[(kl.cap == cap) & (kl.group == "treatment") & (kl.form == form)]
        pb.append(stats.ttest_ind(t.final_B, c.final_B).pvalue)
sel = pd.read_csv(path("selection_strength_check", "replicate_results.csv"))
c06 = sel[(sel.arm == "A") & (sel.group == "control")]
for form in ("superlinear", "saturating"):
    pb.append(stats.ttest_ind(sel[(sel.arm == "A") & (sel.group == "treatment") & (sel.form == form)].final_B, c06.final_B).pvalue)
print(f"  Smallest two-sided p value for a difference in brain size across the 6 conditions: {min(pb):.3f}")
kj = [json.loads(line) for line in open(path("knowledge_limit_check", "out.jsonl"))]
rs = [r for r in kj if r["cap"] == 0.4 and r["form"] == "superlinear" and r["group"] == "treatment"]
print("  0.4B, superlinear: mean care L at generations 50,000, 75,000, 90,000, and 99,000: " + ", ".join(f"{mean_L_at(rs, g):.2f}" for g in (50000, 75000, 90000, 99000)))
