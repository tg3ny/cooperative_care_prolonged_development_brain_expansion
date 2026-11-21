# Knowledge-limit analysis (Supporting Information, text S1(g), Table S6)

Added in release v1.1.0 (November 2026), in response to peer review.

## Question

Knowledge in the model is capped at 0.6 x brain size (K <= 0.6B), an assumption taken from the Cultural Brain
Hypothesis model. A reviewer asked whether the results hold when this limit is tightened or relaxed. This analysis
repeats Simulation 1 with the limit set to 0.4B and to 0.8B.

## What was run

`selection_strength_check/sim_patched.py` with one further edit (`make_patched_cap.py`, `patch_cap.diff`): the
coefficient 0.6 of the storage rule becomes a parameter, `knowledge_cap` (default 0.6). With the default the code is
unchanged: 4 runs of the selection-strength analysis repeated with this code are identical to the stored results,
including their trajectories (`out_validation.jsonl`).

The plan (`PLAN.md`) was written before the runs. It fixes the design, the constants that hold baseline mortality at
about 10% per generation (-0.3 at 0.4B, 7.9 at 0.8B, from the model equations), the outcome measures, and the
predictions. It also discloses a pilot with 5 replicates per condition (saturating return only), kept in
`pilot_2026-11-16/` and not used in the analysis.

Runs (`out.jsonl`, 120 runs; N = 1,000; 100,000 generations; recorded seeds of replicates 1-20 of Simulation 1): for
each limit, 20 control runs, 20 care runs with the superlinear (default) return, and 20 care runs with the saturating
return. The limit applies throughout a run, including the burn-in. Rows for the default limit (0.6B) come from
`selection_strength_check/out_main.jsonl` (same seeds).

## Results

`table_S6.csv` (written by `analyze.py`) holds the values reported in Table S6. Mean ± SEM of 20 replicates at
generation 100,000; differences are care vs control at the same limit; TOST with ±2% bounds.

| Limit | Form | Control L | Care L | ΔL | Control B | Care B | ΔB | TOST p (B) |
|---|---|---|---|---|---|---|---|---|
| 0.4B | Superlinear | 3.58 ± 0.02 | 5.43 ± 0.96 | +51.5% | 38.82 ± 0.14 | 39.16 ± 0.20 | +0.86% | 0.039 |
| 0.4B | Saturating | 3.58 ± 0.02 | 7.28 ± 0.28 | +103.3% | 38.82 ± 0.14 | 38.58 ± 0.32 | -0.62% | 0.065 |
| 0.6B | Superlinear | 5.04 ± 0.02 | 19.60 ± 0.15 | +288.9% | 42.44 ± 0.07 | 42.27 ± 0.19 | -0.41% | 0.0010 |
| 0.6B | Saturating | 5.04 ± 0.02 | 7.27 ± 0.24 | +44.3% | 42.44 ± 0.07 | 42.06 ± 0.22 | -0.89% | 0.024 |
| 0.8B | Superlinear | 6.39 ± 0.02 | 19.80 ± 0.02 | +209.7% | 43.82 ± 0.05 | 43.62 ± 0.11 | -0.46% | 1.9e-06 |
| 0.8B | Saturating | 6.39 ± 0.02 | 8.23 ± 0.24 | +28.6% | 43.82 ± 0.05 | 43.50 ± 0.12 | -0.74% | 5.6e-05 |

Within the range examined (0.4B to 0.8B), brain size under care is within 1% of controls in all 6 conditions;
equivalence is established in 5 and not in the saturating form at 0.4B (p = 0.065). Under the saturating return childhood lengthens at every limit. Under the
superlinear return it approaches the upper limit of L at 0.6B and 0.8B, but at 0.4B, where the baseline duration is
well below the duration (about 5.8) above which this return exceeds the cost of childhood, it lengthens only modestly
in 18 of 20 replicates and rises toward the upper limit in 2; the mean difference is not significant (p = 0.061),
although final L exceeds the largest control value in 16 of 20 replicates. The extension at 0.4B is thus reduced
and inconsistent, not absent.

Relative to the plan: P1 (brain size) held on the criterion |delta B| < 2% in all conditions, with the one equivalence
test noted above falling short; P2 and P4 held; under P3 childhood settled at 7.3 at 0.4B and at 8.2 at 0.8B, above
the predicted 7.25 at the higher limit.

## Files

`PLAN.md` · `make_patched_cap.py`, `patch_cap.diff`, `sim_patched_cap.py` · `worker.py`, `tasks.jsonl`,
`tasks_validation.jsonl` · `out.jsonl`, `out_validation.jsonl` · `analyze.py` -> `table_S6.csv`,
`replicate_results.csv` · `pilot_2026-11-16/`.

## To rerun

```bash
python knowledge_limit_check/make_patched_cap.py selection_strength_check/sim_patched.py knowledge_limit_check/sim_patched_cap.py
python knowledge_limit_check/worker.py knowledge_limit_check/tasks.jsonl knowledge_limit_check/out.jsonl 0 1
python knowledge_limit_check/analyze.py
```

Run with Python 3.13, NumPy 2.5.3, SciPy 1.18.1, and pandas 3.0.5.
