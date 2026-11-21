# Selection-strength analysis (Supporting Information, text S1(f), Figure S7, Table S5)

Added in release v1.1.0 (November 2026), in response to peer review.

## Question

Survival in the model is a logistic function of net fitness with an upper bound of 0.999. Under the default
parameters the evolved populations lie close to that bound (mean survival about 0.996 in control and 0.998 in
care populations at the end of Simulation 1), so selection is weak. In the alternative-model analysis
(`Data_S4_Alternative_Models.csv`), brain size (B) was lower under care when the social return was linear or
saturating. Those forms deliver larger benefits and bring survival to or near the bound. This analysis asks
whether B is still lower under linear and saturating returns when survival is kept below the bound.

## What was run

The deposited simulation (`code/childhood_first_simulation.py`, release v1.0.0) with the edits in `patch.diff`:
a constant (`baseline_offset`) subtracted from net fitness from the start of Phase 2, in control and treatment
alike, and a switch for the linear and saturating benefit forms (expressions copied from
`code/alternative_models_analysis.py`). With the default switches the patched code reproduces the published
replicates (identical up to rounding in the last digit); `analyze.py` checks this for replicates 1-20 of both groups.

The plan (`PLAN.md`) was written before the runs. It fixes the constant (3.5, a baseline mortality of about 10%
per generation), the design, the outcome measures, and the predictions, and it records one amendment: in a
first attempt the constant was applied from generation 0, and some populations collapsed during the burn-in;
those runs are kept in `attempt1_offset_from_generation0.jsonl` and are not used.

Runs (N = 1,000; 100,000 generations; recorded seeds of replicates 1-20 of Simulation 1):

| File | Runs |
|---|---|
| `out_main.jsonl` | Constant 3.5: 20 control, and 20 treatment for each of the superlinear, linear, and saturating forms. No constant: 20 treatment for the linear and for the saturating form. (120 runs) |
| `out_default.jsonl` | No constant, superlinear form: 20 control and 20 treatment. These reproduce the published replicates and supply agent-level survival for the default rows of Table S5. (40 runs, added after the main runs; they change no comparison) |
| `out_validation.jsonl` | 4 runs used to validate the patch before the main runs |

Each line of an output file is one run: final population means, mean survival and the share of agents at the
survival bound in the final generation, and the trajectory sampled every 1,000 generations.

## Results

`table_S5.csv` (written by `analyze.py`) holds the values reported in Table S5. Mean of 20 replicates at
generation 100,000; differences are care vs control within a selection regime; TOST with ±2% bounds.

| Selection | Form | Care L | ΔL | Care B | ΔB | TOST p (B) | Survival under care | Agents at bound |
|---|---|---|---|---|---|---|---|---|
| Default (offset 0) | Superlinear | 8.54 ± 0.41 | +32.0% | 41.05 ± 0.47 | -1.25% | 0.31 | 0.998 | 0% |
| Default (offset 0) | Linear | 7.89 ± 0.44 | +21.9% | 40.89 ± 0.52 | -1.62% | 0.41 | 0.999 | 15% |
| Default (offset 0) | Saturating | 8.39 ± 0.46 | +29.7% | 39.10 ± 1.07 | -5.94% | 0.92 | 0.999 | 79% |
| Baseline mortality ~10% (offset 3.5) | Superlinear | 19.60 ± 0.15 | +288.9% | 42.27 ± 0.19 | -0.41% | 0.0010 | 0.983 | 0% |
| Baseline mortality ~10% (offset 3.5) | Linear | 18.39 ± 0.29 | +264.9% | 42.24 ± 0.19 | -0.47% | 0.0012 | 0.982 | 0% |
| Baseline mortality ~10% (offset 3.5) | Saturating | 7.27 ± 0.24 | +44.3% | 42.06 ± 0.22 | -0.89% | 0.024 | 0.984 | 0% |

Control populations: default, L = 6.47 ± 0.20, B = 41.57 ± 0.41, survival 0.996; baseline mortality ~10%, L = 5.04 ± 0.02,
B = 42.44 ± 0.07, survival 0.911.

With survival below the bound, brain size is statistically equivalent in care and control populations under
all three forms. The form of the return governs how far childhood extends: to the model's upper limit
(L = 20) under the superlinear form, toward it and still rising under the linear form, and to a stable level
of about 7.3 under the saturating form.

## Files

`PLAN.md` plan, predictions, amendment · `patch.diff`, `make_patched.py`, `sim_patched.py` code ·
`worker.py`, `tasks_*.jsonl` run lists · `out_*.jsonl` results · `analyze.py` → `table_S5.csv`,
`replicate_results.csv` · `make_figure_s7.py` → `Figure_S7.png` · `attempt1_*` the discarded first attempt.

## To rerun

```bash
python selection_strength_check/make_patched.py code/childhood_first_simulation.py selection_strength_check/sim_patched.py
python selection_strength_check/worker.py selection_strength_check/tasks_main.jsonl selection_strength_check/out_main.jsonl 0 1
python selection_strength_check/worker.py selection_strength_check/tasks_default.jsonl selection_strength_check/out_default.jsonl 0 1
python selection_strength_check/analyze.py
python selection_strength_check/make_figure_s7.py
```

`worker.py <tasks> <output> <worker index> <number of workers>` skips finished tasks, so several workers can
share one task list. One run takes about 50 seconds. Run with Python 3.13, NumPy 2.5.3, SciPy 1.18.1,
pandas 3.0.5, and Matplotlib 3.11.2.
