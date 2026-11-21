# Reproducibility workflow

## Environment

The package was validated with Python 3.12. Create an isolated environment from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python verify_repository.py
```

On Windows, activate the environment with `.venv\Scripts\activate`.

## Verify deposited results

`verify_repository.py` performs the following checks without rerunning the evolutionary simulations:

1. validates SHA-256 hashes for every deposited code and data file;
2. verifies file dimensions, treatment groups, and replicate counts;
3. verifies master seed and replicate metadata in the JSON outputs;
4. recalculates representative means and independent-samples tests for Simulations 1–3; and
5. checks that the principal results match the manuscript within stated tolerances.

## Regenerate figures from deposited data

```bash
python code/generate_simulation1_figures.py --data-dir . --output-dir figures/generated
python code/generate_figures.py --data-dir . --output-dir figures/generated
python code/visualize_alternative_models.py \
  --input Data_S4_Alternative_Models.csv \
  --output figures/alternative_models
python code/generate_supplementary_figures.py --data-dir . --output-dir figures/generated
python code/make_figure1_model_diagram.py
python code/equilibrium_analysis.py
```

`generate_simulation1_figures.py` regenerates Figures 1 and 2 from `Data_S1_Summary.csv` and `Data_S3_Trajectories.csv.gz`. `visualize_alternative_models.py` regenerates Figure 1—figure supplement 4 (`Figure1_figure_supplement4.pdf`); panels A and B are drawn from the model's benefit functions and default parameters, and panel C is read from `Data_S4_Alternative_Models.csv`. `generate_figures.py` reads the instantaneous and gradual ceiling-lift trajectories from separate condition-specific files and concatenates each control/treatment pair automatically. The split keeps every file below GitHub’s browser-upload limit and ensures that Figure 3—figure supplement 1 is regenerated from the gradual-ramp condition. `generate_supplementary_figures.py` redraws Figure 1—figure supplements 1, 2, 3, and 5 (Figures S1, S2, S3, and S5 of the revised manuscript) from `Data_S3_Trajectories.csv.gz`, `Data_S1_Summary.csv`, and `Data_S2_Sensitivity.csv`; Figure S3 shows replicates 1 to 20 of each condition. `equilibrium_analysis.py` prints auxiliary convergence diagnostics that are not reported in the article.

## Recompute the reported statistics

```bash
python code/reported_statistics.py
```

The script reads only deposited files and prints the statistics reported in the article and its Supporting Information in the order in which they appear. Its conventions are those of the article: sample standard deviations, Cohen's d with the pooled sample SD, two-sided t tests, and equivalence tests with bounds of ±2% of the control mean. The `analysis` blocks of the deposited JSON files use the population SD for Cohen's d, so those values differ from the reported ones in the second decimal place. The two analyses added in v1.1.0 have their own folders, plans, and READMEs (`selection_strength_check/`, `knowledge_limit_check/`); their `analyze.py` scripts print Tables S5 and S6.

## Rerun Simulation 1

```bash
python code/childhood_first_simulation.py \
  --n-replicates 100 --n-jobs -1 --seed 42 \
  --output results/experiment1
```

The principal output mapping is:

- `results/experiment1/Data_S1_Summary.csv` → deposited `Data_S1_Summary.csv`
- `results/experiment1/Data_S3_Trajectories.csv` → deposited `Data_S3_Trajectories.csv.gz`

Run the sensitivity sweep separately:

```bash
python code/childhood_first_simulation.py \
  --sensitivity --n-jobs -1 --seed 42 \
  --output results/sensitivity
```

## Rerun Simulations 2 and 3

```bash
python code/grey_ceiling_simulation.py \
  --n-replicates 100 --n-jobs -1 --mode both --seed 42 \
  --output results/ceiling
```

The `pure` condition corresponds to Simulation 2. The `mixed` run supplies the care-independent energetic comparison used in Simulation 3. Outputs map to `Data_S6_Pure_Summary.csv`, `Data_S7_Mixed_Summary.csv`, `Data_S8_Pure_Results.json`, and `Data_S9_Mixed_Results.json` after applying the manuscript’s stable filenames.

Internal simulation-output names such as `experiment1` and `experiment2_trajectories.csv` are retained in the source code for backward compatibility. In the deposited repository, the two largest trajectory outputs are divided into control and treatment files without changing their contents.

For the gradual ceiling-lift robustness analysis:

```bash
python code/grey_ceiling_simulation.py \
  --n-replicates 100 --n-jobs -1 --mode pure --seed 42 \
  --ceiling-ramp 50000 --output results/gradual_ramp
```

## Rerun alternative models

Run this command from a dedicated results directory because the script writes its result files to the current directory:

```bash
mkdir -p results/alternative_models
cd results/alternative_models
python ../../code/alternative_models_analysis.py --replicates 100 --jobs -1
cd ../..
```

The deposited `Data_S4_Alternative_Models.csv` used 100 replicates per condition per specification. Seeds in the current script are derived deterministically from each model name with `hashlib`; the deposited file was produced by an earlier version that derived seeds from Python's per-process string hash, so a rerun reproduces its results statistically but not bit for bit. The alternative-model analysis is independent of the Simulation 1 run (different seeds), which is why its superlinear row differs slightly from `Data_S1_Summary.csv`.

## Determinism and computing time

Simulations 1 to 3, the sensitivity sweep, and the gradual-ramp analysis use master seed 42 (the alternative-model analysis uses model-name-derived seeds; see above). Replicate-specific seeds are recorded in the summary files of Simulations 1 to 3. In the sensitivity sweep, the 20 seeds of a cell are drawn from a generator seeded with the master seed plus the index of the cell's first run, and the control and treatment runs of a replicate share a seed, so the two runs are identical until care is introduced; cell-level results are deposited (`Data_S2_Sensitivity.csv`), and replicate-level results of the sweep and of the alternative-model analysis are not. The quick heatmap that `childhood_first_simulation.py --sensitivity` writes after a sweep needs seaborn, which is not in `requirements.txt`; that step is skipped if seaborn is absent. Parallel execution can change completion order but not the deterministic result assigned to a recorded replicate seed. Full regeneration is computationally intensive; checking the deposited outputs and regenerating figures is substantially faster.

## Reruns checked for release v1.1.1

Single replicates were rerun with Python 3.12 and the versions pinned in `requirements.txt`, using the deposited code and the recorded seeds: replicates 1 and 57 of Simulation 1, replicate 1 of Simulations 2 and 3 and of the gradual-ramp analysis (control and treatment in each case), 2 runs of each analysis added in v1.1.0, and all 40 runs of one cell of the sensitivity sweep (r = 0.30, e = 1.8). Every rerun reproduced the deposited endpoints, and the deposited trajectories where these exist, to floating-point rounding (relative differences below 1e-12).

## Releases

- `v1.0.0`: snapshot associated with the original submission (Zenodo, https://doi.org/10.5281/zenodo.22694519).
- `v1.1.0`: adds the selection-strength and knowledge-limit analyses and the model-diagram script (Zenodo, https://doi.org/10.5281/zenodo.24408976).
- `v1.1.1`: adds `code/reported_statistics.py` and `code/generate_supplementary_figures.py`, redirects two figure scripts to `figures/generated/`, corrects the printed summary of `code/equilibrium_analysis.py`, and updates the documentation. Simulation code, data, and results are unchanged.

All versions are archived at Zenodo under https://doi.org/10.5281/zenodo.22694518.
