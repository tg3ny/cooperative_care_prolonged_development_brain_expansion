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
python code/generate_figures.py --data-dir . --output-dir figures/generated
python code/visualize_alternative_models.py \
  --input Data_S4_Alternative_Models.csv \
  --output figures/alternative_models
python code/equilibrium_analysis.py
```

`generate_figures.py` reads the instantaneous and gradual ceiling-lift trajectories from separate condition-specific files and concatenates each control/treatment pair automatically. The split keeps every file below GitHub’s browser-upload limit and ensures that Figure 3—figure supplement 1 is regenerated from the gradual-ramp condition.

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
python ../../code/alternative_models_analysis.py --replicates 20 --jobs -1
cd ../..
```

## Determinism and computing time

All headline analyses use master seed 42. Replicate-specific seeds are recorded in the summary files. Parallel execution can change completion order but not the deterministic result assigned to a recorded replicate seed. Full regeneration is computationally intensive; checking the deposited outputs and regenerating figures is substantially faster.

## Versioned public release

For the public eLife release:

1. synchronize these corrected root-layout files with https://github.com/tg3ny/cooperative_care_prolonged_development_brain_expansion;
2. commit the files and tag the exact submitted version as `v1.0.0`;
3. create a GitHub release from that tag;
4. connect the repository to Zenodo and archive the release;
5. add the GitHub release URL, full commit identifier, and Zenodo DOI to the article and cover letter; and
6. retain the immutable repository ZIP supplied in the eLife package as a reviewer-accessible backup.
