# Cooperative care and childhood first simulation

This repository contains the complete Python code and data supporting the manuscript “Cooperative care as a candidate mechanism linking prolonged childhood and brain expansion in human evolution” by Tobias Grossmann. Three agent-based evolutionary simulations test whether a cooperative-care context can extend a developmental-duration trait before brain expansion, accompany later brain expansion when an energetic ceiling is relaxed, and generate a route distinct from care-independent energetic relaxation.

The repository is self-contained for verification and figure regeneration. Trajectory tables are stored as gzip-compressed CSV files. The two largest datasets are divided into control and treatment files so every individual file is below GitHub’s browser-upload limit; the figure-generation script concatenates the pairs automatically.

## Repository contents

```text
code/                         Simulation, analysis, and figure scripts
selection_strength_check/     Selection-strength analysis added in v1.1.0 (plan, code, results, README)
knowledge_limit_check/        Knowledge-limit analysis added in v1.1.0 (plan, code, results, README)
Data_S*.csv, *.json, *.csv.gz Primary results and complete trajectories at repository root
docs/DATA_DICTIONARY.md       File and variable definitions
docs/REPRODUCIBILITY.md       Full regeneration workflow
checksums/SHA256SUMS          Integrity hashes for code and data
verify_repository.py          Automated integrity and result checks
requirements.txt              Tested Python dependencies
CITATION.cff                  Software and dataset citation metadata
LICENSE.md                    Code and data licenses
```

## Quick verification

Python 3.12 is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_repository.py
```

The verification script checks the manifest, data dimensions, replicate counts, random-seed metadata, and the principal numerical results reported in the manuscript.

## Regenerate figures from the deposited data

From the repository root:

```bash
python code/generate_simulation1_figures.py --data-dir . --output-dir figures/generated
python code/generate_figures.py --data-dir . --output-dir figures/generated
python code/visualize_alternative_models.py \
  --input Data_S4_Alternative_Models.csv \
  --output figures/alternative_models
python code/equilibrium_analysis.py
```

The first command regenerates Figures 1 and 2; the second regenerates Figure 3, Figure 4, and Figure 3—figure supplement 1; the third regenerates Figure 1—figure supplement 4. Figure 1—figure supplements 1, 2, 3, and 5 are produced by `childhood_first_simulation.py` as described in `docs/REPRODUCIBILITY.md`.

## Full simulation regeneration

The simulations use 100 independent stochastic replicates per group and master seed 42 (the alternative-model analysis derives its seeds from each model name; see `docs/REPRODUCIBILITY.md`). Full regeneration is computationally intensive:

```bash
python code/childhood_first_simulation.py -n 100 -j -1 -o results/experiment1 --seed 42
python code/grey_ceiling_simulation.py -n 100 -j -1 --mode both -o results/ceiling --seed 42
python code/grey_ceiling_simulation.py -n 100 -j -1 --mode pure \
  --ceiling-ramp 50000 -o results/gradual_ramp --seed 42
```

See `docs/REPRODUCIBILITY.md` for output mappings and the alternative-model workflow.

## Analyses added in v1.1.0

Under the default parameters, survival lies close to its upper bound of 0.999, so selection is weak. The folder `selection_strength_check/` contains an analysis, added during peer review, that repeats Simulation 1 for the superlinear, linear, and saturating forms of the social return with a baseline mortality of about 10% per generation (20 replicates per condition, run with the recorded seeds of replicates 1-20). With survival below the bound, brain size is statistically equivalent in care and control populations under all three forms, and the form of the return governs how far childhood extends. The folder holds the plan written before the runs, the patch to the deposited simulation, all results, and its own README. 

The folder `knowledge_limit_check/` contains a second analysis that repeats Simulation 1 with the limit on knowledge (0.6 x brain size by default) set to 0.4 and 0.8 x brain size. Brain size under care stays within 1% of controls at every limit examined (0.4, 0.6, and 0.8 x brain size); under the default, accelerating return, a large and consistent lengthening of childhood occurs only when the baseline duration, which the limit sets, is not too short. To recompute the summaries and the figure from the deposited results:

```bash
python selection_strength_check/analyze.py
python selection_strength_check/make_figure_s7.py
python knowledge_limit_check/analyze.py
```

## Figure and table numbering in the article

Script comments and docstrings keep the figure names of an earlier manuscript version. They correspond to the revised manuscript as follows.

| Name in this repository | Revised manuscript |
|---|---|
| (model diagram, `code/make_figure1_model_diagram.py`) | Figure 1 |
| Figure 1 | Figure 2 |
| Figure 2 | Figure 3 |
| Figure 3 | Figure 4 |
| Figure 4 | Figure 5 |
| Figure 1-figure supplements 1 to 5 | Figures S1 to S5 |
| Figure 3-figure supplement 1 | Figure S6 |
| (`selection_strength_check/make_figure_s7.py`) | Figure S7 |
| (`selection_strength_check/table_S5.csv`) | Table S5 |
| (`knowledge_limit_check/table_S6.csv`) | Table S6 |

## Key reported results

| Simulation | Result |
|---|---|
| Simulation 1 | Introducing the cooperative-care context increased childhood duration by 32.5% while brain size remained statistically equivalent under the prespecified ±2% bounds. |
| Simulation 2 | With an imposed care-gated ceiling lift, mean treatment brain size increased by approximately 34% relative to its Simulation 1 endpoint. |
| Simulation 3 | Cooperative-care and care-independent energetic conditions reached similar brain sizes but differed in childhood duration. The cooperative-care trajectories are numerically identical to those of Simulation 2 because the raised ceiling never binds. |

## Origin of the code

The simulation code was adapted in Python from the C++ code published with the Cultural Brain Hypothesis model (Muthukrishna M, Doebeli M, Chudek M, Henrich J. The Cultural Brain Hypothesis: how culture drives brain expansion, sociality, and life history. PLoS Comput Biol. 2018;14:e1006504; S1 Code; https://doi.org/10.1371/journal.pcbi.1006504), which is distributed under the Creative Commons Attribution License. Claude Code (Anthropic) was used to make this adaptation and to write the simulation and analysis code to the author's specifications; the author verified the code and its output.

## Model scope

The simulated populations are synthetic. `B`, `L`, `S`, and `K` are dimensionless model constructs rather than direct measurements of brain volume, life-history milestones, social behavior, or empirical knowledge. `Sub_care` is a dimensionless social-context parameter, not a count of helpers, care time, or calories. `B_max` is an imposed energetic constraint in arbitrary units. The model does not include body size, adiposity, explicit helper strategies, or allocation of energy among brain tissue, somatic growth, and reproduction.

The manipulations identify causal consequences within the specified fitness functions; they do not reconstruct historical causation in a fossil population. Stable internal filenames and selected output keys retain the word `experiment` for backward compatibility with the deposited data.

## Data availability and preservation

All data required to reproduce the manuscript’s statistics and figures are included here. The public repository is https://github.com/tg3ny/cooperative_care_prolonged_development_brain_expansion. Version `v1.1.0` identifies the repository snapshot associated with the revised manuscript; it adds the two analyses made in revision and the model-diagram script to `v1.0.0`, the snapshot associated with the original submission (Zenodo, https://doi.org/10.5281/zenodo.22694519). No file of `v1.0.0` other than `README.md`, `CITATION.cff`, and `checksums/SHA256SUMS` was changed. All versions are archived at Zenodo under https://doi.org/10.5281/zenodo.22694518.

## Contact

Tobias Grossmann  
Department of Psychology, University of Virginia  
tg3ny@virginia.edu
