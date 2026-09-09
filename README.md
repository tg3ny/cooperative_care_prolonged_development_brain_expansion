# Cooperative care and childhood first simulation

This repository contains the complete Python code and data supporting the manuscript “Cooperative care as a candidate mechanism linking prolonged childhood and brain expansion in human evolution” by Tobias Grossmann. Three individual-based evolutionary simulations test whether a cooperative-care context can extend a developmental-duration trait before brain expansion, accompany later brain expansion when an energetic ceiling is relaxed, and generate a route distinct from care-independent energetic relaxation.

The repository is self-contained for verification and figure regeneration. Trajectory tables are stored as gzip-compressed CSV files. The two largest datasets are divided into control and treatment files so every individual file is below GitHub’s browser-upload limit; the figure-generation script concatenates the pairs automatically.

## Repository contents

```text
code/                         Simulation, analysis, and figure scripts
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
python code/generate_figures.py --data-dir . --output-dir figures/generated
python code/visualize_alternative_models.py \
  --input Data_S4_Alternative_Models.csv \
  --output figures/alternative_models
python code/equilibrium_analysis.py
```

The first command regenerates Figure 3, Figure 4, and Figure 3—figure supplement 1 from the deposited trajectory and summary files. The Simulation 1 and alternative-model scripts generate the other figures as described in `docs/REPRODUCIBILITY.md`.

## Full simulation regeneration

The simulations use 100 independent stochastic replicates per group and master seed 42. Full regeneration is computationally intensive:

```bash
python code/childhood_first_simulation.py -n 100 -j -1 -o results/experiment1 --seed 42
python code/grey_ceiling_simulation.py -n 100 -j -1 --mode both -o results/ceiling --seed 42
python code/grey_ceiling_simulation.py -n 100 -j -1 --mode pure \
  --ceiling-ramp 50000 -o results/gradual_ramp --seed 42
```

See `docs/REPRODUCIBILITY.md` for output mappings and the alternative-model workflow.

## Key reported results

| Simulation | Result |
|---|---|
| Simulation 1 | Introducing the cooperative-care context increased childhood duration by 32.5% while brain size remained statistically equivalent under the prespecified ±2% bounds. |
| Simulation 2 | With an imposed care-gated ceiling lift, mean treatment brain size increased by approximately 34% relative to its Simulation 1 endpoint. |
| Simulation 3 | Cooperative-care and care-independent energetic conditions reached similar brain sizes but differed in childhood duration. |

## Model scope

The simulated populations are synthetic. `B`, `L`, `S`, and `K` are dimensionless model constructs rather than direct measurements of brain volume, life-history milestones, social behavior, or empirical knowledge. `Sub_care` is a dimensionless social-context parameter, not a count of helpers, care time, or calories. `B_max` is an imposed energetic constraint in arbitrary units. The model does not include body size, adiposity, explicit helper strategies, or allocation of energy among brain tissue, somatic growth, and reproduction.

The manipulations identify causal consequences within the specified fitness functions; they do not reconstruct historical causation in a fossil population. Stable internal filenames and selected output keys retain the word `experiment` for backward compatibility with the deposited data.

## Data availability and preservation

All data required to reproduce the manuscript’s statistics and figures are included here. The public repository is https://github.com/tg3ny/cooperative_care_prolonged_development_brain_expansion. The exact eLife submission version should be tagged as `v1.0.0` and archived with Zenodo; record the release URL, full commit identifier, and Zenodo DOI in the manuscript’s Data Availability Statement.

## Contact

Tobias Grossmann  
Department of Psychology, University of Virginia  
tg3ny@virginia.edu
