# Data dictionary

All deposited CSV, compressed CSV, and JSON data files are stored in the repository root to match the public GitHub layout. CSV files use a header row, comma delimiters, decimal points, and UTF-8 encoding. Compressed files use standard gzip compression and can be opened directly with `pandas.read_csv`. Simulation units are dimensionless model units unless stated otherwise.

## Primary files

| File | Rows | Purpose |
|---|---:|---|
| `Data_S1_Summary.csv` | 200 | Simulation 1 replicate-level endpoints, 100 control and 100 treatment replicates. |
| `Data_S2_Sensitivity.csv` | 81 | Parameter sweep across social benefit rate `r` and convexity exponent `e`. |
| `Data_S3_Trajectories.csv.gz` | 200,000 | Simulation 1 generation-level trajectories. |
| `Data_S4_Alternative_Models.csv` | 5 | Alternative social-benefit model results. |
| `Data_S5_Equilibrium_Statistics.csv` | 8 | Phase-boundary equilibrium diagnostics. |
| `Data_S6_Pure_Summary.csv` | 200 | Simulation 2 replicate endpoints for the care-gated ceiling condition. |
| `Data_S7_Mixed_Summary.csv` | 200 | Simulation 3 replicate endpoints for cooperative-care and care-independent energetic routes. |
| `Data_S8_Pure_Results.json` | — | Simulation 2 parameters and analysis results. |
| `Data_S9_Mixed_Results.json` | — | Simulation 3 parameters and analysis results. |
| `Data_S10_Ramp_Results.json` | — | Gradual ceiling-lift parameters and analysis results. |
| `experiment2_control_trajectories.csv.gz` | 200,000 | Control trajectories for the instantaneous ceiling-lift analysis used in Figure 3. |
| `experiment2_treatment_trajectories.csv.gz` | 200,000 | Treatment trajectories for the instantaneous ceiling-lift analysis used in Figure 3. |
| `gradual_ceiling_control_trajectories.csv.gz` | 200,000 | Control trajectories for the gradual ceiling-lift analysis used in Figure 3—figure supplement 1. |
| `gradual_ceiling_treatment_trajectories.csv.gz` | 200,000 | Treatment trajectories for the gradual ceiling-lift analysis used in Figure 3—figure supplement 1. |

The two 400,000-row trajectory datasets are divided by condition solely to keep every repository file below GitHub’s browser-upload limit. `code/generate_figures.py` concatenates each control/treatment pair automatically; no observations were removed or altered.

## Common identifiers

- `replicate`: integer replicate identifier within a group, 1–100.
- `group`: simulation condition, `control` or `treatment`.
- `seed`: replicate-specific pseudorandom seed derived from master seed 42.
- `generation`: simulated generation.
- `phase`: sequential model phase.

## Trait variables

- `B`, `B_mean`, `B_std`: arbitrary heritable brain-size scalar and its within-population mean and standard deviation; it is not calibrated to brain volume or adjusted for body size.
- `L`, `L_mean`, `L_std`: developmental-duration scalar and its within-population mean and standard deviation; it is not years, growth rate, weaning age, age at maturity, or age at first reproduction.
- `S`, `S_mean`, `S_std`: dimensionless social-learning propensity and its within-population mean and standard deviation; it is not a measured social behavior.
- `K`, `K_mean`, `K_std`: payoff-bearing information derived within each generation; it is not an empirical measure of knowledge.
- `B_max`, `phase3_B_max`: imposed upper energetic constraint on the brain-size scalar in arbitrary units.
- `K_ceiling_frac`: knowledge relative to the brain-dependent storage ceiling.

Prefixes such as `phase2_` and `phase3_` identify the model phase at which an endpoint was measured. Internal source-code names containing `epoch`, `experiment`, or `nutrition_only` are retained for backward compatibility with earlier outputs; the manuscript uses “phase,” “simulation,” and “care-independent energetic route.”

## Derived outcome variables

- `delta_B`, `delta_L`: within-replicate change in brain size or childhood duration across the specified phases.
- `delta_B_percent`, `delta_L_percent`: relative between-condition difference in percent.
- `decoupled` or `decoupling`: whether the prespecified developmental-decoupling criterion was met: childhood extension of at least 10% and absolute brain-size difference below 2%.
- `t_stat_L`, `p_value_L`, `cohens_d`: independent-samples childhood comparison and standardized effect size.
- `p_tost_B`: equivalence-test p-value for brain size under ±2% bounds.

The JSON files contain the complete parameter dictionaries, master seed, replicate count, runtime metadata, and analysis outputs for their respective conditions.
