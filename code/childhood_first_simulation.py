"""
Canonical script for Simulation 1 (developmental decoupling).

Produces the primary simulation outputs for Table 1, Figures 1-2, and related
SI summaries (Tables S2-S3, Figures S1-S3, S5).

Outputs:
  Data_S1_Summary.csv     — Simulation 1 endpoint values (canonical for Table 1)
  Data_S2_Sensitivity.csv — Sensitivity analysis (canonical for Appendix 1-table 3 and Figure 1-figure supplement 5)
  Data_S3_Trajectories.csv — Full trait trajectories (canonical for Figs. 1, 2, S1, S3)

Terminology: Phase 1 (baseline, gen 0-10k), Phase 2 (care divergence, gen 10k-100k).
Master seed: 42. All results are fully deterministic.

Note on internal naming: variable names use 'epoch' (e.g., T_epoch1, epoch2_B)
for backward compatibility with data files. The manuscript uses Phase 1 / Phase 2.
"""
#!/usr/bin/env python3
"""
childhood_first_simulation.py

Evolutionary Simulation Testing the Childhood-First Hypothesis
===============================================================

This simulation tests whether cooperative care can create selection pressure
for extended childhood independent of brain size expansion, as proposed by
the Childhood-First hypothesis (Grossmann, 2026).

THEORETICAL BASIS
-----------------
Some accounts treat extended development primarily as a metabolic consequence
of encephalization. The fossil record instead contains proxy-specific mosaics
that motivate testing whether developmental duration and brain size can vary
partly independently:

  - Dikika (>3 Ma): Australopithecus afarensis combines chimpanzee-like
    dental development with evidence for prolonged brain growth.
  - Dmanisi (~1.77 Ma): limited dental evidence has been interpreted as
    indicating extended growth despite relatively small endocranial volumes.
  - Rising Star (~300 ka): Homo naledi shows some human-like enamel features
    despite a small endocranial volume; inferences about its overall life
    history remain uncertain.

The Childhood-First hypothesis proposes a candidate sequence in which ecological
variation makes cooperative provisioning valuable. When multiple caregivers
provision a dependent juvenile, longer dependence may become selectively
advantageous without requiring simultaneous brain expansion. The model tests
this conditional mechanism; it does not establish the historical origin or
maintenance of cooperation in a particular fossil population.

MODEL MECHANISM
---------------
The model implements a "social fitness" mechanism representing survival
buffering against ecological instability:

    Social_fitness = Care × rate × (L/L_max)^exponent × L_max

The fitness benefit from cooperative care is specified as convex (exponent > 1)
to represent increasing returns to longer care-supported development. This
non-linearity is a model assumption whose consequences are evaluated through
sensitivity analyses:

  - LINEAR returns (exponent = 1.0): Childhood extends but brain SIZE SHRINKS
    as a compensatory trade-off—contradicting the fossil record
  - CONVEX returns (exponent > 1.4): Childhood extends with brain size STASIS
    —matching the Dmanisi and Rising Star pattern

Within the tested parameterizations, sufficient convexity is required to
produce the modeled combination of longer developmental duration and brain-size
stasis. This result is not a direct reconstruction of a fossil taxon's growth.

Without cooperative care (Care = 0), the model-defined social-fitness term is
zero regardless of childhood length. With care, longer childhood yields a
larger value of that imposed fitness component.

SIMULATION DESIGN
-------------------
The simulation runs in two phases across 100,000 generations:

  Phase 1 (generations 0-10,000): Baseline
    - Both groups evolve without cooperative care
    - Brain size ceiling at B_max = 50 (metabolic constraint)
    - Establishes common evolutionary starting point

  Phase 2 (generations 10,000-100,000): Treatment divergence
    - Treatment group receives cooperative care (Care = 0.6)
    - Control group remains without care (Care = 0.0)
    - Brain ceiling maintained—tests whether childhood can extend
      without brain expansion (Decoupling Hypothesis)

HYPOTHESIS TESTED
-----------------
Decoupling Hypothesis: Cooperative care enables childhood extension without
    requiring brain expansion. Operationalized as:
    - Treatment and control L differ significantly (p < 0.005, two-sided)
    - Treatment B equivalent to Control B (|ΔB| < 2%)

USAGE
-----
    # Quick test (10 replicates per group)
    python childhood_first_simulation.py -n 10 -o ./test_results

    # Full simulation (100 replicates per group)
    python childhood_first_simulation.py -n 100 -o ./full_results

    # SENSITIVITY ANALYSIS (Figure 1-figure supplement 5)
    python childhood_first_simulation.py --sensitivity -o ./results -j -1

OUTPUT
------
    experiment_summary.csv  : Per-replicate summary statistics
    full_trajectories.csv   : Generation-by-generation evolution data
    experiment_results.json : Complete analysis with hypothesis tests
    sensitivity_analysis.csv: Grid search results (if --sensitivity used)

REQUIREMENTS
------------
    numpy, pandas, scipy, joblib (for parallel execution)

    Install with: pip install numpy pandas scipy joblib

CITATION
--------
Author: Tobias Grossmann (tg3ny@virginia.edu)
"""

import os
import sys
import json
import time
import argparse
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from dataclasses import dataclass, asdict, replace
from typing import List, Dict, Tuple, Optional

try:
    from joblib import Parallel, delayed
    PARALLEL_AVAILABLE = True
except ImportError:
    PARALLEL_AVAILABLE = False
    print("Warning: joblib not available. Running in serial mode.")


# =============================================================================
# MODEL PARAMETERS
# =============================================================================

@dataclass
class SimulationParameters:
    """
    Parameters for the Childhood-First evolutionary simulation.
    
    Parameters are expressed in arbitrary model units and selected to produce
    stable dynamics while allowing clear tests of the stated predictions.
    """
    
    # --- Population ---
    N: int = 1000                    # Population size
    T_total: int = 100_000           # Total generations
    T_epoch1: int = 10_000           # End of baseline phase
    T_epoch2: int = 100_000          # End of treatment phase (runs to completion)
    
    # --- Metabolic Costs ---
    alpha: float = 0.1               # Linear brain metabolic cost
    gamma: float = 0.001             # Quadratic brain metabolic cost
    beta: float = 0.12               # Temporal cost of childhood per unit L
    k_ceiling: float = 0.5           # Steepness of brain size ceiling
    
    # --- Knowledge/Learning ---
    lambda_: float = 0.5             # Benefit coefficient for knowledge
    learning_exponent: float = 1.15  # Superlinearity of childhood learning
    asocial_rate: float = 0.10       # Base rate of asocial learning
    
    # --- Social Fitness (Key Mechanism) ---
    # These parameters control the social benefit from cooperative care.
    # The social benefit creates selection for longer childhood when care
    # is present, independent of brain size.
    social_rate: float = 0.30        # Strength of social fitness benefit
    social_exponent: float = 1.8     # How superlinearly L affects social benefit
    
    # --- Cooperative Care ---
    Sub_care_control: float = 0.0    # Care level for control group
    Sub_care_treatment: float = 0.6  # Care level for treatment group
    
    # --- Evolution ---
    mu: float = 0.01                 # Mutation probability per individual
    sigma_mut: float = 0.1           # Mutation effect size (scaled per trait)
    
    # --- Trait Bounds ---
    B_max: float = 50.0              # Brain ceiling (metabolic constraint)
    L_max: float = 20.0              # Maximum childhood duration
    
    # --- Statistical Thresholds ---
    rel_dB_threshold: float = 0.02   # Equivalence threshold for brain size (±2%)
    rel_dL_threshold: float = 0.10   # Minimum childhood extension (≥10%)
    alpha_level: float = 0.005       # Significance level for hypothesis tests


# =============================================================================
# POPULATION CLASS
# =============================================================================

class Population:
    """
    Represents an evolving population with heritable traits for brain size (B),
    childhood duration (L), and social learning propensity (S).
    
    The fitness function implements the social fitness mechanism that creates
    selection for extended childhood under cooperative care conditions.
    """
    
    def __init__(self, params: SimulationParameters, seed: int):
        """
        Initialize population with random trait values.
        
        Args:
            params: Simulation parameters
            seed: Random seed for reproducibility
        """
        self.params = params
        self.rng = np.random.default_rng(seed)
        self.N = params.N
        
        # Environmental state (set by simulation)
        self.B_max = params.B_max
        self.Sub_care = 0.0
        
        # Initialize traits with variation
        self.B = self.rng.uniform(5, 15, self.N)      # Brain size
        self.L = self.rng.uniform(1, 3, self.N)       # Childhood duration
        self.S = self.rng.uniform(0.3, 0.7, self.N)   # Social learning propensity
        self.K = np.zeros(self.N)                      # Knowledge (derived)
        
        # For social learning across generations
        self.prev_K = None
        self.prev_fitness = None
    
    def step(self) -> Dict[str, float]:
        """
        Execute one generation: learning, selection, reproduction.
        
        Returns:
            Dictionary of population statistics
        """
        p = self.params
        
        # =====================================================================
        # PHASE 1: LEARNING
        # =====================================================================
        
        # Asocial learning: individual exploration, scales with brain and childhood
        # Superlinear with L to capture developmental advantage of longer childhood
        L_factor = np.power(self.L, p.learning_exponent)
        K_asocial = self.B * L_factor * p.asocial_rate
        
        # Social learning: acquiring knowledge from previous generation
        K_social = np.zeros(self.N)
        if self.prev_K is not None:
            # Fitness-biased transmission: learn from successful individuals
            probs = self.prev_fitness / np.sum(self.prev_fitness)
            models = self.rng.choice(len(self.prev_K), size=(self.N, 3), p=probs)
            K_social = np.mean(self.prev_K[models], axis=1)
        
        # Combine learning types, weighted by social learning propensity
        K_from_social = self.S * K_social * (self.L / p.L_max)
        K_from_asocial = (1 - self.S) * K_asocial
        K_total = K_from_social + K_from_asocial
        
        # Brain capacity constraint on knowledge
        self.K = np.minimum(K_total, self.B * 0.6)
        
        # =====================================================================
        # PHASE 2: FITNESS CALCULATION
        # =====================================================================
        
        # --- Costs ---
        
        # Structural/metabolic cost of brain (linear + quadratic)
        C_struct = p.alpha * self.B + p.gamma * (self.B ** 2)
        
        # Soft ceiling on brain size (sigmoid penalty)
        ceiling_penalty = 1.0 / (1.0 + np.exp(-p.k_ceiling * (self.B - self.B_max)))
        C_struct += 10.0 * ceiling_penalty
        
        # Temporal cost of childhood (fixed, not reduced by care in this model)
        C_time = p.beta * self.L
        
        total_cost = C_struct + C_time
        
        # --- Benefits ---
        
        # Knowledge benefit (survival/foraging advantage)
        benefit_knowledge = p.lambda_ * self.K
        
        # Social fitness benefit (KEY MECHANISM)
        #
        # This represents fitness advantages from being a well-encultured
        # social partner:
        #   - Preferred as cooperation partner (partner choice)
        #   - Stronger social relationships (built over extended childhood)
        #   - Better cultural integration (skills, norms, conventions)
        #
        # CRITICAL: This benefit REQUIRES cooperative care context.
        # Without caregivers (Sub_care = 0), there is no social environment
        # to provide these benefits, regardless of childhood length.
        #
        # With caregivers (Sub_care > 0), the benefit scales SUPERLINEARLY
        # with childhood duration because:
        #   - Relationships deepen over time
        #   - Social skills compound
        #   - Cultural knowledge integrates
        
        L_normalized = self.L / p.L_max
        social_benefit = (self.Sub_care * p.social_rate * np.power(L_normalized, p.social_exponent) * p.L_max)
        
        total_benefit = benefit_knowledge + social_benefit
        
        # --- Net Fitness ---
        net_fitness = total_benefit - total_cost
        
        # Convert to survival probability via logistic function
        P_survival = 1.0 / (1.0 + np.exp(-net_fitness))
        P_survival = np.clip(P_survival, 0.001, 0.999)
        
        # Store for next generation's social learning
        self.prev_K = self.K.copy()
        self.prev_fitness = P_survival.copy()
        
        # =====================================================================
        # PHASE 3: SELECTION
        # =====================================================================
        
        survived = self.rng.random(self.N) < P_survival
        survivor_idx = np.where(survived)[0]
        
        # Ensure minimum viable population
        if len(survivor_idx) < 2:
            survivor_idx = np.arange(self.N)
        
        # =====================================================================
        # PHASE 4: REPRODUCTION
        # =====================================================================
        
        # Sexual reproduction with random mating among survivors
        parents = self.rng.choice(survivor_idx, size=(self.N, 2))
        p1, p2 = parents[:, 0], parents[:, 1]
        
        # Mendelian inheritance (random parent for each trait)
        new_B = np.where(self.rng.random(self.N) < 0.5, self.B[p1], self.B[p2])
        new_L = np.where(self.rng.random(self.N) < 0.5, self.L[p1], self.L[p2])
        new_S = np.where(self.rng.random(self.N) < 0.5, self.S[p1], self.S[p2])
        
        # Mutation
        mutants = self.rng.random(self.N) < p.mu
        n_mutants = np.sum(mutants)
        new_B[mutants] += self.rng.normal(0, p.sigma_mut * 5.0, n_mutants)
        new_L[mutants] += self.rng.normal(0, p.sigma_mut * 1.0, n_mutants)
        new_S[mutants] += self.rng.normal(0, p.sigma_mut * 0.1, n_mutants)
        
        # Apply bounds
        self.B = np.clip(new_B, 0.1, self.B_max)
        self.L = np.clip(new_L, 0.1, p.L_max)
        self.S = np.clip(new_S, 0.0, 1.0)
        
        # Return population statistics
        return {
            'B_mean': float(np.mean(self.B)),
            'B_std': float(np.std(self.B)),
            'L_mean': float(np.mean(self.L)),
            'L_std': float(np.std(self.L)),
            'K_mean': float(np.mean(self.K)),
            'S_mean': float(np.mean(self.S)),
        }


# =============================================================================
# SIMULATION FUNCTIONS
# =============================================================================

def run_replicate(
    replicate_id: int,
    group: str,
    seed: int,
    params: SimulationParameters
) -> Dict:
    """
    Run a single replicate of the simulation.
    
    Args:
        replicate_id: Identifier for this replicate
        group: 'control' or 'treatment'
        seed: Random seed
        params: Simulation parameters
    
    Returns:
        Dictionary containing full results for this replicate
    """
    pop = Population(params, seed)
    history = []
    
    # -------------------------------------------------------------------------
    # EPOCH 1: Baseline (both groups identical)
    # -------------------------------------------------------------------------
    pop.B_max = params.B_max
    pop.Sub_care = 0.0
    
    for gen in range(params.T_epoch1):
        stats = pop.step()
        if gen % 100 == 0:  # Sample every 100 generations
            history.append({
                'generation': gen,
                'phase': 1,
                **stats
            })
    
    # -------------------------------------------------------------------------
    # EPOCH 2: Treatment divergence
    # -------------------------------------------------------------------------
    pop.Sub_care = params.Sub_care_treatment if group == 'treatment' else 0.0
    
    for gen in range(params.T_epoch1, params.T_epoch2):
        stats = pop.step()
        if gen % 100 == 0:
            history.append({
                'generation': gen,
                'phase': 2,
                **stats
            })
    
    # Record end-of-phase-2 values (primary outcome)
    epoch2_B = stats['B_mean']
    epoch2_L = stats['L_mean']
    
    return {
        'replicate_id': replicate_id,
        'group': group,
        'seed': seed,
        'history': history,
        'epoch2_B': epoch2_B,
        'epoch2_L': epoch2_L,
        'final_B': epoch2_B,
        'final_L': epoch2_L,
        'final_S': stats['S_mean'],
        'final_K': stats['K_mean'],
    }


def analyze_results(
    results: List[Dict],
    params: SimulationParameters
) -> Dict:
    """
    Perform statistical analysis of simulation results.
    
    Tests:
        Decoupling Hypothesis: L_treatment > L_control AND B_treatment ≈ B_control
    
    Args:
        results: List of replicate result dictionaries
        params: Simulation parameters
    
    Returns:
        Dictionary containing all analysis results
    """
    # Separate groups
    ctrl = [r for r in results if r['group'] == 'control']
    treat = [r for r in results if r['group'] == 'treatment']
    
    # Extract primary outcomes
    ctrl_B = np.array([r['epoch2_B'] for r in ctrl])
    ctrl_L = np.array([r['epoch2_L'] for r in ctrl])
    treat_B = np.array([r['epoch2_B'] for r in treat])
    treat_L = np.array([r['epoch2_L'] for r in treat])
    
    # Descriptive statistics
    n_ctrl, n_treat = len(ctrl), len(treat)
    
    L_ctrl_mean, L_ctrl_std = np.mean(ctrl_L), np.std(ctrl_L)
    L_treat_mean, L_treat_std = np.mean(treat_L), np.std(treat_L)
    B_ctrl_mean, B_ctrl_std = np.mean(ctrl_B), np.std(ctrl_B)
    B_treat_mean, B_treat_std = np.mean(treat_B), np.std(treat_B)
    
    L_ctrl_se = L_ctrl_std / np.sqrt(n_ctrl)
    L_treat_se = L_treat_std / np.sqrt(n_treat)
    B_ctrl_se = B_ctrl_std / np.sqrt(n_ctrl)
    B_treat_se = B_treat_std / np.sqrt(n_treat)
    
    # Effect sizes
    rel_dB = (B_treat_mean - B_ctrl_mean) / B_ctrl_mean
    rel_dL = (L_treat_mean - L_ctrl_mean) / L_ctrl_mean
    
    # Cohen's d
    pooled_std_L = np.sqrt(((n_ctrl-1)*L_ctrl_std**2 + (n_treat-1)*L_treat_std**2) / 
                           (n_ctrl + n_treat - 2))
    pooled_std_B = np.sqrt(((n_ctrl-1)*B_ctrl_std**2 + (n_treat-1)*B_treat_std**2) / 
                           (n_ctrl + n_treat - 2))
    
    cohen_d_L = (L_treat_mean - L_ctrl_mean) / pooled_std_L
    cohen_d_B = (B_treat_mean - B_ctrl_mean) / pooled_std_B
    
    # -------------------------------------------------------------------------
    # Decoupling Hypothesis Test
    # -------------------------------------------------------------------------
    # Two-sided treatment-control comparison, matching the manuscript tables.
    L_ttest = scipy_stats.ttest_ind(treat_L, ctrl_L)
    L_significant = L_ttest.pvalue < params.alpha_level
    
    # B should be equivalent (within threshold)
    B_equivalent = abs(rel_dB) <= params.rel_dB_threshold
    
    # Decoupling hypothesis supported if both conditions met
    decoupling_supported = L_significant and B_equivalent
    
    # -------------------------------------------------------------------------
    # Decoupling criterion (effect size based)
    # -------------------------------------------------------------------------
    decoupled = (abs(rel_dB) <= params.rel_dB_threshold and 
                 rel_dL >= params.rel_dL_threshold)
    
    return {
        # Sample sizes
        'n_control': n_ctrl,
        'n_treatment': n_treat,
        
        # Childhood duration
        'L_ctrl_mean': float(L_ctrl_mean),
        'L_ctrl_std': float(L_ctrl_std),
        'L_ctrl_se': float(L_ctrl_se),
        'L_treat_mean': float(L_treat_mean),
        'L_treat_std': float(L_treat_std),
        'L_treat_se': float(L_treat_se),
        
        # Brain size
        'B_ctrl_mean': float(B_ctrl_mean),
        'B_ctrl_std': float(B_ctrl_std),
        'B_ctrl_se': float(B_ctrl_se),
        'B_treat_mean': float(B_treat_mean),
        'B_treat_std': float(B_treat_std),
        'B_treat_se': float(B_treat_se),
        
        # Effect sizes
        'rel_dL': float(rel_dL),
        'rel_dB': float(rel_dB),
        'cohen_d_L': float(cohen_d_L),
        'cohen_d_B': float(cohen_d_B),
        
        # Decoupling hypothesis test
        'L_tstatistic': float(L_ttest.statistic),
        'L_pvalue': float(L_ttest.pvalue),
        'L_significant': bool(L_significant),
        'B_equivalent': bool(B_equivalent),
        'decoupling_supported': bool(decoupling_supported),
        
        # Overall
        'decoupled': bool(decoupled),
    }


# =============================================================================
# SENSITIVITY ANALYSIS (GRID SEARCH)
# =============================================================================

def run_sensitivity_analysis(args, base_params):
    """
    Performs a grid search over social_rate and social_exponent to test
    the robustness of the decoupling effect.
    
    Generates data for Figure 1-figure supplement 5 (heatmap). 20 replicates per cell per condition.
    """
    print("\n" + "="*70)
    print(" RUNNING SENSITIVITY ANALYSIS (GRID SEARCH)")
    print("="*70)
    
    # Define the parameter space
    # Grid centered on the default values (r = 0.30, e = 1.8)
    rates = np.linspace(0.10, 0.50, 9)      # 0.10, 0.15, ..., 0.50
    exponents = np.linspace(1.0, 2.6, 9)    # 1.0, 1.2, ..., 2.6
    
    # Create grid of tasks
    tasks = []
    task_id = 0
    reps_per_cell = 20  # Lower reps for grid search to save time (sufficient for heatmaps)
    
    print(f" Grid: {len(rates)} rates x {len(exponents)} exponents")
    print(f" Replicates per cell: {reps_per_cell}")
    print(f" Total simulations: {len(rates) * len(exponents) * reps_per_cell * 2} (Treatment + Control)")
    
    for r in rates:
        for e in exponents:
            # Create specific params for this cell
            cell_params = replace(base_params)
            cell_params.social_rate = r
            cell_params.social_exponent = e
            
            # Generate seeds for this cell
            cell_rng = np.random.default_rng(args.seed + task_id)
            seeds = cell_rng.integers(0, int(1e9), reps_per_cell)
            
            for i in range(reps_per_cell):
                # We need both Control and Treatment for each cell to calculate Delta L
                # Control run (always needed as baseline, though technically stable across params)
                tasks.append((task_id, 'control', seeds[i], cell_params, r, e))
                # Treatment run
                tasks.append((task_id, 'treatment', seeds[i], cell_params, r, e))
                task_id += 1

    # Wrapper for parallel execution that returns params too
    def run_grid_replicate(tid, grp, sd, p, rate, exp):
        res = run_replicate(tid, grp, sd, p)
        return {
            'rate': rate,
            'exponent': exp,
            'group': grp,
            'epoch2_L': res['epoch2_L'],
            'epoch2_B': res['epoch2_B']
        }

    # Execute
    if PARALLEL_AVAILABLE and args.n_jobs != 1:
        results = Parallel(n_jobs=args.n_jobs, verbose=5)(
            delayed(run_grid_replicate)(*t) for t in tasks
        )
    else:
        results = [run_grid_replicate(*t) for t in tasks]

    # Process Results into a DataFrame
    df = pd.DataFrame(results)
    
    # Aggregate by cell
    summary = df.groupby(['rate', 'exponent', 'group']).agg({
        'epoch2_L': 'mean',
        'epoch2_B': 'mean'
    }).reset_index()
    
    # Pivot to compare Treatment vs Control
    pivot = summary.pivot(index=['rate', 'exponent'], columns='group')
    
    # Calculate % Change
    # (Treatment - Control) / Control
    analysis = pd.DataFrame()
    analysis['rate'] = pivot.index.get_level_values('rate')
    analysis['exponent'] = pivot.index.get_level_values('exponent')
    
    # L change
    L_treat = pivot[('epoch2_L', 'treatment')].values
    L_ctrl = pivot[('epoch2_L', 'control')].values
    analysis['pct_change_L'] = (L_treat - L_ctrl) / L_ctrl * 100.0
    
    # B change
    B_treat = pivot[('epoch2_B', 'treatment')].values
    B_ctrl = pivot[('epoch2_B', 'control')].values
    analysis['pct_change_B'] = (B_treat - B_ctrl) / B_ctrl * 100.0
    
    # Define "Decoupling" boolean
    # L extends >= 10% AND B changes < 2%
    analysis['is_decoupled'] = (analysis['pct_change_L'] >= 10.0) & (abs(analysis['pct_change_B']) < 2.0)
    
    # Save
    os.makedirs(args.output, exist_ok=True)
    analysis.to_csv(f"{args.output}/sensitivity_analysis.csv", index=False)
    analysis.to_csv(f"{args.output}/Data_S2_Sensitivity.csv", index=False)
    print(f"\nSensitivity data saved to: {args.output}/sensitivity_analysis.csv")
    
    # Optional: Generate Quick Plot (if matplotlib available)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # Matrix for Heatmap
        heatmap_data = analysis.pivot(index='exponent', columns='rate', values='pct_change_L')
        
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(heatmap_data, annot=True, fmt=".1f", cmap="viridis", cbar_kws={'label': '% Change in Childhood Duration'})
        plt.title('Sensitivity Analysis: Effect of Social Parameters on Childhood Extension')
        plt.gca().invert_yaxis() # Put low exponents at bottom
        
        # Highlight the decoupling zone
        # We can overlay stippling or a contour for 'is_decoupled'
        # For now, just saving the heatmap
        plt.savefig(f"{args.output}/sensitivity_heatmap_preliminary.pdf")
        print(f"Preliminary plot saved to: {args.output}/sensitivity_heatmap_preliminary.pdf")
        
    except ImportError:
        print("Matplotlib/Seaborn not found. Skipping plot generation.")


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main entry point for the simulation."""
    
    parser = argparse.ArgumentParser(
        description='Childhood-First Hypothesis: Evolutionary Simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python childhood_first_simulation.py -n 10 -o ./test
  python childhood_first_simulation.py -n 100 -o ./results
  python childhood_first_simulation.py -n 100 --social-rate 0.25 -o ./results
  
  # Run Sensitivity Analysis
  python childhood_first_simulation.py --sensitivity -o ./results
        """
    )
    
    parser.add_argument('-n', '--n-replicates', type=int, default=100,
                        help='Number of replicates per group (default: 100)')
    parser.add_argument('-o', '--output', type=str, default='./results',
                        help='Output directory (default: ./results)')
    parser.add_argument('-j', '--n-jobs', type=int, default=-1,
                        help='Number of parallel jobs (-1 for all cores)')
    parser.add_argument('--social-rate', type=float, default=0.30,
                        help='Social fitness rate parameter (default: 0.30)')
    parser.add_argument('--social-exp', type=float, default=1.8,
                        help='Social fitness exponent (default: 1.8)')
    parser.add_argument('--sub-care', type=float, default=0.6,
                        help='Treatment cooperative care level (default: 0.6)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Master random seed (default: 42)')
    
    # Sensitivity-analysis flag
    parser.add_argument('--sensitivity', action='store_true',
                        help='Run the parameter sweep for Figure 1-figure supplement 5 instead of a single simulation')
    
    args = parser.parse_args()
    
    # Create parameters (base)
    params = SimulationParameters(
        social_rate=args.social_rate,
        social_exponent=args.social_exp,
        Sub_care_treatment=args.sub_care,
    )
    
    # BRANCH: Run sensitivity analysis or the standard simulation
    if args.sensitivity:
        run_sensitivity_analysis(args, params)
        return
    
    # ... [Existing standard simulation logic] ...
    
    N_REPS = args.n_replicates
    
    # Print header
    print("=" * 70)
    print(" CHILDHOOD-FIRST HYPOTHESIS: EVOLUTIONARY SIMULATION")
    print("=" * 70)
    print()
    print(" Testing whether cooperative care creates selection for extended")
    print(" childhood independent of brain size expansion.")
    print()
    print(" Author: Tobias Grossmann (tg3ny@virginia.edu)")
    print()
    print(f" Replicates per group: {N_REPS}")
    print(f" Total simulations: {N_REPS * 2}")
    print(f" Generations per simulation: {params.T_total:,}")
    print()
    print(" Key parameters:")
    print(f"   Social fitness rate: {params.social_rate}")
    print(f"   Social fitness exponent: {params.social_exponent}")
    print(f"   Treatment cooperative care: {params.Sub_care_treatment}")
    print(f"   Knowledge benefit (lambda): {params.lambda_}")
    print()
    print("=" * 70)
    print()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Generate seeds
    master_rng = np.random.default_rng(args.seed)
    seeds_ctrl = master_rng.integers(0, int(1e9), N_REPS)
    seeds_treat = master_rng.integers(0, int(1e9), N_REPS)
    
    # Prepare tasks
    tasks = []
    for i in range(N_REPS):
        tasks.append((i + 1, 'control', int(seeds_ctrl[i]), params))
    for i in range(N_REPS):
        tasks.append((i + 1, 'treatment', int(seeds_treat[i]), params))
    
    # Run simulations
    print(f"Running {len(tasks)} simulations...")
    start_time = time.time()
    
    if PARALLEL_AVAILABLE and args.n_jobs != 1:
        results = Parallel(n_jobs=args.n_jobs, verbose=10)(
            delayed(run_replicate)(rid, grp, seed, p) 
            for rid, grp, seed, p in tasks
        )
    else:
        results = []
        for i, (rid, grp, seed, p) in enumerate(tasks):
            if (i + 1) % 10 == 0:
                print(f"  Completed {i + 1}/{len(tasks)}")
            results.append(run_replicate(rid, grp, seed, p))
    
    elapsed = (time.time() - start_time) / 60.0
    print(f"\nCompleted in {elapsed:.1f} minutes")
    print()
    
    # -------------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------------
    
    # Summary CSV
    summary_data = [{
        'replicate': r['replicate_id'],
        'group': r['group'],
        'seed': r['seed'],
        'epoch2_B': r['epoch2_B'],
        'epoch2_L': r['epoch2_L'],
        'final_B': r['final_B'],
        'final_L': r['final_L'],
        'final_S': r['final_S'],
        'final_K': r['final_K'],
    } for r in results]
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(f"{args.output}/experiment_summary.csv", index=False)
    summary_df.to_csv(f"{args.output}/Data_S1_Summary.csv", index=False)
    
    # Full trajectories CSV
    trajectory_data = []
    for r in results:
        for h in r['history']:
            trajectory_data.append({
                'replicate': r['replicate_id'],
                'group': r['group'],
                **h
            })
    
    traj_df = pd.DataFrame(trajectory_data)
    traj_df.to_csv(f"{args.output}/full_trajectories.csv", index=False)
    traj_df.to_csv(f"{args.output}/Data_S3_Trajectories.csv", index=False)
    
    # -------------------------------------------------------------------------
    # Analysis
    # -------------------------------------------------------------------------
    
    print("=" * 70)
    print(" RESULTS")
    print("=" * 70)
    
    analysis = analyze_results(results, params)
    
    print()
    print("-" * 70)
    print(" DECOUPLING HYPOTHESIS TEST")
    print(" Does cooperative care enable childhood extension without brain expansion?")
    print("-" * 70)
    print()
    print(f"  Childhood Duration (L):")
    print(f"    Control:   {analysis['L_ctrl_mean']:.2f} +/- {analysis['L_ctrl_se']:.2f}")
    print(f"    Treatment: {analysis['L_treat_mean']:.2f} +/- {analysis['L_treat_se']:.2f}")
    print(f"    Delta = {analysis['rel_dL']*100:+.1f}%, Cohen's d = {analysis['cohen_d_L']:.2f}")
    print(f"    t = {analysis['L_tstatistic']:.2f}, p = {analysis['L_pvalue']:.6f}")
    sig_marker = "***" if analysis['L_pvalue'] < 0.001 else "**" if analysis['L_pvalue'] < 0.01 else "*" if analysis['L_pvalue'] < 0.05 else ""
    print(f"    Significant (p < 0.005): {'YES ' + sig_marker if analysis['L_significant'] else 'NO'}")
    print()
    print(f"  Brain Size (B):")
    print(f"    Control:   {analysis['B_ctrl_mean']:.2f} +/- {analysis['B_ctrl_se']:.2f}")
    print(f"    Treatment: {analysis['B_treat_mean']:.2f} +/- {analysis['B_treat_se']:.2f}")
    print(f"    Delta = {analysis['rel_dB']*100:+.1f}%, Cohen's d = {analysis['cohen_d_B']:.2f}")
    print(f"    Equivalent (|Delta| < 2%): {'YES' if analysis['B_equivalent'] else 'NO'}")
    print()
    print(f"  >>> Decoupling Hypothesis: {'*** SUPPORTED ***' if analysis['decoupling_supported'] else 'NOT SUPPORTED'}")
    
    print()
    print("-" * 70)
    print(" DECOUPLING CRITERION (EFFECT SIZE)")
    print("-" * 70)
    print()
    print(f"  |Delta B| < 2%:  {abs(analysis['rel_dB']*100):.1f}% {'[PASS]' if abs(analysis['rel_dB']) < 0.02 else '[FAIL]'}")
    print(f"  Delta L >= 10%:  {analysis['rel_dL']*100:.1f}% {'[PASS]' if analysis['rel_dL'] >= 0.10 else '[FAIL]'}")
    print()
    print(f"  >>> Decoupling: {'*** YES ***' if analysis['decoupled'] else 'NO'}")
    
    print()
    print("=" * 70)
    if analysis['decoupling_supported']:
        print(" CHILDHOOD-FIRST HYPOTHESIS: SUPPORTED")
        print(" Cooperative care enables childhood extension without requiring")
        print(" brain expansion, matching the fossil record at Dmanisi and Rising Star.")
    elif analysis['decoupled']:
        print(" DECOUPLING OBSERVED")
        print(" Brain-childhood decoupling achieved (effect size criterion met).")
    else:
        print(" CHILDHOOD-FIRST HYPOTHESIS: NOT SUPPORTED")
    print("=" * 70)
    
    # -------------------------------------------------------------------------
    # Save analysis JSON
    # -------------------------------------------------------------------------
    
    output_json = {
        'model': 'Childhood-First Evolutionary Simulation',
        'version': '2.1',
        'reference': 'Grossmann, T. (tg3ny@virginia.edu)',
        'parameters': asdict(params),
        'n_replicates_per_group': N_REPS,
        'master_seed': args.seed,
        'runtime_minutes': elapsed,
        'analysis': analysis,
    }
    
    with open(f"{args.output}/experiment_results.json", 'w') as f:
        json.dump(output_json, f, indent=2)
    
    print()
    print(f" Results saved to: {args.output}/")
    print("   - experiment_summary.csv")
    print("   - full_trajectories.csv")
    print("   - experiment_results.json")
    print()


if __name__ == "__main__":
    main()
