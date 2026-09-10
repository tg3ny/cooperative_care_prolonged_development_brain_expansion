"""
Canonical script for Simulations 2 and 3 (care-gated grey-ceiling analyses).

Produces the authoritative outputs for Tables 2-3, Figures 3-4, and the
gradual-ramp robustness analysis (Figure 3-figure supplement 1).

Outputs:
  Simulation 2 (care-gated ceiling condition):
    Data_S6_Pure_Summary.csv   — endpoint values (canonical for Table 2)
    Data_S8_Pure_Results.json  — full analysis with parameters
    experiment2_trajectories.csv — generation-level trajectories (for Fig. 3)

  Simulation 3 (cooperative-care vs care-independent energetic routes):
    Data_S7_Mixed_Summary.csv  — endpoint values (canonical for Table 3)
    Data_S9_Mixed_Results.json — full analysis with parameters

  Robustness:
    Data_S10_Ramp_Results.json — gradual ceiling ramp (for Figure 3-figure supplement 1)

Terminology: Phase 1 (baseline), Phase 2 (care divergence), Phase 3 (ceiling lift).
Master seed: 42. All results are fully deterministic.
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
    Parameters for the Care-Gated Grey Ceiling simulation.

    Phases 1-2 use identical parameters to Simulation 1 (Grossmann, 2026),
    which itself builds on the Cultural Brain Hypothesis model
    (Muthukrishna et al., 2018). Phase 3 introduces the care-gated
    ceiling-lift mechanism.
    """

    # --- Population ---
    N: int = 1000                    # Population size
    T_phase1: int = 10_000           # End of baseline phase
    T_phase2: int = 100_000          # End of treatment phase
    T_phase3: int = 200_000          # End of ceiling-lift phase

    # --- Metabolic Costs ---
    # (Identical to Simulation 1)
    alpha: float = 0.1               # Linear brain metabolic cost
    gamma: float = 0.001             # Quadratic brain metabolic cost
    beta: float = 0.12               # Temporal cost of childhood per unit L
    k_ceiling: float = 0.5           # Steepness of brain size ceiling

    # --- Knowledge/Learning ---
    # (Identical to Simulation 1)
    lambda_: float = 0.5             # Benefit coefficient for knowledge
    learning_exponent: float = 1.15  # Superlinearity of childhood learning
    asocial_rate: float = 0.10       # Base rate of asocial learning

    # --- Social Fitness (Key Mechanism) ---
    # (Identical to Simulation 1)
    social_rate: float = 0.30        # Strength of social fitness benefit
    social_exponent: float = 1.8     # How superlinearly L affects social benefit

    # --- Cooperative Care ---
    # (Identical to Simulation 1)
    Sub_care_control: float = 0.0    # Care level for control group
    Sub_care_treatment: float = 0.6  # Care level for treatment group

    # --- Evolution ---
    # (Identical to Simulation 1)
    mu: float = 0.01                 # Mutation probability per individual
    sigma_mut: float = 0.1           # Mutation effect size (scaled per trait)

    # --- Trait Bounds ---
    B_max_initial: float = 50.0      # Brain ceiling during Phases 1-2
    L_max: float = 20.0              # Maximum childhood duration

    # --- Care-Gated Ceiling Lift (Phase 3) ---
    #
    # The grey ceiling on brain size is lifted as a function of
    # cooperative care:
    #
    #   B_max = B_max_initial + exogenous_lift + Sub_care * care_lift
    #
    # care_lift: Imposed care-gated energetic headroom. It formalizes the
    #   hypothesis that pooled provisioning can relax an energetic constraint;
    #   it is not derived from measured caloric transfers.
    #
    # exogenous_lift: Care-independent energetic headroom, applied equally
    #   to both groups without assigning it to a specific historical cause.
    #
    # PURE condition:  exogenous_lift=0,  care_lift=100
    #   Treatment B_max = 50 + 0 + 0.6*100 = 110
    #   Control   B_max = 50 + 0 + 0.0*100 = 50  (unchanged)
    #
    # MIXED condition: exogenous_lift=20, care_lift=80
    #   Treatment B_max = 50 + 20 + 0.6*80 = 118
    #   Control   B_max = 50 + 20 + 0.0*80 = 70
    #
    care_lift: float = 100.0         # Care-dependent ceiling lift
    exogenous_lift: float = 0.0      # Care-independent ceiling lift

    # --- Ceiling Lift Timing ---
    ceiling_ramp: int = 0            # Generations over which to raise ceiling
                                     # 0 = instantaneous at Phase 3 start
                                     # >0 = linear ramp over this many gens

    # --- Statistical Thresholds ---
    alpha_level: float = 0.005       # Significance level for hypothesis tests


# =============================================================================
# POPULATION CLASS
# =============================================================================

class Population:
    """
    Represents an evolving population with heritable traits for brain size (B),
    childhood duration (L), and social learning propensity (S).

    This class is identical to the Population class in Simulation 1
    (childhood_first_simulation.py; Grossmann, 2026), which was adapted
    from the Cultural Brain Hypothesis model (Muthukrishna et al., 2018).

    Modifications from Muthukrishna et al. (2018):
      - Structural costs (alpha*B + gamma*B^2) separated from temporal
        costs (beta*L)
      - Social fitness mechanism gated by cooperative care (Sub_care)
      - Group size (N) is not a co-evolving trait; population size is fixed
      - Knowledge ceiling (K <= 0.6*B) retained from CBH
    """

    def __init__(self, params: SimulationParameters, seed: int):
        self.params = params
        self.rng = np.random.default_rng(seed)
        self.N = params.N

        # Environmental state (set by simulation phases)
        self.B_max = params.B_max_initial
        self.Sub_care = 0.0

        # Initialize traits with variation
        self.B = self.rng.uniform(5, 15, self.N)
        self.L = self.rng.uniform(1, 3, self.N)
        self.S = self.rng.uniform(0.3, 0.7, self.N)
        self.K = np.zeros(self.N)

        # For social learning across generations
        self.prev_K = None
        self.prev_fitness = None

    def step(self) -> Dict[str, float]:
        """
        Execute one generation: learning, selection, reproduction.
        Identical to Simulation 1. Phase 3 operates through self.B_max.
        """
        p = self.params

        # === LEARNING ===
        L_factor = np.power(self.L, p.learning_exponent)
        K_asocial = self.B * L_factor * p.asocial_rate

        K_social = np.zeros(self.N)
        if self.prev_K is not None:
            probs = self.prev_fitness / np.sum(self.prev_fitness)
            models = self.rng.choice(len(self.prev_K), size=(self.N, 3), p=probs)
            K_social = np.mean(self.prev_K[models], axis=1)

        K_from_social = self.S * K_social * (self.L / p.L_max)
        K_from_asocial = (1 - self.S) * K_asocial
        K_total = K_from_social + K_from_asocial

        # Brain capacity constraint on knowledge (shared with CBH)
        self.K = np.minimum(K_total, self.B * 0.6)

        # === FITNESS ===
        # Costs
        C_struct = p.alpha * self.B + p.gamma * (self.B ** 2)
        ceiling_penalty = 1.0 / (1.0 + np.exp(-p.k_ceiling * (self.B - self.B_max)))
        C_struct += 10.0 * ceiling_penalty
        C_time = p.beta * self.L
        total_cost = C_struct + C_time

        # Benefits
        benefit_knowledge = p.lambda_ * self.K
        L_normalized = self.L / p.L_max
        social_benefit = (self.Sub_care * p.social_rate
                          * np.power(L_normalized, p.social_exponent)
                          * p.L_max)
        total_benefit = benefit_knowledge + social_benefit

        net_fitness = total_benefit - total_cost
        P_survival = 1.0 / (1.0 + np.exp(-net_fitness))
        P_survival = np.clip(P_survival, 0.001, 0.999)

        self.prev_K = self.K.copy()
        self.prev_fitness = P_survival.copy()

        # === SELECTION ===
        survived = self.rng.random(self.N) < P_survival
        survivor_idx = np.where(survived)[0]
        if len(survivor_idx) < 2:
            survivor_idx = np.arange(self.N)

        # === REPRODUCTION ===
        parents = self.rng.choice(survivor_idx, size=(self.N, 2))
        p1, p2 = parents[:, 0], parents[:, 1]

        new_B = np.where(self.rng.random(self.N) < 0.5, self.B[p1], self.B[p2])
        new_L = np.where(self.rng.random(self.N) < 0.5, self.L[p1], self.L[p2])
        new_S = np.where(self.rng.random(self.N) < 0.5, self.S[p1], self.S[p2])

        mutants = self.rng.random(self.N) < p.mu
        n_mutants = np.sum(mutants)
        new_B[mutants] += self.rng.normal(0, p.sigma_mut * 5.0, n_mutants)
        new_L[mutants] += self.rng.normal(0, p.sigma_mut * 1.0, n_mutants)
        new_S[mutants] += self.rng.normal(0, p.sigma_mut * 0.1, n_mutants)

        self.B = np.clip(new_B, 0.1, self.B_max)
        self.L = np.clip(new_L, 0.1, p.L_max)
        self.S = np.clip(new_S, 0.0, 1.0)

        return {
            'B_mean': float(np.mean(self.B)),
            'B_std': float(np.std(self.B)),
            'L_mean': float(np.mean(self.L)),
            'L_std': float(np.std(self.L)),
            'K_mean': float(np.mean(self.K)),
            'K_std': float(np.std(self.K)),
            'S_mean': float(np.mean(self.S)),
            'S_std': float(np.std(self.S)),
            'K_ceiling_frac': float(np.mean(
                self.K >= (self.B * 0.6 - 0.01)
            )),
        }


# =============================================================================
# SIMULATION
# =============================================================================

def compute_phase3_bmax(
    gen: int,
    group: str,
    params: SimulationParameters
) -> float:
    """
    Compute B_max for a given generation in Phase 3.

    The ceiling lift is gated by the cooperative-care parameter:
      B_max = B_max_initial + exogenous_lift + Sub_care * care_lift

    For the treatment group (Sub_care = 0.6), the care-gated term implements
    the hypothesis that pooled provisioning can relax an energetic constraint.
    For the control group (Sub_care = 0.0), only the care-independent term
    contributes. Both terms are imposed model manipulations.

    The lift can be instantaneous or gradual (linear ramp).
    """
    sub_care = (params.Sub_care_treatment if group == 'treatment'
                else params.Sub_care_control)

    # Target B_max for this group
    target = (params.B_max_initial
              + params.exogenous_lift
              + sub_care * params.care_lift)

    if params.ceiling_ramp <= 0:
        return target

    # Gradual ramp from initial to target
    progress = min(1.0, (gen - params.T_phase2) / params.ceiling_ramp)
    return params.B_max_initial + progress * (target - params.B_max_initial)


def run_replicate(
    replicate_id: int,
    group: str,
    seed: int,
    params: SimulationParameters
) -> Dict:
    """
    Run a single replicate through all three phases.

    Phases 1-2 reproduce Simulation 1 exactly (same code, same seeds).
    Phase 3 applies the care-gated ceiling lift.
    """
    pop = Population(params, seed)
    history = []

    # -----------------------------------------------------------------
    # PHASE 1: Baseline (identical to Simulation 1, Phase 1)
    # -----------------------------------------------------------------
    pop.B_max = params.B_max_initial
    pop.Sub_care = 0.0

    for gen in range(params.T_phase1):
        stats = pop.step()
        if gen % 100 == 0:
            history.append({
                'generation': gen,
                'phase': 1,
                'B_max': pop.B_max,
                **stats
            })

    # -----------------------------------------------------------------
    # PHASE 2: Treatment divergence (identical to Simulation 1, Phase 2)
    # -----------------------------------------------------------------
    pop.Sub_care = params.Sub_care_treatment if group == 'treatment' else 0.0

    for gen in range(params.T_phase1, params.T_phase2):
        stats = pop.step()
        if gen % 100 == 0:
            history.append({
                'generation': gen,
                'phase': 2,
                'B_max': pop.B_max,
                **stats
            })

    phase2_B = stats['B_mean']
    phase2_L = stats['L_mean']
    phase2_S = stats['S_mean']
    phase2_K = stats['K_mean']

    # -----------------------------------------------------------------
    # PHASE 3: Care-gated ceiling lift (NEW)
    #
    # Cooperative care serves a DUAL ROLE:
    #   (1) Social fitness: creates selection for extended childhood
    #       and enhanced social learning (established in Phase 2)
    #   (2) Metabolic subsidy: allomaternal provisioning lifts the
    #       grey ceiling, funding brain expansion (this phase)
    #
    # B_max = B_max_initial + exogenous_lift + Sub_care * care_lift
    #
    # The treatment group's elevated social learning from Phase 2
    # means their knowledge is at the storage ceiling (K = 0.6B),
    # so the marginal benefit of brain expansion is immediately
    # realized. The control group (if it gets any ceiling lift via
    # exogenous_lift) benefits less because social learning is lower
    # and knowledge is less constrained.
    # -----------------------------------------------------------------

    for gen in range(params.T_phase2, params.T_phase3):
        pop.B_max = compute_phase3_bmax(gen, group, params)

        stats = pop.step()
        if gen % 100 == 0:
            history.append({
                'generation': gen,
                'phase': 3,
                'B_max': pop.B_max,
                **stats
            })

    phase3_B = stats['B_mean']
    phase3_L = stats['L_mean']
    phase3_S = stats['S_mean']
    phase3_K = stats['K_mean']

    return {
        'replicate_id': replicate_id,
        'group': group,
        'seed': seed,
        'history': history,
        'phase2_B': phase2_B,
        'phase2_L': phase2_L,
        'phase2_S': phase2_S,
        'phase2_K': phase2_K,
        'phase3_B': phase3_B,
        'phase3_L': phase3_L,
        'phase3_S': phase3_S,
        'phase3_K': phase3_K,
        'phase3_B_max': pop.B_max,
        'delta_B': phase3_B - phase2_B,
        'delta_L': phase3_L - phase2_L,
    }


# =============================================================================
# ANALYSIS
# =============================================================================

def analyze_results(
    results: List[Dict],
    params: SimulationParameters
) -> Dict:
    """
    Statistical analysis of the care-gated grey-ceiling simulation.
    """
    ctrl = [r for r in results if r['group'] == 'control']
    treat = [r for r in results if r['group'] == 'treatment']
    n_ctrl, n_treat = len(ctrl), len(treat)

    def desc(arr):
        return {
            'mean': float(np.mean(arr)),
            'std': float(np.std(arr)),
            'se': float(np.std(arr) / np.sqrt(len(arr))),
        }

    # --- Phase 2 Validation ---
    p2_ctrl_B = np.array([r['phase2_B'] for r in ctrl])
    p2_ctrl_L = np.array([r['phase2_L'] for r in ctrl])
    p2_ctrl_S = np.array([r['phase2_S'] for r in ctrl])
    p2_treat_B = np.array([r['phase2_B'] for r in treat])
    p2_treat_L = np.array([r['phase2_L'] for r in treat])
    p2_treat_S = np.array([r['phase2_S'] for r in treat])

    # --- Phase 3 Outcomes ---
    p3_ctrl_B = np.array([r['phase3_B'] for r in ctrl])
    p3_ctrl_L = np.array([r['phase3_L'] for r in ctrl])
    p3_ctrl_S = np.array([r['phase3_S'] for r in ctrl])
    p3_ctrl_K = np.array([r['phase3_K'] for r in ctrl])
    p3_treat_B = np.array([r['phase3_B'] for r in treat])
    p3_treat_L = np.array([r['phase3_L'] for r in treat])
    p3_treat_S = np.array([r['phase3_S'] for r in treat])
    p3_treat_K = np.array([r['phase3_K'] for r in treat])

    ctrl_delta_B = np.array([r['delta_B'] for r in ctrl])
    treat_delta_B = np.array([r['delta_B'] for r in treat])
    ctrl_delta_L = np.array([r['delta_L'] for r in ctrl])
    treat_delta_L = np.array([r['delta_L'] for r in treat])

    # --- Ceiling values ---
    ctrl_bmax = params.B_max_initial + params.exogenous_lift
    treat_bmax = (params.B_max_initial + params.exogenous_lift
                  + params.Sub_care_treatment * params.care_lift)

    # === Hypothesis Tests ===

    # PRIMARY: Two-sided treatment-control comparison of brain expansion.
    dB_ttest = scipy_stats.ttest_ind(treat_delta_B, ctrl_delta_B)
    pooled_std_dB = np.sqrt(
        ((n_ctrl - 1) * np.std(ctrl_delta_B)**2
         + (n_treat - 1) * np.std(treat_delta_B)**2)
        / (n_ctrl + n_treat - 2)
    )
    cohen_d_dB = ((np.mean(treat_delta_B) - np.mean(ctrl_delta_B))
                  / pooled_std_dB) if pooled_std_dB > 0 else 0.0

    # SECONDARY: Two-sided Phase 3 brain comparison.
    p3_B_ttest = scipy_stats.ttest_ind(p3_treat_B, p3_ctrl_B)
    pooled_std_p3B = np.sqrt(
        ((n_ctrl - 1) * np.std(p3_ctrl_B)**2
         + (n_treat - 1) * np.std(p3_treat_B)**2)
        / (n_ctrl + n_treat - 2)
    )
    cohen_d_p3B = ((np.mean(p3_treat_B) - np.mean(p3_ctrl_B))
                   / pooled_std_p3B) if pooled_std_p3B > 0 else 0.0

    # VALIDATION: Did treatment actually expand?
    treat_expanded = scipy_stats.ttest_1samp(treat_delta_B, 0)
    ctrl_expanded = scipy_stats.ttest_1samp(ctrl_delta_B, 0)

    # EXPLORATORY: Childhood
    p3_L_ttest = scipy_stats.ttest_ind(p3_treat_L, p3_ctrl_L)
    # EXPLORATORY: Knowledge
    p3_K_ttest = scipy_stats.ttest_ind(p3_treat_K, p3_ctrl_K)

    primary_supported = (dB_ttest.pvalue < params.alpha_level
                         and np.mean(treat_delta_B) > np.mean(ctrl_delta_B))

    return {
        'n_control': n_ctrl,
        'n_treatment': n_treat,

        # Ceiling configuration
        'ctrl_B_max_phase3': ctrl_bmax,
        'treat_B_max_phase3': treat_bmax,

        # Phase 2 Validation
        'phase2_ctrl_B': desc(p2_ctrl_B),
        'phase2_ctrl_L': desc(p2_ctrl_L),
        'phase2_ctrl_S': desc(p2_ctrl_S),
        'phase2_treat_B': desc(p2_treat_B),
        'phase2_treat_L': desc(p2_treat_L),
        'phase2_treat_S': desc(p2_treat_S),
        'phase2_L_extension_pct': float(
            (np.mean(p2_treat_L) - np.mean(p2_ctrl_L))
            / np.mean(p2_ctrl_L) * 100
        ),
        'phase2_B_change_pct': float(
            (np.mean(p2_treat_B) - np.mean(p2_ctrl_B))
            / np.mean(p2_ctrl_B) * 100
        ),
        'phase2_S_enhancement_pct': float(
            (np.mean(p2_treat_S) - np.mean(p2_ctrl_S))
            / np.mean(p2_ctrl_S) * 100
        ),

        # Phase 3 Endpoints
        'phase3_ctrl_B': desc(p3_ctrl_B),
        'phase3_ctrl_L': desc(p3_ctrl_L),
        'phase3_ctrl_S': desc(p3_ctrl_S),
        'phase3_ctrl_K': desc(p3_ctrl_K),
        'phase3_treat_B': desc(p3_treat_B),
        'phase3_treat_L': desc(p3_treat_L),
        'phase3_treat_S': desc(p3_treat_S),
        'phase3_treat_K': desc(p3_treat_K),

        # Brain Expansion
        'ctrl_delta_B': desc(ctrl_delta_B),
        'treat_delta_B': desc(treat_delta_B),
        'differential_expansion': float(
            np.mean(treat_delta_B) - np.mean(ctrl_delta_B)
        ),

        # Primary Hypothesis Test
        'dB_t_statistic': float(dB_ttest.statistic),
        'dB_p_value': float(dB_ttest.pvalue),
        'dB_cohen_d': float(cohen_d_dB),
        'dB_significant': bool(dB_ttest.pvalue < params.alpha_level),

        # Secondary
        'p3B_t_statistic': float(p3_B_ttest.statistic),
        'p3B_p_value': float(p3_B_ttest.pvalue),
        'p3B_cohen_d': float(cohen_d_p3B),

        # Validation
        'ctrl_expanded_p': float(ctrl_expanded.pvalue),
        'treat_expanded_p': float(treat_expanded.pvalue),
        'treat_expanded': bool(treat_expanded.pvalue < 0.05),

        # Exploratory
        'ctrl_delta_L': desc(ctrl_delta_L),
        'treat_delta_L': desc(treat_delta_L),
        'p3L_t_statistic': float(p3_L_ttest.statistic),
        'p3L_p_value': float(p3_L_ttest.pvalue),
        'p3K_t_statistic': float(p3_K_ttest.statistic),
        'p3K_p_value': float(p3_K_ttest.pvalue),

        # Overall
        'grey_ceiling_hypothesis_supported': bool(primary_supported),
    }


# =============================================================================
# RUN CONDITION
# =============================================================================

def run_condition(
    condition_name: str,
    params: SimulationParameters,
    n_reps: int,
    master_seed: int,
    output_dir: str,
    n_jobs: int,
):
    """Run a single simulation condition (pure or mixed)."""

    # Effective ceilings
    ctrl_bmax = params.B_max_initial + params.exogenous_lift
    treat_bmax = (params.B_max_initial + params.exogenous_lift
                  + params.Sub_care_treatment * params.care_lift)

    print()
    print("=" * 70)
    print(f" CONDITION: {condition_name.upper()}")
    print("=" * 70)
    print()
    print(f" Care-gated ceiling lift (Hrdy-Isler-van Schaik mechanism):")
    print(f"   B_max = B_max_initial + exogenous_lift + Sub_care * care_lift")
    print(f"   B_max = {params.B_max_initial} + {params.exogenous_lift}"
          f" + Sub_care * {params.care_lift}")
    print()
    print(f" Effective ceilings in Phase 3:")
    print(f"   Control   (Sub_care={params.Sub_care_control}): "
          f"B_max = {ctrl_bmax:.0f}")
    print(f"   Treatment (Sub_care={params.Sub_care_treatment}): "
          f"B_max = {treat_bmax:.0f}")
    print()
    if params.ceiling_ramp > 0:
        print(f" Gradual ramp over {params.ceiling_ramp:,} generations")
    else:
        print(f" Instantaneous lift at generation {params.T_phase2:,}")
    print()
    print(f" Replicates per group: {n_reps}")
    print()

    # Generate seeds (same as Simulation 1)
    master_rng = np.random.default_rng(master_seed)
    seeds_ctrl = master_rng.integers(0, int(1e9), n_reps)
    seeds_treat = master_rng.integers(0, int(1e9), n_reps)

    tasks = []
    for i in range(n_reps):
        tasks.append((i + 1, 'control', int(seeds_ctrl[i]), params))
    for i in range(n_reps):
        tasks.append((i + 1, 'treatment', int(seeds_treat[i]), params))

    print(f"Running {len(tasks)} simulations...")
    start_time = time.time()

    if PARALLEL_AVAILABLE and n_jobs != 1:
        results = Parallel(n_jobs=n_jobs, verbose=10)(
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

    # --- Save ---
    cond_dir = os.path.join(output_dir, condition_name)
    os.makedirs(cond_dir, exist_ok=True)

    summary_data = [{
        'replicate': r['replicate_id'],
        'group': r['group'],
        'seed': r['seed'],
        'phase2_B': r['phase2_B'],
        'phase2_L': r['phase2_L'],
        'phase2_S': r['phase2_S'],
        'phase2_K': r['phase2_K'],
        'phase3_B': r['phase3_B'],
        'phase3_L': r['phase3_L'],
        'phase3_S': r['phase3_S'],
        'phase3_K': r['phase3_K'],
        'phase3_B_max': r['phase3_B_max'],
        'delta_B': r['delta_B'],
        'delta_L': r['delta_L'],
    } for r in results]
    pd.DataFrame(summary_data).to_csv(
        f"{cond_dir}/experiment2_summary.csv", index=False)

    trajectory_data = []
    for r in results:
        for h in r['history']:
            trajectory_data.append({
                'replicate': r['replicate_id'],
                'group': r['group'],
                **h
            })
    pd.DataFrame(trajectory_data).to_csv(
        f"{cond_dir}/experiment2_trajectories.csv", index=False)

    # --- Analysis ---
    analysis = analyze_results(results, params)

    print()
    print("-" * 70)
    print(" PHASE 2 VALIDATION (should match Simulation 1)")
    print("-" * 70)
    print(f"  Childhood Extension: "
          f"{analysis['phase2_L_extension_pct']:+.1f}%")
    print(f"  Brain Size Change:   "
          f"{analysis['phase2_B_change_pct']:+.1f}%")
    print(f"  Social Learning Enhancement: "
          f"{analysis['phase2_S_enhancement_pct']:+.1f}%")

    print()
    print("-" * 70)
    print(" CARE-GATED GREY CEILING HYPOTHESIS TEST")
    print("-" * 70)
    print()
    print(f"  Effective B_max:  Control = {ctrl_bmax:.0f},  "
          f"Treatment = {treat_bmax:.0f}")
    print()
    print(f"  Brain Expansion (Phase 3 - Phase 2):")
    print(f"    Control:   {analysis['ctrl_delta_B']['mean']:+.2f} "
          f"+/- {analysis['ctrl_delta_B']['se']:.2f}")
    print(f"    Treatment: {analysis['treat_delta_B']['mean']:+.2f} "
          f"+/- {analysis['treat_delta_B']['se']:.2f}")
    print(f"    Differential: "
          f"{analysis['differential_expansion']:+.2f}")
    print(f"    t = {analysis['dB_t_statistic']:.2f}, "
          f"p = {analysis['dB_p_value']:.6f}, "
          f"Cohen's d = {analysis['dB_cohen_d']:.2f}")
    sig = ("***" if analysis['dB_p_value'] < 0.001
           else "**" if analysis['dB_p_value'] < 0.01
           else "*" if analysis['dB_p_value'] < 0.05 else "")
    print(f"    Significant (p < 0.005): "
          f"{'YES ' + sig if analysis['dB_significant'] else 'NO'}")
    print()
    print(f"  Phase 3 Brain Size:")
    print(f"    Control:   {analysis['phase3_ctrl_B']['mean']:.2f} "
          f"+/- {analysis['phase3_ctrl_B']['se']:.2f}")
    print(f"    Treatment: {analysis['phase3_treat_B']['mean']:.2f} "
          f"+/- {analysis['phase3_treat_B']['se']:.2f}")
    print()
    print(f"  Treatment expanded: "
          f"p = {analysis['treat_expanded_p']:.2e} "
          f"{'[YES]' if analysis['treat_expanded'] else '[NO]'}")
    print()

    print(f"  Exploratory - Childhood Duration (Phase 3):")
    print(f"    Control:   {analysis['phase3_ctrl_L']['mean']:.2f} "
          f"+/- {analysis['phase3_ctrl_L']['se']:.2f}")
    print(f"    Treatment: {analysis['phase3_treat_L']['mean']:.2f} "
          f"+/- {analysis['phase3_treat_L']['se']:.2f}")
    print(f"    t = {analysis['p3L_t_statistic']:.2f}, "
          f"p = {analysis['p3L_p_value']:.6f}")
    print()
    print(f"  Exploratory - Knowledge (Phase 3):")
    print(f"    Control:   {analysis['phase3_ctrl_K']['mean']:.2f} "
          f"+/- {analysis['phase3_ctrl_K']['se']:.2f}")
    print(f"    Treatment: {analysis['phase3_treat_K']['mean']:.2f} "
          f"+/- {analysis['phase3_treat_K']['se']:.2f}")
    print(f"    t = {analysis['p3K_t_statistic']:.2f}, "
          f"p = {analysis['p3K_p_value']:.6f}")

    print()
    print("=" * 70)
    if analysis['grey_ceiling_hypothesis_supported']:
        print(f" {condition_name.upper()}: *** HYPOTHESIS SUPPORTED ***")
        print()
        print(" Within this parameterization, the care condition produces")
        print(" brain expansion through two specified model components:")
        print("   (1) a model-defined learning-storage gap")
        print("   (2) an imposed care-gated energetic ceiling lift")
        print()
        print(" This is a within-model counterfactual result, not evidence")
        print(" that the same causal sequence occurred in a fossil population.")
    else:
        print(f" {condition_name.upper()}: HYPOTHESIS NOT SUPPORTED")
    print("=" * 70)

    # Save JSON
    output_json = {
        'model': 'Care-Gated Grey Ceiling (Simulation 2)',
        'version': '2.0',
        'condition': condition_name,
        'reference': 'Grossmann, T. (tg3ny@virginia.edu)',
        'relationship_to_prior_work': {
            'experiment_1': 'Grossmann (2026) Childhood-First Simulation 1',
            'CBH': 'Muthukrishna et al. (2018) PLoS Comput. Biol.',
            'grey_ceiling': 'Isler & van Schaik (2012) Curr. Anthropol.',
            'cooperative_breeding': 'Hrdy (2009) Mothers and Others',
            'note': ('Phases 1-2 reproduce Simulation 1 exactly. '
                     'Phase 3 imposes care-gated and, where specified, '
                     'care-independent energetic ceiling lifts.'),
        },
        'parameters': asdict(params),
        'effective_ceilings': {
            'control': ctrl_bmax,
            'treatment': treat_bmax,
        },
        'n_replicates_per_group': n_reps,
        'master_seed': master_seed,
        'runtime_minutes': elapsed,
        'analysis': analysis,
    }

    with open(f"{cond_dir}/experiment2_results.json", 'w') as f:
        json.dump(output_json, f, indent=2)

    print()
    print(f" Results saved to: {cond_dir}/")
    print()

    return analysis


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Care-Gated Grey Ceiling: Brain Expansion via '
                    'Cooperative Care',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python grey_ceiling_simulation.py -n 100 --mode pure -o ./results
  python grey_ceiling_simulation.py -n 100 --mode mixed -o ./results
  python grey_ceiling_simulation.py -n 100 --mode both -o ./results
  python grey_ceiling_simulation.py -n 100 --ceiling-ramp 50000
        """
    )

    parser.add_argument('-n', '--n-replicates', type=int, default=100)
    parser.add_argument('-o', '--output', type=str, default='./results_exp2')
    parser.add_argument('-j', '--n-jobs', type=int, default=-1)
    parser.add_argument('--mode', type=str, default='both',
                        choices=['pure', 'mixed', 'both'],
                        help='pure: ceiling fully care-gated; '
                             'mixed: small exogenous + care-gated; '
                             'both: run both conditions (default)')
    parser.add_argument('--care-lift', type=float, default=None,
                        help='Care-dependent ceiling lift (overrides mode)')
    parser.add_argument('--exogenous-lift', type=float, default=None,
                        help='Exogenous ceiling lift (overrides mode)')
    parser.add_argument('--ceiling-ramp', type=int, default=0,
                        help='Generations for gradual ceiling lift; '
                             '0 = instantaneous (default: 0)')
    parser.add_argument('--T-phase3', type=int, default=200_000)
    parser.add_argument('--seed', type=int, default=42)

    args = parser.parse_args()

    # Print header
    print("=" * 70)
    print(" CARE-GATED GREY-CEILING SIMULATIONS")
    print(" Conditional routes to brain-size expansion")
    print("=" * 70)
    print()
    print(" Testing consequences of a model-defined learning-storage gap")
    print(" and imposed care-gated or care-independent ceiling lifts.")
    print()
    print(" Simulations 2 and 3, following:")
    print("   Grossmann (2026) Childhood-First Simulation 1")
    print()
    print(" Implementing the mechanism proposed by:")
    print("   Hrdy (2009) - cooperative breeding and food sharing")
    print("   Isler & van Schaik (2012) - breaking through the grey ceiling")
    print()
    print(" Building on:")
    print("   Muthukrishna et al. (2018) Cultural Brain Hypothesis model")
    print()
    print(f" Author: Tobias Grossmann (tg3ny@virginia.edu)")
    print()

    os.makedirs(args.output, exist_ok=True)

    # Determine which conditions to run
    conditions = []

    if args.care_lift is not None or args.exogenous_lift is not None:
        # Custom parameters override mode
        p = SimulationParameters(
            care_lift=args.care_lift if args.care_lift is not None else 100.0,
            exogenous_lift=(args.exogenous_lift
                            if args.exogenous_lift is not None else 0.0),
            ceiling_ramp=args.ceiling_ramp,
            T_phase3=args.T_phase3,
        )
        conditions.append(('custom', p))
    else:
        if args.mode in ('pure', 'both'):
            conditions.append(('pure', SimulationParameters(
                care_lift=100.0,
                exogenous_lift=0.0,
                ceiling_ramp=args.ceiling_ramp,
                T_phase3=args.T_phase3,
            )))
        if args.mode in ('mixed', 'both'):
            conditions.append(('mixed', SimulationParameters(
                care_lift=80.0,
                exogenous_lift=20.0,
                ceiling_ramp=args.ceiling_ramp,
                T_phase3=args.T_phase3,
            )))

    # Run conditions
    all_results = {}
    for cond_name, params in conditions:
        analysis = run_condition(
            condition_name=cond_name,
            params=params,
            n_reps=args.n_replicates,
            master_seed=args.seed,
            output_dir=args.output,
            n_jobs=args.n_jobs,
        )
        all_results[cond_name] = analysis

    # Summary across conditions
    if len(all_results) > 1:
        print()
        print("=" * 70)
        print(" SUMMARY ACROSS CONDITIONS")
        print("=" * 70)
        for name, a in all_results.items():
            supported = "SUPPORTED" if a['grey_ceiling_hypothesis_supported'] else "NOT SUPPORTED"
            print(f"  {name:8s}: delta_B treat={a['treat_delta_B']['mean']:+.1f}, "
                  f"ctrl={a['ctrl_delta_B']['mean']:+.1f}, "
                  f"diff={a['differential_expansion']:+.1f}, "
                  f"p={a['dB_p_value']:.4f} -> {supported}")
        print("=" * 70)


if __name__ == "__main__":
    main()
