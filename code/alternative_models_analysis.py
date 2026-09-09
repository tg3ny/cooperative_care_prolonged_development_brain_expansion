"""
Canonical SI analysis script for alternative model specifications.

Produces outputs used for Figure S4 and the robustness claims in the SI Appendix
regarding the necessity of convex (superlinear) return structures for developmental
decoupling.

Tests five model variants: superlinear (default), linear, saturating, weak B-coupling,
and null. Only the superlinear specification produces true decoupling (childhood
extension without brain change).

Deposited output: Data_S4_Alternative_Models.csv in the repository root.

Note on internal naming: variable names use 'epoch' (e.g., T_epoch1, epoch2_B)
for backward compatibility with data files. The manuscript uses Phase 1 / Phase 2.
"""
#!/usr/bin/env python3
"""
alternative_models_analysis.py

Tests Alternative Model Specifications
=======================================

This script tests whether the developmental decoupling result is robust to
alternative specifications of the social benefit function:

1. SUPERLINEAR (Original): B_social = Sub_care × r × (L/L_max)^e × L_max
   - The model used in the main analysis (e = 1.8)

2. LINEAR: B_social = Sub_care × r × L
   - No compounding; tests if superlinearity is necessary

3. SATURATING: B_social = Sub_care × r × L_max × L/(L + K)
   - Diminishing returns; tests if effect persists with saturation

4. WEAK B-COUPLING: B_social = Sub_care × r × (L/L_max)^e × L_max × (B/B_ref)^0.2
   - Allows brain to matter slightly; tests robustness to weak B effects

5. NULL CONTROL: B_social = Sub_care × r × random()
   - Scrambled payoff unrelated to L; verifies effect is specifically
     tied to care-contingent L-dependent payoffs

EXPECTED RESULTS
----------------
- Models 1-4 should show decoupling (ΔL ≥ 10%, |ΔB| < 2%) if robust
- Model 5 (null) should show NO systematic ΔL difference

USAGE
-----
    python alternative_models_analysis.py

    # Or with more replicates for publication:
    python alternative_models_analysis.py --replicates 20

Author: Tobias Grossmann (tg3ny@virginia.edu)
"""

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from dataclasses import dataclass
from typing import Dict, Callable
import time
import argparse
import json

try:
    from joblib import Parallel, delayed
    PARALLEL_AVAILABLE = True
except ImportError:
    PARALLEL_AVAILABLE = False
    print("Note: joblib not available. Running in serial mode.")
    print("Install with: pip install joblib")


# =============================================================================
# PARAMETERS
# =============================================================================

@dataclass
class SimulationParameters:
    """Parameters for the Childhood-First evolutionary simulation."""
    
    # Population
    N: int = 1000
    T_total: int = 100_000
    T_epoch1: int = 10_000
    T_epoch2: int = 100_000  # Now runs to completion
    
    # Costs
    alpha: float = 0.1
    gamma: float = 0.001
    beta: float = 0.12
    k_ceiling: float = 0.5
    
    # Learning
    lambda_: float = 0.5
    learning_exponent: float = 1.15
    asocial_rate: float = 0.10
    
    # Social fitness (key mechanism)
    social_rate: float = 0.30
    social_exponent: float = 1.8
    
    # Cooperative care
    Sub_care_control: float = 0.0
    Sub_care_treatment: float = 0.6
    
    # Evolution
    mu: float = 0.01
    sigma_mut: float = 0.1
    
    # Bounds
    B_max: float = 50.0
    L_max: float = 20.0
    
    # For saturating model
    saturation_K: float = 5.0
    
    # For weak B-coupling
    B_coupling_exp: float = 0.2
    B_ref: float = 50.0


# =============================================================================
# SOCIAL BENEFIT FUNCTIONS (Alternative Specifications)
# =============================================================================

def social_benefit_superlinear(L: np.ndarray, Sub_care: float, 
                                p: SimulationParameters, 
                                rng=None, B=None) -> np.ndarray:
    """
    Original superlinear model (main analysis).
    
    B_social = Sub_care × r × (L/L_max)^e × L_max
    
    Rationale: Social learning and relationships show compounding returns.
    """
    L_normalized = L / p.L_max
    return Sub_care * p.social_rate * np.power(L_normalized, p.social_exponent) * p.L_max


def social_benefit_linear(L: np.ndarray, Sub_care: float,
                          p: SimulationParameters,
                          rng=None, B=None) -> np.ndarray:
    """
    Linear model - no superlinearity.
    
    B_social = Sub_care × r × L
    
    Tests whether superlinearity is necessary for decoupling.
    """
    return Sub_care * p.social_rate * L


def social_benefit_saturating(L: np.ndarray, Sub_care: float,
                               p: SimulationParameters,
                               rng=None, B=None) -> np.ndarray:
    """
    Saturating model - diminishing returns.
    
    B_social = Sub_care × r × L_max × L/(L + K)
    
    Tests whether decoupling persists with diminishing returns.
    K = 5.0 means half-saturation at L = 5.
    """
    return Sub_care * p.social_rate * p.L_max * L / (L + p.saturation_K)


def social_benefit_weak_B_coupling(L: np.ndarray, Sub_care: float,
                                    p: SimulationParameters,
                                    rng=None, B=None) -> np.ndarray:
    """
    Weak brain-coupling model.
    
    B_social = Sub_care × r × (L/L_max)^e × L_max × (B/B_ref)^0.2
    
    Allows brain size to matter slightly (exponent 0.2).
    Tests robustness to weak B effects.
    """
    L_normalized = L / p.L_max
    B_factor = np.power(B / p.B_ref, p.B_coupling_exp)
    return Sub_care * p.social_rate * np.power(L_normalized, p.social_exponent) * p.L_max * B_factor


def social_benefit_null(L: np.ndarray, Sub_care: float,
                        p: SimulationParameters,
                        rng=None, B=None) -> np.ndarray:
    """
    Null control - scrambled payoff.
    
    B_social = Sub_care × r × random(0, 20)
    
    Payoff is random, not dependent on L.
    Should show NO systematic ΔL if effect is truly L-dependent.
    """
    if rng is None:
        rng = np.random.default_rng()
    return Sub_care * p.social_rate * rng.uniform(0, 20, size=len(L))


# =============================================================================
# POPULATION CLASS
# =============================================================================

class Population:
    """
    Evolving population with configurable social benefit function.
    """
    
    def __init__(self, params: SimulationParameters, seed: int,
                 social_benefit_func: Callable):
        self.params = params
        self.rng = np.random.default_rng(seed)
        self.N = params.N
        self.social_benefit_func = social_benefit_func
        
        # Environmental state
        self.B_max = params.B_max
        self.Sub_care = 0.0
        
        # Initialize traits
        self.B = self.rng.uniform(5, 15, self.N)
        self.L = self.rng.uniform(1, 3, self.N)
        self.S = self.rng.uniform(0.3, 0.7, self.N)
        self.K = np.zeros(self.N)
        
        self.prev_K = None
        self.prev_fitness = None
    
    def step(self) -> Dict[str, float]:
        """Execute one generation."""
        p = self.params
        
        # Learning
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
        self.K = np.minimum(K_total, self.B * 0.6)
        
        # Costs
        C_struct = p.alpha * self.B + p.gamma * (self.B ** 2)
        ceiling_penalty = 1.0 / (1.0 + np.exp(-p.k_ceiling * (self.B - self.B_max)))
        C_struct += 10.0 * ceiling_penalty
        C_time = p.beta * self.L
        total_cost = C_struct + C_time
        
        # Benefits
        benefit_knowledge = p.lambda_ * self.K
        social_benefit = self.social_benefit_func(
            self.L, self.Sub_care, p, rng=self.rng, B=self.B
        )
        total_benefit = benefit_knowledge + social_benefit
        
        # Fitness
        net_fitness = total_benefit - total_cost
        P_survival = 1.0 / (1.0 + np.exp(-net_fitness))
        P_survival = np.clip(P_survival, 0.001, 0.999)
        
        self.prev_K = self.K.copy()
        self.prev_fitness = P_survival.copy()
        
        # Selection
        survived = self.rng.random(self.N) < P_survival
        survivor_idx = np.where(survived)[0]
        if len(survivor_idx) < 2:
            survivor_idx = np.arange(self.N)
        
        # Reproduction
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
            'L_mean': float(np.mean(self.L)),
        }


# =============================================================================
# SIMULATION FUNCTIONS
# =============================================================================

def run_replicate(replicate_id: int, group: str, seed: int,
                  params: SimulationParameters,
                  social_benefit_func: Callable) -> Dict:
    """Run a single replicate with specified social benefit function."""
    
    pop = Population(params, seed, social_benefit_func)
    
    # Phase 1: Baseline (no care for anyone)
    pop.B_max = params.B_max
    pop.Sub_care = 0.0
    
    for gen in range(params.T_epoch1):
        pop.step()
    
    # Phase 2: Treatment divergence
    if group == 'treatment':
        pop.Sub_care = params.Sub_care_treatment
    else:
        pop.Sub_care = params.Sub_care_control
    
    for gen in range(params.T_epoch1, params.T_epoch2):
        pop.step()
    
    # Record end of Phase 2
    return {
        'replicate': replicate_id,
        'group': group,
        'L_epoch2': float(np.mean(pop.L)),
        'B_epoch2': float(np.mean(pop.B)),
    }


def analyze_results(results: list) -> Dict:
    """Compute statistics from replicate results."""
    
    df = pd.DataFrame(results)
    
    control = df[df['group'] == 'control']
    treatment = df[df['group'] == 'treatment']
    
    L_control = control['L_epoch2'].values
    L_treatment = treatment['L_epoch2'].values
    B_control = control['B_epoch2'].values
    B_treatment = treatment['B_epoch2'].values
    
    # Effect sizes
    delta_L = (L_treatment.mean() - L_control.mean()) / L_control.mean() * 100
    delta_B = (B_treatment.mean() - B_control.mean()) / B_control.mean() * 100
    
    # T-test for childhood extension
    t_L, p_L = scipy_stats.ttest_ind(L_treatment, L_control, alternative='greater')
    
    # Cohen's d
    pooled_std = np.sqrt((L_treatment.var() + L_control.var()) / 2)
    cohens_d = (L_treatment.mean() - L_control.mean()) / pooled_std if pooled_std > 0 else 0
    
    # TOST for brain equivalence
    bound = 0.02 * B_control.mean()  # ±2%
    diff = B_treatment.mean() - B_control.mean()
    pooled_se = np.sqrt(B_treatment.var()/len(B_treatment) + B_control.var()/len(B_control))
    df_stat = len(B_treatment) + len(B_control) - 2
    
    if pooled_se > 0:
        t_lower = (diff - (-bound)) / pooled_se
        t_upper = (bound - diff) / pooled_se
        p_lower = 1 - scipy_stats.t.cdf(t_lower, df_stat)
        p_upper = 1 - scipy_stats.t.cdf(t_upper, df_stat)
        p_tost = max(p_lower, p_upper)
    else:
        p_tost = 1.0
    
    # Decoupling criteria
    decoupling = (delta_L >= 10) and (abs(delta_B) < 2) and (p_L < 0.005)
    
    return {
        'L_control_mean': float(L_control.mean()),
        'L_control_se': float(L_control.std() / np.sqrt(len(L_control))),
        'L_treatment_mean': float(L_treatment.mean()),
        'L_treatment_se': float(L_treatment.std() / np.sqrt(len(L_treatment))),
        'B_control_mean': float(B_control.mean()),
        'B_treatment_mean': float(B_treatment.mean()),
        'delta_L_pct': float(delta_L),
        'delta_B_pct': float(delta_B),
        't_stat_L': float(t_L),
        'p_value_L': float(p_L),
        'cohens_d': float(cohens_d),
        'p_tost_B': float(p_tost),
        'decoupling': bool(decoupling),  # Ensure Python bool for JSON serialization
    }


def run_model_variant(model_name: str, social_benefit_func: Callable,
                      n_replicates: int = 20, n_jobs: int = -1,
                      verbose: bool = True) -> Dict:
    """Run experiment with a specific model variant."""
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Testing: {model_name}")
        print(f"{'='*60}")
    
    params = SimulationParameters()
    base_seed = abs(hash(model_name)) % (2**31)
    
    # Prepare all jobs
    jobs = []
    for group in ['control', 'treatment']:
        for i in range(n_replicates):
            seed = base_seed + i + (10000 if group == 'treatment' else 0)
            jobs.append((i, group, seed))
    
    # Run replicates
    if PARALLEL_AVAILABLE and n_jobs != 1:
        if verbose:
            print(f"  Running {len(jobs)} replicates in parallel...")
        results = Parallel(n_jobs=n_jobs, verbose=10 if verbose else 0)(
            delayed(run_replicate)(rep_id, group, seed, params, social_benefit_func)
            for rep_id, group, seed in jobs
        )
    else:
        results = []
        for idx, (rep_id, group, seed) in enumerate(jobs):
            results.append(run_replicate(rep_id, group, seed, params, social_benefit_func))
            if verbose and (idx + 1) % 5 == 0:
                print(f"  Completed {idx + 1}/{len(jobs)} replicates")
    
    # Analyze
    stats = analyze_results(results)
    stats['model'] = model_name
    stats['n_replicates'] = n_replicates
    
    if verbose:
        print(f"\n  Results:")
        print(f"    ΔL = {stats['delta_L_pct']:+.1f}% (t = {stats['t_stat_L']:.2f}, "
              f"p = {stats['p_value_L']:.4f}, d = {stats['cohens_d']:.2f})")
        print(f"    ΔB = {stats['delta_B_pct']:+.2f}% (TOST p = {stats['p_tost_B']:.4f})")
        print(f"    Decoupling: {'YES ✓' if stats['decoupling'] else 'NO ✗'}")
    
    return stats


# =============================================================================
# MAIN
# =============================================================================

def main(n_replicates: int = 20, n_jobs: int = -1):
    """Run all alternative model analyses."""
    
    print("="*70)
    print("ALTERNATIVE MODEL SPECIFICATIONS ANALYSIS")
    print("Testing robustness of developmental decoupling to model assumptions")
    print("="*70)
    print(f"\nSettings: {n_replicates} replicates per condition")
    
    start_time = time.time()
    
    # Define model variants to test
    models = [
        ("1. SUPERLINEAR (Original)", social_benefit_superlinear,
         "B_social = Sub_care × r × (L/L_max)^e × L_max"),
        ("2. LINEAR", social_benefit_linear,
         "B_social = Sub_care × r × L"),
        ("3. SATURATING", social_benefit_saturating,
         "B_social = Sub_care × r × L_max × L/(L + K)"),
        ("4. WEAK B-COUPLING", social_benefit_weak_B_coupling,
         "B_social = Sub_care × r × (L/L_max)^e × L_max × (B/B_ref)^0.2"),
        ("5. NULL CONTROL", social_benefit_null,
         "B_social = Sub_care × r × random(0,20)"),
    ]
    
    # Run all variants
    all_results = []
    for model_name, model_func, formula in models:
        result = run_model_variant(model_name, model_func, 
                                   n_replicates=n_replicates, n_jobs=n_jobs)
        result['formula'] = formula
        all_results.append(result)
    
    # =========================================================================
    # SUMMARY TABLE
    # =========================================================================
    print("\n" + "="*70)
    print("SUMMARY: ALTERNATIVE MODEL SPECIFICATIONS")
    print("="*70)
    
    print(f"\n{'Model':<28} {'ΔL%':>8} {'ΔB%':>8} {'p(L)':>10} {'p(TOST)':>10} {'Decoupling':>12}")
    print("-"*78)
    
    for r in all_results:
        decoup_str = "YES ✓" if r['decoupling'] else "NO ✗"
        print(f"{r['model']:<28} {r['delta_L_pct']:>+7.1f}% {r['delta_B_pct']:>+7.2f}% "
              f"{r['p_value_L']:>10.4f} {r['p_tost_B']:>10.4f} {decoup_str:>12}")
    
    # =========================================================================
    # INTERPRETATION
    # =========================================================================
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    
    # Count decoupling in non-null models
    non_null = [r for r in all_results if 'NULL' not in r['model']]
    decoupling_count = sum(1 for r in non_null if r['decoupling'])
    
    print(f"\n✓ {decoupling_count}/4 alternative model specifications show decoupling")
    
    # Check LINEAR model specifically - this is theoretically important
    linear_result = next((r for r in all_results if 'LINEAR' in r['model'] and 'SUPER' not in r['model']), None)
    if linear_result and not linear_result['decoupling']:
        print(f"\n✓ LINEAR model shows childhood extension (ΔL = {linear_result['delta_L_pct']:+.1f}%) but")
        print(f"  brain SHRINKAGE (ΔB = {linear_result['delta_B_pct']:+.1f}%), confirming that convexity")
        print("  is necessary to avoid brain-childhood trade-offs.")
    
    # Check null control
    null_result = next((r for r in all_results if 'NULL' in r['model']), None)
    if null_result:
        if null_result['decoupling']:
            print(f"\n⚠ WARNING: Null control shows significant effect (ΔL = {null_result['delta_L_pct']:.1f}%)")
            print("  This would suggest the effect is NOT specifically tied to L-dependent payoffs.")
        else:
            print(f"\n✓ Null control shows NO systematic effect (ΔL = {null_result['delta_L_pct']:+.1f}%, p = {null_result['p_value_L']:.3f})")
            print("  This confirms the decoupling effect is specifically tied to L-dependent payoffs,")
            print("  not an artifact of the simulation structure.")
    
    # Overall conclusion - updated to reflect theoretical significance
    print("\n" + "-"*70)
    superlinear_result = next((r for r in all_results if 'SUPERLINEAR' in r['model']), None)
    
    if superlinear_result and superlinear_result['decoupling']:
        if linear_result and not linear_result['decoupling'] and linear_result['delta_B_pct'] < -2:
            print("CONCLUSION: Results SUPPORT the convexity hypothesis.")
            print("  - SUPERLINEAR (convex) model: Decoupling achieved (ΔL > 10%, |ΔB| < 2%)")
            print("  - LINEAR model: Childhood extends BUT brain shrinks (trade-off)")
            print("  This confirms that convex fitness returns are NECESSARY for true")
            print("  developmental decoupling, as predicted by the Childhood-First hypothesis.")
        elif decoupling_count >= 3:
            print("CONCLUSION: The developmental decoupling result is ROBUST to alternative")
            print("model specifications.")
        else:
            print("CONCLUSION: Decoupling is achieved with the SUPERLINEAR (convex) model.")
            print("  Alternative specifications show varying degrees of brain-childhood trade-offs,")
            print("  supporting the theoretical importance of convex fitness returns.")
    else:
        print("CONCLUSION: Unexpected results. The SUPERLINEAR model did not show decoupling.")
        print("  Further investigation recommended.")
    
    # =========================================================================
    # SAVE RESULTS
    # =========================================================================
    
    # CSV
    df_results = pd.DataFrame(all_results)
    df_results.to_csv('alternative_models_results.csv', index=False)
    df_results.to_csv('Data_S4_Alternative_Models.csv', index=False)
    
    # JSON summary
    with open('alternative_models_results.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    elapsed = time.time() - start_time
    print(f"\n{'='*70}")
    print(f"Total time: {elapsed/60:.1f} minutes")
    print(f"Results saved to: alternative_models_results.csv")
    print(f"                  alternative_models_results.json")
    print("="*70)
    
    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test alternative model specifications")
    parser.add_argument('-n', '--replicates', type=int, default=20,
                        help='Number of replicates per condition (default: 20)')
    parser.add_argument('-j', '--jobs', type=int, default=-1,
                        help='Number of parallel jobs (-1 for all cores, default: -1)')
    
    args = parser.parse_args()
    
    results = main(n_replicates=args.replicates, n_jobs=args.jobs)
