"""
Auxiliary convergence diagnostic script.

Reports the dispersion of replicate means near the end of each phase of Simulation 1 and the change in the
mean trajectories over the final 2,000 generations; output: Data_S5_Equilibrium_Statistics.csv.
Not a source of the statistics reported in the article.

The script does not establish that populations are at equilibrium. As the Supporting Information states
(text S1(c)), trajectories are still changing slowly at the end of each phase, endpoint values characterize a
fixed computational horizon, and generations are computational iterations, not a calibrated chronology.
The printed summary was corrected accordingly in release v1.1.1; the computed statistics are unchanged.

Note on internal naming: variable names use 'epoch' for backward
compatibility with data files. The manuscript uses Phase 1 / Phase 2.
"""
import pandas as pd
import numpy as np
from pathlib import Path

def load_trajectory_data(filepath: str) -> pd.DataFrame:
    """Load trajectory data from CSV."""
    return pd.read_csv(filepath)

def analyze_epoch_equilibrium(df: pd.DataFrame, group: str, 
                               gen_start: int, gen_end: int,
                               trait: str = 'L_mean') -> dict:
    """
    Analyze equilibrium statistics for a specific phase window.
    
    Parameters:
    -----------
    df : DataFrame
        Trajectory data
    group : str
        'treatment' or 'control'
    gen_start : int
        Start generation for analysis window
    gen_end : int
        End generation for analysis window
    trait : str
        Trait column to analyze (default: 'L_mean')
    
    Returns:
    --------
    dict : Statistics including mean, std, CV, and range
    """
    mask = (df['group'] == group) & \
           (df['generation'] >= gen_start) & \
           (df['generation'] <= gen_end)
    
    subset = df.loc[mask, trait]
    
    if len(subset) == 0:
        return None
    
    mean_val = subset.mean()
    std_val = subset.std()
    cv = (std_val / mean_val) * 100 if mean_val != 0 else np.nan
    
    return {
        'mean': mean_val,
        'std': std_val,
        'cv_percent': cv,
        'min': subset.min(),
        'max': subset.max(),
        'n_samples': len(subset)
    }

def calculate_rate_of_change(df: pd.DataFrame, group: str,
                              gen1: int, gen2: int,
                              trait: str = 'L_mean') -> dict:
    """
    Calculate rate of change between two generation points.
    
    Returns:
    --------
    dict : Values at each point and percent change
    """
    # Get mean across replicates at each generation
    mask1 = (df['group'] == group) & (df['generation'] == gen1)
    mask2 = (df['group'] == group) & (df['generation'] == gen2)
    
    val1 = df.loc[mask1, trait].mean()
    val2 = df.loc[mask2, trait].mean()
    
    pct_change = ((val2 - val1) / val1) * 100 if val1 != 0 else np.nan
    
    return {
        f'value_at_{gen1}': val1,
        f'value_at_{gen2}': val2,
        'absolute_change': val2 - val1,
        'percent_change': pct_change
    }

def main():
    """Run equilibrium convergence analysis."""
    
    # Find data file - searches common locations automatically
    data_paths = [
        # Public GitHub repository root layout
        Path('./Data_S3_Trajectories.csv.gz'),
        Path('./Data_S3_Trajectories.csv'),
        # Submission package naming convention
        Path('./data/Data_S3_Trajectories.csv.gz'),
        Path('./data/Data_S3_Trajectories.csv'),
        Path('./results/Data_S3_Trajectories.csv'),
        Path('./Data_S3_Trajectories.csv'),
        # Default output from childhood_first_simulation.py
        Path('./results/full_trajectories.csv'),
        Path('../results/full_trajectories.csv'),
        # Current directory
        Path('./full_trajectories.csv'),
        # Common alternative locations
        Path('./data/full_trajectories.csv'),
        Path('../data/full_trajectories.csv'),
        # Uploaded files location
        Path('data/full_trajectories.csv'),
    ]
    
    data_file = None
    for p in data_paths:
        if p.exists():
            data_file = p
            break
    
    if data_file is None:
        print("ERROR: Could not find trajectory data file")
        print("\nSearched locations:")
        for p in data_paths:
            print(f"  - {p}")
        print("\nPlease ensure you have run childhood_first_simulation.py first:")
        print("  python childhood_first_simulation.py -n 100 -o ./results -j -1")
        print("\nOr specify the path by editing data_paths in this script.")
        return
    
    print("=" * 70)
    print("EQUILIBRIUM CONVERGENCE ANALYSIS")
    print("Childhood-First Hypothesis Simulation")
    print("=" * 70)
    print(f"\nData source: {data_file}")
    
    # Load data
    df = load_trajectory_data(data_file)
    print(f"Loaded {len(df)} trajectory records")
    print(f"Groups: {df['group'].unique()}")
    print(f"Generation range: {df['generation'].min()} - {df['generation'].max()}")
    
    # Define phase boundaries (now only 2 phases)
    phases = {
        'phase 1 (Baseline)': {'start': 0, 'end': 10000, 'window_start': 9500},
        'phase 2 (Cooperative Care)': {'start': 10000, 'end': 100000, 'window_start': 99500},
    }
    
    print("\n" + "=" * 70)
    print("PART 1: EQUILIBRIUM AT EPOCH BOUNDARIES")
    print("=" * 70)
    print("\nAnalyzing final 500 generations of each epoch...")
    print("(Coefficient of variation of replicate means; larger values indicate greater dispersion among replicates)\n")
    
    for epoch_name, bounds in phases.items():
        print(f"\n{epoch_name}")
        print("-" * 50)
        
        for group in ['treatment', 'control']:
            stats = analyze_epoch_equilibrium(
                df, group, 
                bounds['window_start'], bounds['end'],
                'L_mean'
            )
            
            if stats:
                print(f"  {group.capitalize()} group (L):")
                print(f"    Mean: {stats['mean']:.3f} ± {stats['std']:.3f}")
                print(f"    CV: {stats['cv_percent']:.2f}%")
                print(f"    Range: [{stats['min']:.3f}, {stats['max']:.3f}]")
                print(f"    N samples: {stats['n_samples']}")
    
    print("\n" + "=" * 70)
    print("PART 2: RATE OF CHANGE IN FINAL GENERATIONS")
    print("=" * 70)
    print("\nMeasuring trait change over final 2,000 generations...")
    print("(Change in the mean trajectory; a nonzero change means that the trait is still moving at the horizon)\n")
    
    for trait, trait_name in [('L_mean', 'Childhood Duration (L)'), 
                               ('B_mean', 'Brain Size (B)')]:
        print(f"\n{trait_name}:")
        print("-" * 50)
        
        for group in ['treatment', 'control']:
            change = calculate_rate_of_change(
                df, group, 98000, 99900, trait
            )
            
            print(f"  {group.capitalize()} group:")
            print(f"    Value at gen 98,000: {change['value_at_98000']:.3f}")
            print(f"    Value at gen 99,900: {change['value_at_99900']:.3f}")
            print(f"    Absolute change: {change['absolute_change']:+.3f}")
            print(f"    Percent change: {change['percent_change']:+.2f}%")
    
    print("\n" + "=" * 70)
    print("PART 3: CROSS-EPOCH COMPARISON")
    print("=" * 70)
    print("\nComparing treatment group L values across phases...\n")
    
    epoch_means = {}
    for epoch_name, bounds in phases.items():
        stats = analyze_epoch_equilibrium(
            df, 'treatment',
            bounds['window_start'], bounds['end'],
            'L_mean'
        )
        if stats:
            epoch_means[epoch_name] = stats['mean']
            print(f"  {epoch_name}: L = {stats['mean']:.3f}")
    
    # Calculate changes between phases
    if len(epoch_means) == 2:
        epoch_list = list(epoch_means.keys())
        vals = list(epoch_means.values())
        
        print(f"\n  Change Phase 1 → Phase 2: {((vals[1]-vals[0])/vals[0])*100:+.1f}%")
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
These diagnostics describe the state of Simulation 1 at fixed computational
horizons. They do not show that populations are at equilibrium: replicate
means remain dispersed, and mean childhood duration in the treatment group
is still increasing at generation 100,000. Endpoint values should therefore
be read as horizon-dependent (Supporting Information, text S1(c)).
Generations are computational iterations and are not converted to years.
""")
    
    # Save summary statistics to CSV
    output_data = []
    for epoch_name, bounds in phases.items():
        for group in ['treatment', 'control']:
            for trait in ['L_mean', 'B_mean']:
                stats = analyze_epoch_equilibrium(
                    df, group,
                    bounds['window_start'], bounds['end'],
                    trait
                )
                if stats:
                    output_data.append({
                        'phase': epoch_name,
                        'group': group,
                        'trait': trait,
                        'mean': stats['mean'],
                        'std': stats['std'],
                        'cv_percent': stats['cv_percent'],
                        'n_samples': stats['n_samples']
                    })
    
    output_df = pd.DataFrame(output_data)
    output_file = Path('equilibrium_statistics.csv')
    output_df.to_csv(output_file, index=False)
    print(f"\nDetailed statistics saved to: {output_file}")

if __name__ == '__main__':
    main()
