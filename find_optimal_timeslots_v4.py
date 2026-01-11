"""
find_optimal_timeslots_v4.py - Grid Search for Optimal Timeslots
==============================================================
Searches over timeslot combinations to find the best for each policy.

Usage:
    python find_optimal_timeslots.py

Output:
    - Prints optimal timeslots for dedicated and pooled policies
    - Saves grid search results to CSV files
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

from simulation_core_v4 import (
    MRISimulationDES,
    TYPE1_TIMESLOT_RANGE, TYPE2_TIMESLOT_RANGE,
    DEFAULT_SIMULATION_DAYS, DEFAULT_NUM_REPLICATIONS,
    TYPE1_DURATION_MEAN, TYPE1_DURATION_STD,
    TYPE2_DURATION_MEAN, TYPE2_DURATION_STD,
    TYPE1_ARRIVALS_PER_DAY, TYPE2_ARRIVALS_PER_DAY
)


def run_monte_carlo(timeslot_type1: int, timeslot_type2: int, policy: str,
                    num_replications: int = 30, num_days: int = DEFAULT_SIMULATION_DAYS) -> Dict:
    """
    Run Monte Carlo simulation with multiple replications.
    
    Returns:
        Dictionary with mean and std of each metric
    """
    results = []
    
    for rep in range(num_replications):
        sim = MRISimulationDES(
            timeslot_type1=timeslot_type1,
            timeslot_type2=timeslot_type2,
            policy=policy,
            num_days=num_days,
            seed=42 + rep
        )
        sim.run()
        results.append(sim.calculate_statistics())
    
    # Aggregate
    aggregated = {}
    for metric in results[0].keys():
        values = [r[metric] for r in results]
        aggregated[f'{metric}_mean'] = np.mean(values)
        aggregated[f'{metric}_std'] = np.std(values)
    
    return aggregated


def calculate_score(results: Dict) -> float:
    """
    Calculate objective score (lower is better).
    
    Balances wait time, overtime, and utilization.
    """
    wait_score = results['avg_wait_days_mean'] * 10                              # Heavily weight wait time
    overtime_score = results['avg_overtime_minutes_mean'] * 0.1                  # Penalize overtime
    util_penalty = (1 - results['avg_utilization_mean']) * 5                     # Penalize low utilization
    
    return wait_score + overtime_score + util_penalty


def grid_search(policy: str,
                type1_range: range = TYPE1_TIMESLOT_RANGE,
                type2_range: range = TYPE2_TIMESLOT_RANGE,
                num_replications: int = 30) -> pd.DataFrame:
    """
    Perform grid search over all timeslot combinations.
    
    Args:
        policy: 'dedicated' or 'pooled'
        type1_range: Range of Type 1 timeslots to test
        type2_range: Range of Type 2 timeslots to test
        num_replications: Replications per combination
    
    Returns:
        DataFrame with results for each combination
    """
    results = []
    total = len(type1_range) * len(type2_range)
    count = 0
    
    print(f"\nGrid search for {policy.upper()} policy:")
    print(f"  Type 1 range: {list(type1_range)} minutes")
    print(f"  Type 2 range: {list(type2_range)} minutes")
    print(f"  Total combinations: {total}")
    print(f"  Replications per combination: {num_replications}")
    print()
    
    for ts1 in type1_range:
        for ts2 in type2_range:
            count += 1
            print(f"  [{count:2d}/{total}] Type1={ts1}min, Type2={ts2}min...", end='\r')
            
            mc_results = run_monte_carlo(ts1, ts2, policy, num_replications)
            score = calculate_score(mc_results)
            
            results.append({
                'timeslot_type1': ts1,
                'timeslot_type2': ts2,
                'policy': policy,
                'score': score,
                'avg_wait_days': mc_results['avg_wait_days_mean'],
                'avg_wait_std': mc_results['avg_wait_days_std'],
                'avg_overtime_min': mc_results['avg_overtime_minutes_mean'],
                'overtime_pct': mc_results['proportion_overtime_days_mean'] * 100,
                'avg_utilization': mc_results['avg_utilization_mean'],
                'avg_daily_patients': mc_results['avg_daily_patients_mean']
            })
    
    print()
    return pd.DataFrame(results)


def find_optimal(policy: str, num_replications: int = 30) -> Tuple[int, int, pd.DataFrame]:
    """
    Find optimal timeslots for a policy.
    
    Returns:
        (optimal_ts1, optimal_ts2, full_grid_results)
    """
    grid_df = grid_search(policy, num_replications=num_replications)
    
    # Find best (lowest score)
    best_idx = grid_df['score'].idxmin()
    best = grid_df.loc[best_idx]
    
    ts1, ts2 = int(best['timeslot_type1']), int(best['timeslot_type2'])
    
    print(f"\n  OPTIMAL for {policy.upper()}:")
    print(f"    Type 1 timeslot: {ts1} min")
    print(f"    Type 2 timeslot: {ts2} min")
    print(f"    Wait time:       {best['avg_wait_days']:.2f} days (±{best['avg_wait_std']:.2f})")
    print(f"    Overtime:        {best['overtime_pct']:.1f}% of days")
    print(f"    Utilization:     {best['avg_utilization']*100:.1f}%")
    print(f"    Score:           {best['score']:.2f}")
    
    return ts1, ts2, grid_df


def main():
    """Main function: find optimal timeslots for both policies."""
    
    print("=" * 70)
    print("MRI SCHEDULING - OPTIMAL TIMESLOT SEARCH")
    print("=" * 70)
    print(f"\nParameters from Part 1 Bootstrap Analysis:")
    print(f"  Type 1: Mean={TYPE1_DURATION_MEAN:.2f}min, SD={TYPE1_DURATION_STD:.2f}min, "
          f"Arrivals={TYPE1_ARRIVALS_PER_DAY:.1f}/day")
    print(f"  Type 2: Mean={TYPE2_DURATION_MEAN:.2f}min, SD={TYPE2_DURATION_STD:.2f}min, "
          f"Arrivals={TYPE2_ARRIVALS_PER_DAY:.1f}/day")
    
    # Find optimal for dedicated policy
    print("\n" + "-" * 70)
    ded_ts1, ded_ts2, ded_grid = find_optimal('dedicated', num_replications=30)
    
    # Find optimal for pooled policy
    print("\n" + "-" * 70)
    pool_ts1, pool_ts2, pool_grid = find_optimal('pooled', num_replications=30)
    
    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY OF OPTIMAL TIMESLOTS")
    print("=" * 70)
    print(f"\n  DEDICATED (Old System): Type1={ded_ts1}min, Type2={ded_ts2}min")
    print(f"  POOLED (New System):    Type1={pool_ts1}min, Type2={pool_ts2}min")
    
    # Save results
    ded_grid.to_csv('grid_search_dedicated.csv', index=False)
    pool_grid.to_csv('grid_search_pooled.csv', index=False)
    
    # Save optimal timeslots to file for compare_policies.py to read
    optimal_df = pd.DataFrame([
        {'policy': 'dedicated', 'timeslot_type1': ded_ts1, 'timeslot_type2': ded_ts2},
        {'policy': 'pooled', 'timeslot_type1': pool_ts1, 'timeslot_type2': pool_ts2}
    ])
    optimal_df.to_csv('optimal_timeslots.csv', index=False)
    
    print(f"\nResults saved:")
    print(f"  - grid_search_dedicated.csv")
    print(f"  - grid_search_pooled.csv")
    print(f"  - optimal_timeslots.csv")
    print("=" * 70)
    
    return {
        'dedicated': {'ts1': ded_ts1, 'ts2': ded_ts2, 'grid': ded_grid},
        'pooled': {'ts1': pool_ts1, 'ts2': pool_ts2, 'grid': pool_grid}
    }


if __name__ == "__main__":
    results = main()
