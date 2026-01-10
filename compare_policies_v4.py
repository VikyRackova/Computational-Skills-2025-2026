"""
compare_policies.py - Compare Policies at Optimal Timeslots
============================================================
Runs detailed comparison of dedicated vs pooled policies using their
respective optimal timeslots (from find_optimal_timeslots.py).

Usage:
    python find_optimal_timeslots.py   # Run first to find optimal timeslots
    python compare_policies.py         # Then run this to compare
"""

import numpy as np
import pandas as pd
import os
from typing import Dict

from simulation_core import (
    MRISimulationDES,
    DEFAULT_SIMULATION_DAYS
)


def run_detailed_simulation(timeslot_type1: int, timeslot_type2: int, 
                            policy: str, num_replications: int = 100,
                            num_days: int = DEFAULT_SIMULATION_DAYS) -> Dict:
    """Run Monte Carlo simulation and return aggregated results."""
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
    
    aggregated = {}
    for metric in results[0].keys():
        values = [r[metric] for r in results]
        aggregated[f'{metric}_mean'] = np.mean(values)
        aggregated[f'{metric}_std'] = np.std(values)
    
    return aggregated


def load_optimal_timeslots() -> Dict:
    """Load optimal timeslots from CSV (created by find_optimal_timeslots.py)."""
    if os.path.exists('optimal_timeslots.csv'):
        df = pd.read_csv('optimal_timeslots.csv')
        result = {}
        for _, row in df.iterrows():
            result[row['policy']] = {
                'ts1': int(row['timeslot_type1']),
                'ts2': int(row['timeslot_type2'])
            }
        return result
    else:
        print("Warning: optimal_timeslots.csv not found. Using default values.")
        return {
            'dedicated': {'ts1': 28, 'ts2': 45},
            'pooled': {'ts1': 28, 'ts2': 45}
        }


def create_results_table(ded_results: Dict, pool_results: Dict,
                         ded_ts: Dict, pool_ts: Dict) -> pd.DataFrame:
    """Create summary results table."""
    rows = []
    
    rows.append({
        'Metric': 'Timeslot Type 1 (min)',
        'Dedicated': ded_ts['ts1'],
        'Pooled': pool_ts['ts1'],
        'Difference': '-'
    })
    rows.append({
        'Metric': 'Timeslot Type 2 (min)',
        'Dedicated': ded_ts['ts2'],
        'Pooled': pool_ts['ts2'],
        'Difference': '-'
    })
    
    ded_wait = ded_results['avg_wait_days_mean']
    pool_wait = pool_results['avg_wait_days_mean']
    diff = ((ded_wait - pool_wait) / ded_wait * 100) if ded_wait > 0 else 0
    rows.append({
        'Metric': 'Avg Wait Time (days)',
        'Dedicated': f"{ded_wait:.2f} (±{ded_results['avg_wait_days_std']:.2f})",
        'Pooled': f"{pool_wait:.2f} (±{pool_results['avg_wait_days_std']:.2f})",
        'Difference': f"{diff:+.1f}%"
    })
    
    ded_util = ded_results['avg_utilization_mean'] * 100
    pool_util = pool_results['avg_utilization_mean'] * 100
    diff = pool_util - ded_util
    rows.append({
        'Metric': 'Avg Utilization (%)',
        'Dedicated': f"{ded_util:.1f}",
        'Pooled': f"{pool_util:.1f}",
        'Difference': f"{diff:+.1f}pp"
    })
    
    ded_ot = ded_results['proportion_overtime_days_mean'] * 100
    pool_ot = pool_results['proportion_overtime_days_mean'] * 100
    diff = pool_ot - ded_ot
    rows.append({
        'Metric': 'Overtime Frequency (%)',
        'Dedicated': f"{ded_ot:.1f}",
        'Pooled': f"{pool_ot:.1f}",
        'Difference': f"{diff:+.1f}pp"
    })
    
    ded_tp = ded_results['avg_daily_patients_mean']
    pool_tp = pool_results['avg_daily_patients_mean']
    diff = ((pool_tp - ded_tp) / ded_tp * 100) if ded_tp > 0 else 0
    rows.append({
        'Metric': 'Avg Daily Patients',
        'Dedicated': f"{ded_tp:.1f}",
        'Pooled': f"{pool_tp:.1f}",
        'Difference': f"{diff:+.1f}%"
    })
    
    return pd.DataFrame(rows)


def main():
    """Main function: compare policies at optimal timeslots."""
    
    # Load optimal timeslots
    optimal = load_optimal_timeslots()
    ded_ts = optimal['dedicated']
    pool_ts = optimal['pooled']
    
    # Run simulations
    ded_results = run_detailed_simulation(
        ded_ts['ts1'], ded_ts['ts2'], 'dedicated', num_replications=100
    )
    pool_results = run_detailed_simulation(
        pool_ts['ts1'], pool_ts['ts2'], 'pooled', num_replications=100
    )
    
    # Create and save results table
    results_table = create_results_table(ded_results, pool_results, ded_ts, pool_ts)
    results_table.to_csv('comparison_results.csv', index=False)
    
    # Print summary
    print(f"Dedicated: Type1={ded_ts['ts1']}min, Type2={ded_ts['ts2']}min")
    print(f"  Wait: {ded_results['avg_wait_days_mean']:.2f}d, Util: {ded_results['avg_utilization_mean']*100:.1f}%")
    
    print(f"Pooled: Type1={pool_ts['ts1']}min, Type2={pool_ts['ts2']}min")
    print(f"  Wait: {pool_results['avg_wait_days_mean']:.2f}d, Util: {pool_results['avg_utilization_mean']*100:.1f}%")
    
    # Recommendation
    ded_wait = ded_results['avg_wait_days_mean']
    pool_wait = pool_results['avg_wait_days_mean']
    
    if pool_wait < ded_wait:
        improvement = ((ded_wait - pool_wait) / ded_wait * 100)
        print(f"\nRecommendation: POOLED (wait time reduced by {improvement:.1f}%)")
    else:
        print(f"\nRecommendation: DEDICATED")
    
    return {
        'dedicated': {'timeslots': ded_ts, 'results': ded_results},
        'pooled': {'timeslots': pool_ts, 'results': pool_results}
    }


if __name__ == "__main__":
    comparison = main()
