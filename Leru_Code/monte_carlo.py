import numpy as np
import pandas as pd
import copy
from typing import Dict
import os

# IMPORTANT:
# Change "simulation_core" to the filename where your classes/functions live (without .py)
from mriSim import (
    DEFAULT_SIMULATION_DAYS,
    DEFAULT_NUM_REPLICATIONS,
    WORK_MINUTES_PER_DAY,
    extract_type2_empirical_data,
    RandomGeneratorPatients,
    MRISimulationDedicated,
    MRISimulationFlexible,
    avg_extra_wait_min,
    avg_extra_wait_days,
)


def total_overtime(overtime_by_day):
    m0 = sum(day[0] for day in overtime_by_day)
    m1 = sum(day[1] for day in overtime_by_day)
    return float(m0), float(m1), float(m0 + m1)


def avg_wait_days(patients, patient_type: int) -> float:
    waits = [(p.scheduled_time - p.call_time) for p in patients if p.patient_type == patient_type]
    if len(waits) == 0:
        return 0.0
    return float(np.mean(waits) / WORK_MINUTES_PER_DAY)


def run_one_replication(
    rep: int,
    seed: int,
    num_days: int,
    timeslot_type1: int,
    timeslot_type2: int,
    type2_durations_min: np.ndarray,
    type2_interarrival_min: np.ndarray,
) -> Dict[str, float]:

    rng_patients = RandomGeneratorPatients(
        emp_type2_data=type2_durations_min,
        interarrival_type2=type2_interarrival_min,
        seed=seed
    )

    patients_base = rng_patients.generate_patients(nr_days_limit=num_days)

    patients_dedicated = copy.deepcopy(patients_base)
    patients_flexible = copy.deepcopy(patients_base)

    sim_ded = MRISimulationDedicated(
        timeslot_type1=timeslot_type1,
        timeslot_type2=timeslot_type2,
        patients_to_schedule=patients_dedicated,
        num_days=num_days,
        seed=seed
    )
    ded_patients, ded_overtime = sim_ded.run()

    sim_flex = MRISimulationFlexible(
        timeslot_type1=timeslot_type1,
        timeslot_type2=timeslot_type2,
        patients_to_schedule=patients_flexible,
        num_days=num_days,
        seed=seed
    )
    flex_patients, flex_overtime = sim_flex.run()

    d_m0, d_m1, d_tot = total_overtime(ded_overtime)
    f_m0, f_m1, f_tot = total_overtime(flex_overtime)

    n_total = len(patients_base)
    n_type1 = sum(1 for p in patients_base if p.patient_type == 1)
    n_type2 = n_total - n_type1

    row = {
        "replication": rep,
        "seed": seed,
        "num_days": num_days,
        "timeslot_type1": timeslot_type1,
        "timeslot_type2": timeslot_type2,
        "n_patients": n_total,
        "n_type1": n_type1,
        "n_type2": n_type2,

        # Dedicated KPIs
        "ded_overtime_m0_min": d_m0,
        "ded_overtime_m1_min": d_m1,
        "ded_overtime_total_min": d_tot,
        "ded_avg_wait_type1_days": avg_wait_days(ded_patients, 1),
        "ded_avg_wait_type2_days": avg_wait_days(ded_patients, 2),
        "ded_avg_extra_wait_type1_min": avg_extra_wait_min(ded_patients, 1),
        "ded_avg_extra_wait_type2_min": avg_extra_wait_min(ded_patients, 2),
        "ded_avg_extra_wait_type1_days": avg_extra_wait_days(ded_patients, 1),
        "ded_avg_extra_wait_type2_days": avg_extra_wait_days(ded_patients, 2),

        # Flexible KPIs
        "flex_overtime_m0_min": f_m0,
        "flex_overtime_m1_min": f_m1,
        "flex_overtime_total_min": f_tot,
        "flex_avg_wait_type1_days": avg_wait_days(flex_patients, 1),
        "flex_avg_wait_type2_days": avg_wait_days(flex_patients, 2),
        "flex_avg_extra_wait_type1_min": avg_extra_wait_min(flex_patients, 1),
        "flex_avg_extra_wait_type2_min": avg_extra_wait_min(flex_patients, 2),
        "flex_avg_extra_wait_type1_days": avg_extra_wait_days(flex_patients, 1),
        "flex_avg_extra_wait_type2_days": avg_extra_wait_days(flex_patients, 2),

        # Differences (Flexible - Dedicated)
        "diff_overtime_total_min": f_tot - d_tot,
        "diff_avg_wait_type1_days": avg_wait_days(flex_patients, 1) - avg_wait_days(ded_patients, 1),
        "diff_avg_wait_type2_days": avg_wait_days(flex_patients, 2) - avg_wait_days(ded_patients, 2),
        "diff_avg_extra_wait_type1_min": avg_extra_wait_min(flex_patients, 1) - avg_extra_wait_min(ded_patients, 1),
        "diff_avg_extra_wait_type2_min": avg_extra_wait_min(flex_patients, 2) - avg_extra_wait_min(ded_patients, 2),
    }

    return row

if __name__ == "__main__":
    # -----------------------------
    # Settings
    # -----------------------------
    csv_path = "ScanRecords_Copy.csv"   # put this next to monte_carlo.py (or change path)
    num_days = DEFAULT_SIMULATION_DAYS
    num_replications = 10000

    # fixed timeslot lengths (minutes)
    timeslot_type1 = 32
    timeslot_type2 = 64

    base_seed = 42

    # -----------------------------
    # Load empirical Type 2 once
    # -----------------------------
    type2_durations_min, _, type2_interarrival_min = extract_type2_empirical_data(csv_path)
    type2_interarrival_min = type2_interarrival_min[type2_interarrival_min > 0]

    # -----------------------------
    # Monte Carlo loop
    # -----------------------------
    rows = []
    for rep in range(num_replications):
        seed = base_seed + rep
        rows.append(
            run_one_replication(
                rep=rep,
                seed=seed,
                num_days=num_days,
                timeslot_type1=timeslot_type1,
                timeslot_type2=timeslot_type2,
                type2_durations_min=type2_durations_min,
                type2_interarrival_min=type2_interarrival_min,
            )
        )

    results = pd.DataFrame(rows)
    results.to_csv("monte_carlo_results.csv", index=False)

    # Quick summary prints
    print("\nSaved: monte_carlo_results.csv")
    print("Replications:", num_replications)

    print("\nMean overtime total (min):")
    print("Dedicated:", round(float(results["ded_overtime_total_min"].mean()), 2))
    print("Flexible :", round(float(results["flex_overtime_total_min"].mean()), 2))
    print("Diff (Flex - Ded):", round(float(results["diff_overtime_total_min"].mean()), 2))

    print("\nMean avg wait Type 1 (days): Ded =", round(float(results["ded_avg_wait_type1_days"].mean()), 3),
          " Flex =", round(float(results["flex_avg_wait_type1_days"].mean()), 3))
    print("Mean avg wait Type 2 (days): Ded =", round(float(results["ded_avg_wait_type2_days"].mean()), 3),
          " Flex =", round(float(results["flex_avg_wait_type2_days"].mean()), 3))

    print("\nMean extra wait Type 1 (min): Ded =", round(float(results["ded_avg_extra_wait_type1_min"].mean()), 2),
          " Flex =", round(float(results["flex_avg_extra_wait_type1_min"].mean()), 2))
    print("Mean extra wait Type 2 (min): Ded =", round(float(results["ded_avg_extra_wait_type2_min"].mean()), 2),
          " Flex =", round(float(results["flex_avg_extra_wait_type2_min"].mean()), 2))
