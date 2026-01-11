import numpy as np
import pandas as pd
import copy

# Change "simulation_core" to your core file name (without .py)
from mriSim import (
    DEFAULT_SIMULATION_DAYS,
    DEFAULT_NUM_REPLICATIONS,
    WORK_MINUTES_PER_DAY,
    extract_type2_empirical_data,
    RandomGeneratorPatients,
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


def pre_generate_patient_lists(
    num_replications: int,
    base_seed: int,
    num_days: int,
    type2_durations_min: np.ndarray,
    type2_interarrival_min: np.ndarray
):
    patient_lists = []

    for rep in range(num_replications):
        seed = base_seed + rep

        rng_patients = RandomGeneratorPatients(
            emp_type2_data=type2_durations_min,
            interarrival_type2=type2_interarrival_min,
            seed=seed
        )

        patients = rng_patients.generate_patients(nr_days_limit=num_days)
        patient_lists.append(patients)

    return patient_lists


if __name__ == "__main__":
    # -----------------------------
    # Settings
    # -----------------------------
    csv_path = "ScanRecords_Copy.csv"
    num_days = DEFAULT_SIMULATION_DAYS
    num_replications = DEFAULT_NUM_REPLICATIONS
    base_seed = 42

    # Timeslot ranges (inclusive)
    type1_range = range(24, 34, 1)   # 20..30
    type2_range = range(57, 68, 1)   # 35..45

    # -----------------------------
    # Empirical Type 2 inputs
    # -----------------------------
    type2_durations_min, _, type2_interarrival_min = extract_type2_empirical_data(csv_path)
    type2_interarrival_min = type2_interarrival_min[type2_interarrival_min > 0]

    # -----------------------------
    # Pre-generate patient lists (CRN)
    # -----------------------------
    patient_lists = pre_generate_patient_lists(
        num_replications=num_replications,
        base_seed=base_seed,
        num_days=num_days,
        type2_durations_min=type2_durations_min,
        type2_interarrival_min=type2_interarrival_min
    )

    # -----------------------------
    # Grid search (Flexible)
    # -----------------------------
    rows = []

    for ts1 in type1_range:
        for ts2 in type2_range:
            rep_overtime_total = []
            rep_overtime_m0 = []
            rep_overtime_m1 = []

            rep_wait_t1 = []
            rep_wait_t2 = []

            rep_extra_t1_min = []
            rep_extra_t2_min = []

            rep_n_total = []
            rep_n_t1 = []
            rep_n_t2 = []

            for rep in range(num_replications):
                patients_rep = copy.deepcopy(patient_lists[rep])

                sim = MRISimulationFlexible(
                    timeslot_type1=ts1,
                    timeslot_type2=ts2,
                    patients_to_schedule=patients_rep,
                    num_days=num_days,
                    seed=base_seed + rep
                )

                patients_done, overtime_by_day = sim.run()

                m0, m1, tot = total_overtime(overtime_by_day)

                rep_overtime_m0.append(m0)
                rep_overtime_m1.append(m1)
                rep_overtime_total.append(tot)

                rep_wait_t1.append(avg_wait_days(patients_done, 1))
                rep_wait_t2.append(avg_wait_days(patients_done, 2))

                rep_extra_t1_min.append(avg_extra_wait_min(patients_done, 1))
                rep_extra_t2_min.append(avg_extra_wait_min(patients_done, 2))

                n_total = len(patients_done)
                n_t1 = sum(1 for p in patients_done if p.patient_type == 1)
                n_t2 = n_total - n_t1

                rep_n_total.append(n_total)
                rep_n_t1.append(n_t1)
                rep_n_t2.append(n_t2)

            row = {
                "timeslot_type1": ts1,
                "timeslot_type2": ts2,

                "mean_overtime_total_min": float(np.mean(rep_overtime_total)),
                "std_overtime_total_min": float(np.std(rep_overtime_total, ddof=1)),

                "mean_avg_wait_type1_days": float(np.mean(rep_wait_t1)),
                "std_avg_wait_type1_days": float(np.std(rep_wait_t1, ddof=1)),
                "mean_avg_wait_type2_days": float(np.mean(rep_wait_t2)),
                "std_avg_wait_type2_days": float(np.std(rep_wait_t2, ddof=1)),

                "mean_avg_extra_wait_type1_min": float(np.mean(rep_extra_t1_min)),
                "std_avg_extra_wait_type1_min": float(np.std(rep_extra_t1_min, ddof=1)),
                "mean_avg_extra_wait_type2_min": float(np.mean(rep_extra_t2_min)),
                "std_avg_extra_wait_type2_min": float(np.std(rep_extra_t2_min, ddof=1)),

                "mean_avg_extra_wait_type1_days": float(np.mean(np.array(rep_extra_t1_min) / WORK_MINUTES_PER_DAY)),
                "mean_avg_extra_wait_type2_days": float(np.mean(np.array(rep_extra_t2_min) / WORK_MINUTES_PER_DAY)),
            }

            rows.append(row)

    results = pd.DataFrame(rows)

    # -----------------------------
    # Identify "best" (flexible)
    # -----------------------------
    results_sorted = results.sort_values(
        by=[
            "mean_overtime_total_min",
            "mean_avg_wait_type1_days",
            "mean_avg_wait_type2_days",
            "mean_avg_extra_wait_type1_min",
            "mean_avg_extra_wait_type2_min",
        ],
        ascending=True
    ).reset_index(drop=True)

    results_sorted.to_csv("flexible_timeslot_grid_results3.csv", index=False)
    results_sorted.head(10).to_csv("flexible_timeslot_grid_top103.csv", index=False)

    best = results_sorted.iloc[0]
    print("\nBEST (Flexible) by: overtime_total -> avg_waits -> extra_waits")
    print("timeslot_type1 =", int(best["timeslot_type1"]),
          "timeslot_type2 =", int(best["timeslot_type2"]))
    print("mean_overtime_total_min =", round(float(best["mean_overtime_total_min"]), 2))
    print("mean_avg_wait_type1_days =", round(float(best["mean_avg_wait_type1_days"]), 3),
          "mean_avg_wait_type2_days =", round(float(best["mean_avg_wait_type2_days"]), 3))
    print("mean_avg_extra_wait_type1_min =", round(float(best["mean_avg_extra_wait_type1_min"]), 2),
          "mean_avg_extra_wait_type2_min =", round(float(best["mean_avg_extra_wait_type2_min"]), 2))
    print("\nSaved: flexible_timeslot_grid_results.csv and flexible_timeslot_grid_top10.csv")
