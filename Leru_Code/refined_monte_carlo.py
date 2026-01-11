import numpy as np
import pandas as pd
import copy

from mriSim import (
    DEFAULT_SIMULATION_DAYS,
    WORK_MINUTES_PER_DAY,
    extract_type2_empirical_data,
    RandomGeneratorPatients,
    MRISimulationDedicated,
    MRISimulationFlexible,
    avg_extra_wait_min,
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


def pre_generate_patient_lists(num_replications, base_seed, num_days, type2_durations_min, type2_interarrival_min):
    patient_lists = []
    for rep in range(num_replications):
        seed = base_seed + rep
        rng_patients = RandomGeneratorPatients(
            emp_type2_data=type2_durations_min,
            interarrival_type2=type2_interarrival_min,
            seed=seed
        )
        patient_lists.append(rng_patients.generate_patients(nr_days_limit=num_days))
    return patient_lists


if __name__ == "__main__":
    # -----------------------------
    # CHOOSE POLICY HERE
    # -----------------------------
    policy = "flexible"   # "dedicated" or "flexible"

    # -----------------------------
    # Settings
    # -----------------------------
    csv_path = "ScanRecords_Copy.csv"
    num_days = DEFAULT_SIMULATION_DAYS
    base_seed = 42
    refined_replications = 1000

    if policy == "dedicated":
        top10_path = "dedicated_timeslot_grid_top101.csv"
        out_path = "dedicated_top10_refined_1000.csv"
        SimClass = MRISimulationDedicated
    else:
        top10_path = "flexible_timeslot_grid_top103.csv"
        out_path = "flexible_top10_refined_1000.csv"
        SimClass = MRISimulationFlexible

    # -----------------------------
    # Load top 10 combos for that policy
    # -----------------------------
    top10 = pd.read_csv(top10_path)
    combos = list(zip(top10["timeslot_type1"].astype(int), top10["timeslot_type2"].astype(int)))

    # -----------------------------
    # Empirical Type 2 inputs
    # -----------------------------
    type2_durations_min, _, type2_interarrival_min = extract_type2_empirical_data(csv_path)
    type2_interarrival_min = type2_interarrival_min[type2_interarrival_min > 0]

    # -----------------------------
    # Pre-generate patient lists (CRN)
    # -----------------------------
    patient_lists = pre_generate_patient_lists(
        num_replications=refined_replications,
        base_seed=base_seed,
        num_days=num_days,
        type2_durations_min=type2_durations_min,
        type2_interarrival_min=type2_interarrival_min
    )

    # -----------------------------
    # Run refined MC for top 10
    # -----------------------------
    rows = []

    for ts1, ts2 in combos:
        rep_overtime_m0 = []
        rep_overtime_m1 = []
        rep_overtime_tot = []

        rep_wait_t1 = []
        rep_wait_t2 = []

        rep_extra_t1 = []
        rep_extra_t2 = []

        for rep in range(refined_replications):
            patients_rep = copy.deepcopy(patient_lists[rep])

            sim = SimClass(
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
            rep_overtime_tot.append(tot)

            rep_wait_t1.append(avg_wait_days(patients_done, 1))
            rep_wait_t2.append(avg_wait_days(patients_done, 2))

            rep_extra_t1.append(avg_extra_wait_min(patients_done, 1))
            rep_extra_t2.append(avg_extra_wait_min(patients_done, 2))

        rows.append({
            "policy": policy,
            "timeslot_type1": ts1,
            "timeslot_type2": ts2,
            "replications": refined_replications,

            "mean_overtime_total_min": float(np.mean(rep_overtime_tot)),
            "std_overtime_total_min": float(np.std(rep_overtime_tot, ddof=1)),

            "mean_avg_wait_type1_days": float(np.mean(rep_wait_t1)),
            "std_avg_wait_type1_days": float(np.std(rep_wait_t1, ddof=1)),
            "mean_avg_wait_type2_days": float(np.mean(rep_wait_t2)),
            "std_avg_wait_type2_days": float(np.std(rep_wait_t2, ddof=1)),

            "mean_avg_extra_wait_type1_min": float(np.mean(rep_extra_t1)),
            "std_avg_extra_wait_type1_min": float(np.std(rep_extra_t1, ddof=1)),
            "mean_avg_extra_wait_type2_min": float(np.mean(rep_extra_t2)),
            "std_avg_extra_wait_type2_min": float(np.std(rep_extra_t2, ddof=1)),
        })

    refined = pd.DataFrame(rows).sort_values(
        by=[
            "mean_overtime_total_min",
            "mean_avg_wait_type1_days",
            "mean_avg_wait_type2_days",
            "mean_avg_extra_wait_type1_min",
            "mean_avg_extra_wait_type2_min",
        ],
        ascending=True
    ).reset_index(drop=True)

    refined.to_csv(out_path, index=False)

    print("Policy:", policy)
    print("Saved:", out_path)
    print("\nTop result:")
    print(refined.head(1))
