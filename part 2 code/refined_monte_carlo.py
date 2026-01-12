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
)


def avg_overtime_per_day(overtime_by_day) -> float:
    if len(overtime_by_day) == 0:
        return 0.0
    total = sum(day[0] + day[1] for day in overtime_by_day)
    return float(total / len(overtime_by_day))


def avg_wait_days(patients) -> float:
    waits = [(p.scheduled_time - p.call_time) for p in patients]
    if len(waits) == 0:
        return 0.0
    return float(np.mean(waits) / WORK_MINUTES_PER_DAY)


def avg_extra_wait_min(patients) -> float:
    extra = [(p.scan_start_time - p.scheduled_time) for p in patients]
    if len(extra) == 0:
        return 0.0
    return float(np.mean(extra))


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

    # choose the policy here
    policy = "flexible"   # "dedicated" or "flexible"


    csv_path = "ScanRecords_Copy.csv"
    num_days = DEFAULT_SIMULATION_DAYS
    base_seed = 42
    refined_replications = 1000

    if policy == "dedicated":
        top10_path = "dedicated_timeslot_grid_top103.csv"
        out_path = "dedicated_top10_refined_1000.csv"
        SimClass = MRISimulationDedicated
    else:
        top10_path = "flexible_timeslot_grid_top103.csv"
        out_path = "flexible_top10_refined_1000.csv"
        SimClass = MRISimulationFlexible


    top10 = pd.read_csv(top10_path)
    combos = list(zip(top10["timeslot_type1"].astype(int), top10["timeslot_type2"].astype(int)))


    type2_durations_min, _, type2_interarrival_min = extract_type2_empirical_data(csv_path)
    type2_interarrival_min = type2_interarrival_min[type2_interarrival_min > 0]


    patient_lists = pre_generate_patient_lists(
        num_replications=refined_replications,
        base_seed=base_seed,
        num_days=num_days,
        type2_durations_min=type2_durations_min,
        type2_interarrival_min=type2_interarrival_min
    )

    rows = []

    for ts1, ts2 in combos:
        rep_avg_overtime_per_day = []
        rep_avg_wait_days = []
        rep_avg_extra_wait_min = []

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

            rep_avg_overtime_per_day.append(avg_overtime_per_day(overtime_by_day))
            rep_avg_wait_days.append(avg_wait_days(patients_done))
            rep_avg_extra_wait_min.append(avg_extra_wait_min(patients_done))

        rows.append({
            "policy": policy,
            "timeslot_type1": ts1,
            "timeslot_type2": ts2,
            "replications": refined_replications,

            "mean_avg_overtime_per_day_min": float(np.mean(rep_avg_overtime_per_day)),
            "std_avg_overtime_per_day_min": float(np.std(rep_avg_overtime_per_day, ddof=1)),

            "mean_avg_wait_days": float(np.mean(rep_avg_wait_days)),
            "std_avg_wait_days": float(np.std(rep_avg_wait_days, ddof=1)),

            "mean_avg_extra_wait_min": float(np.mean(rep_avg_extra_wait_min)),
            "std_avg_extra_wait_min": float(np.std(rep_avg_extra_wait_min, ddof=1)),
        })

    refined = pd.DataFrame(rows).sort_values(
        by=[
            "mean_avg_overtime_per_day_min",
            "mean_avg_wait_days",
            "mean_avg_extra_wait_min",
        ],
        ascending=True
    ).reset_index(drop=True)

    refined.to_csv(out_path, index=False)

    print("Policy:", policy)
    print("Saved:", out_path)
    print("\nTop result:")
    print(refined.head(1))
