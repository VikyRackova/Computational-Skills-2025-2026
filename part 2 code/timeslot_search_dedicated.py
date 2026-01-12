import numpy as np
import pandas as pd
import copy

from mriSim import (
    DEFAULT_SIMULATION_DAYS,
    DEFAULT_NUM_REPLICATIONS,
    WORK_MINUTES_PER_DAY,
    extract_type2_empirical_data,
    RandomGeneratorPatients,
    MRISimulationDedicated,
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

        # for each replication generate patients
        rng_patients = RandomGeneratorPatients(
            emp_type2_data=type2_durations_min,
            interarrival_type2=type2_interarrival_min,
            seed=seed
        )
        patients = rng_patients.generate_patients(nr_days_limit=num_days)

        patient_lists.append(patients)

    return patient_lists


if __name__ == "__main__":

    csv_path = "ScanRecords_Copy.csv"
    num_days = DEFAULT_SIMULATION_DAYS
    num_replications = DEFAULT_NUM_REPLICATIONS
    base_seed = 42

    # timeslot ranges
    type1_range = range(24, 35, 2)
    type2_range = range(40, 68, 2)

    # the empirical inputs for type 2
    type2_durations_min, _, type2_interarrival_min = extract_type2_empirical_data(csv_path)
    type2_interarrival_min = type2_interarrival_min[type2_interarrival_min > 0]

    # pregenerate the list of aptients
    patient_lists = pre_generate_patient_lists(
        num_replications=num_replications,
        base_seed=base_seed,
        num_days=num_days,
        type2_durations_min=type2_durations_min,
        type2_interarrival_min=type2_interarrival_min
    )

    # now we search for through all the combination of timeslots
    rows = []

    for ts1 in type1_range:
        for ts2 in type2_range:

            rep_avg_overtime_per_day = []
            rep_avg_wait_days = []
            rep_avg_extra_wait_min = []

            for rep in range(num_replications):
                patients_rep = copy.deepcopy(patient_lists[rep])

                sim = MRISimulationDedicated(
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

            row = {
                "timeslot_type1": ts1,
                "timeslot_type2": ts2,

                "mean_avg_overtime_per_day_min": float(np.mean(rep_avg_overtime_per_day)),
                "std_avg_overtime_per_day_min": float(np.std(rep_avg_overtime_per_day, ddof=1)),

                "mean_avg_wait_days": float(np.mean(rep_avg_wait_days)),
                "std_avg_wait_days": float(np.std(rep_avg_wait_days, ddof=1)),

                "mean_avg_extra_wait_min": float(np.mean(rep_avg_extra_wait_min)),
                "std_avg_extra_wait_min": float(np.std(rep_avg_extra_wait_min, ddof=1)),
            }

            rows.append(row)

    results = pd.DataFrame(rows)

    # and now we choose the best one
    results_sorted = results.sort_values(
        by=[
            "mean_avg_overtime_per_day_min",
            "mean_avg_wait_days",
            "mean_avg_extra_wait_min",
        ],
        ascending=True
    ).reset_index(drop=True)

    # save full table and the top 10 based on overtime
    results_sorted.to_csv("dedicated_timeslot_grid_results_final.csv", index=False)
    results_sorted.head(10).to_csv("dedicated_timeslot_grid_top10_final.csv", index=False)

    # print the best on overtime
    best = results_sorted.iloc[0]
    print("\nBEST (Dedicated) by: avg_overtime/day -> avg_wait_days -> avg_extra_wait_min")
    print("timeslot_type1 =", int(best["timeslot_type1"]),
          "timeslot_type2 =", int(best["timeslot_type2"]))
    print("mean_avg_overtime_per_day_min =", round(float(best["mean_avg_overtime_per_day_min"]), 3))
    print("mean_avg_wait_days =", round(float(best["mean_avg_wait_days"]), 4))
    print("mean_avg_extra_wait_min =", round(float(best["mean_avg_extra_wait_min"]), 3))
    print("\nSaved: dedicated_timeslot_grid_results3.csv and dedicated_timeslot_grid_top103.csv")
