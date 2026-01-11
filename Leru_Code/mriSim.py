import numpy as np
import pandas as pd
import copy
import heapq
import os
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from enum import Enum

# CONFIGURATION PARAMETERS (from Part 1 Bootstrap Analysis)

# Scan duration parameters
TYPE1_DURATION_MEAN = 25.71  # Mean scan duration Type 1 (minutes)
TYPE1_DURATION_STD = 5.84  # Std deviation Type 1 (minutes)
TYPE2_DURATION_MEAN = 40.71  # Mean scan duration Type 2 (minutes)
TYPE2_DURATION_STD = 12.07  # Std deviation Type 2 (minutes)

# Arrival parameters
TYPE1_ARRIVALS_PER_DAY = 16.95  # Mean daily arrivals Type 1 (Poisson λ)
TYPE2_ARRIVALS_PER_DAY = 9.95  # Mean daily arrivals Type 2 (Poisson λ)

# Operating hours
WORK_START_HOUR = 8.0  # Facility opens at 8:00
WORK_END_HOUR = 17.0  # Facility closes at 17:00
WORK_HOURS_PER_DAY = WORK_END_HOUR - WORK_START_HOUR  # 9 hours per day
WORK_MINUTES_PER_DAY = WORK_HOURS_PER_DAY * 60  # 540 minutes per day

# Time constants
MINUTES_PER_HOUR = 60
MINUTES_PER_DAY = 24 * 60  # 1440 minutes

# Simulation defaults
DEFAULT_SIMULATION_DAYS = 21
DEFAULT_NUM_REPLICATIONS = 100

# Timeslot search ranges
TYPE1_TIMESLOT_RANGE = range(24, 36, 3)
TYPE2_TIMESLOT_RANGE = range(7, 13, 3)

@dataclass
class Patient:
    """Represents a patient's journey through the system."""
    patient_id: int  # Unique identifier
    patient_type: int  # 1 or 2
    call_time: float  # When they called (absolute minutes)
    call_day: int  # Day number of call
    scheduled_day: int = -1  # Day of appointment
    scheduled_time: float = -1  # Scheduled start time
    machine_id: int = -1  # Assigned machine (0 or 1)
    scan_duration: float = 0.0  # Actual scan duration
    scan_start_time: float = -1  # Actual start time
    scan_end_time: float = -1  # Actual end time

    def wait_time(self) -> float:
        """Calculate waiting time from call to appointment in days."""
        if self.scheduled_day >= 0:
            return self.scheduled_time - self.call_time
        return 0.0


# RANDOM NUMBER GENERATOR
class RandomGeneratorPatients:
    """Generates random variates for the simulation."""

    def __init__(self, emp_type2_data, interarrival_type2, seed: int = 42):
        self.rng = np.random.default_rng(seed)
        self.type2_empirical = emp_type2_data
        self.interarrival_type2 = interarrival_type2

    def generate_interarrival_time(self, patient_type: int) -> float:
        """Generate exponential interarrival time (minutes)."""
        if patient_type == 1:
            rate = TYPE1_ARRIVALS_PER_DAY
            minute_rate = rate / WORK_HOURS_PER_DAY / MINUTES_PER_HOUR
            return self.rng.exponential(1.0 / minute_rate)
        else:
            return self.rng.choice(self.interarrival_type2)

    def generate_scan_duration(self, patient_type: int) -> float:
        """Generate scan duration (minutes)."""
        if patient_type == 1:
            duration = self.rng.normal(TYPE1_DURATION_MEAN, TYPE1_DURATION_STD)
            return max(8.0, duration)
        else:
            return self.rng.choice(self.type2_empirical)

    def generate_patients(self, nr_days_limit: float) -> List[Patient]:
        """
        Pre-generate patients up to time_limit (in WORKING minutes).
        We generate Type 1 and Type 2 arrival streams independently, then merge and sort.

        Notes:
        - call_time is measured in 'working minutes' since day 0 at 08:00.
        - call_day is 0-based: 0, 1, 2, ...
        """
        patients: List[Patient] = []

        time_limit = nr_days_limit * WORK_MINUTES_PER_DAY
        # Generate Type 1 stream
        t1 = 0.0
        while True:
            t1 += self.generate_interarrival_time(patient_type=1)
            if t1 > time_limit:
                break

            call_day = int(t1 // WORK_MINUTES_PER_DAY)
            scan_duration = self.generate_scan_duration(patient_type=1)

            patients.append(Patient(
                patient_id=-1,  # assigned after sorting
                patient_type=1,
                call_time=t1,
                call_day=call_day,
                scan_duration=scan_duration
            ))

        # Generate Type 2 stream
        t2 = 0.0
        while True:
            t2 += self.generate_interarrival_time(patient_type=2)
            if t2 > time_limit:
                break

            call_day = int(t2 // WORK_MINUTES_PER_DAY)
            scan_duration = self.generate_scan_duration(patient_type=2)

            patients.append(Patient(
                patient_id=-1,  # assigned after sorting
                patient_type=2,
                call_time=t2,
                call_day=call_day,
                scan_duration=scan_duration
            ))

        # Merge streams: sort by call_time, then assign chronological ids
        patients.sort(key=lambda p: p.call_time)
        for i, p in enumerate(patients):
            p.patient_id = i

        return patients


# DISCRETE EVENT SIMULATION ENGINE

class MRISimulationDedicated:

    def __init__(self, timeslot_type1: int, timeslot_type2: int, patients_to_schedule: List[Patient],
                 num_days: int = DEFAULT_SIMULATION_DAYS, seed: int = 42):
        # Configuration
        self.timeslot_type1 = timeslot_type1
        self.timeslot_type2 = timeslot_type2
        self.num_days = num_days
        self.patients_to_schedule = patients_to_schedule

        # Simulation clock
        self.current_time = 0.0
        self.current_day = 0

        # Overtime per day per machine: overtime_by_day[d] = [m0_overtime, m1_overtime]
        self.overtime_by_day: List[List[float]] = [[0.0, 0.0] for _ in range(self.num_days)]

        # Two pointers = next available slot start time per machine (absolute working minutes)
        self.next_slot_start = [0.0, 0.0]

    def schedule_all_patients(self) -> None:
        # Make sure we schedule in call order
        patients = sorted(self.patients_to_schedule, key=lambda p: p.call_time)

        for p in patients:
            # Dedicated assignment
            m = 0 if p.patient_type == 1 else 1
            L = self.timeslot_type1 if p.patient_type == 1 else self.timeslot_type2

            # Earliest allowed is next working day (cannot schedule same day as call)
            earliest_day_start = (p.call_day + 1) * WORK_MINUTES_PER_DAY

            # Move pointer forward if it's behind earliest allowed
            if self.next_slot_start[m] < earliest_day_start:
                self.next_slot_start[m] = float(earliest_day_start)

            # Ensure the slot fits in the current day. If not, jump to next day start.
            day = int(self.next_slot_start[m] // WORK_MINUTES_PER_DAY)
            day_end = (day + 1) * WORK_MINUTES_PER_DAY
            if self.next_slot_start[m] + L > day_end:
                self.next_slot_start[m] = float((day + 1) * WORK_MINUTES_PER_DAY)
                day = day + 1  # (optional, for clarity)

            # Assign scheduled info
            p.machine_id = m
            p.scheduled_time = self.next_slot_start[m]
            p.scheduled_day = int(p.scheduled_time // WORK_MINUTES_PER_DAY)

            # Reserve the slot (pointer moves exactly once)
            self.next_slot_start[m] += L

    def run(self) -> Tuple[List[Patient], List[List[float]]]:
        self.schedule_all_patients()

        # Extend num_days if we scheduled beyond initial horizon
        max_day = 0
        for p in self.patients_to_schedule:
            if p.scheduled_day > max_day:
                max_day = p.scheduled_day
        if max_day + 1 > self.num_days:
            extra = (max_day + 1) - self.num_days
            for _ in range(extra):
                self.overtime_by_day.append([0.0, 0.0])
            self.num_days = max_day + 1

        # Group patients by day and machine (simple and fast)
        patients_by_day_machine: List[List[List[Patient]]] = [
            [[], []] for _ in range(self.num_days)
        ]
        for p in self.patients_to_schedule:
            if 0 <= p.scheduled_day < self.num_days:
                patients_by_day_machine[p.scheduled_day][p.machine_id].append(p)

        # Execute day-by-day and record overtime
        for day in range(self.num_days):
            self.current_day = day
            day_start = day * WORK_MINUTES_PER_DAY
            day_end = day_start + WORK_MINUTES_PER_DAY

            for m in [0, 1]:
                todays = sorted(patients_by_day_machine[day][m], key=lambda p: p.scheduled_time)

                machine_time = float(day_start)  # machine opens at 08:00 each day

                for p in todays:
                    p.scan_start_time = max(p.scheduled_time, machine_time)
                    p.scan_end_time = p.scan_start_time + p.scan_duration
                    machine_time = p.scan_end_time

                self.overtime_by_day[day][m] = max(0.0, machine_time - day_end)

        return self.patients_to_schedule, self.overtime_by_day


class MRISimulationFlexible:

    def __init__(self, timeslot_type1: int, timeslot_type2: int, patients_to_schedule: List[Patient],
                 num_days: int = DEFAULT_SIMULATION_DAYS, seed: int = 42):
        # Configuration
        self.timeslot_type1 = timeslot_type1
        self.timeslot_type2 = timeslot_type2
        self.num_days = num_days
        self.patients_to_schedule = patients_to_schedule

        # Simulation clock
        self.current_time = 0.0
        self.current_day = 0

        # Overtime per day per machine: overtime_by_day[d] = [m0_overtime, m1_overtime]
        self.overtime_by_day: List[List[float]] = [[0.0, 0.0] for _ in range(self.num_days)]

        # Two pointers = next available slot start time per machine (absolute working minutes)
        self.next_slot_start = [0.0, 0.0]

    def schedule_all_patients(self) -> None:
        # Schedule in call order
        patients = sorted(self.patients_to_schedule, key=lambda p: p.call_time)

        for p in patients:
            # Slot length depends on patient type
            L = self.timeslot_type1 if p.patient_type == 1 else self.timeslot_type2

            # Earliest allowed is next working day (cannot schedule same day as call)
            earliest_day_start = (p.call_day + 1) * WORK_MINUTES_PER_DAY

            # Compute candidate scheduled time on each machine
            candidate_times = [0.0, 0.0]
            for m in [0, 1]:
                t = self.next_slot_start[m]

                # must be at/after earliest allowed
                if t < earliest_day_start:
                    t = float(earliest_day_start)

                # ensure slot fits in the day; if not, jump to next day start
                day = int(t // WORK_MINUTES_PER_DAY)
                day_end = (day + 1) * WORK_MINUTES_PER_DAY
                if t + L > day_end:
                    t = float((day + 1) * WORK_MINUTES_PER_DAY)

                candidate_times[m] = t

            # Choose the machine with earliest candidate time (tie -> machine 0)
            if candidate_times[0] <= candidate_times[1]:
                chosen_machine = 0
            else:
                chosen_machine = 1

            scheduled_time = candidate_times[chosen_machine]

            # Assign scheduled info
            p.machine_id = chosen_machine
            p.scheduled_time = scheduled_time
            p.scheduled_day = int(scheduled_time // WORK_MINUTES_PER_DAY)

            # Reserve the slot (advance only the chosen machine pointer)
            self.next_slot_start[chosen_machine] = scheduled_time + L

    def run(self) -> Tuple[List[Patient], List[List[float]]]:
        self.schedule_all_patients()

        # Extend num_days if we scheduled beyond initial horizon
        max_day = 0
        for p in self.patients_to_schedule:
            if p.scheduled_day > max_day:
                max_day = p.scheduled_day

        if max_day + 1 > self.num_days:
            extra = (max_day + 1) - self.num_days
            for _ in range(extra):
                self.overtime_by_day.append([0.0, 0.0])
            self.num_days = max_day + 1

        # Group patients by day and machine
        patients_by_day_machine: List[List[List[Patient]]] = [
            [[], []] for _ in range(self.num_days)
        ]
        for p in self.patients_to_schedule:
            if 0 <= p.scheduled_day < self.num_days:
                patients_by_day_machine[p.scheduled_day][p.machine_id].append(p)

        # Execute day-by-day and record overtime
        for day in range(self.num_days):
            self.current_day = day
            day_start = day * WORK_MINUTES_PER_DAY
            day_end = day_start + WORK_MINUTES_PER_DAY

            for m in [0, 1]:
                todays = sorted(patients_by_day_machine[day][m], key=lambda p: p.scheduled_time)

                machine_time = float(day_start)  # machine opens at 08:00 each day

                for p in todays:
                    p.scan_start_time = max(p.scheduled_time, machine_time)
                    p.scan_end_time = p.scan_start_time + p.scan_duration
                    machine_time = p.scan_end_time

                self.overtime_by_day[day][m] = max(0.0, machine_time - day_end)

        return self.patients_to_schedule, self.overtime_by_day


def extract_type2_empirical_data(csv_path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns:
        type2_durations_min: np.ndarray of scan durations (minutes)
        type2_call_times_workmin: np.ndarray of call times in 'working minutes' timeline
        type2_interarrival_min: np.ndarray of interarrival times (minutes) between consecutive Type 2 calls
    """
    df = pd.read_csv(csv_path)
    df["Date"] = pd.to_datetime(df["Date"])

    # Filter Type 2
    type2 = df[df["PatientType"].str.strip() == "Type 2"].copy()

    # Duration: hours -> minutes
    type2_durations_min = (type2["Duration"].astype(float) * 60.0).to_numpy()

    # Build call time in a working-minutes timeline:
    # - each business day counts as WORK_MINUTES_PER_DAY
    # - within-day time is (Time - WORK_START_HOUR) * 60
    start_date = type2["Date"].min()

    type2["workday_index"] = type2["Date"].apply(lambda d: len(pd.bdate_range(start_date, d)) - 1)
    type2["within_day_min"] = (type2["Time"].astype(float) - WORK_START_HOUR) * 60.0
    type2["call_time_workmin"] = type2["workday_index"] * WORK_MINUTES_PER_DAY + type2["within_day_min"]

    # Sort by call time and compute interarrivals
    type2 = type2.sort_values("call_time_workmin").reset_index(drop=True)
    type2_call_times_workmin = type2["call_time_workmin"].to_numpy()

    type2_interarrival_min = np.diff(type2_call_times_workmin)  # already in minutes

    return type2_durations_min, type2_call_times_workmin, type2_interarrival_min

def avg_extra_wait_min(patients, patient_type: int) -> float:
    extra = [(p.scan_start_time - p.scheduled_time) for p in patients if p.patient_type == patient_type]
    if len(extra) == 0:
        return 0.0
    return float(np.mean(extra))

def avg_extra_wait_days(patients, patient_type: int) -> float:
    return avg_extra_wait_min(patients, patient_type) / WORK_MINUTES_PER_DAY


'''
if __name__ == "__main__":
    # -----------------------------
    # 1) Empirical Type 2 inputs
    # -----------------------------
    csv_path = "ScanRecords_Copy.csv"
    type2_durations_min, _, type2_interarrival_min = extract_type2_empirical_data(csv_path)
    type2_interarrival_min = type2_interarrival_min[type2_interarrival_min > 0]

    # -----------------------------
    # 2) Generate ONE patient list
    # -----------------------------
    seed = 42
    num_days = DEFAULT_SIMULATION_DAYS

    rng_patients = RandomGeneratorPatients(
        emp_type2_data=type2_durations_min,
        interarrival_type2=type2_interarrival_min,
        seed=seed
    )

    patients_base = rng_patients.generate_patients(nr_days_limit=num_days)

    # Make two identical copies (because the sims update patient fields)
    patients_dedicated = copy.deepcopy(patients_base)
    patients_flexible = copy.deepcopy(patients_base)

    # -----------------------------
    # 3) Run both simulations
    # -----------------------------
    timeslot_type1 = 24
    timeslot_type2 = 40

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

    # -----------------------------
    # 4) Compare outputs (simple KPIs)
    # -----------------------------
    def total_overtime(overtime_by_day):
        m0 = sum(day[0] for day in overtime_by_day)
        m1 = sum(day[1] for day in overtime_by_day)
        return m0, m1, m0 + m1

    def avg_wait_days(patients, patient_type):
        waits = [(p.scheduled_time - p.call_time) for p in patients if p.patient_type == patient_type]
        if len(waits) == 0:
            return 0.0
        return float(np.mean(waits) / WORK_MINUTES_PER_DAY)

    d_m0, d_m1, d_tot = total_overtime(ded_overtime)
    f_m0, f_m1, f_tot = total_overtime(flex_overtime)

    print("\n=== Dedicated ===")
    print("Overtime m0:", round(d_m0, 2), "min | m1:", round(d_m1, 2), "min | total:", round(d_tot, 2), "min")
    print("Avg wait Type 1:", round(avg_wait_days(ded_patients, 1), 3), "working-days")
    print("Avg wait Type 2:", round(avg_wait_days(ded_patients, 2), 3), "working-days")

    print("\n=== Flexible ===")
    print("Overtime m0:", round(f_m0, 2), "min | m1:", round(f_m1, 2), "min | total:", round(f_tot, 2), "min")
    print("Avg wait Type 1:", round(avg_wait_days(flex_patients, 1), 3), "working-days")
    print("Avg wait Type 2:", round(avg_wait_days(flex_patients, 2), 3), "working-days")

    # -----------------------------
    # 5) Optional: save both patient result tables
    # -----------------------------
    pd.DataFrame([p.__dict__ for p in ded_patients]).to_csv("patients_dedicated_results.csv", index=False)
    pd.DataFrame([p.__dict__ for p in flex_patients]).to_csv("patients_flexible_results.csv", index=False)
    pd.DataFrame(ded_overtime, columns=["overtime_m0", "overtime_m1"]).to_csv("overtime_dedicated.csv", index=False)
    pd.DataFrame(flex_overtime, columns=["overtime_m0", "overtime_m1"]).to_csv("overtime_flexible.csv", index=False)
    print("\n=== Dedicated ===")
    print("Avg extra wait Type 1:",
          round(avg_extra_wait_days(ded_patients, 1), 4), "working-days |",
          round(avg_extra_wait_min(ded_patients, 1), 2), "min")
    print("Avg extra wait Type 2:",
          round(avg_extra_wait_days(ded_patients, 2), 4), "working-days |",
          round(avg_extra_wait_min(ded_patients, 2), 2), "min")

    print("\n=== Flexible ===")
    print("Avg extra wait Type 1:",
          round(avg_extra_wait_days(flex_patients, 1), 4), "working-days |",
          round(avg_extra_wait_min(flex_patients, 1), 2), "min")
    print("Avg extra wait Type 2:",
          round(avg_extra_wait_days(flex_patients, 2), 4), "working-days |",
          round(avg_extra_wait_min(flex_patients, 2), 2), "min")
'''