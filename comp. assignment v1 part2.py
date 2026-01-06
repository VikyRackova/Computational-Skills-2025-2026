import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple
import heapq
from collections import defaultdict

TYPE1_DURATION_MEAN = 25.71  # minutes - from Part 1 bootstrap
TYPE1_DURATION_SD = 5.84
TYPE1_ARRIVALS_PER_DAY = 16.95
TYPE2_DURATION_MEAN = 40.71
TYPE2_DURATION_SD = 12.07
TYPE2_ARRIVALS_PER_DAY = 9.95
WORK_START = 8.0  # 8:00 AM
WORK_END = 17.0   # 5:00 PM
WORK_HOURS = WORK_END - WORK_START
LAMBDA_TYPE1 = TYPE1_ARRIVALS_PER_DAY / WORK_HOURS  # arrivals per hour
LAMBDA_TYPE2 = TYPE2_ARRIVALS_PER_DAY / WORK_HOURS

class EventType:
    ARRIVAL = "arrival"
    SCAN_START = "scan_start"
    SCAN_END = "scan_end"

@dataclass
class Patient:
    id: int
    patient_type: int
    call_time: float  # hours from simulation start
    call_day: int
    scan_duration: float  # minutes
    scheduled_time: float = None
    scheduled_day: int = None
    scan_start_time: float = None
    scan_end_time: float = None
    assigned_machine: int = None

@dataclass
class Event:
    time: float
    event_type: str
    patient: Patient = None
    machine: int = None
    def __lt__(self, other):
        return self.time < other.time

@dataclass
class MachineState:
    busy: bool = False
    current_patient: Patient = None
    schedule: List[Tuple[float, Patient]] = None
    def __post_init__(self):
        if self.schedule is None:
            self.schedule = []

class MRISimulation:
    def __init__(self, timeslot_type1, timeslot_type2, policy="separate", simulation_days=20, random_seed=None):
        self.timeslot_type1 = timeslot_type1
        self.timeslot_type2 = timeslot_type2
        self.policy = policy
        self.simulation_days = simulation_days
        if random_seed is not None:
            np.random.seed(random_seed)
        self.current_time = 0.0
        self.current_day = 0
        self.patient_counter = 0
        self.event_list = []
        self.machines = {1: MachineState(), 2: MachineState()}
        self.patients = []
        self.stats = defaultdict(list)
        
    def generate_scan_duration(self, patient_type):
        if patient_type == 1:
            return max(5.0, np.random.normal(TYPE1_DURATION_MEAN, TYPE1_DURATION_SD))
        else:
            return max(10.0, np.random.normal(TYPE2_DURATION_MEAN, TYPE2_DURATION_SD))
    
    def generate_interarrival_time(self, patient_type):
        if patient_type == 1:
            return np.random.exponential(1.0 / LAMBDA_TYPE1)
        else:
            return np.random.exponential(1.0 / LAMBDA_TYPE2)
    
    def schedule_event(self, event):
        heapq.heappush(self.event_list, event)
    
    def get_timeslot_length(self, patient_type):
        return self.timeslot_type1 if patient_type == 1 else self.timeslot_type2
    
    def time_to_day_hour(self, time):
        day = int(time / WORK_HOURS)
        hour = WORK_START + (time % WORK_HOURS)
        return day, hour
    
    def find_next_available_slot(self, patient):
        call_day, call_hour = self.time_to_day_hour(patient.call_time)
        timeslot_hours = self.get_timeslot_length(patient.patient_type) / 60.0
        search_day = call_day + 1  # earliest is next day
        machines_to_check = [patient.patient_type] if self.policy == "separate" else [1, 2]
        
        for day in range(search_day, search_day + 30):  # search up to 30 days
            for machine_id in machines_to_check:
                day_schedule = [(t, p) for t, p in self.machines[machine_id].schedule if self.time_to_day_hour(t)[0] == day]
                day_schedule.sort()
                candidate_time = day * WORK_HOURS
                
                if not day_schedule:  # no appointments, schedule at start
                    return day, candidate_time, machine_id
                
                for i, (start_time, p) in enumerate(day_schedule):
                    slot_end = start_time + self.get_timeslot_length(p.patient_type) / 60.0
                    if i == 0 and candidate_time + timeslot_hours <= start_time:  # fits before first
                        return day, candidate_time, machine_id
                    if i < len(day_schedule) - 1:
                        next_start = day_schedule[i + 1][0]
                        if slot_end + timeslot_hours <= next_start:  # fits between
                            return day, slot_end, machine_id
                    else:  # after last appointment
                        day_end = day * WORK_HOURS + WORK_HOURS
                        if slot_end + timeslot_hours <= day_end:
                            return day, slot_end, machine_id
        return None, None, None  # no slot found
    
    def handle_arrival(self, event):
        patient = event.patient
        sched_day, sched_time, machine_id = self.find_next_available_slot(patient)
        if sched_time is None:
            return  # skip if no slot found
        patient.scheduled_time = sched_time
        patient.scheduled_day = sched_day
        patient.assigned_machine = machine_id
        self.machines[machine_id].schedule.append((sched_time, patient))
        self.schedule_event(Event(time=sched_time, event_type=EventType.SCAN_START, patient=patient, machine=machine_id))
        next_arrival = event.time + self.generate_interarrival_time(patient.patient_type)
        next_day, _ = self.time_to_day_hour(next_arrival)
        if next_day < self.simulation_days:
            new_patient = Patient(id=self.patient_counter, patient_type=patient.patient_type, call_time=next_arrival, call_day=next_day, scan_duration=self.generate_scan_duration(patient.patient_type))
            self.patient_counter += 1
            self.patients.append(new_patient)
            self.schedule_event(Event(time=next_arrival, event_type=EventType.ARRIVAL, patient=new_patient))
    
    def handle_scan_start(self, event):
        machine = self.machines[event.machine]
        machine.busy = True
        machine.current_patient = event.patient
        event.patient.scan_start_time = event.time
        scan_end_time = event.time + event.patient.scan_duration / 60.0
        event.patient.scan_end_time = scan_end_time
        self.schedule_event(Event(time=scan_end_time, event_type=EventType.SCAN_END, patient=event.patient, machine=event.machine))
    
    def handle_scan_end(self, event):
        machine = self.machines[event.machine]
        machine.busy = False
        machine.current_patient = None
    
    def run(self):
        for patient_type in [1, 2]:  # initialize first arrivals
            first_arrival = self.generate_interarrival_time(patient_type)
            patient = Patient(id=self.patient_counter, patient_type=patient_type, call_time=first_arrival, call_day=0, scan_duration=self.generate_scan_duration(patient_type))
            self.patient_counter += 1
            self.patients.append(patient)
            self.schedule_event(Event(time=first_arrival, event_type=EventType.ARRIVAL, patient=patient))
        
        while self.event_list:
            event = heapq.heappop(self.event_list)
            self.current_time = event.time
            day, _ = self.time_to_day_hour(event.time)
            if day >= self.simulation_days:
                break
            if event.event_type == EventType.ARRIVAL:
                self.handle_arrival(event)
            elif event.event_type == EventType.SCAN_START:
                self.handle_scan_start(event)
            elif event.event_type == EventType.SCAN_END:
                self.handle_scan_end(event)
        self.calculate_statistics()
    
    def calculate_statistics(self):
        scheduled = [p for p in self.patients if p.scheduled_time is not None]
        if not scheduled:
            return
        waiting_days = [(p.scheduled_time - p.call_time) * WORK_HOURS * 60 / 60.0 / 24.0 for p in scheduled]  # convert to days
        self.stats['mean_waiting_days'] = np.mean(waiting_days)
        overtime_minutes = []
        for p in scheduled:
            if p.scan_end_time is not None:
                day, end_hour = self.time_to_day_hour(p.scan_end_time)
                if end_hour > WORK_END:
                    overtime_minutes.append((end_hour - WORK_END) * 60)
        self.stats['total_overtime_hours'] = sum(overtime_minutes) / 60.0
        self.stats['pct_overtime_scans'] = len(overtime_minutes) / len(scheduled) * 100
        total_scan = sum([p.scan_duration for p in scheduled if p.scan_end_time is not None])
        available = self.simulation_days * WORK_HOURS * 60 * 2  # 2 machines
        self.stats['utilization_pct'] = (total_scan / available) * 100
        self.stats['patients_served'] = len(scheduled)

def run_experiment(ts1, ts2, policy, n_reps=50):
    results = []
    for i in range(n_reps):
        sim = MRISimulation(ts1, ts2, policy, simulation_days=20, random_seed=i)
        sim.run()
        results.append(sim.stats)
    agg = {}
    for key in results[0].keys():
        vals = [r[key] for r in results]
        agg[f'{key}_mean'] = np.mean(vals)
        agg[f'{key}_std'] = np.std(vals)
    return agg

def objective_score(results):  # lower is better
    return (results['mean_waiting_days_mean'] * 2 + results['pct_overtime_scans_mean'] * 1 + (100 - results['utilization_pct_mean']))

def grid_search(type1_range, type2_range, policy, n_reps=30):
    results = []
    total = len(type1_range) * len(type2_range)
    num = 0
    for ts1 in type1_range:
        for ts2 in type2_range:
            num += 1
            print(f"  Testing {num}/{total}: Type1={ts1}min, Type2={ts2}min", end="\r")
            sim_res = run_experiment(ts1, ts2, policy, n_reps)
            results.append({'ts1': ts1, 'ts2': ts2, 'policy': policy, 'score': objective_score(sim_res), 'waiting_days': sim_res['mean_waiting_days_mean'], 'overtime_pct': sim_res['pct_overtime_scans_mean'], 'utilization_pct': sim_res['utilization_pct_mean']})
    print()
    return pd.DataFrame(results)

if __name__ == "__main__":
    print("=" * 70)
    print("FINDING OPTIMAL TIMESLOTS")
    print("=" * 70)
    
    print("\nSEPARATE policy:")
    sep_grid = grid_search(range(28, 43, 2), range(58, 73, 2), "separate", n_reps=30)  # test Type1: 28-42, Type2: 58-72
    best_sep = sep_grid.loc[sep_grid['score'].idxmin()]
    
    print("\nPOOLED policy:")
    pool_grid = grid_search(range(28, 43, 2), range(58, 73, 2), "pooled", n_reps=30)
    best_pool = pool_grid.loc[pool_grid['score'].idxmin()]
    
    print("\n" + "=" * 70)
    print("OPTIMAL TIMESLOTS")
    print("=" * 70)
    print(f"\nSEPARATE: Type1={best_sep['ts1']:.0f}min, Type2={best_sep['ts2']:.0f}min | Wait={best_sep['waiting_days']:.2f}d, Overtime={best_sep['overtime_pct']:.1f}%, Util={best_sep['utilization_pct']:.1f}%, Score={best_sep['score']:.2f}")
    print(f"POOLED:   Type1={best_pool['ts1']:.0f}min, Type2={best_pool['ts2']:.0f}min | Wait={best_pool['waiting_days']:.2f}d, Overtime={best_pool['overtime_pct']:.1f}%, Util={best_pool['utilization_pct']:.1f}%, Score={best_pool['score']:.2f}")
    
    improvement = ((best_sep['score'] - best_pool['score']) / best_sep['score']) * 100
    print(f"\n{'POOLED RECOMMENDED' if best_pool['score'] < best_sep['score'] else 'SEPARATE RECOMMENDED'} (Score improvement: {abs(improvement):.1f}%)")
    print("=" * 70)
    
    sep_grid.to_csv('optimal_timeslots_separate.csv', index=False)  # save detailed results
    pool_grid.to_csv('optimal_timeslots_pooled.csv', index=False)
    print("\nResults saved to CSV files.")