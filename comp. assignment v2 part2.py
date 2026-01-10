import numpy as np                                                              # Numerical computing and random number generation
import pandas as pd                                                             # Data manipulation and CSV export
from dataclasses import dataclass                                               # Decorator for creating data classes
from typing import List, Tuple                                                  # Type hints for better code clarity
import heapq                                                                    # Priority queue for Future Event List (FEL)
from collections import defaultdict                                             # Dictionary with default values for stats

# --- CONFIGURATION (from Part 1 bootstrap analysis) ---
TYPE1_DURATION_MEAN = 25.71                                                     # Mean scan duration Type 1 (minutes)
TYPE1_DURATION_SD = 5.84                                                        # Standard deviation Type 1 (minutes)
TYPE1_ARRIVALS_PER_DAY = 16.95                                                  # Expected daily arrivals Type 1

TYPE2_DURATION_MEAN = 40.71                                                     # Mean scan duration Type 2 (minutes)
TYPE2_DURATION_SD = 12.07                                                       # Standard deviation Type 2 (minutes)
TYPE2_ARRIVALS_PER_DAY = 9.95                                                   # Expected daily arrivals Type 2

WORK_START = 8.0                                                                # Facility opens at 8:00 AM
WORK_END = 17.0                                                                 # Facility closes at 5:00 PM
WORK_HOURS = WORK_END - WORK_START                                              # 9 working hours per day
LAMBDA_TYPE1 = TYPE1_ARRIVALS_PER_DAY / WORK_HOURS                              # Arrival rate Type 1 (per hour) for Poisson process
LAMBDA_TYPE2 = TYPE2_ARRIVALS_PER_DAY / WORK_HOURS                              # Arrival rate Type 2 (per hour) for Poisson process

class EventType:                                                                # Enumeration of discrete event types
    ARRIVAL = "arrival"                                                         # Patient calls to schedule appointment
    SCAN_START = "scan_start"                                                   # Patient's MRI scan begins
    SCAN_END = "scan_end"                                                       # Patient's MRI scan completes

@dataclass
class Patient:                                                                  # Data structure representing a patient entity
    id: int                                                                     # Unique patient identifier
    patient_type: int                                                           # 1 or 2 (determines scan characteristics)
    call_time: float                                                            # Simulation time when patient called
    call_day: int                                                               # Day number when patient called
    scan_duration: float                                                        # Actual scan duration in minutes (randomly generated)
    scheduled_time: float = None                                                # Simulation time of scheduled appointment
    scheduled_day: int = None                                                   # Day number of scheduled appointment
    scan_start_time: float = None                                               # Actual start time of scan
    scan_end_time: float = None                                                 # Actual end time of scan
    assigned_machine: int = None                                                # Machine ID (1 or 2) assigned to patient

@dataclass
class Event:                                                                    # Data structure representing a simulation event
    time: float                                                                 # Simulation time when event occurs
    event_type: str                                                             # Type of event (ARRIVAL, SCAN_START, SCAN_END)
    patient: Patient = None                                                     # Patient associated with this event
    machine: int = None                                                         # Machine associated with this event
    def __lt__(self, other):                                                    # Comparison method for priority queue ordering
        return self.time < other.time                                           # Earlier events have higher priority

@dataclass
class MachineState:                                                             # Data structure tracking MRI machine state
    busy: bool = False                                                          # True if machine is currently scanning
    current_patient: Patient = None                                             # Patient currently being scanned
    schedule: List[Tuple[float, Patient]] = None                                # List of scheduled appointments (time, patient)
    def __post_init__(self):                                                    # Called after dataclass initialization
        if self.schedule is None:                                               # Initialize empty schedule if not provided
            self.schedule = []                                                  # Mutable default requires this pattern

class MRISimulation:                                                            # Main discrete event simulation class
    def __init__(self, timeslot_type1, timeslot_type2, policy="separate", simulation_days=20, random_seed=None):
        self.timeslot_type1 = timeslot_type1                                    # Scheduled timeslot length for Type 1 (minutes)
        self.timeslot_type2 = timeslot_type2                                    # Scheduled timeslot length for Type 2 (minutes)
        self.policy = policy                                                    # "separate" or "pooled" scheduling policy
        self.simulation_days = simulation_days                                  # Number of working days to simulate
        if random_seed is not None:                                             # Set seed for reproducibility
            np.random.seed(random_seed)                                         # Enables replication of results
        self.current_time = 0.0                                                 # Current simulation clock
        self.current_day = 0                                                    # Current day in simulation
        self.patient_counter = 0                                                # Counter for generating unique patient IDs
        self.event_list = []                                                    # Future Event List (FEL) as min-heap
        self.machines = {1: MachineState(), 2: MachineState()}                  # Two MRI machines
        self.patients = []                                                      # List of all patients in simulation
        self.stats = defaultdict(list)                                          # Dictionary to store performance statistics
        
    def generate_scan_duration(self, patient_type):                             # Generate random scan duration from normal distribution
        if patient_type == 1:                                                   # Type 1 patient
            return max(5.0, np.random.normal(TYPE1_DURATION_MEAN, TYPE1_DURATION_SD))   # Minimum 5 minutes
        else:                                                                   # Type 2 patient
            return max(10.0, np.random.normal(TYPE2_DURATION_MEAN, TYPE2_DURATION_SD))  # Minimum 10 minutes
    
    def generate_interarrival_time(self, patient_type):                         # Generate time until next arrival (exponential)
        if patient_type == 1:                                                   # Type 1 patient
            return np.random.exponential(1.0 / LAMBDA_TYPE1)                    # Exponential with mean = 1/lambda hours
        else:                                                                   # Type 2 patient
            return np.random.exponential(1.0 / LAMBDA_TYPE2)                    # Exponential with mean = 1/lambda hours
    
    def schedule_event(self, event):                                            # Add event to Future Event List
        heapq.heappush(self.event_list, event)                                  # Maintains heap property (sorted by time)
    
    def get_timeslot_length(self, patient_type):                                # Get timeslot duration for patient type
        return self.timeslot_type1 if patient_type == 1 else self.timeslot_type2  # Return appropriate timeslot
    
    def time_to_day_hour(self, time):                                           # Convert simulation time to (day, hour)
        day = int(time / WORK_HOURS)                                            # Day number (0, 1, 2, ...)
        hour = WORK_START + (time % WORK_HOURS)                                 # Clock hour (8.0 to 17.0)
        return day, hour                                                        # Return tuple of day and hour
    
    def find_next_available_slot(self, patient):                                # Find earliest appointment slot for patient
        call_day, call_hour = self.time_to_day_hour(patient.call_time)          # Get day and hour of patient's call
        timeslot_hours = self.get_timeslot_length(patient.patient_type) / 60.0  # Convert timeslot from minutes to hours
        search_day = call_day + 1                                               # Cannot schedule same day (case rule)
        machines_to_check = [patient.patient_type] if self.policy == "separate" else [1, 2]  # Which machines to check
        
        for day in range(search_day, search_day + 30):                          # Search up to 30 days ahead
            for machine_id in machines_to_check:                                # Check each relevant machine
                day_schedule = [(t, p) for t, p in self.machines[machine_id].schedule if self.time_to_day_hour(t)[0] == day]  # Get this day's appointments
                day_schedule.sort()                                             # Sort by appointment time
                candidate_time = day * WORK_HOURS                               # Start of day in simulation time
                
                if not day_schedule:                                            # No appointments on this day
                    return day, candidate_time, machine_id                      # Schedule at start of day
                
                for i, (start_time, p) in enumerate(day_schedule):              # Check each existing appointment
                    slot_end = start_time + self.get_timeslot_length(p.patient_type) / 60.0  # When this slot ends
                    if i == 0 and candidate_time + timeslot_hours <= start_time:  # Gap before first appointment?
                        return day, candidate_time, machine_id                  # Schedule at start of day
                    if i < len(day_schedule) - 1:                               # Not the last appointment
                        next_start = day_schedule[i + 1][0]                     # Next appointment's start time
                        if slot_end + timeslot_hours <= next_start:             # Gap between appointments?
                            return day, slot_end, machine_id                    # Schedule in the gap
                    else:                                                       # After last appointment
                        day_end = day * WORK_HOURS + WORK_HOURS                 # End of working day
                        if slot_end + timeslot_hours <= day_end:                # Fits before closing?
                            return day, slot_end, machine_id                    # Schedule after last appointment
        return None, None, None                                                 # No slot found (shouldn't happen)
    
    def handle_arrival(self, event):                                            # Process patient arrival (call) event
        patient = event.patient                                                 # Get patient from event
        sched_day, sched_time, machine_id = self.find_next_available_slot(patient)  # Find appointment slot
        if sched_time is None:                                                  # No slot found (edge case)
            return                                                              # Skip this patient
        patient.scheduled_time = sched_time                                     # Record scheduled time
        patient.scheduled_day = sched_day                                       # Record scheduled day
        patient.assigned_machine = machine_id                                   # Record assigned machine
        self.machines[machine_id].schedule.append((sched_time, patient))        # Add to machine's schedule
        self.schedule_event(Event(time=sched_time, event_type=EventType.SCAN_START, patient=patient, machine=machine_id))  # Schedule scan start
        next_arrival = event.time + self.generate_interarrival_time(patient.patient_type)  # Generate next arrival time
        next_day, _ = self.time_to_day_hour(next_arrival)                       # Get day of next arrival
        if next_day < self.simulation_days:                                     # If within simulation period
            new_patient = Patient(id=self.patient_counter, patient_type=patient.patient_type, call_time=next_arrival, call_day=next_day, scan_duration=self.generate_scan_duration(patient.patient_type))  # Create new patient
            self.patient_counter += 1                                           # Increment patient counter
            self.patients.append(new_patient)                                   # Add to patient list
            self.schedule_event(Event(time=next_arrival, event_type=EventType.ARRIVAL, patient=new_patient))  # Schedule their arrival
    
    def handle_scan_start(self, event):                                         # Process scan start event
        machine = self.machines[event.machine]                                  # Get machine state
        machine.busy = True                                                     # Mark machine as busy
        machine.current_patient = event.patient                                 # Record current patient
        event.patient.scan_start_time = event.time                              # Record actual start time
        scan_end_time = event.time + event.patient.scan_duration / 60.0         # Calculate end time (duration in hours)
        event.patient.scan_end_time = scan_end_time                             # Record end time
        self.schedule_event(Event(time=scan_end_time, event_type=EventType.SCAN_END, patient=event.patient, machine=event.machine))  # Schedule scan end
    
    def handle_scan_end(self, event):                                           # Process scan end event
        machine = self.machines[event.machine]                                  # Get machine state
        machine.busy = False                                                    # Mark machine as available
        machine.current_patient = None                                          # Clear current patient
    
    def run(self):                                                              # Main simulation loop
        for patient_type in [1, 2]:                                             # Initialize first arrival for each type
            first_arrival = self.generate_interarrival_time(patient_type)       # Generate first arrival time
            patient = Patient(id=self.patient_counter, patient_type=patient_type, call_time=first_arrival, call_day=0, scan_duration=self.generate_scan_duration(patient_type))  # Create first patient
            self.patient_counter += 1                                           # Increment counter
            self.patients.append(patient)                                       # Add to patient list
            self.schedule_event(Event(time=first_arrival, event_type=EventType.ARRIVAL, patient=patient))  # Schedule arrival
        
        while self.event_list:                                                  # Process events until FEL is empty
            event = heapq.heappop(self.event_list)                              # Get next event (earliest)
            self.current_time = event.time                                      # Advance simulation clock
            day, _ = self.time_to_day_hour(event.time)                          # Get current day
            if day >= self.simulation_days:                                     # Past simulation horizon
                break                                                           # End simulation
            if event.event_type == EventType.ARRIVAL:                           # Arrival event
                self.handle_arrival(event)                                      # Handle it
            elif event.event_type == EventType.SCAN_START:                      # Scan start event
                self.handle_scan_start(event)                                   # Handle it
            elif event.event_type == EventType.SCAN_END:                        # Scan end event
                self.handle_scan_end(event)                                     # Handle it
        self.calculate_statistics()                                             # Compute performance metrics
    
    def calculate_statistics(self):                                             # Calculate performance metrics
        scheduled = [p for p in self.patients if p.scheduled_time is not None]  # Filter to scheduled patients
        if not scheduled:                                                       # No patients scheduled
            return                                                              # Nothing to calculate
        
        # --- WAITING TIME (in working days) ---
        waiting_days = [(p.scheduled_time - p.call_time) / WORK_HOURS for p in scheduled]  # Convert sim time difference to days
        self.stats['mean_waiting_days'] = np.mean(waiting_days)                 # Average waiting time
        
        # --- OVERTIME (in minutes past 17:00) ---
        overtime_minutes = []                                                   # List to store overtime occurrences
        for p in scheduled:                                                     # For each scheduled patient
            if p.scan_end_time is not None:                                     # If scan completed
                day, end_hour = self.time_to_day_hour(p.scan_end_time)          # Get end day and hour
                if end_hour > (WORK_END + 0.001):                               # Past 17:00 (with tolerance)?
                    overtime_minutes.append((end_hour - WORK_END) * 60)         # Calculate overtime in minutes
        
        total_ot_min = sum(overtime_minutes)                                    # Total overtime minutes
        self.stats['avg_daily_overtime_min'] = total_ot_min / self.simulation_days  # Average overtime per day
        self.stats['pct_overtime_scans'] = (len(overtime_minutes) / len(scheduled)) * 100 if scheduled else 0  # % scans with overtime
        
        # --- UTILIZATION ---
        total_scan = sum([p.scan_duration for p in scheduled if p.scan_end_time is not None])  # Total scan minutes
        available = self.simulation_days * WORK_HOURS * 60 * 2                  # Total available minutes (2 machines)
        self.stats['utilization_pct'] = (total_scan / available) * 100          # Utilization percentage
        self.stats['patients_served'] = len(scheduled)                          # Total patients served

def run_experiment(ts1, ts2, policy, n_reps=50):                                # Run multiple replications
    results = []                                                                # Store results from each replication
    for i in range(n_reps):                                                     # For each replication
        sim = MRISimulation(ts1, ts2, policy, simulation_days=20, random_seed=i)  # Create simulation with unique seed
        sim.run()                                                               # Run simulation
        results.append(sim.stats)                                               # Collect statistics
    agg = {}                                                                    # Aggregated results
    for key in results[0].keys():                                               # For each metric
        vals = [r[key] for r in results]                                        # Get values from all replications
        agg[f'{key}_mean'] = np.mean(vals)                                      # Calculate mean
        agg[f'{key}_std'] = np.std(vals)                                        # Calculate standard deviation
    return agg                                                                  # Return aggregated statistics

def calculate_cost_score(results):                                              # Calculate weighted cost score
    w_wait = 10.0                                                               # Weight: cost per day of waiting
    w_overtime = 5.0                                                            # Weight: cost per minute of overtime
    w_idle = 1.0                                                                # Weight: cost per % idle time

    wait_val = results['mean_waiting_days_mean']                                # Extract mean waiting time
    ot_val = results['avg_daily_overtime_min_mean']                             # Extract mean overtime
    util_val = results['utilization_pct_mean']                                  # Extract mean utilization
    
    score_wait = w_wait * wait_val                                              # Waiting time cost component
    score_ot = w_overtime * ot_val                                              # Overtime cost component
    score_idle = w_idle * (100 - util_val)                                      # Idle time cost component

    total_score = score_wait + score_ot + score_idle                            # Total weighted cost
    return total_score                                                          # Lower is better

def grid_search(type1_range, type2_range, policy, n_reps=30):                   # Search over timeslot combinations
    results = []                                                                # Store all results
    total = len(type1_range) * len(type2_range)                                 # Total combinations to test
    num = 0                                                                     # Counter for progress display
    for ts1 in type1_range:                                                     # For each Type 1 timeslot
        for ts2 in type2_range:                                                 # For each Type 2 timeslot
            num += 1                                                            # Increment counter
            print(f"  Testing {num}/{total}: Type1={ts1}min, Type2={ts2}min", end="\r")  # Progress indicator
            
            sim_res = run_experiment(ts1, ts2, policy, n_reps)                  # Run experiment
            final_score = calculate_cost_score(sim_res)                         # Calculate cost score
            
            results.append({                                                    # Store results
                'ts1': ts1,                                                     # Type 1 timeslot
                'ts2': ts2,                                                     # Type 2 timeslot
                'policy': policy,                                               # Policy used
                'score': final_score,                                           # Cost score
                'waiting_days': sim_res['mean_waiting_days_mean'],              # Mean waiting days
                'overtime_min': sim_res['avg_daily_overtime_min_mean'],         # Mean overtime min/day
                'utilization_pct': sim_res['utilization_pct_mean']              # Mean utilization %
            })
    print()                                                                     # Newline after progress
    return pd.DataFrame(results)                                                # Return as DataFrame

if __name__ == "__main__":                                                      # Main execution block
    print("=" * 70)
    print("FINDING OPTIMAL TIMESLOTS (Weighted Cost Function)")
    print("=" * 70)
    
    print("\nSEPARATE policy:")
    sep_grid = grid_search(range(28, 30, 1), range(45, 46, 1), "separate", n_reps=30)  # Search Type1: 26-31, Type2: 41-46
    best_sep = sep_grid.loc[sep_grid['score'].idxmin()]                         # Find best (lowest score)
    
    print("\nPOOLED policy:")
    pool_grid = grid_search(range(28, 30, 1), range(45, 46, 1), "pooled", n_reps=30)  # Same range for pooled
    best_pool = pool_grid.loc[pool_grid['score'].idxmin()]                      # Find best (lowest score)
    
    print("\n" + "=" * 70)
    print("OPTIMAL TIMESLOTS")
    print("=" * 70)
    print(f"\nSEPARATE: Type1={best_sep['ts1']:.0f}min, Type2={best_sep['ts2']:.0f}min")
    print(f"Metrics:  Wait={best_sep['waiting_days']:.2f} days")
    print(f"          Overtime={best_sep['overtime_min']:.1f} min/day")
    print(f"          Util={best_sep['utilization_pct']:.1f}%")
    print(f"          Total Cost Score={best_sep['score']:.2f}")
    
    print(f"\nPOOLED:   Type1={best_pool['ts1']:.0f}min, Type2={best_pool['ts2']:.0f}min")
    print(f"Metrics:  Wait={best_pool['waiting_days']:.2f} days")
    print(f"          Overtime={best_pool['overtime_min']:.1f} min/day")
    print(f"          Util={best_pool['utilization_pct']:.1f}%")
    print(f"          Total Cost Score={best_pool['score']:.2f}")
    
    improvement = ((best_sep['score'] - best_pool['score']) / best_sep['score']) * 100  # Calculate % improvement
    print(f"\n{'POOLED RECOMMENDED' if best_pool['score'] < best_sep['score'] else 'SEPARATE RECOMMENDED'} (Cost reduction: {abs(improvement):.1f}%)")
    print("=" * 70)
    
    sep_grid.to_csv('optimal_timeslots_separate.csv', index=False)              # Save separate results
    pool_grid.to_csv('optimal_timeslots_pooled.csv', index=False)               # Save pooled results
    print("\nResults saved to CSV files.")