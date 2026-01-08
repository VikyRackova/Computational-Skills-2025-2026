"""
MRI Facility Scheduling Simulation - Corrected Version
========================================================
This simulation models MRI facility scheduling for two patient types.
Key improvements over original:
1. Fixed time unit handling throughout
2. Proper queue management with cascading delays
3. Correct waiting time calculation
4. True pooled policy implementation (earliest slot across both machines)
5. Comprehensive performance metrics
"""

import numpy as np                                                              # Numerical operations and random number generation
import pandas as pd                                                             # Data handling and CSV export
from dataclasses import dataclass, field                                        # For clean data structures
from typing import List, Tuple, Optional, Dict                                  # Type hints for clarity
import heapq                                                                    # Priority queue for event list
from collections import defaultdict                                             # For collecting statistics

# === CONFIGURATION (from Part 1 statistical analysis) ===
TYPE1_DURATION_MEAN = 25.71                                                     # Mean scan duration Type 1 (minutes)
TYPE1_DURATION_SD = 5.84                                                        # Std dev scan duration Type 1 (minutes)
TYPE1_ARRIVALS_PER_DAY = 16.95                                                  # Expected daily arrivals Type 1

TYPE2_DURATION_MEAN = 40.71                                                     # Mean scan duration Type 2 (minutes)
TYPE2_DURATION_SD = 12.07                                                       # Std dev scan duration Type 2 (minutes)
TYPE2_ARRIVALS_PER_DAY = 9.95                                                   # Expected daily arrivals Type 2

WORK_START = 8.0                                                                # Facility opens at 8:00 (hours)
WORK_END = 17.0                                                                 # Facility closes at 17:00 (hours)
WORK_HOURS = WORK_END - WORK_START                                              # 9 working hours per day

# Arrival rates per hour (for exponential inter-arrival times)
LAMBDA_TYPE1 = TYPE1_ARRIVALS_PER_DAY / WORK_HOURS                              # ~1.88 arrivals/hour for Type 1
LAMBDA_TYPE2 = TYPE2_ARRIVALS_PER_DAY / WORK_HOURS                              # ~1.11 arrivals/hour for Type 2


class EventType:                                                                # Enum-like class for event types
    """Event types that can occur in the simulation."""
    ARRIVAL = "arrival"                                                         # Patient calls to make appointment
    SCAN_START = "scan_start"                                                   # Patient's scan begins
    SCAN_END = "scan_end"                                                       # Patient's scan completes


@dataclass
class Patient:                                                                  # Data class to store patient info
    """Represents a patient in the system."""
    id: int                                                                     # Unique patient identifier
    patient_type: int                                                           # 1 or 2 (determines scan characteristics)
    call_time: float                                                            # When patient called (simulation time units)
    call_day: int                                                               # Which day patient called
    scan_duration: float                                                        # Actual scan duration (minutes)
    scheduled_time: Optional[float] = None                                      # Scheduled appointment time (sim units)
    scheduled_day: Optional[int] = None                                         # Scheduled appointment day
    scan_start_time: Optional[float] = None                                     # Actual scan start (may differ due to delays)
    scan_end_time: Optional[float] = None                                       # Actual scan end time
    assigned_machine: Optional[int] = None                                      # Which MRI machine (1 or 2)


@dataclass
class Event:                                                                    # Data class for simulation events
    """Represents an event in the discrete event simulation."""
    time: float                                                                 # When event occurs (simulation time)
    event_type: str                                                             # Type of event (ARRIVAL, SCAN_START, SCAN_END)
    patient: Optional[Patient] = None                                           # Associated patient (if any)
    machine: Optional[int] = None                                               # Associated machine (if any)
    
    def __lt__(self, other):                                                    # Required for heapq ordering
        """Allow events to be sorted by time in priority queue."""
        return self.time < other.time                                           # Earlier events have priority


@dataclass
class MachineState:                                                             # Tracks state of each MRI machine
    """Tracks the state of an MRI machine."""
    busy: bool = False                                                          # Is machine currently scanning?
    current_patient: Optional[Patient] = None                                   # Patient being scanned (if any)
    busy_until: float = 0.0                                                     # Time when machine becomes free
    schedule: List[Tuple[float, Patient]] = field(default_factory=list)         # List of (time, patient) appointments
    total_scan_time: float = 0.0                                                # Accumulated scan time (for utilization)


class MRISimulation:                                                            # Main simulation class
    """
    Discrete Event Simulation for MRI facility scheduling.
    
    Implements two policies:
    - "separate": Machine 1 serves Type 1 patients, Machine 2 serves Type 2
    - "pooled": Both machines can serve both patient types
    
    Key modeling assumptions:
    - Patients call during working hours (8:00-17:00)
    - Earliest appointment is next working day
    - Timeslots are fixed length per patient type
    - Actual scan durations vary (normally distributed)
    - If scan exceeds timeslot, subsequent patients must wait (cascading delay)
    - Overtime occurs when scans extend past 17:00
    """
    
    def __init__(self, 
                 timeslot_type1: float,                                         # Timeslot length for Type 1 (minutes)
                 timeslot_type2: float,                                         # Timeslot length for Type 2 (minutes)
                 policy: str = "separate",                                      # "separate" or "pooled"
                 simulation_days: int = 20,                                     # Number of days to simulate
                 random_seed: Optional[int] = None):                            # For reproducibility
        
        self.timeslot_type1 = timeslot_type1                                    # Store timeslot for Type 1
        self.timeslot_type2 = timeslot_type2                                    # Store timeslot for Type 2
        self.policy = policy                                                    # Store scheduling policy
        self.simulation_days = simulation_days                                  # Store simulation length
        
        if random_seed is not None:                                             # Set random seed if provided
            np.random.seed(random_seed)                                         # Ensures reproducible results
        
        # Initialize simulation state
        self.current_time = 0.0                                                 # Current simulation time
        self.current_day = 0                                                    # Current simulation day
        self.patient_counter = 0                                                # For unique patient IDs
        self.event_list: List[Event] = []                                       # Future Event List (FEL)
        
        # Initialize two MRI machines
        self.machines: Dict[int, MachineState] = {                              # Dictionary of machine states
            1: MachineState(),                                                  # Machine 1
            2: MachineState()                                                   # Machine 2
        }
        
        self.patients: List[Patient] = []                                       # All patients in simulation
        self.stats: Dict = defaultdict(list)                                    # Collected statistics
        
    def generate_scan_duration(self, patient_type: int) -> float:               # Generate random scan duration
        """
        Generate scan duration from normal distribution.
        Applies minimum bounds to prevent unrealistic values.
        """
        if patient_type == 1:                                                   # Type 1 patient
            duration = np.random.normal(TYPE1_DURATION_MEAN, TYPE1_DURATION_SD) # Draw from N(25.71, 5.84)
            return max(5.0, duration)                                           # Minimum 5 minutes
        else:                                                                   # Type 2 patient
            duration = np.random.normal(TYPE2_DURATION_MEAN, TYPE2_DURATION_SD) # Draw from N(40.71, 12.07)
            return max(10.0, duration)                                          # Minimum 10 minutes
    
    def generate_interarrival_time(self, patient_type: int) -> float:           # Generate time between arrivals
        """
        Generate inter-arrival time from exponential distribution.
        Per case description: Poisson arrivals -> exponential inter-arrivals.
        """
        if patient_type == 1:                                                   # Type 1 patient
            return np.random.exponential(1.0 / LAMBDA_TYPE1)                    # Mean = 1/lambda hours
        else:                                                                   # Type 2 patient
            return np.random.exponential(1.0 / LAMBDA_TYPE2)                    # Mean = 1/lambda hours
    
    def schedule_event(self, event: Event) -> None:                             # Add event to FEL
        """Add event to the Future Event List (maintains heap property)."""
        heapq.heappush(self.event_list, event)                                  # Push onto min-heap by time
    
    def get_timeslot_length(self, patient_type: int) -> float:                  # Get timeslot for patient type
        """Return the timeslot length in minutes for given patient type."""
        return self.timeslot_type1 if patient_type == 1 else self.timeslot_type2
    
    def sim_time_to_day_hour(self, sim_time: float) -> Tuple[int, float]:       # Convert sim time to day/hour
        """
        Convert simulation time to (day, hour_of_day).
        
        Simulation time is measured in hours from start.
        Each working day has WORK_HOURS (9) hours.
        
        Example: sim_time = 10.5 -> day 1, hour 9.5 (1.5 hours into day 1)
                 But we add WORK_START, so hour = 8 + 1.5 = 9.5
        """
        day = int(sim_time / WORK_HOURS)                                        # Integer division gives day number
        hour_within_day = sim_time % WORK_HOURS                                 # Remainder is hours into that day
        actual_hour = WORK_START + hour_within_day                              # Add start time (8:00)
        return day, actual_hour                                                 # Return (day, clock hour)
    
    def day_hour_to_sim_time(self, day: int, hour: float) -> float:             # Convert day/hour to sim time
        """Convert (day, hour_of_day) to simulation time."""
        return day * WORK_HOURS + (hour - WORK_START)                           # day * 9 + (hour - 8)
    
    def get_day_end_sim_time(self, day: int) -> float:                          # Get sim time for end of day
        """Get simulation time for end of working day (17:00)."""
        return day * WORK_HOURS + WORK_HOURS                                    # day * 9 + 9 = next day start
    
    def find_next_available_slot(self, patient: Patient) -> Tuple[Optional[int], Optional[float], Optional[int]]:
        """
        Find the earliest available timeslot for a patient.
        
        Rules from case description:
        - Cannot schedule same day as call
        - Should ideally finish by 17:00 (but overtime is allowed)
        
        For SEPARATE policy: only check the machine for patient's type
        For POOLED policy: check both machines and pick earliest slot
        
        Returns: (scheduled_day, scheduled_sim_time, machine_id) or (None, None, None)
        """
        call_day, _ = self.sim_time_to_day_hour(patient.call_time)              # Get day patient called
        timeslot_hours = self.get_timeslot_length(patient.patient_type) / 60.0  # Convert minutes to hours
        search_start_day = call_day + 1                                         # Earliest is next day
        
        # Determine which machines to check based on policy
        if self.policy == "separate":                                           # Separate policy
            machines_to_check = [patient.patient_type]                          # Type 1 -> Machine 1, Type 2 -> Machine 2
        else:                                                                   # Pooled policy
            machines_to_check = [1, 2]                                          # Check both machines
        
        best_slot = None                                                        # Track best (earliest) slot found
        best_machine = None                                                     # Track which machine has best slot
        best_day = None                                                         # Track which day has best slot
        
        # Search up to 60 days ahead (should always find a slot)
        for day in range(search_start_day, search_start_day + 60):              # Loop through future days
            day_start = self.day_hour_to_sim_time(day, WORK_START)              # 8:00 in sim time
            day_end = self.day_hour_to_sim_time(day, WORK_END)                  # 17:00 in sim time
            
            for machine_id in machines_to_check:                                # Check each relevant machine
                # Get this machine's schedule for this day
                day_schedule = [                                                # Filter appointments for this day
                    (t, p) for t, p in self.machines[machine_id].schedule 
                    if self.sim_time_to_day_hour(t)[0] == day                   # Only appointments on this day
                ]
                day_schedule.sort(key=lambda x: x[0])                           # Sort by appointment time
                
                # Find gaps in schedule where new appointment could fit
                available_slots = []                                            # List of (start_time, end_time) gaps
                
                if not day_schedule:                                            # No appointments this day
                    available_slots.append((day_start, day_end))                # Entire day is available
                else:                                                           # Some appointments exist
                    # Check gap before first appointment
                    first_appt_time = day_schedule[0][0]                        # First appointment start
                    if day_start + timeslot_hours <= first_appt_time:           # Can fit before first?
                        available_slots.append((day_start, first_appt_time))    # Add this gap
                    
                    # Check gaps between appointments
                    for i in range(len(day_schedule) - 1):                      # Loop through adjacent pairs
                        current_end = day_schedule[i][0] + self.get_timeslot_length(day_schedule[i][1].patient_type) / 60.0
                        next_start = day_schedule[i + 1][0]                     # Next appointment start
                        if current_end + timeslot_hours <= next_start:          # Gap big enough?
                            available_slots.append((current_end, next_start))   # Add this gap
                    
                    # Check gap after last appointment
                    last_appt = day_schedule[-1]                                # Last appointment
                    last_end = last_appt[0] + self.get_timeslot_length(last_appt[1].patient_type) / 60.0
                    if last_end + timeslot_hours <= day_end:                    # Can fit after last?
                        available_slots.append((last_end, day_end))             # Add this gap
                    elif last_end < day_end:                                    # Partial slot (may cause overtime)
                        available_slots.append((last_end, day_end + 2))         # Allow overtime (up to 2 hours)
                
                # Find earliest slot that fits
                for slot_start, slot_end in available_slots:                    # Check each gap
                    if slot_start + timeslot_hours <= slot_end:                 # Timeslot fits in gap
                        # Check if this is better than current best
                        if best_slot is None or slot_start < best_slot:         # Found earlier slot
                            best_slot = slot_start                              # Update best slot
                            best_machine = machine_id                           # Update best machine
                            best_day = day                                      # Update best day
                
                # If found a slot today, can stop searching (for this machine)
                if best_slot is not None and best_day == day:                   # Found slot today
                    break                                                       # Move to next machine
            
            # If found slots on this day, don't search further days
            if best_day == day and best_slot is not None:                       # Found slot today
                break                                                           # Stop searching
        
        return best_day, best_slot, best_machine                                # Return best slot found
    
    def handle_arrival(self, event: Event) -> None:                             # Process patient arrival (call)
        """
        Handle patient arrival (call for appointment).
        
        1. Find earliest available slot
        2. Schedule the appointment
        3. Create scan_start event for appointment time
        4. Generate next arrival of same type
        """
        patient = event.patient                                                 # Get patient from event
        
        # Find and book appointment
        sched_day, sched_time, machine_id = self.find_next_available_slot(patient)
        
        if sched_time is None:                                                  # No slot found (shouldn't happen)
            return                                                              # Skip this patient
        
        # Record appointment details
        patient.scheduled_time = sched_time                                     # Store scheduled time
        patient.scheduled_day = sched_day                                       # Store scheduled day
        patient.assigned_machine = machine_id                                   # Store assigned machine
        
        # Add to machine's schedule
        self.machines[machine_id].schedule.append((sched_time, patient))        # Add to schedule list
        
        # Schedule scan start event
        self.schedule_event(Event(                                              # Create scan start event
            time=sched_time,                                                    # At scheduled appointment time
            event_type=EventType.SCAN_START,                                    # Type: scan start
            patient=patient,                                                    # For this patient
            machine=machine_id                                                  # On this machine
        ))
        
        # Generate next arrival of same patient type
        interarrival = self.generate_interarrival_time(patient.patient_type)    # Random inter-arrival time
        next_arrival_time = event.time + interarrival                           # Add to current time
        next_day, _ = self.sim_time_to_day_hour(next_arrival_time)              # Get day of next arrival
        
        # Handle day boundary crossing
        while next_arrival_time % WORK_HOURS > WORK_HOURS:                      # If past end of day
            next_arrival_time = (next_day + 1) * WORK_HOURS                     # Start of next day
            next_day += 1                                                       # Increment day
        
        if next_day < self.simulation_days:                                     # Still within simulation period
            new_patient = Patient(                                              # Create new patient
                id=self.patient_counter,                                        # Unique ID
                patient_type=patient.patient_type,                              # Same type as current
                call_time=next_arrival_time,                                    # When they call
                call_day=next_day,                                              # Which day they call
                scan_duration=self.generate_scan_duration(patient.patient_type) # Random scan duration
            )
            self.patient_counter += 1                                           # Increment counter
            self.patients.append(new_patient)                                   # Add to patient list
            
            self.schedule_event(Event(                                          # Create arrival event
                time=next_arrival_time,                                         # When patient calls
                event_type=EventType.ARRIVAL,                                   # Type: arrival
                patient=new_patient                                             # For new patient
            ))
    
    def handle_scan_start(self, event: Event) -> None:                          # Process scan start
        """
        Handle scan start event.
        
        Key insight: Actual start may be delayed if previous scan ran over.
        This models cascading delays from duration variability.
        """
        machine = self.machines[event.machine]                                  # Get machine state
        patient = event.patient                                                 # Get patient
        
        # Determine actual start time (may be delayed if machine still busy)
        actual_start = max(event.time, machine.busy_until)                      # Can't start until machine free
        
        # Update machine state
        machine.busy = True                                                     # Machine is now busy
        machine.current_patient = patient                                       # Store current patient
        
        # Record actual start time
        patient.scan_start_time = actual_start                                  # May differ from scheduled
        
        # Calculate scan end time using ACTUAL duration (not timeslot)
        scan_duration_hours = patient.scan_duration / 60.0                      # Convert minutes to hours
        scan_end_time = actual_start + scan_duration_hours                      # When scan will finish
        
        # Update machine availability
        machine.busy_until = scan_end_time                                      # Machine busy until scan ends
        patient.scan_end_time = scan_end_time                                   # Store end time
        
        # Track scan time for utilization
        machine.total_scan_time += patient.scan_duration                        # Add to machine's total
        
        # Schedule scan end event
        self.schedule_event(Event(                                              # Create scan end event
            time=scan_end_time,                                                 # When scan finishes
            event_type=EventType.SCAN_END,                                      # Type: scan end
            patient=patient,                                                    # For this patient
            machine=event.machine                                               # On this machine
        ))
    
    def handle_scan_end(self, event: Event) -> None:                            # Process scan end
        """Handle scan end event - machine becomes available."""
        machine = self.machines[event.machine]                                  # Get machine state
        machine.busy = False                                                    # Machine now idle
        machine.current_patient = None                                          # No patient being scanned
    
    def run(self) -> None:                                                      # Main simulation loop
        """
        Run the discrete event simulation.
        
        Algorithm:
        1. Initialize with first arrivals of each type
        2. While events remain and within time limit:
           a. Pop earliest event from FEL
           b. Update simulation clock
           c. Handle event (triggers new events)
        3. Calculate statistics
        """
        # Initialize: create first arrival of each patient type
        for patient_type in [1, 2]:                                             # Both patient types
            first_arrival = self.generate_interarrival_time(patient_type)       # Random first arrival
            
            patient = Patient(                                                  # Create first patient
                id=self.patient_counter,                                        # Unique ID
                patient_type=patient_type,                                      # Type 1 or 2
                call_time=first_arrival,                                        # When they call
                call_day=0,                                                     # Day 0
                scan_duration=self.generate_scan_duration(patient_type)         # Random duration
            )
            self.patient_counter += 1                                           # Increment counter
            self.patients.append(patient)                                       # Add to list
            
            self.schedule_event(Event(                                          # Create arrival event
                time=first_arrival,                                             # When patient calls
                event_type=EventType.ARRIVAL,                                   # Type: arrival
                patient=patient                                                 # For this patient
            ))
        
        # Main simulation loop
        while self.event_list:                                                  # While events remain
            event = heapq.heappop(self.event_list)                              # Get earliest event
            self.current_time = event.time                                      # Update clock
            
            day, _ = self.sim_time_to_day_hour(event.time)                      # Get current day
            
            # Stop if past simulation period (but allow ongoing scans to complete)
            if day >= self.simulation_days and event.event_type == EventType.ARRIVAL:
                continue                                                        # Skip new arrivals past end
            
            # Dispatch to appropriate handler
            if event.event_type == EventType.ARRIVAL:                           # Arrival event
                self.handle_arrival(event)                                      # Handle it
            elif event.event_type == EventType.SCAN_START:                      # Scan start event
                self.handle_scan_start(event)                                   # Handle it
            elif event.event_type == EventType.SCAN_END:                        # Scan end event
                self.handle_scan_end(event)                                     # Handle it
        
        self.calculate_statistics()                                             # Compute performance metrics
    
    def calculate_statistics(self) -> None:                                     # Compute performance metrics
        """
        Calculate comprehensive performance statistics.
        
        Metrics:
        1. Waiting time: Days between call and appointment
        2. Overtime: Minutes worked past 17:00
        3. Utilization: % of available time used for scanning
        4. Patient throughput
        """
        # Filter to completed patients
        scheduled = [p for p in self.patients if p.scheduled_time is not None]  # Have appointment
        completed = [p for p in scheduled if p.scan_end_time is not None]       # Completed scan
        
        if not completed:                                                       # No completed patients
            return                                                              # Nothing to compute
        
        # === WAITING TIME (in working days) ===
        waiting_days = []                                                       # Store waiting times
        for p in scheduled:                                                     # For each scheduled patient
            wait = (p.scheduled_time - p.call_time) / WORK_HOURS                # Convert sim time diff to days
            waiting_days.append(wait)                                           # Add to list
        
        self.stats['mean_waiting_days'] = np.mean(waiting_days)                 # Average wait
        self.stats['max_waiting_days'] = np.max(waiting_days)                   # Maximum wait
        self.stats['std_waiting_days'] = np.std(waiting_days)                   # Std dev of wait
        
        # === OVERTIME (in minutes) ===
        overtime_by_day = defaultdict(float)                                    # Track overtime per day
        for p in completed:                                                     # For each completed patient
            day, end_hour = self.sim_time_to_day_hour(p.scan_end_time)          # Get day and hour of end
            if end_hour > WORK_END:                                             # Past 17:00?
                overtime_minutes = (end_hour - WORK_END) * 60                   # Minutes past closing
                overtime_by_day[day] = max(overtime_by_day[day], overtime_minutes)  # Track max per day
        
        total_overtime = sum(overtime_by_day.values())                          # Total overtime
        days_with_overtime = len(overtime_by_day)                               # Count days with overtime
        
        self.stats['total_overtime_min'] = total_overtime                       # Total overtime minutes
        self.stats['avg_daily_overtime_min'] = total_overtime / self.simulation_days  # Average per day
        self.stats['days_with_overtime'] = days_with_overtime                   # Number of days with OT
        self.stats['pct_days_overtime'] = (days_with_overtime / self.simulation_days) * 100  # % of days
        
        # === UTILIZATION ===
        # For separate policy, each machine serves specific type
        # For pooled policy, both machines serve both types
        if self.policy == "separate":                                           # Separate policy
            # Machine 1 serves Type 1 only
            machine1_available = self.simulation_days * WORK_HOURS * 60         # Minutes available
            machine1_used = self.machines[1].total_scan_time                    # Minutes used
            
            # Machine 2 serves Type 2 only
            machine2_available = self.simulation_days * WORK_HOURS * 60         # Minutes available
            machine2_used = self.machines[2].total_scan_time                    # Minutes used
            
            self.stats['utilization_machine1_pct'] = (machine1_used / machine1_available) * 100
            self.stats['utilization_machine2_pct'] = (machine2_used / machine2_available) * 100
            self.stats['utilization_pct'] = ((machine1_used + machine2_used) / (2 * machine1_available)) * 100
            
        else:                                                                   # Pooled policy
            total_available = 2 * self.simulation_days * WORK_HOURS * 60        # Both machines
            total_used = sum(m.total_scan_time for m in self.machines.values()) # Total scan time
            self.stats['utilization_pct'] = (total_used / total_available) * 100
            self.stats['utilization_machine1_pct'] = (self.machines[1].total_scan_time / 
                                                       (self.simulation_days * WORK_HOURS * 60)) * 100
            self.stats['utilization_machine2_pct'] = (self.machines[2].total_scan_time / 
                                                       (self.simulation_days * WORK_HOURS * 60)) * 100
        
        # === PATIENT THROUGHPUT ===
        self.stats['patients_served'] = len(completed)                          # Total completed
        self.stats['patients_type1'] = len([p for p in completed if p.patient_type == 1])
        self.stats['patients_type2'] = len([p for p in completed if p.patient_type == 2])
        
        # === DELAY FROM SCHEDULE ===
        delays = []                                                             # Track actual vs scheduled start
        for p in completed:                                                     # For completed patients
            delay = (p.scan_start_time - p.scheduled_time) * 60                 # Delay in minutes
            delays.append(max(0, delay))                                        # Only positive delays
        
        self.stats['mean_delay_min'] = np.mean(delays)                          # Average delay
        self.stats['max_delay_min'] = np.max(delays)                            # Maximum delay


def run_experiment(ts1: float, ts2: float, policy: str, n_reps: int = 50) -> dict:
    """
    Run multiple replications of the simulation and aggregate results.
    
    Args:
        ts1: Timeslot length for Type 1 patients (minutes)
        ts2: Timeslot length for Type 2 patients (minutes)
        policy: "separate" or "pooled"
        n_reps: Number of replications
    
    Returns:
        Dictionary with mean and std of each metric
    """
    results = []                                                                # Store results from each rep
    
    for i in range(n_reps):                                                     # For each replication
        sim = MRISimulation(                                                    # Create simulation
            timeslot_type1=ts1,                                                 # Type 1 timeslot
            timeslot_type2=ts2,                                                 # Type 2 timeslot
            policy=policy,                                                      # Policy choice
            simulation_days=20,                                                 # Simulate 20 days
            random_seed=i                                                       # Different seed each rep
        )
        sim.run()                                                               # Run simulation
        results.append(sim.stats)                                               # Store statistics
    
    # Aggregate results
    agg = {}                                                                    # Aggregated statistics
    for key in results[0].keys():                                               # For each metric
        vals = [r[key] for r in results]                                        # Get all values
        agg[f'{key}_mean'] = np.mean(vals)                                      # Compute mean
        agg[f'{key}_std'] = np.std(vals)                                        # Compute std dev
    
    return agg                                                                  # Return aggregated results


def calculate_cost_score(results: dict) -> float:                               # Compute weighted cost
    """
    Calculate weighted cost score for optimization.
    
    The score balances three objectives:
    1. Minimize waiting time (patient satisfaction)
    2. Minimize overtime (staff cost, work-life balance)
    3. Maximize utilization (efficiency)
    
    Lower score is better.
    """
    # Weights for each objective (tune based on hospital priorities)
    w_wait = 10.0                                                               # Cost per day of waiting
    w_overtime = 5.0                                                            # Cost per minute of overtime
    w_idle = 1.0                                                                # Cost per % idle time
    
    # Extract metrics
    wait_val = results['mean_waiting_days_mean']                                # Average waiting time
    ot_val = results['avg_daily_overtime_min_mean']                             # Average daily overtime
    util_val = results['utilization_pct_mean']                                  # Utilization percentage
    
    # Calculate component scores
    score_wait = w_wait * wait_val                                              # Waiting time cost
    score_ot = w_overtime * ot_val                                              # Overtime cost
    score_idle = w_idle * (100 - util_val)                                      # Idle time cost
    
    return score_wait + score_ot + score_idle                                   # Total cost


def grid_search(type1_range: range, type2_range: range, policy: str, n_reps: int = 30) -> pd.DataFrame:
    """
    Search over grid of timeslot combinations to find optimal settings.
    
    Args:
        type1_range: Range of timeslots to try for Type 1
        type2_range: Range of timeslots to try for Type 2
        policy: "separate" or "pooled"
        n_reps: Replications per combination
    
    Returns:
        DataFrame with results for each combination
    """
    results = []                                                                # Store all results
    total = len(type1_range) * len(type2_range)                                 # Total combinations
    num = 0                                                                     # Counter for progress
    
    for ts1 in type1_range:                                                     # For each Type 1 timeslot
        for ts2 in type2_range:                                                 # For each Type 2 timeslot
            num += 1                                                            # Increment counter
            print(f"  Testing {num}/{total}: Type1={ts1}min, Type2={ts2}min", end="\r")
            
            sim_res = run_experiment(ts1, ts2, policy, n_reps)                  # Run experiment
            final_score = calculate_cost_score(sim_res)                         # Calculate cost
            
            results.append({                                                    # Store results
                'ts1': ts1,                                                     # Type 1 timeslot
                'ts2': ts2,                                                     # Type 2 timeslot
                'policy': policy,                                               # Policy
                'score': final_score,                                           # Cost score
                'waiting_days': sim_res['mean_waiting_days_mean'],              # Avg wait
                'max_waiting_days': sim_res['max_waiting_days_mean'],           # Max wait
                'overtime_min': sim_res['avg_daily_overtime_min_mean'],         # Avg overtime
                'pct_days_overtime': sim_res['pct_days_overtime_mean'],         # % days with OT
                'utilization_pct': sim_res['utilization_pct_mean'],             # Utilization
                'mean_delay_min': sim_res['mean_delay_min_mean'],               # Avg delay
                'patients_served': sim_res['patients_served_mean']              # Throughput
            })
    
    print()                                                                     # Newline after progress
    return pd.DataFrame(results)                                                # Return as DataFrame


if __name__ == "__main__":                                                      # Main execution block
    print("=" * 80)
    print("MRI FACILITY SCHEDULING OPTIMIZATION (Corrected Simulation)")
    print("=" * 80)
    
    # Search over reasonable range of timeslots
    # Type 1: mean=25.71, std=5.84 -> try 26-32 (mean to ~mean+1std)
    # Type 2: mean=40.71, std=12.07 -> try 42-50 (mean to ~mean+0.8std)
    
    print("\nTesting SEPARATE policy...")
    sep_grid = grid_search(range(25, 30, 1), range(40, 44, 1), "separate", n_reps=30)
    best_sep = sep_grid.loc[sep_grid['score'].idxmin()]
    
    print("\nTesting POOLED policy...")
    pool_grid = grid_search(range(25, 30, 1), range(40, 44, 1), "pooled", n_reps=30)
    best_pool = pool_grid.loc[pool_grid['score'].idxmin()]
    
    print("\n" + "=" * 80)
    print("OPTIMAL RESULTS")
    print("=" * 80)
    
    print(f"\n{'SEPARATE POLICY':^40}")
    print("-" * 40)
    print(f"Optimal timeslots: Type1={best_sep['ts1']:.0f}min, Type2={best_sep['ts2']:.0f}min")
    print(f"Mean waiting time: {best_sep['waiting_days']:.2f} days")
    print(f"Max waiting time:  {best_sep['max_waiting_days']:.2f} days")
    print(f"Daily overtime:    {best_sep['overtime_min']:.1f} min/day")
    print(f"Days with OT:      {best_sep['pct_days_overtime']:.1f}%")
    print(f"Utilization:       {best_sep['utilization_pct']:.1f}%")
    print(f"Mean delay:        {best_sep['mean_delay_min']:.1f} min")
    print(f"Cost score:        {best_sep['score']:.2f}")
    
    print(f"\n{'POOLED POLICY':^40}")
    print("-" * 40)
    print(f"Optimal timeslots: Type1={best_pool['ts1']:.0f}min, Type2={best_pool['ts2']:.0f}min")
    print(f"Mean waiting time: {best_pool['waiting_days']:.2f} days")
    print(f"Max waiting time:  {best_pool['max_waiting_days']:.2f} days")
    print(f"Daily overtime:    {best_pool['overtime_min']:.1f} min/day")
    print(f"Days with OT:      {best_pool['pct_days_overtime']:.1f}%")
    print(f"Utilization:       {best_pool['utilization_pct']:.1f}%")
    print(f"Mean delay:        {best_pool['mean_delay_min']:.1f} min")
    print(f"Cost score:        {best_pool['score']:.2f}")
    
    # Compare policies
    print(f"\n{'COMPARISON':^40}")
    print("-" * 40)
    improvement = ((best_sep['score'] - best_pool['score']) / best_sep['score']) * 100
    
    if best_pool['score'] < best_sep['score']:
        print(f"RECOMMENDATION: POOLED policy")
        print(f"Cost reduction vs separate: {improvement:.1f}%")
    else:
        print(f"RECOMMENDATION: SEPARATE policy")
        print(f"Cost reduction vs pooled: {-improvement:.1f}%")
    
    # Save results
    sep_grid.to_csv('optimal_timeslots_separate_corrected.csv', index=False)
    pool_grid.to_csv('optimal_timeslots_pooled_corrected.csv', index=False)
    print("\nResults saved to CSV files.")
    print("=" * 80)
