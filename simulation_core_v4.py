"""
simulation_core.py - Core DES Simulation Components
====================================================
Contains all building blocks for the MRI scheduling simulation:
- Configuration parameters (from Part 1 analysis)
- Data structures (Event, Patient, Machine)
- Random number generation
- DES simulation engine
"""

import numpy as np
import pandas as pd
import heapq
import os
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from enum import Enum

# CONFIGURATION PARAMETERS (from Part 1 Bootstrap Analysis)

# Scan duration parameters
TYPE1_DURATION_MEAN = 25.71                                                     # Mean scan duration Type 1 (minutes)
TYPE1_DURATION_STD = 5.84                                                       # Std deviation Type 1 (minutes)
TYPE2_DURATION_MEAN = 40.71                                                     # Mean scan duration Type 2 (minutes)
TYPE2_DURATION_STD = 12.07                                                      # Std deviation Type 2 (minutes)

# Arrival parameters
TYPE1_ARRIVALS_PER_DAY = 16.95                                                  # Mean daily arrivals Type 1 (Poisson λ)
TYPE2_ARRIVALS_PER_DAY = 9.95                                                   # Mean daily arrivals Type 2 (Poisson λ)

# Operating hours
WORK_START_HOUR = 8.0                                                           # Facility opens at 8:00
WORK_END_HOUR = 17.0                                                            # Facility closes at 17:00
WORK_HOURS_PER_DAY = WORK_END_HOUR - WORK_START_HOUR                            # 9 hours per day
WORK_MINUTES_PER_DAY = WORK_HOURS_PER_DAY * 60                                  # 540 minutes per day

# Time constants
MINUTES_PER_HOUR = 60
MINUTES_PER_DAY = 24 * 60                                                       # 1440 minutes

# Simulation defaults
DEFAULT_SIMULATION_DAYS = 30
DEFAULT_NUM_REPLICATIONS = 50

# Timeslot search ranges
TYPE1_TIMESLOT_RANGE = range(24, 36, 2)                                         # 24, 26, 28, 30, 32, 34
TYPE2_TIMESLOT_RANGE = range(40, 56, 2)                                         # 40, 42, 44, 46, 48, 50, 52, 54


# DATA STRUCTURES

class EventType(Enum):
    """Types of discrete events in the simulation."""
    PATIENT_CALL = 1                                                            # Patient calls to schedule
    DAY_START = 2                                                               # Working day begins
    SCAN_START = 3                                                              # MRI scan begins
    SCAN_END = 4                                                                # MRI scan completes
    DAY_END = 5                                                                 # Working day ends


@dataclass
class Patient:
    """Represents a patient's journey through the system."""
    patient_id: int                                                             # Unique identifier
    patient_type: int                                                           # 1 or 2
    call_time: float                                                            # When they called (absolute minutes)
    call_day: int                                                               # Day number of call
    scheduled_day: int = -1                                                     # Day of appointment
    scheduled_time: float = -1                                                  # Scheduled start time
    machine_id: int = -1                                                        # Assigned machine (0 or 1)
    scan_duration: float = 0.0                                                  # Actual scan duration
    scan_start_time: float = -1                                                 # Actual start time
    scan_end_time: float = -1                                                   # Actual end time
    
    def wait_time_days(self) -> float:
        """Calculate waiting time from call to appointment in days."""
        if self.scheduled_day >= 0:
            return self.scheduled_day - self.call_day
        return 0.0


@dataclass(order=True)
class Event:
    """A discrete event, ordered by time for the priority queue."""
    time: float                                                                 # Event time (sort key)
    event_type: EventType = field(compare=False)                                # Type of event
    patient: Optional[Patient] = field(compare=False, default=None)             # Associated patient
    machine_id: int = field(compare=False, default=-1)                          # Associated machine


@dataclass
class Machine:
    """Represents an MRI machine's state."""
    machine_id: int                                                             # 0 or 1
    is_busy: bool = False                                                       # Currently scanning?
    current_patient: Optional[Patient] = None                                   # Patient being scanned
    schedule: Dict[int, List[Patient]] = field(default_factory=dict)            # day -> patients
    queue: List[Patient] = field(default_factory=list)                          # Today's waiting patients
    total_scan_time_today: float = 0.0                                          # Minutes scanned today
    patients_served_today: int = 0                                              # Patients completed today
    
    def reset_daily_stats(self):
        """Reset at start of each day."""
        self.is_busy = False
        self.current_patient = None
        self.total_scan_time_today = 0.0
        self.patients_served_today = 0
        self.queue = []
    
    def get_scheduled_time_for_day(self, day: int, get_timeslot_func) -> float:
        """Get total scheduled minutes for a day."""
        if day not in self.schedule:
            return 0.0
        return sum(get_timeslot_func(p.patient_type) for p in self.schedule[day])
    
    def count_slots_for_day(self, day: int) -> int:
        """Count scheduled slots for a day."""
        return len(self.schedule.get(day, []))
    
    def add_to_schedule(self, day: int, patient: Patient):
        """Add patient to schedule."""
        if day not in self.schedule:
            self.schedule[day] = []
        self.schedule[day].append(patient)
    
    def load_queue_for_day(self, day: int):
        """Load scheduled patients into queue."""
        self.queue = list(self.schedule.get(day, []))


# RANDOM NUMBER GENERATOR

class RandomGenerator:
    """Generates random variates for the simulation."""
    
    def __init__(self, seed: int = 42):
        self.rng = np.random.default_rng(seed)
        self.type2_empirical = None
        self._load_empirical_data()
    
    def _load_empirical_data(self):
        """Try to load empirical Type 2 data from CSV."""
        paths = ['ScanRecords__1_.csv', '../ScanRecords__1_.csv', 
                 '/mnt/user-data/uploads/ScanRecords__1_.csv']
        for path in paths:
            try:
                if os.path.exists(path):
                    df = pd.read_csv(path)
                    type2 = df[df['PatientType'] == 'Type 2']['Duration']
                    self.type2_empirical = type2.values * 60                    # Hours to minutes
                    return
            except:
                continue
    
    def generate_interarrival_time(self, patient_type: int) -> float:
        """Generate exponential interarrival time (minutes)."""
        rate = TYPE1_ARRIVALS_PER_DAY if patient_type == 1 else TYPE2_ARRIVALS_PER_DAY
        minute_rate = rate / WORK_HOURS_PER_DAY / MINUTES_PER_HOUR
        return self.rng.exponential(1.0 / minute_rate)
    
    def generate_scan_duration(self, patient_type: int) -> float:
        """Generate scan duration (minutes)."""
        if patient_type == 1:
            duration = self.rng.normal(TYPE1_DURATION_MEAN, TYPE1_DURATION_STD)
            return max(5.0, duration)
        else:
            if self.type2_empirical is not None and len(self.type2_empirical) > 0:
                return self.rng.choice(self.type2_empirical)
            else:
                duration = self.rng.normal(TYPE2_DURATION_MEAN, TYPE2_DURATION_STD)
                return max(10.0, duration)


# DISCRETE EVENT SIMULATION ENGINE

class MRISimulationDES:
    """
    Discrete Event Simulation for MRI Scheduling.
    
    Uses a Future Event List (priority queue) to process events chronologically.
    Supports 'dedicated' (old) and 'pooled' (new) scheduling policies.
    """
    
    def __init__(self, timeslot_type1: int, timeslot_type2: int,
                 policy: str = 'dedicated', num_days: int = DEFAULT_SIMULATION_DAYS,
                 seed: int = 42):
        # Configuration
        self.timeslot_type1 = timeslot_type1
        self.timeslot_type2 = timeslot_type2
        self.policy = policy
        self.num_days = num_days
        
        # Random generator
        self.rng = RandomGenerator(seed)
        
        # Simulation clock
        self.current_time = 0.0
        self.current_day = 0
        
        # Future Event List (min-heap)
        self.FEL: List[Event] = []
        heapq.heapify(self.FEL)
        
        # System state
        self.machines = [Machine(machine_id=0), Machine(machine_id=1)]
        self.patient_counter = 0
        self.all_patients: List[Patient] = []
        self.completed_patients: List[Patient] = []
        
        # Statistics
        self.daily_overtime: List[float] = []
        self.daily_utilization: Dict[int, List[float]] = {0: [], 1: []}
        self.daily_patients_served: Dict[int, List[int]] = {0: [], 1: []}
    
    # Time utilities
    
    def get_timeslot(self, patient_type: int) -> int:
        return self.timeslot_type1 if patient_type == 1 else self.timeslot_type2
    
    def time_to_day(self, time: float) -> int:
        return int(time // MINUTES_PER_DAY)
    
    def day_to_time(self, day: int, minutes_from_midnight: float = 0) -> float:
        return day * MINUTES_PER_DAY + minutes_from_midnight
    
    def get_work_start(self, day: int) -> float:
        return self.day_to_time(day, WORK_START_HOUR * MINUTES_PER_HOUR)
    
    def get_work_end(self, day: int) -> float:
        return self.day_to_time(day, WORK_END_HOUR * MINUTES_PER_HOUR)
    
    # --- Event list management ---
    
    def schedule_event(self, event: Event):
        heapq.heappush(self.FEL, event)
    
    def get_next_event(self) -> Optional[Event]:
        return heapq.heappop(self.FEL) if self.FEL else None
    
    # --- Scheduling logic ---
    
    def find_earliest_slot_dedicated(self, patient: Patient) -> Tuple[int, int, float]:
        """DEDICATED: Type 1 → Machine 0, Type 2 → Machine 1."""
        machine_id = 0 if patient.patient_type == 1 else 1
        machine = self.machines[machine_id]
        timeslot = self.get_timeslot(patient.patient_type)
        slots_per_day = int(WORK_MINUTES_PER_DAY // timeslot)
        
        search_day = patient.call_day + 1
        while search_day < patient.call_day + 100:
            if machine.count_slots_for_day(search_day) < slots_per_day:
                slot_index = machine.count_slots_for_day(search_day)
                start_min = WORK_START_HOUR * MINUTES_PER_HOUR + slot_index * timeslot
                return (machine_id, search_day, self.day_to_time(search_day, start_min))
            search_day += 1
        return (machine_id, search_day, self.get_work_start(search_day))
    
    def find_earliest_slot_pooled(self, patient: Patient) -> Tuple[int, int, float]:
        """POOLED: Find earliest slot across both machines with load balancing."""
        timeslot = self.get_timeslot(patient.patient_type)
        search_day = patient.call_day + 1
        
        while search_day < patient.call_day + 100:
            best_machine, best_start = -1, float('inf')
            
            for mid in [0, 1]:
                machine = self.machines[mid]
                scheduled = machine.get_scheduled_time_for_day(search_day, self.get_timeslot)
                
                if scheduled + timeslot <= WORK_MINUTES_PER_DAY:
                    start = WORK_START_HOUR * MINUTES_PER_HOUR + scheduled
                    if start < best_start or (start == best_start and 
                        scheduled < self.machines[best_machine].get_scheduled_time_for_day(
                            search_day, self.get_timeslot)):
                        best_start, best_machine = start, mid
            
            if best_machine != -1:
                return (best_machine, search_day, self.day_to_time(search_day, best_start))
            search_day += 1
        
        return (0, search_day, self.get_work_start(search_day))
    
    def schedule_patient(self, patient: Patient):
        """Schedule patient based on policy."""
        if self.policy == 'pooled':
            mid, day, time = self.find_earliest_slot_pooled(patient)
        else:
            mid, day, time = self.find_earliest_slot_dedicated(patient)
        
        patient.machine_id, patient.scheduled_day, patient.scheduled_time = mid, day, time
        self.machines[mid].add_to_schedule(day, patient)
    
    # --- Event handlers ---
    
    def handle_patient_call(self, event: Event):
        """Handle patient call: schedule appointment and generate next arrival."""
        patient = event.patient
        self.schedule_patient(patient)
        
        # Generate next call
        interarrival = self.rng.generate_interarrival_time(patient.patient_type)
        next_time = self.current_time + interarrival
        next_day = self.time_to_day(next_time)
        time_in_day = next_time % MINUTES_PER_DAY
        
        # Adjust to working hours
        if time_in_day < WORK_START_HOUR * MINUTES_PER_HOUR:
            next_time = self.get_work_start(next_day)
        elif time_in_day >= WORK_END_HOUR * MINUTES_PER_HOUR:
            next_day += 1
            next_time = self.get_work_start(next_day)
        
        if next_day < self.num_days:
            self.patient_counter += 1
            new_patient = Patient(self.patient_counter, patient.patient_type, next_time, next_day)
            self.all_patients.append(new_patient)
            self.schedule_event(Event(next_time, EventType.PATIENT_CALL, new_patient))
    
    def handle_day_start(self, event: Event):
        """Handle day start: load queues, start scans, schedule day end."""
        day = self.time_to_day(self.current_time)
        
        for machine in self.machines:
            machine.reset_daily_stats()
            machine.load_queue_for_day(day)
            if machine.queue:
                self.start_next_scan(machine)
        
        self.schedule_event(Event(self.get_work_end(day), EventType.DAY_END))
    
    def start_next_scan(self, machine: Machine):
        """Start next scan on machine."""
        if not machine.queue:
            machine.is_busy = False
            return
        
        patient = machine.queue.pop(0)
        machine.is_busy = True
        machine.current_patient = patient
        patient.scan_duration = self.rng.generate_scan_duration(patient.patient_type)
        patient.scan_start_time = self.current_time
        
        self.schedule_event(Event(
            self.current_time + patient.scan_duration,
            EventType.SCAN_END, patient, machine.machine_id
        ))
    
    def handle_scan_end(self, event: Event):
        """Handle scan completion: update stats, start next patient."""
        patient = event.patient
        machine = self.machines[event.machine_id]
        
        patient.scan_end_time = self.current_time
        self.completed_patients.append(patient)
        machine.total_scan_time_today += patient.scan_duration
        machine.patients_served_today += 1
        machine.is_busy = False
        machine.current_patient = None
        
        if machine.queue:
            self.start_next_scan(machine)
    
    def handle_day_end(self, event: Event):
        """Handle day end: record stats, schedule next day."""
        day = self.time_to_day(self.current_time)
        total_overtime = 0.0
        
        for machine in self.machines:
            overtime = max(0, machine.total_scan_time_today - WORK_MINUTES_PER_DAY)
            total_overtime += overtime
            util = min(machine.total_scan_time_today / WORK_MINUTES_PER_DAY, 1.0)
            self.daily_utilization[machine.machine_id].append(util)
            self.daily_patients_served[machine.machine_id].append(machine.patients_served_today)
        
        self.daily_overtime.append(total_overtime)
        
        # Schedule next day
        next_day = day + 1
        if next_day < self.num_days + 10:
            has_patients = any(next_day in m.schedule and m.schedule[next_day] for m in self.machines)
            if has_patients or next_day < self.num_days:
                self.schedule_event(Event(self.get_work_start(next_day), EventType.DAY_START))
    
    # --- Main simulation ---
    
    def initialize(self):
        """Set up initial events."""
        for ptype in [1, 2]:
            self.patient_counter += 1
            patient = Patient(self.patient_counter, ptype, self.get_work_start(0), 0)
            self.all_patients.append(patient)
            self.schedule_event(Event(self.get_work_start(0), EventType.PATIENT_CALL, patient))
        
        self.schedule_event(Event(self.get_work_start(0), EventType.DAY_START))
    
    def run(self):
        """Execute the simulation."""
        self.initialize()
        
        while self.FEL:
            event = self.get_next_event()
            if event is None:
                break
            
            self.current_time = event.time
            self.current_day = self.time_to_day(self.current_time)
            
            if event.event_type == EventType.PATIENT_CALL:
                self.handle_patient_call(event)
            elif event.event_type == EventType.DAY_START:
                self.handle_day_start(event)
            elif event.event_type == EventType.SCAN_END:
                self.handle_scan_end(event)
            elif event.event_type == EventType.DAY_END:
                self.handle_day_end(event)
    
    def calculate_statistics(self) -> Dict:
        """Calculate summary statistics."""
        stats = {}
        
        # Overtime
        if self.daily_overtime:
            stats['avg_overtime_minutes'] = np.mean(self.daily_overtime)
            stats['proportion_overtime_days'] = np.mean([o > 0 for o in self.daily_overtime])
        else:
            stats['avg_overtime_minutes'] = 0
            stats['proportion_overtime_days'] = 0
        
        # Wait time
        wait_times = [p.wait_time_days() for p in self.completed_patients]
        stats['avg_wait_days'] = np.mean(wait_times) if wait_times else 0
        
        # Utilization
        all_util = self.daily_utilization[0] + self.daily_utilization[1]
        stats['avg_utilization'] = np.mean(all_util) if all_util else 0
        
        # Throughput
        if self.daily_patients_served[0] and self.daily_patients_served[1]:
            totals = [a + b for a, b in zip(self.daily_patients_served[0], self.daily_patients_served[1])]
            stats['avg_daily_patients'] = np.mean(totals)
        else:
            stats['avg_daily_patients'] = 0
        
        return stats
