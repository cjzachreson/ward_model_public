

module Entities

export Bed,
       BedTable,
       BedBlock,
       Room,
       HCW, # this is weird - should be Nurse... wonder how it compiled?
       Roster,
       Shift,
       Patient

using DataFrames
using Main.EventEntityInterface 

import Main.Events

abstract type Bed_T end

abstract type BedTable_T end 

#TODO: change "bedblock" this means something else in wards. 
abstract type BedBlock_T end 

abstract type Room_T end

abstract type Roster_T end

abstract type Shift_T end 

#abstract type Patient_T end #now defined in event-entity interface module  

#subtypes: 
abstract type Nurse_T <: HCW_T end
#TODO: add concrete types, initialisers, functions etc.
# for the agent types below: 
abstract type Doctor_T <: HCW_T end
abstract type AlliedHealth_T <: HCW_T end 
abstract type SpecialNurse_T <: HCW_T end
abstract type Cleaner_T <: HCW_T end
abstract type WardNurse_T <: HCW_T end 

# custom types: 
# (day, shift) -> list of events 
#TODO: include and enforce scheduled event times (events can occur later but not earlier than scheduled)
ServiceSchedule_T = Dict{Tuple{Int64, Int64}, Array{Events.Event, 1}}

#concrete types: 
mutable struct Bed <: Bed_T
    id::Int64
    room_id::Int64
    block_id::Int64
    is_occupied::Bool
    is_assigned::Bool
    patient_id::Int64
    hcw_ids::Set{Int64}
    #random constructor. default parameters assigned in separate function
    Bed() = new() 
end

function default_init_Bed!(bed::Bed_T)
    bed.id = -1
    bed.room_id = -1
    bed.block_id = -1
    bed.is_occupied =  false
    bed.is_assigned = false
    bed.patient_id = -1
    bed.hcw_ids = Set{Int64}()
end


# bed table is a data frame for keeping track of bed characteristics. 
mutable struct BedTable <: BedTable_T
    data::DataFrame
    # empty constructor
    BedTable() = new(DataFrame(bed_id = Int64[], 
                               room_id = Int64[], 
                               block_id = Int64[]))
end

# bed block (group of beds assigned together)
mutable struct BedBlock <: BedBlock_T
    id::Int64
    bed_ids::Array{Int64, 1}
    #random constructor 
    BedBlock() = new()
end

#BedBlock default initialiser 
function default_init_BedBlock!(bedblock::BedBlock_T)
    bedblock.id = -1
    bedblock.bed_ids = []
end

# Ward room
mutable struct Room <: Room_T
    id::Int64
    bed_ids::Array{Int64, 1}
    volume::Float64
    #random constructor
    Room() = new()
end

function default_init_Room!(room::Room_T)
    room.id = -1
    room.bed_ids = []
    room.volume = -1.0
end

# Healthcare worker Nurse
mutable struct Nurse <: Nurse_T
    id::Int64
    #TODO: update to container, may need multiple assignments on shift. 
    bedblock_id::Int64 
    worked_last_shift::Bool
    worked_shift_before_last::Bool
    n_consecutive_shifts::Int64 #number of consecutive days worked
    n_shifts_this_roster::Int64
    is_eligible::Bool

    # filled on the fly at start of each shift 
    assigned_tasks::Array{Events.Event, 1}

    #break schedule
    break_schedule::Array{Events.Event, 1}

    # dynamically-updated during shift: 
    is_busy::Bool 

    #random constructor
    Nurse() = new()
end

function default_init_Nurse!(hcw::Nurse_T)
    hcw.id = -1
    #TODO: may need to be vector if same blocks are used for day and night shift. 
    hcw.bedblock_id = -1 
    hcw.worked_last_shift = false
    hcw.worked_shift_before_last = false 
    hcw.n_consecutive_shifts = 0
    hcw.n_shifts_this_roster = 0
    hcw.is_eligible = true
    hcw.assigned_tasks = Array{Events.Event, 1}()
    hcw.break_schedule = Array{Events.Event, 1}()
    hcw.is_busy = false 
end


#Roster
mutable struct Roster <: Roster_T
    #assignments links the shift index (key) to the set of HCW 
    # ids working that shift (values)
    assignments::Dict{Int64, Array{Int64, 1}}
    Roster() = new()
end

function default_init_Roster!(roster::Roster_T)
    roster.assignments = Dict{Int64, Array{Int64, 1}}()
end


#Shift
mutable struct Shift <: Shift_T
    id::Int64
    time_of_day::String
    day_index::Int64
    roster_day::Int64
    # this will contain details of which hcw 
    # is assigned which bed blocks. 
    # want 2-way map:
    bb_to_nurses::Dict{Int64, Array{Int64, 1}}
    nurse_to_bbs::Dict{Int64, Array{Int64, 1}}

    Shift() = new()
end

function default_init_Shift!(shift::Shift_T)
    shift.id = -1
    shift.time_of_day = "None"
    shift.day_index = -1
    shift.roster_day = -1
    
    shift.bb_to_nurses = Dict{Int64, Array{Int64, 1}}()
    shift.nurse_to_bbs = Dict{Int64, Array{Int64, 1}}()
end



# Patient: 

mutable struct Patient <: Patient_T
    id::Int64
    bed_id::Int64
    length_of_stay::Int64 #length of stay in shifts. 

    # keeps track of where they 
    # are in the service schedule. 
    next_event_index::Int64

    # more agent variables can be included
    # such as needs level, meds needed, etc. if relevant 

    # the service schedule indicates the patients current 
    # and future needs. 
    # maps (day, shift) => array of events 
    service_schedule::ServiceSchedule_T
    #NOTE: events have names and durations 

    #random constructor. default parameters assigned in separate function
    Patient() = new() 
end

#custom type: 
PatientList_T = Dict{Int64, Entities.Patient}

function default_init_Patient!(patient::Patient_T)
    patient.id = -1
    patient.bed_id = -1
    patient.length_of_stay = -1
    patient.next_event_index = 0
    patient.service_schedule = ServiceSchedule_T()
end




#### end module
end