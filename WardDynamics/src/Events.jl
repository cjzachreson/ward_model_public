module Events

#setup entity type structure (these are specified in Entities modle)

# abstract type Entity_T end

# abstract type HCW_T <: Entity_T end
# abstract type Patient_T <: Entity_T end

#NOTE: all instances of Events.Entity_T now 
# refer to EventEntityInterface.Entity_T 
#TODO: disambiguate? 
using Main.EventEntityInterface 

abstract type Event_T end

#TODO add two quantities to Event struct:
# 1) time elapsed since event was queued (useful for setting priorities)
# 2) time elapsed since Event was initiated (for temporal implementation) 
 

mutable struct Event <: Event_T 
    name::String #name of the event 
    duration::Float64 #duration in minutes
    priority::String # a label corresponding to priority level
    #links the types of entities to the number required.  
    entities_required::Dict{String, Int64} 
    # records the set of entities that are actually
    # involved with the event. 
    entities::Dict{String, Array{Entity_T, 1}}

    # more properties can be included here
    # like the time since the event was queued 
    # or the time since the event was initiated.
    # this should be in range: 1:n_days_in_roster.  
    day::Int64
    # shift: This should have value 1:n_shifts_per_day
    shift::Int64
    start_time::Int64 #TODO: check data type. 
    finish_time::Int64 
    elapsed_time::Int64 

    scheduled_start_time::Int64 # for scheduled events only.
    
    requeue_flag::Bool # if this is true, the event will get 
    # requeued if it's left unfinished by the previous shift. 

    #random initialiser
    Event() = new()
end

#default initialiser 
function default_init_Event!(event::Event, 
                             day::Int64, 
                             shift::Int64)
    event.name = "None"
    event.duration = -1.0
    event.priority = "None"
    event.entities_required = Dict{String, Int64}()
    event.entities = Dict{String, Array{Entity_T, 1}}()

    event.day = day
    event.shift = shift
    event.start_time = -1
    event.finish_time = -1
    event.elapsed_time = -1

    # priority is just a string, so I'll needed
    #  a function that parses priority. 

    event.scheduled_start_time = -1

    event.requeue_flag = false #default to false. 
end




# using ints for time.... NOTE: may want to 
# switch to floats later for generality. 
function start_Event!(event::Event, t_start::Int64)
    # NOTE: time is relative to shift start time.
    # currently cannot account for events lasting 
    # longer than one shift - this could be handled
    # by passing an event to the next shift.

    event.start_time = t_start 
    # NOTE: does not update event. 
    # this must be called separately. 
    
    # all nurses involved are now busy: 
    # TODO: will need to extend this to all hcws or staff involved. 
    for nurse_i in event.entities["nurse"]
        nurse_i.is_busy = true 
    end

end

function finish_Event!(event::Event, t_finish::Int64)::Bool
    event.finish_time = t_finish
    # free-up any HCWs that were occupied by the event just finished: 

    # if event.name == "Handover"
    #     println("handing over")
    # end


    for nurse_i in event.entities["nurse"]
        nurse_i.is_busy = false 
    end

    is_finished = true 

    return is_finished

end

function step_Event!(event::Event, dt::Int64, current_time::Int64)::Bool

    is_finished = false 

    if event.start_time < 0 #event has not started yet
        println("WARNING: stepping an event ($(event.name)) that has not yet started")
    end
    
    event.elapsed_time += dt 

    if event.elapsed_time >= event.duration
       is_finished =  finish_Event!(event, current_time)
    end

    return is_finished

end




#note: these could be read and initialised from an 
# external table. TODO: once general structure of events
# is finalised, build initialisation from input table. 

# Patient events 


function init_homogeneous_mixing!(event::Event, 
                                  n_hcw::Int64, 
                                  patients::Array{<:Entity_T, 1})

    
    n_patients = size(patients, 1)

    event.name = "Homogeneous Mixing"
    event.duration = 480 #shift duration in minutes
    event.priority = "High"
    event.entities_required["patient"] = n_patients
    event.entities_required["nurse"] = n_hcw

    event.entities["patient"] = patients
    event.entities["nurse"] = Entity_T[]

    event.requeue_flag = true #not sure this will ever arise... 

end

function init_admission!(event::Event, patient::Entity_T)
    event.name = "Admission"
    event.duration = 60
    event.priority = "High"

    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["patient"] = 1
    event.entities_required["nurse"] = 1

    event.entities["patient"] = [patient]
    event.entities["nurse"] = Entity_T[]


    # entity set still incomplete at initialisation 
    # it will be filled when the event is 
    # triggered to occur. 
    event.requeue_flag = true #not sure this will ever arise... 

end

# added start time 
function init_vitals!(event::Event, patient::Entity_T, start_time::Int64)
    event.name = "Vitals"
    event.duration = 10
    event.priority = "Required"
    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["patient"] = 1
    event.entities_required["nurse"] = 1

    event.entities["patient"] = [patient]
    event.entities["nurse"] = Entity_T[]

    event.scheduled_start_time = start_time #start time in minutes from start of shift. 

end

# added start time. 
function init_meds!(event::Event, patient::Entity_T, start_time::Int64)
    event.name = "Medications"
    event.duration = 10
    event.priority = "High"
    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["patient"] = 1
    event.entities_required["nurse"] = 1

    event.entities["patient"] = [patient]
    event.entities["nurse"] = Entity_T[]

    event.scheduled_start_time = start_time # start time in minutes from start of shift 
end

function init_bath!(event::Event, patient::Entity_T)
    event.name = "Bath"
    event.duration = 30
    event.priority = "Required"
    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["patient"] = 1
    event.entities_required["nurse"] = 1

    event.entities["patient"] = [patient]
    event.entities["nurse"] = Entity_T[]

    event.requeue_flag = true #maybe? It's only once per day... 
end

function init_discharge!(event::Event, patient::Entity_T, start_time::Int64)
    event.name = "Discharge"
    event.duration = 60
    event.priority = "High"

    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["patient"] = 1
    event.entities_required["nurse"] = 1

    event.entities["patient"] = [patient]
    event.entities["nurse"] = Entity_T[]

    event.scheduled_start_time = start_time 

    event.requeue_flag = true #a must. 

end


# HCW events 
# TODO: HCW-initiated events. 

function init_long_break!(event::Event, 
                           nurse::Entity_T,
                           scheduled_start_time::Int64 )
    event.name = "Long Break"
    event.duration = 60
    event.priority = "Required"
    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["nurse"] = 1

    event.entities["nurse"] = [nurse]

    event.scheduled_start_time = scheduled_start_time

end

function init_short_break!(event::Event, 
                           nurse::Entity_T,
                           scheduled_start_time::Int64)
    event.name = "Short Break"
    event.duration = 10
    event.priority = "Required"
    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["nurse"] = 1

    event.entities["nurse"] = [nurse]

    event.scheduled_start_time = scheduled_start_time
end

function init_handover!(event::Event, 
                        patient::Patient_T)
    event.name = "Handover"
    event.duration = 5
    event.priority = "High"
    # types of entities involved and the 
    # number of each that are required:
    event.entities_required["patient"] = 1
    event.entities_required["nurse"] = 2

    event.entities["patient"] = [patient]
    # NOTE: order signifies which nurse is out and which is in 
    # [last_shift, this_shift]
    event.entities["nurse"] = []
end






### end module
end