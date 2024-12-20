module NurseFunctions

import Main.Entities 
import Main.Events



function patient_load(nurse::Entities.Nurse_T, 
                      bedblocks::Dict{Int64, Entities.BedBlock}, 
                      beds::Dict{Int64, Entities.Bed})::Int64
    bed_ids = bedblocks[nurse.bedblock_id].bed_ids
    n_assigned = 0
    for i in bed_ids
        bed_i = beds[i]
        if bed_i.is_occupied
            if in(nurse.id, bed_i.hcw_ids)
                n_assigned += 1
            end
        end
    end
    return n_assigned 
end


function assign_nurse_to_event!(nurse::Entities.Nurse_T, event::Events.Event_T)
    # add the hcw to the event's entities. 

    # sometimes events can involve multiple patients with the same hcw
    # in these cases we only want to add the nurse once. 
    if !in(nurse, event.entities["nurse"]) #doesn't happen in most cases
        push!(event.entities["nurse"], nurse)
    end
end

function assign_event_to_nurse!(nurse::Entities.Nurse_T, event::Events.Event_T)
    # add the hcw to the event's entities. 
    push!(nurse.assigned_tasks, event)
end


function assign_bedblocks!(nurse::Entities.Nurse_T, shift::Entities.Shift)
    nurse.bedblock_id = shift.nurse_to_bbs[nurse.id][1]
end

function remove_task!(nurse::Entities.Nurse_T, task::Events.Event_T)
    if in(task, nurse.assigned_tasks)
        #locate index of task
        #TODO: this is broken - deleting by value may delete tasks that are not unique. 
        #deleteat!(nurse.assigned_tasks, find(nurse.assigned_tasks .== task))
    else
        println("task not found")
    end

end

function clear_assigned_tasks!(nurse::Entities.Nurse_T)
    empty!(nurse.assigned_tasks)
end

# TODO: check break schedule make sure break schedule is initialised and cleared properly. 
#add breaks to schedule: 
function schedule_breaks!(nurse::Entities.Nurse_T, 
                          shift_duration::Int64, 
                          day::Int64, 
                          shift::Int64)

    # this is hard-coded for now, but can be 
    # extended and generalised as per requirements. 
    t_break_1 = Int64(shift_duration * (1/4))
    short_break_1 = Events.Event()
    Events.default_init_Event!(short_break_1, day, shift)
    Events.init_short_break!(short_break_1, nurse, t_break_1)
    
    t_break_2 = Int64(shift_duration * (2/4))
    long_break_1 = Events.Event()
    Events.default_init_Event!(long_break_1, day, shift)
    Events.init_long_break!(long_break_1, nurse, t_break_2)
    
    t_break_3 = Int64(shift_duration * (3/4))
    short_break_2 = Events.Event()
    Events.default_init_Event!(short_break_2, day, shift)
    Events.init_short_break!(short_break_2, nurse, t_break_3)
    
    push!(nurse.break_schedule, short_break_1)
    push!(nurse.break_schedule, long_break_1)
    push!(nurse.break_schedule, short_break_2)

end

function clear_break_schedule!(nurse::Entities.Nurse_T)
    empty!(nurse.break_schedule)
end




### end module 
end