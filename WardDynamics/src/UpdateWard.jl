module UpdateWard

#using Main.InitialiseWard
import Main.Configuration
import Main.EventEntityInterface
import Main.Events
import Main.Entities
import Main.PatientFunctions
import Main.NurseFunctions 
import Main.InitialiseWard 


# simple objects
#unfinished_tasks = Array{Events.Event_T, 1}()


function find_available_beds(beds::Dict{Int64, Entities.Bed})::Array{Int64, 1} # returns ids of bed
    available_bed_ids = Array{Int64, 1}()
    # sort for reproducibility
    for i in sort(collect(keys(beds)))
        bed_i = beds[i]
        if !bed_i.is_occupied
            push!(available_bed_ids, bed_i.id)
        end
    end

    # sort ids to control random assignment: 
    #sort!(available_bed_ids)

    return available_bed_ids
end

function admit_next_patient!(w::InitialiseWard.Ward,
                             c::Configuration.Config,
                             day_index::Int64,
                             shift_index::Int64)::Int64

    patient_id = 0
    available_bed_ids = find_available_beds(w.beds)
    if !isempty(available_bed_ids)
        # find next patient id. #NOTE: this may become redundant, but for now useful. 
        ids = collect(keys(w.admitted_patients)) #NOTE: slightly unstable if ids are not the same as keys. 
        if isempty(ids)
            patient_id = 1
        else
            patient_id = maximum(ids) + 1
        end 

        p = PatientFunctions.create_patient(patient_id)

        # find an available bed 

        # pick a random bed from those available. 
        #TODO this may need to be updated if there is a selection process,
        # but I think we'll be mostly dealing with full wards so the options
        # will be limited. 

        # the lines below admit a patient and initialise their service schedule
        # these can go in 'UpdateWard' 
    
        #rand_index = rand(c.rng_trajectories, 1:size(available_bed_ids, 1))

        #bed_id = available_bed_ids[rand_index]

        bed_id = rand(c.rng_trajectories, available_bed_ids)

        # admit the patient 
        #NOTE: this adds the patient to the admitted_patients.all dict with ID as key. 
        LoS_i = Int64(ceil(rand(c.rng_LoS, c.length_of_stay)))
        println(LoS_i)
        PatientFunctions.admit_patient!(p, 
                                        w.admitted_patients, 
                                        bed_id,
                                        day_index,
                                        shift_index,
                                        LoS_i)

        w.beds[bed_id].is_occupied = true

        # fill patient's service schedule. 
        PatientFunctions.fill_service_schedule!(p, 
                                                c.shifts_per_day,
                                                day_index,
                                                shift_index)
    end

    return patient_id

end


function queue_handover!(admitted_patients::Entities.PatientList_T, 
                        day::Int64, 
                        shift::Int64)

    for (id, p) in admitted_patients
        if !isempty(p.service_schedule)
            if p.service_schedule[(day, shift)][1].name != "Handover"
                PatientFunctions.add_handover_to_schedule!(p, day, shift)
            end
        end
    end

end


# compile a set of events that should happen. 

# cycle through the patients
# and create a queue of events that should happen. 
# NOTE: might be more efficient to pre-allocate and modify empty input. 

function queue_shift_duties(admitted_patients::Entities.PatientList_T,
                            day_index::Int64, 
                            shift_index::Int64)::Array{Events.Event, 1}


    # add handover to top of each patient's service schedule
    if (day_index, shift_index) != (1, 1)
        queue_handover!(admitted_patients, day_index, shift_index)
    end

    # sort for reproducibility
    patient_ids = sort(collect(keys(admitted_patients)))

    shift_duties = Array{Events.Event, 1}()

    n_queued = 1

    while n_queued > 0

        n_queued = 0

        # iterates over admitted patients in random (but fixed) order
        #for (id, p) in admitted_patients
        
        # to ensure handover happens before admission, 
        # iterate from smallest index to largest. 
        #NOTE: this introduces a repeated order of service. 
        for id in patient_ids
            p = admitted_patients[id]

            if !haskey(p.service_schedule, (day_index, shift_index))
                println("WARNING: patient $(p.id) has no service schedule.")
            end

            schedule = p.service_schedule[(day_index, shift_index)]

            if (!isempty(schedule))
                push!(shift_duties, popfirst!(schedule))
                n_queued += 1
            else
                n_queued += 0
            end

        end

    end

    return shift_duties
end

function queue_shift_duties_homogeneous(admitted_patients::Entities.PatientList_T,
                                        n_hcw::Int64,
                                        day_index::Int64, 
                                        shift_index::Int64)::Array{Events.Event, 1}

    patients = collect(values(admitted_patients))

    shift_duties = Array{Events.Event, 1}()

    homogeneous_mixing = Events.Event()
    Events.default_init_Event!(homogeneous_mixing, day_index, shift_index)
    Events.init_homogeneous_mixing!(homogeneous_mixing, n_hcw, patients)

    push!(shift_duties, homogeneous_mixing)

    return shift_duties # should return a one-element vector. 


end

#this gets invoked when a handover event finishes. 
function finish_handover!(finished_events::Array{Events.Event_T, 1},
                          beds::Dict{Int64, Entities.Bed})
    for e in finished_events
        if e.name == "Handover"
            # transfer assignment from last shift nurse to new nurse. 
            bed_id = e.entities["patient"][1].bed_id 
            bed = beds[bed_id]
            last_nurse_id = e.entities["nurse"][1].id
            new_nurse_id = e.entities["nurse"][2].id
            delete!(bed.hcw_ids, last_nurse_id)
            push!(bed.hcw_ids, new_nurse_id)
            # NOTE:  a patient_load state variable for hcws could be 
            # modified here, speeding up the clockout checks. 
        end
    end
end

function clockout_nurses!(nurses::Dict{Int64, <:Entities.HCW_T},
                          nurses_on_ward::Set{Int64}, 
                          ids_to_clockout::Array{Int64, 1})
    for i in ids_to_clockout
        NurseFunctions.clear_assigned_tasks!(nurses[i])
        delete!(nurses_on_ward, i)
    end
end

# this is a bit greedy (will check p_load after every bed)
# but it is nice structurally because it only needs the 
# finished events and beds as input. This could be made 
# faster by adding state variables to the nurse entities 
# and updating them after each handover to keep track of p_load.  
function queue_clockouts(finished_events::Array{Events.Event_T, 1},
                         bedblocks::Dict{Int64, Entities.BedBlock}, 
                         beds::Dict{Int64, Entities.Bed})::Array{Int64, 1}

    ids_to_clockout = Array{Int64, 1}()

    for e in finished_events 

        if e.name == "Handover"
            last_nurse = e.entities["nurse"][1]
            # check if nurse has handed over all patients: 
            p_load = NurseFunctions.patient_load(last_nurse,
                                                 bedblocks,
                                                 beds)
            if p_load == 0
                push!(ids_to_clockout, last_nurse.id)
            end
        end

    end

    return ids_to_clockout

end

###
# TODO: consider what to do with scheduled start times for requeued tasks 
function requeue_unfinished_tasks!(shift_queue::Array{<:Events.Event_T, 1},
                                   w::InitialiseWard.Ward,
                                   c::Configuration.Config,
                                   day_index::Int64,
                                   shift_index::Int64
                                   )::Array{<:Events.Event_T, 1}

    unfinished_events = Array{Events.Event_T, 1}()

    if shift_index == c.shifts_per_day
        day_index += 1
        shift_index = 1
    else
        shift_index += 1
    end

    # iterate through unfinished events (still in shift queue)
    reverse!(shift_queue)

    for e in shift_queue
        

        pushfirst!(unfinished_events, e)

        # put unfinished tasks back into service schedule
        # for the associated patient. 

        if e.requeue_flag 

            if haskey(e.entities,"patient")
                for i in e.entities_required["patient"] #should only have one entry. 
                    p_id = e.entities["patient"][i].id
                    if !haskey(w.admitted_patients, p_id)
                        println(e.name)
                    end
                    sched = w.admitted_patients[p_id].service_schedule
                    if haskey(sched, (day_index, shift_index))
                        pushfirst!(sched[(day_index, shift_index)], e)
                    else
                        sched[(day_index, shift_index)] = [e]
                    end
                end
            end

            #remove hcw assignments from unfinished tasks: 
            #TODO: generalise this to any hcw
            # and extend to any other information added to 
            # an event during allocation or dynamics. 
            if haskey(e.entities, "nurse")
                empty!(e.entities["nurse"])
            end

        end
            
        
    end

    empty!(shift_queue)

    return unfinished_events

end


function display_events(q::Array{<:Events.Event_T, 1})
    for event_i in q
        if !isempty(event_i.entities["nurse"])
            if haskey(event_i.entities, "patient")
                event_str = 
                    [rpad("day: $(event_i.day)", 10, " ") * 
                    rpad("shift: $(event_i.shift)", 15, " ") * 
                    rpad("patient: $(event_i.entities["patient"][1].id)", 15, " ") * 
                    rpad("task: $(event_i.name)", 20, " ") * 
                    rpad("nurse: $(getproperty.(event_i.entities["nurse"], :id))", 15, " ") *
                    rpad("started at: $(event_i.start_time)", 20, " ")*
                    rpad("duration: $(event_i.duration)", 20, " ")*  
                    rpad("priority: $(event_i.priority)", 20, " ")]
            else
                event_str = 
                    [rpad("day: $(event_i.day)", 10, " ") * 
                    rpad("shift: $(event_i.shift)", 15, " ") * 
                    rpad("patient: NA", 15, " ") * 
                    rpad("task: $(event_i.name)", 20, " ") * 
                    rpad("nurse: $(getproperty.(event_i.entities["nurse"], :id))", 15, " ") *
                    rpad("started at: $(event_i.start_time)", 20, " ")*
                    rpad("duration: $(event_i.duration)", 20, " ")*  
                    rpad("priority: $(event_i.priority)", 20, " ")]
            end

        else 
            event_str = 
                [rpad("day: $(event_i.day)", 10, " ") * 
                rpad("shift: $(event_i.shift)", 15, " ") * 
                rpad("patient: $(event_i.entities["patient"][1].id)", 15, " ") * 
                rpad("task: $(event_i.name)", 20, " ") * 
                rpad("started at: $(event_i.start_time)", 20, " ")*
                rpad("duration: $(event_i.duration)", 20, " ")*  
                rpad("priority: $(event_i.priority)", 20, " ")]
        end

        println(event_str)
    end
end


function assign_bedblocks_to_hcws!(hcws::Dict{Int64, <:Entities.HCW_T},
                                   shift::Entities.Shift,
                                   hcw_ids::Array{Int64, 1})

    for hcw_id in hcw_ids 
        hcw = hcws[hcw_id]
        NurseFunctions.assign_bedblocks!(hcw, shift)
    end


end

#### 2 methods:
function assign_HCWs_to_tasks!(q::Array{<:Events.Event_T, 1},
                               hcws::Dict{Int64, <:Entities.HCW_T},
                               shift::Entities.Shift,
                               beds::Dict{Int64, Entities.Bed})
    for event_i in q
        # identify patient(s)
        p_list_i = event_i.entities["patient"]
        for p in p_list_i # should only have one element [TODO: ensure each hcw id only gets added once if there are multiple patient ids associated. ]
            bed_p = beds[p.bed_id]
            bb_id_p = bed_p.block_id
            nurse_ids_bb = shift.bb_to_nurses[bb_id_p] 
            for nurse_id in nurse_ids_bb
                nurse = hcws[nurse_id]
                NurseFunctions.assign_nurse_to_event!(nurse, event_i)
            end
        end
    end

end

# multiple dispatch (this method uses two shift inputs and 
# is applicable to handover, or any other events requiring 
# hcws from different shifts.)
function assign_HCWs_to_tasks!(q::Array{<:Events.Event_T, 1},
                                hcws::Dict{Int64, <:Entities.HCW_T},
                                shift::Entities.Shift,
                                last_shift::Entities.Shift,
                                beds::Dict{Int64, Entities.Bed})
    for event_i in q
        
        # identify patient(s)

        p_list_i = event_i.entities["patient"]

        if event_i.name == "Handover"

            
            for p in p_list_i # should only have one element 
                bed_p = beds[p.bed_id]
                bb_id_p = bed_p.block_id
    
                nurse_ids_bb_s1 = shift.bb_to_nurses[bb_id_p] 
                nurse_ids_bb_s2 = last_shift.bb_to_nurses[bb_id_p]
    
                # NOTE: ORDER MATTERS [last_shift, this_shift]
                #last shift
                for nurse_id in nurse_ids_bb_s2
                    nurse = hcws[nurse_id]
                    NurseFunctions.assign_nurse_to_event!(nurse, event_i)
                end

                #this shift 
                for nurse_id in nurse_ids_bb_s1
                    nurse = hcws[nurse_id]
                    NurseFunctions.assign_nurse_to_event!(nurse, event_i)
                end
    
            end

            #println("$(event_i.name)")


        else

            for p in p_list_i # should only have one element 
                bed_p = beds[p.bed_id]
                bb_id_p = bed_p.block_id

                nurse_ids_bb = shift.bb_to_nurses[bb_id_p] 

                for nurse_id in nurse_ids_bb
                    nurse = hcws[nurse_id]
                    NurseFunctions.assign_nurse_to_event!(nurse, event_i)
                end
            end

        end
    end

end
####


function display_HCW_tasks(ids::Set{Int64}, 
                           hcws::Dict{Int64, <:Entities.HCW_T})
    # 
    println("*****")
    println("all tasks for HCWs scheduled for this shift:")
    for id in ids
        tasks = hcws[id].assigned_tasks
        for event_i in tasks 
            event_str = 
            [rpad("day: $(event_i.day)", 10, " ") * 
            rpad("shift: $(event_i.shift)", 15, " ") * 
            rpad("patient: $(event_i.entities["patient"][1].id)", 15, " ") * 
            rpad("nurse: $id", 15, " ") *
            rpad("task: $(event_i.name)", 20, " ") * 
            rpad("duration: $(event_i.duration)", 20, " ")*  
            rpad("priority: $(event_i.priority)", 20, " ")]

            println(event_str)
        end
    end
end

# produce the list of nurses that are not occupied. 
function find_available_hcw_ids(ids::Array{Int64, 1},
                                hcws::Dict{Int64, <:Entities.HCW_T})::Set{Int64}
    ids_available = Set{Int64}()
    for (id) in ids
        if !(hcws[id].is_busy)
            push!(ids_available, id)
        end
    end
    return ids_available 
end


function add_breaks_to_shift_queue!(t::Int64, 
                                   hcw_ids::Set{Int64},
                                   hcws::Dict{Int64, <:Entities.HCW_T},
                                   shift_queue::Array{<:Events.Event_T, 1})
    for id in hcw_ids 
        hcw = hcws[id]
        
        if !isempty(hcw.break_schedule)
            #println("breaks?")
            t_next_break = hcw.break_schedule[1].scheduled_start_time
            if t >= t_next_break 
                next_break = popfirst!(hcw.break_schedule)
                pushfirst!(shift_queue, next_break)
                
            end
        end
    end

end

# TODO: assesses each task to see if the conditions required for its initiation are 
# present. This can be extended to other types of events or tasks. 
function queue_next_events!(ids::Set{Int64},
                            shift_queue::Array{<:Events.Event_T, 1},
                            t_remaining::Int64,
                            t::Int64)::Array{Events.Event_T, 1}

    next_events = Array{Events.Event_T, 1}()

    indices_to_remove = Array{Int64, 1}()

    termination_flag = false 

    # iterate through tasks and try to assign them 
    #NOTE: termination flags can be used to stop the loop
    # eg. if all patients are being serviced or all staff are busy. 
    for i in 1:size(shift_queue, 1)

        next_task = shift_queue[i] #iterating through in priority order
    
        # flags below decide if the event can be queued:
        # 1) is there enough time left?
        time_flag = false # default 
        # 2) are the require hcws available?
        availability_flag = true # default 
        # 3) has the event's start time been reached?
        start_time_flag = false # default

        if next_task.duration <= t_remaining
            time_flag = true
            #println("not enough time")
        end

        # has the scheduled time been reached? 
        t_next_task = next_task.scheduled_start_time
        if t >= t_next_task
            start_time_flag = true 
        end


        required_hcws = next_task.entities["nurse"]

        for required_hcw in required_hcws
            if required_hcw.is_busy
                availability_flag = false
                #println("required hcws are not available")
            end
        end

        if (time_flag && availability_flag && start_time_flag)
            # assign task to available hcws 
            # flag hcws as busy 
            # remove hcw ids from set of available ids. 
            for hcw in required_hcws
                hcw.is_busy = true
                delete!(ids, hcw.id)
                NurseFunctions.assign_event_to_nurse!(hcw, next_task)
            end
            
            push!(next_events, next_task)
            push!(indices_to_remove, i)
        end

        if isempty(ids)
            termination_flag = true
            #println("no more available workers")
        end

        if termination_flag
            break
        end

    end

    deleteat!(shift_queue, indices_to_remove)

    return next_events
end


function queue_discharges(finished_events::Array{Events.Event_T, 1})::Array{Int64, 1}
    
    ids_to_discharge = Array{Int64, 1}()
    # produces array of patient ids to discharge from ward: 
    for e in finished_events
        if e.name == "Discharge"
            for p in e.entities["patient"]
                push!(ids_to_discharge, p.id)
            end
        end
    end

    return ids_to_discharge

end

function discharge_patients!(p_list::Entities.PatientList_T, 
                             ids_to_discharge::Array{Int64, 1},
                             beds::Dict{Int64, Entities.Bed})
    # remove patient from ward

        # make bed available


    # delete patient objects and clear records. 

    if !isempty(ids_to_discharge)

        for id in ids_to_discharge  
            beds[p_list[id].bed_id].is_occupied = false
            beds[p_list[id].bed_id].patient_id = -1 #as in default init. 

            PatientFunctions.discharge_patient!(id, p_list)
        end

    end

    # TODO: check that all references to patient are gone so it
    # can be completely deleted. 



end




### end module
end