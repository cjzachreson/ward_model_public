# contains functions related to patients. 

module PatientFunctions 

import Main.Entities 
import Main.Events


#test
function create_patient(id::Int64)::Entities.Patient
    #minimal property assignment (no bed assignment etc.)
    p = Entities.Patient()
    Entities.default_init_Patient!(p)
    p.id = id
    return p
end

# TODO: include and enforce scheduled event times (to prevent events happening too early)
# TODO add a function that checks to see if all scheduled events are scheduled in the correct order. (i.e., scheduled start times increase)
function add_event_to_service_schedule!(i::Tuple{Int64, Int64},
                                        event::Events.Event,
                                        schedule::Entities.ServiceSchedule_T)
    if haskey(schedule, i)
        push!(schedule[i], event)
    else
        schedule[i] = [event]
    end

end


function add_event_to_top_of_schedule!(i::Tuple{Int64, Int64},
                                       event::Events.Event,
                                       schedule::Entities.ServiceSchedule_T)
    if haskey(schedule, i)
        pushfirst!(schedule[i], event)
    else
        schedule[i] = [event]
    end

end

function admit_patient!(p::Entities.Patient, 
                       p_list::Entities.PatientList_T, 
                       bed_id::Int64,
                       day::Int64,
                       shift::Int64,
                       length_of_stay::Int64)

    p_list[p.id] = p

    #compute length of stay: 
    p.length_of_stay = length_of_stay #an integer number of shifts e.g., staying for 9 shifts (3 days)
    # TODO : the above is a placeholder, 
    # this will be replaced with a stochastic 
    # implementation requiring an rng input arg. 

    #link to assigned bed (bed assignment will be done elsewhere)
    p.bed_id = bed_id 
    #initialise event schedule with admission event at high priority. 
    #create admission event
    #initialise as blank
    admission = Events.Event()
    #update to default init
    Events.default_init_Event!(admission, day, shift)
    #update to admission event specifiers. 
    Events.init_admission!(admission, p)


    add_event_to_service_schedule!((day, shift),
                                    admission,
                                    p.service_schedule)

end

# TODO: include scheduled times for vitals. 
# TODO: schedules should be accessible as input tables, not hard-coded. 
function add_vitals_to_schedule!(p::Entities.Patient,
                                 day::Int64,
                                 shift::Int64,
                                 start_time::Int64) #start time in minutes. 
    
    #vitals
    vitals = Events.Event()
    Events.default_init_Event!(vitals, day, shift)
    Events.init_vitals!(vitals, p, start_time)
    add_event_to_service_schedule!((day, shift),
                                   vitals,
                                   p.service_schedule)
    
end


function add_meds_to_schedule!(p::Entities.Patient,
                               day::Int64,
                               shift::Int64,
                               start_time::Int64)
    
    #medications
    meds = Events.Event()
    Events.default_init_Event!(meds, day, shift)
    Events.init_meds!(meds, p, start_time)
    add_event_to_service_schedule!((day, shift),
                                   meds,
                                   p.service_schedule)
    
end


function add_bath_to_schedule!(p::Entities.Patient,
                               day::Int64,
                               shift::Int64)
    
    #bath
    bath = Events.Event()
    Events.default_init_Event!(bath, day, shift)
    Events.init_bath!(bath, p)
    add_event_to_service_schedule!((day, shift),
                                   bath,
                                   p.service_schedule)
    
end

function add_discharge_to_schedule!(p::Entities.Patient,
    day::Int64,
    shift::Int64,
    start_time::Int64)
    
    #discharge 
    discharge = Events.Event()
    Events.default_init_Event!(discharge, day, shift)
    Events.init_discharge!(discharge, p, start_time)
    add_event_to_service_schedule!((day, shift),
                                   discharge,
                                   p.service_schedule)
    
end

#NOTE: handover gets added at the beginning of each shift,
# not at admission. 
#TODO: consider scheduling all tasks day-wise or shift-wise 
function add_handover_to_schedule!(p::Entities.Patient,
                                   day::Int64,
                                   shift::Int64)

    #handover 
    handover = Events.Event()
    Events.default_init_Event!(handover, day, shift)
    Events.init_handover!(handover, p)
    add_event_to_top_of_schedule!((day, shift),
                                   handover,
                                   p.service_schedule)
    

end

#TODO: implement handover as a nurse-initiated task
# otherwise it gets superceded by unfinished tasks. 
function fill_service_schedule!(p::Entities.Patient, 
                                shifts_per_day::Int64,
                                day_index::Int64,
                                shift_index::Int64)
                                # remember integer args are not 
                                # mutable, so we don't need to
                                # worry about increasing them
                                #inside the function. 

    #for each day, iterate through shifts
    # for each shift, schedule vitals, and meds
    # ensure bath happens once per day 
    # we'll go from there...

    #determine length of stay (days)


    
    #n_days = ceil((p.length_of_stay + (shift_index - 1))/shifts_per_day)

    n_shifts = p.length_of_stay 

    #iterate through shifts 

    shift_index -= 1

    for i in 1:n_shifts


        if shift_index < shifts_per_day
            shift_index += 1
        else
            shift_index = 1
            day_index += 1
        end


        #if i != 1
         #   add_handover_to_schedule!(p, day_index, shift_index)
        #end

        # TODO: add and enforce scheduled start times 
        # TODO: these should be specified in input file, not hardcoded in module. 
        start_time_vitals_1 = 60
        start_time_meds = 120
        start_time_vitals_2 = 300

        if i < n_shifts

            #vitals 
            add_vitals_to_schedule!(p, day_index, shift_index, start_time_vitals_1)

            #meds
            add_meds_to_schedule!(p, day_index, shift_index, start_time_meds)

            #vitals 
            add_vitals_to_schedule!(p, day_index, shift_index, start_time_vitals_2)


            if shift_index == shifts_per_day
                # TODO: decide when baths should be scheduled, if at all
                add_bath_to_schedule!(p, day_index, shift_index) #for now, baths are remaining unscheduled 
            end

        

        elseif i == n_shifts 
            # at end of LoS, discharge patient
            #vitals 
            start_time_discharge = start_time_vitals_1 + 5
            add_vitals_to_schedule!(p, day_index, shift_index, start_time_vitals_1)
            add_discharge_to_schedule!(p, day_index, shift_index, start_time_discharge)
        end

    end


end


function discharge_patient!(id::Int64, 
                            p_list::Entities.PatientList_T)

    delete!(p_list, id)


end



function display_service_schedule(p::Entities.Patient)

    println("service schedule for patient: $(p.id)")
    for i in sort(collect(keys(p.service_schedule)))

        events_i = p.service_schedule[i]

        for event_i in events_i
            event_str = 
                [rpad("day: $(event_i.day)", 10, " ") * 
                rpad("shift: $(event_i.shift)", 15, " ") *  
                rpad("task: $(event_i.name)", 20, " ") * 
                rpad("duration: $(event_i.duration)", 20, " ")*  
                rpad("priority: $(event_i.priority)", 20, " ")]
            println(event_str)
        end
    end

end



##### end module 
end