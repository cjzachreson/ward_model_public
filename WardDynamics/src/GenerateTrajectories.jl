module GenerateTrajectories 

import Main.Configuration
import Main.Utilities
import Main.EventEntityInterface
import Main.Events 
import Main.Entities
import Main.WardModelOutput
import Main.InitialiseWard
import Main.PatientFunctions
import Main.NurseFunctions
import Main.UpdateWard


function generate_trajectories!(w::InitialiseWard.Ward, 
                                c::Configuration.Config,
                                nt::WardModelOutput.Agent_Trajectories,
                                pt::WardModelOutput.Agent_Trajectories)

    
    # TODO: patient trajectories too of course (but test on nurses first)

    WardModelOutput.init_agent_trajectories!(nt,
                                             w.nurses)

    # nurses on ward will need to be updated
    # dynamically to simulate shift handover. 
    nurses_on_ward = Set{Int64}()

    n_days_burn_in = Int64(c.n_roster_periods_burn_in * c.days_per_roster)
    days_post_burn_in = 0 



    for d in 1:w.n_days

        if d > n_days_burn_in 
            days_post_burn_in += 1
        end


        for s in 1:c.shifts_per_day
    
            println("****** day $d, shift $s ******")
    
            # discrete time dynamics for the shift are implemented here: 
            # this will result in the set of tasks actually performed, 
            # and the time required for each. NOTE: make sure to provide output 
            # that can trace each agent (eg. hcw and patient) to the room they occupy 
            # at each timepoint (critical for wells-riley risk calcs.)
    
            # initiate the discrete timeseries: 
            # duration of shift in minutes: 
            tf = Int64(8 * 60)
            t0 = Int64(0)
            #5-minute timestep, timestep should be smaller than all 
            # expected event durations in the shift. NOTE: this is a soft
            # requirement, and must be so because instantiated events 
            # will have random durations that may not be evenly divisible by
            # dt. 
            
            #set of timesteps: 
            t_steps = range(t0, tf, step = c.dt)
    
    
    
            # NOTE: shift_index is used to locate the shift in the roster
            # the index s is used to index events scheduled on shift s of day d. 
            w.shift_index += 1
            
            shift = w.shifts[w.shift_index]
            # sorting for reproducibility 
            nurses_this_shift = sort(collect(keys(shift.nurse_to_bbs)))
            UpdateWard.assign_bedblocks_to_hcws!(w.nurses,
                                                 shift,
                                                 nurses_this_shift)
            
            #schedule breaks 
            # TODO: move thisto UpdateWard module. 
            for id in nurses_this_shift
                NurseFunctions.schedule_breaks!(w.nurses[id], tf, d, s)
            end
    
    
            # all nurses arriving in the new shift are added to the list
            # attending to the ward. nurses will be removed from the 
            # set after handover of patients to the next shift. 
            # TODO: check - at the start of each shift, 
            #   the nurses present should be all those from the 
            #   previous shift and all those from the current shift,
            #   if there are others present then something has gone wrong
            #   and they haven't been clocked-out at handover. 
            Utilities.push_vector_to_set!(nurses_on_ward, nurses_this_shift)
    
            #TODO: draw from poisson dist with admissions 
            # per shift as rate. 
            admissions_this_shift = Int64(ceil(w.admissions_per_shift)) 
    
            for p in 1:admissions_this_shift
                #note: the function name below is a misnomer - 
                # this function does not actually 'admit' the patient
                # in the sense that the 'admission' 
                # evnet has not actually occurred yet.
                
                # need to check if ward is full, otherwise 
                # the unfinished task queue will back up. 
                p_id = UpdateWard.admit_next_patient!(w, c, d, s)
                pt_entry = WardModelOutput.init_agent_trajectory(p_id)
                pt[p_id] = pt_entry
                #NOTE: if the function returns 0, we're "bedblocked" i.e., full. 
                #PatientFunctions.display_service_schedule(admitted_patients[p_id])
            end
    
            # now generate event list from patients on ward. 
            # this should pick up unfinished tasks from last shift. 
            if c.homogeneous
                shift_queue = UpdateWard.queue_shift_duties_homogeneous(w.admitted_patients,
                                                                        size(nurses_this_shift, 1),
                                                                        d,
                                                                        s)
            else 
                shift_queue = UpdateWard.queue_shift_duties(w.admitted_patients, d, s)
            end
    
            #println("all patient service tasks (unassigned)")
            #UpdateWard.display_events(shift_queue)
    
            #UpdateWard.display_events(shift_queue)
    
            # once full shfit queue is created, iterate through and
            # add the necessary HCWs for each task 
            #NOTE: this does not assign the task to the HCW assigned_tasks vec
            # task assignment is done during the time loop.
    
            if (d, s) == (1, 1)
                UpdateWard.assign_HCWs_to_tasks!(shift_queue,
                                                 w.nurses,
                                                 shift,
                                                 w.beds)
            else 
                last_shift = w.shifts[w.shift_index - 1]
                #using multiple dispatch. 
                UpdateWard.assign_HCWs_to_tasks!(shift_queue,
                                                 w.nurses,
                                                 shift,
                                                 last_shift,# additional argument for handover
                                                 w.beds)
            end
    
    
            all_past_events = Array{Events.Event_T, 1}()
            # consider datatype for current_events
            # Dict is good for adding and deleting elements
            # but not sure what the most useful key values would be. 
    
            current_events = Dict{Int64, Events.Event_T}()
    
            event_index = 0
    
            past_events_shift = Array{Events.Event_T, 1}()
    
            finished_events = Array{Events.Event_T, 1}()  
    
            for t in t_steps 
    
                # iniitate new events - all required agents must be available. 
                # for now this is pretty easy - the tasks are one-to-one nurse to patient, 
                # so just iterate through nurses who are not busy and initiate their 
                # next tasks. Later there will be events involving different kinds of 
                # staff members, and even combinations (i.e., doctor + nurse for discharge)
                # 2 nurses for handover, etc. 
    
                ids_not_busy = UpdateWard.find_available_hcw_ids(nurses_this_shift, 
                                                                 w.nurses)
    
                #NOTE: this model will do weird things if the ward is not very full... 
                # i.e., nurses will just do all the tasks immediately, even though these
                # are supposed to be spaced out over the shift... 
                
                # initiate events from hcw assignments: 
                
                t_remaining = tf - t 
                
    
                # add HCW-initiated events (i.e., breaks - these cannot be allocated apriori)
                UpdateWard.add_breaks_to_shift_queue!(t, 
                                                      ids_not_busy, 
                                                      w.nurses, 
                                                      shift_queue)
    
                #TODO: what about long events that span multiple shifts? I'll need to consider this
                # once the event durations are drawn from random distributions that could possibly
                # select durations longer than a shift length. 
                events_to_initiate = UpdateWard.queue_next_events!(ids_not_busy, 
                                                                   shift_queue, # removes next events
                                                                   t_remaining,
                                                                   t)
    
                for e in events_to_initiate
                    event_index += 1
                    current_events[event_index] = e
                    Events.start_Event!(e, t)
                end
    
                # step active events
                for (e_i, e) in current_events 
                    # stepping finishes any events that have reached their endpoint
                    is_finished = Events.step_Event!(e, c.dt, t)
                    if is_finished 
                        push!(finished_events, e)
                        delete!(current_events, e_i)
                    end
                end
    
                # update agent trajectories from finished events: 
                # TODO: put this in separate function. 

                ### TODO: TEST:only record after burn-in is done. 
                if d > n_days_burn_in 

                    for fe in finished_events 
        
                        #location name
                        
                        #location id 
                        loc_id = -1
                        if in(fe.name, ["Long Break", "Short Break"])
                            loc_id = -10
                            loc_name = "break room"
                        else
                            bed_id = fe.entities["patient"][1].bed_id
                            loc_id = w.beds[bed_id].room_id
                            loc_name = "room $loc_id"
                        end
        
                        es = WardModelOutput.init_event_summary(fe, 
                                                                loc_name, 
                                                                loc_id, 
                                                                s, #values in [1, 2, 3]
                                                                days_post_burn_in)
                        #add event summary to hcw trajectories 
                        if haskey(fe.entities, "nurse")
                            for hcw in fe.entities["nurse"]
                                nurse_id = hcw.id 
                                push!(nt[nurse_id].event_sequence, es)
                            end
                        end

                        #add event summary to patient trajectories 
                        if haskey(fe.entities, "patient")
                            for patient in fe.entities["patient"]
                                patient_id = patient.id
                                push!(pt[patient_id].event_sequence, es)
                            end
                        end

                    end
                end
    
    
    
                # transfer finished events from current_events to past_events 
                append!(all_past_events, finished_events)
                append!(past_events_shift, finished_events)
                # clear list of finished events for next step. 
    
                #check for discharge events in set of finished events. 
                ids_to_discharge = UpdateWard.queue_discharges(finished_events)
    
                UpdateWard.finish_handover!(finished_events, w.beds)
    
                # check completed handovers to see if nurses are ready to leave
                # if a nurse from last shift has no more patients assigned 
                # then handover is done and they can clock out. 
                ids_to_clockout = UpdateWard.queue_clockouts(finished_events,
                                                             w.bedblocks, 
                                                             w.beds)
    
                # #display finished events: 
                # if !isempty(finished_events)
                #     println("*** finished events: ")
                #     UpdateWard.display_events(finished_events)
                #     #println("*** current events: ")
                #     #UpdateWard.display_events(collect(values(current_events)))
                #     #println("****")
                # end
    
                empty!(finished_events)  
    
                # may need to copy by value for past events, to remove 
                # references to ids of discharged patients.
                
    
                UpdateWard.discharge_patients!(w.admitted_patients,
                                               ids_to_discharge,
                                               w.beds)
    


                # NOTE: unplanned overtime could be tabulated here.
                # would need to add hcw state variable for time on shift. 
                UpdateWard.clockout_nurses!(w.nurses, #clears assigned tasks 
                                            nurses_on_ward, #removes id from ward
                                            ids_to_clockout) 
    
                #TODO: ensure tasks in-progress are completed by end of shift or
                #  passed to the next shift. 
                # tasks may need to be pushed to next day or next shift. 
                # May need to update length-of-stay if discharge is not possible. 
                
            end
    
            #display tasks finished this shift: 
                        #display finished events: 
            if !isempty(past_events_shift)
                println("*** finished tasks: ")
                #UpdateWard.display_events(past_events_shift)
            end
    
            # clear break schedule 
            #TODO: test this 
            for id in nurses_this_shift
                NurseFunctions.clear_break_schedule!(w.nurses[id])
            end

            # removes any assigned tasks from hcws on duty and transfers them
            # to the task queue for the next shift. 
            
            #UpdateWard.queue_unfinished_tasks!(nurses, UpdateWard.unfinished_tasks)
    
            # move unfinished tasks from nurse assignments back to patient's 
            # service schedule. 
            #TODO: update with new task assignment specs. 
            # check that priority order is not getting flipped. 
            # empties shift queue. 
            unfinished_tasks = UpdateWard.requeue_unfinished_tasks!(shift_queue, #gets emptied
                                                                    w, 
                                                                    c,
                                                                    d, 
                                                                    s)
    
            println("**** unfinished tasks:")
            #UpdateWard.display_events(unfinished_tasks)
    
            # transfer uncompleted tasks back to shift queue for next shift. 
            # incomplete tasks are contained in the list of assigned tasks for each
            # hcw. Note that some event properties will need to be cleared. 
    
     
            #TODO: transitions between tasks must be implemented too. 
    
        end
    end



    

end


function insert_patient_wait_times!(pt::WardModelOutput.Agent_Trajectories,
                                   c::Configuration.Config)

    dt = c.dt
    t_start_shift = 0
    t_end_shift = 1440 / c.shifts_per_day # minuts per shift 

    for (id, traj) in pt # agent id -> event sequences 

        # copy event sequence with no waiting: 
        e_seq_original = traj.event_sequence
        e_seq_with_waiting = WardModelOutput.Event_Sequence()

        for es_i in 1:(length(e_seq_original) - 1) # array of event summaries

            es_this = e_seq_original[es_i]

            es_next = e_seq_original[es_i + 1]

            t_this_finish = es_this.finish_time # recall these are timestamps
            t_next_start = es_next.start_time

            # check if they're in the same shift: 
            if t_this_finish.shift == t_next_start.shift
                t_wait_start = t_this_finish.time + dt
                t_wait_finish = t_next_start.time - dt

                # make one wait event 
                #most traits are the same as the event before the wait
                wait_1 = deepcopy(es_this) 
                wait_1.event_name = "Wait"
                wait_1.start_time.time = t_wait_start
                wait_1.finish_time.time = t_wait_finish

                push!(e_seq_with_waiting, es_this)
                push!(e_seq_with_waiting, wait_1)
            else

                t_wait_1_start = t_this_finish.time + dt
                t_wait_1_finish = t_end_shift 

                t_wait_2_start = t_start_shift
                t_wait_2_finish = t_next_start.time - dt 
               
                push!(e_seq_with_waiting, es_this)

                # make two wait events 
                wait_1 = deepcopy(es_this) 
                wait_1.event_name = "Wait"
                wait_1.start_time.time = t_wait_1_start
                wait_1.finish_time.time = t_wait_1_finish

                push!(e_seq_with_waiting, wait_1)

                if t_next_start.time != t_start_shift
                    wait_2 = deepcopy(es_next)
                    wait_2.event_name = "Wait"
                    wait_2.start_time.time = t_wait_2_start
                    wait_2.finish_time.time = t_wait_2_finish
                    push!(e_seq_with_waiting, wait_2)
                end

            end
        end

        if !isempty(e_seq_original)
            push!(e_seq_with_waiting, last(e_seq_original))
            traj.event_sequence = e_seq_with_waiting
        else
            println("empty event sequence for agent $id")
            traj.event_sequence = e_seq_original
        end

        

    end


end





# end module 
end