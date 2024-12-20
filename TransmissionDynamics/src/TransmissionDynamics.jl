module TransmissionDynamics

using Random
using DataFrames
using Main.AirflowTransmissionInterface

import Main.EntitiesTD 
import Main.AirflowDynamics
import Main.Utilities
import Main.OutputStructures



abstract type Agent_Lists_T end
#abstract type Config_T end
abstract type Shift_State_T end
abstract type Experiment_Constants_T end

mutable struct Experiment_Constants<:Config_T
    trajectories_hcw::DataFrame
    trajectories_patient::DataFrame  

    flow_network::EntitiesTD.Flow_Network
    rooms::EntitiesTD.Rooms_T

    output_dirname::String

    model_type::String

    index_case_type::String

    experiment_type::String

    detailed_output_flag::Bool

    Experiment_Constants() = new() 

end

function init_Experiment_Constants!(ex::Experiment_Constants,
                                   fname_trajectories_hcw::String, 
                                   fname_trajectories_patient::String, 
                                   fname_ward::String,
                                   output_dirname::String,
                                   model_type::String,
                                   index_case_type::String,
                                   experiment_type::String,
                                   detailed_output_flag::Bool)

    ex.trajectories_hcw = Utilities.parse_trajectories_from_file(fname_trajectories_hcw)
    ex.trajectories_patient = Utilities.parse_trajectories_from_file(fname_trajectories_patient)

    #network structure of ward
    ex.flow_network = EntitiesTD.read_flow_network_from_file(fname_ward)
    ex.rooms = EntitiesTD.read_rooms_from_file(fname_ward)
    
    #EntitiesTD.add_flows_to_rooms!(ex.rooms, ex.flow_network)


    ex.output_dirname = output_dirname 

    ex.model_type = model_type
    ex.index_case_type = index_case_type
    ex.experiment_type = experiment_type

    # if true, simulations produce detailed output for each
    # run including case linelists and timeseries of cases
    ex.detailed_output_flag = detailed_output_flag

end


mutable struct Config<:Config_T
    #TODO: add comments 
    # TODO: ensure all DataFrames are used as static objects
    # these might be OK (they hold the full list of trajectories)
    trajectories_hcw::DataFrame
    trajectories_patient::DataFrame

    hcw_ids::Array{Int64, 1}
    patient_ids::Array{Int64, 1}

    shifts_per_day::Int64
    n_days::Int64
    pathogen::EntitiesTD.Pathogen
    t0::Float64
    minutes_per_shift::Int64
    minutes_per_day::Int64
    dt_days::Float64
    dt_minutes::Int64
    timesteps_g::Array{Float64, 1}

    flow_network::EntitiesTD.Flow_Network
    rooms::EntitiesTD.Rooms_T

    output_structure::OutputStructures.Output_Structure_T 

    output_dirname::String

    rng_infections::MersenneTwister

    rng_index_cases::MersenneTwister

    R0_calib_flag::Bool

    index_case_type::String

    n_index_cases::Int64

    index_case_ids_hcws::Array{Int64, 1}
    index_case_ids_patients::Array{Int64, 1}

    q_vs_t_flag::Bool # tells whether or not to record timeseries of quanta concentrations 

    mask_effect::Float64 

    Config() = new()
end

#TODO: include masking parameters 
function init_config!(c::Config, 
                      ex::Experiment_Constants,
                      output_dirname::String,
                      seed_infections::Int64,
                      seed_index_cases::Int64,
                      beta::Float64,
                      q_vs_t_flag::Bool,
                      n_index_cases::Int64,
                      mask_effect::Float64)
    # trajectories and temporal structure: 
    c.trajectories_hcw = ex.trajectories_hcw#Utilities.parse_trajectories_from_file(fname_trajectories_hcw)
    c.trajectories_patient = ex.trajectories_patient#Utilities.parse_trajectories_from_file(fname_trajectories_patient)

    c.hcw_ids = unique(c.trajectories_hcw.agent_id)
    c.patient_ids = unique(c.trajectories_patient.agent_id)

    c.shifts_per_day = 3
    c.n_days = maximum(c.trajectories_hcw.finish_day)

    c.pathogen = EntitiesTD.Pathogen()
    EntitiesTD.init_SARS_CoV_2!(c.pathogen, beta)

    c.t0 = 0.0
    c.minutes_per_shift = 8*60
    c.minutes_per_day = 24*60
    c.dt_minutes = 5

    c.dt_days = round(Float64(c.dt_minutes)/Float64(c.minutes_per_day), digits = 10)

    c.timesteps_g = collect(0.0:c.dt_days:Float64(c.n_days + 0.1)) #adding on a buffer 0.1 days to prevent bad indexing at the end of the timeseries. 

    #network structure of ward
    c.flow_network = ex.flow_network#EntitiesTD.read_flow_network_from_file(fname_ward)

    #broken - rooms are not constant. 
    c.rooms = deepcopy(EntitiesTD.clear_rooms(ex.rooms))#EntitiesTD.read_rooms_from_file(fname_ward)

    c.q_vs_t_flag = q_vs_t_flag


    #null initialiser 
    c.output_structure = OutputStructures.Output_Structure()
    #default (empty) initialisation 
    OutputStructures.default_init_Output_Structure!(c.output_structure)

    if q_vs_t_flag 
        # populate room_id values for viral load timeseries output: 
        OutputStructures.add_room_ids_to_q_vs_t!(c.output_structure, c.rooms)
    end

    c.output_dirname = output_dirname  

    c.rng_infections = MersenneTwister(seed_infections)

    c.rng_index_cases = MersenneTwister(seed_index_cases)

    c.R0_calib_flag = false # default to false. 

    c.index_case_type = ex.index_case_type

    c.n_index_cases = n_index_cases

    c.index_case_ids_hcws = Array{Int64, 1}()
    c.index_case_ids_patients = Array{Int64, 1}()

    c.mask_effect = mask_effect 

end


mutable struct SIR_Agents <: Agent_Lists_T
    all::EntitiesTD.Agents_T
    infected::EntitiesTD.Agents_T
    susceptible::EntitiesTD.Agents_T
    recovered::EntitiesTD.Agents_T
    SIR_Agents() = new()
end

function default_init_SIR_Agents!(a::SIR_Agents)
    a.all = EntitiesTD.Agents_T()
    a.infected = EntitiesTD.Agents_T()
    a.susceptible = EntitiesTD.Agents_T()
    a.recovered = EntitiesTD.Agents_T()
end

function init_SIR_Agents!(a::SIR_Agents, ids::Array{Int64, 1}, label::String)
    EntitiesTD.add_new_agents!(ids, a.all, label)
end

#TODO: refactor to avoid dynamic allocation of dataframes. 
function init_Patient_Locations!(a::SIR_Agents, s::DataFrame)

    #println("initialising patient locations")

    for (id, agent) in a.all
        first_location = Utilities.extract_first_location_id(s, id)
        #NOTE: this produces the id of the first location the agent visits, 
        # not their initial location (which might not exist e.g., for a patient not yet admitted)
        # this is fine because the update alglorithm only iterates through agents 
        # who are present on a given shift. 
        
        agent.current_room_id = first_location 
    end

    #println("done initialising patient locations")
end


function infect_index_cases!(c::Config, a::SIR_Agents, ids_to_infect::Array{Int64, 1})

    all_ids = collect(keys(a.all))

    for id in all_ids
        if in(id, ids_to_infect)
            #initialise infections 

            #println("infected index case: $id")


            EntitiesTD.new_agent_infection!(a.all[id], c.pathogen, c.t0, c.rng_infections)
            a.all[id].is_index_case = true 
            a.infected[id] = a.all[id]
            #add to linelist 
            OutputStructures.fill_Infection_Data!(c.output_structure, #includes dict of infections 
                                                  id, #agent_id,
                                                  a.infected[id].infections[end],#infection,
                                                  true,#is_index,
                                                  c.t0,#t,
                                                  1,#d,
                                                  1,#s,
                                                  a.infected[id].label)#type)
        else
            a.susceptible[id] = a.all[id]
        end
    end

end


#TODO: make sure this is stable. 
function initialise_infections_hcw!(c::Config, a::SIR_Agents, n::Int64)
    first_traj = EntitiesTD.extract_shift_trajectories(c.trajectories_hcw, 1, 1)

    # this is too constrained - I think it's producing odd fluctuations
    # select from anyone present in the first day 
    # TODO: test what happens if this pool is expanded. 
    ids_first_shift = collect(keys(first_traj.agent_trajectories))
    
    
    ids_to_infect = Array{Int64, 1}()

    shuffle!(c.rng_index_cases, ids_first_shift)

    #println("ids on first shift: $ids_first_shift")

    for i in 1:n
        if size(ids_first_shift, 1) >= n
            push!(ids_to_infect, ids_first_shift[i]) # TODO: apply as argument
        else 
            println("index case error: not enough agents to infect")
        end
    end

    c.index_case_ids_hcws = ids_to_infect

    infect_index_cases!(c, a, ids_to_infect)

end

function initialise_infections_patient!(c::Config, a::SIR_Agents, n::Int64)
    first_traj = EntitiesTD.extract_shift_trajectories(c.trajectories_patient, 1, 1)
    ids_first_shift = collect(keys(first_traj.agent_trajectories))
    ids_to_infect = Array{Int64, 1}()
    
    # randomising index case: 
    shuffle!(c.rng_infections, ids_first_shift)

    for i in 1:n
        if size(ids_first_shift, 1) >= n
            push!(ids_to_infect, ids_first_shift[i]) # TODO: apply as argument
        else 
            println("index case error: not enough agents to infect")
        end
    end

    c.index_case_ids_patients = ids_to_infect

    infect_index_cases!(c, a, ids_to_infect)

end

function initialise_infections_by_ID!(c::Config, 
                                      a::SIR_Agents, 
                                      IDs::Array{Int64, 1})
    ids_to_infect = IDs
    infect_index_cases!(c, a, ids_to_infect)
end

#shift-wise update functions: 

mutable struct Shift_State<:Shift_State_T

    #TODO: refactor dataframes 
   trajectories_hcw::EntitiesTD.Shift_Trajectories
   trajectories_patient::EntitiesTD.Shift_Trajectories

   ids_present_hcw::Array{Int64, 1}
   ids_present_patient::Array{Int64, 1}

   t_i::Float64
   t_f::Float64

   day_index::Int64
   shift_index::Int64 
   
   Shift_State() = new()
end

#NOTE: c.trajectories and s.trajectories refer to different things -this is confusing
#TODO: rename one of the properties (maybe trajectories_df for config.)
function init_Shift_State!(s::Shift_State, c::Config, d_i::Int64, s_i::Int64)

     
    s.trajectories_hcw = EntitiesTD.extract_shift_trajectories(c.trajectories_hcw, d_i, s_i)
    s.trajectories_patient = EntitiesTD.extract_shift_trajectories(c.trajectories_patient, d_i, s_i)

    s.ids_present_hcw = collect(keys(s.trajectories_hcw.agent_trajectories))
    s.ids_present_patient = collect(keys(s.trajectories_patient.agent_trajectories))

    t_minutes = (d_i - 1) * c.minutes_per_day +  (s_i - 1) * c.minutes_per_shift

    s.t_i = Float64(t_minutes) / Float64(c.minutes_per_day) 

    #CHECK that time is being indexed properly. 
    s.t_f = s.t_i + (1.0 / Float64(c.shifts_per_day)) - c.dt_days

    s.day_index = d_i
    s.shift_index = s_i
end

# TODO: check that quanta computation is correct (verify equations)
function compute_quanta(infection::EntitiesTD.Infection, c::Config, a::EntitiesTD.AgentTD)::Float64
    # NOTE: any front-line NPIs could be implemented here 
    # beta_t corresponds to quanta per day
    if (c.R0_calib_flag && (!a.is_index_case))
        quanta = 0.0
    else
        #quanta = infection.beta_t * c.dt #NOTE: removing c.dt term here (the rate is used in forard-euler)
        quanta = infection.beta_t
    end

    return quanta
end


function recover_agents!(a::SIR_Agents, ids::Array{Int64, 1})

    # keep track of which agents are recovered. 
    for id_recovered in ids
        delete!(a.infected, id_recovered) 
        a.recovered[id_recovered] = a.all[id_recovered]
    end

end

# evaluation of mask wearing
function evaluate_masking_effect(c::Config, room::EntitiesTD.RoomTD, agent::EntitiesTD.AgentTD)::Float64
    # see if the agent is a nurse: 
    mask_effect = 0.0
    #ensure the nurse label aligns with the label used at initialisation of SIR agents 
    if agent.label == "nurse"
        if contains(room.name, "patient room")
            mask_effect =  c.mask_effect

        end
    end
    return mask_effect
end

#step active infections 
function update_infections_hcw!(a::SIR_Agents, c::Config, s::Shift_State) 
    ids_recovering = Array{Int64, 1}()
 
    for (id, agent) in a.infected
        if in(id, s.ids_present_hcw)

            infection = last(agent.infections)
            EntitiesTD.step_infection!(infection, s.t_i)

            if infection.recovered
                push!(ids_recovering, id)
            end

            EntitiesTD.compute_FoI!(infection)
            
            quanta = compute_quanta(infection, c, agent)

            #agent_trajectory is a dataframe 
            agent_trajectory = s.trajectories_hcw.agent_trajectories[id]

            n_events = size(agent_trajectory, 1)
            
            for e_i in 1:n_events
                event = agent_trajectory[e_i, :]
                room_id = event.location_id



                if haskey(c.rooms, room_id)

                    t_ind_1 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                                s.t_i,
                                                                c.minutes_per_day,
                                                                event.start_time)

                    t_ind_2 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                                s.t_i,
                                                                c.minutes_per_day,
                                                                event.finish_time)

                    if c.timesteps_g[t_ind_2] > s.t_f #avoids double-counting shedding between shifts. 
                        t_ind_2 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                                    s.t_i,
                                                                    c.minutes_per_day,
                                                                    c.minutes_per_shift - c.dt_minutes)
                    end


                    # TODO: add masking effect here (modifies quanta shed into room)
                    mask_effect = evaluate_masking_effect(c, c.rooms[room_id], agent)

                    quanta_masked = quanta * (1.0 - mask_effect)

                    for t_i in t_ind_1:1:t_ind_2

                        t_step = c.timesteps_g[t_i]

                        EntitiesTD.add_shedding_to_room!(c.rooms[room_id], quanta_masked, t_step)
                    end
                
                else 
                    println("warning: location $room_id not found in config")
                
                end  



            end
        end
    end
    recover_agents!(a, ids_recovering)
end


function update_infections_patient!(a::SIR_Agents, c::Config, s::Shift_State) 
    ids_recovering = Array{Int64, 1}()

    for (id, agent) in a.infected

        if in(id, s.ids_present_patient)
            
            infection = last(agent.infections)
            EntitiesTD.step_infection!(infection, s.t_i)

            if infection.recovered
                push!(ids_recovering, id)
            end

            EntitiesTD.compute_FoI!(infection)

            quanta = compute_quanta(infection, c, agent)
     
            agent_trajectory = s.trajectories_patient.agent_trajectories[id]

            n_events = size(agent_trajectory, 1)
            
            for e_i in 1:n_events
                event = agent_trajectory[e_i, :]
                room_id = event.location_id

                if haskey(c.rooms, room_id)

                    t_ind_1 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                                s.t_i,
                                                                c.minutes_per_day,
                                                                event.start_time)
                    t_ind_2 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                                s.t_i,
                                                                c.minutes_per_day,
                                                                event.finish_time)
                        # avoids double-counting shedding. 
                        if c.timesteps_g[t_ind_2] > s.t_f 
                            t_ind_2 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                                    s.t_i,
                                                                    c.minutes_per_day,
                                                                    c.minutes_per_shift - c.dt_minutes)
                        end

                    mask_effect = evaluate_masking_effect(c, c.rooms[room_id], agent)

                    quanta_masked = quanta * (1.0 - mask_effect)

                    for t_i in t_ind_1:1:t_ind_2

                        t_step = c.timesteps_g[t_i]

                        EntitiesTD.add_shedding_to_room!(c.rooms[room_id], quanta_masked, t_step)
                    end

                else 
                    println("warning: location not found in config")
                end
            end


        end
    end
    # keep track of which agents are recovered. 
    recover_agents!(a, ids_recovering)
end



function update_quanta_forward_euler!(s::Shift_State, c::Config_T)



    t_ind_1 = Utilities.shift_time_2_tstep_index(c.timesteps_g, 
                                                 s.t_i, 
                                                 c.minutes_per_day, 
                                                 0)

    t_ind_2 = Utilities.shift_time_2_tstep_index(c.timesteps_g, 
                                                 s.t_i, 
                                                 c.minutes_per_day, 
                                                 c.minutes_per_shift)

    time_indices = (t_ind_1 + 1):1:t_ind_2

    
    for i in time_indices

        t_i = c.timesteps_g[i - 1]
    
        t_f = c.timesteps_g[i]
    
        # TODO: make the fine-scale factor a config property 
        q_out = AirflowDynamics.quanta_flow_forward_euler_optimised!(c, t_i, t_f, 5.0)
        #q_out = AirflowDynamics.quanta_flow_forward_euler_optimised_no_q_vs_t!(c, t_i, t_f, 5.0)

    end

end


function update_quanta_forward_euler_no_q_vs_t!(s::Shift_State, c::Config_T)



    t_ind_1 = Utilities.shift_time_2_tstep_index(c.timesteps_g, 
                                                 s.t_i, 
                                                 c.minutes_per_day, 
                                                 0)

    t_ind_2 = Utilities.shift_time_2_tstep_index(c.timesteps_g, 
                                                 s.t_i, 
                                                 c.minutes_per_day, 
                                                 c.minutes_per_shift)

    time_indices = (t_ind_1 + 1):1:t_ind_2

    
    for i in time_indices

        t_i = c.timesteps_g[i - 1]
    
        t_f = c.timesteps_g[i]
    
        #q_out = AirflowDynamics.quanta_flow_forward_euler_optimised!(c, t_i, t_f, 5.0)
        q_out = AirflowDynamics.quanta_flow_forward_euler_optimised_no_q_vs_t!(c, t_i, t_f, 5.0)

    end

end



# function update_quanta_super_sink!(s::Shift_State, c::Config_T )


    #     for step in 1:size(s.timesteps, 1)
    #         t = s.timesteps[step]
    #         if step > 1
    #             t_previous = s.timesteps[step-1]
    #         else
    #             t_previous = -1.0*c.dt
    #         end

    #         for (room_id, room) in c.rooms

    #             if haskey(room.viral_quanta, t_previous)

    #                 #note that -1 is the id for the super sink. 
    #                 # decay per minute
    #                 decay_rate = c.flow_network.st_links_2_lambdas[(room_id, -1)]
    #                 #println("decay rate: $decay_rate")

    #                 dt_minutes = c.dt * c.minutes_per_day
    #                 #decay per timestep: 
    #                 #decay_rate *= dt_minutes

    #                 room.viral_quanta[t] = room.viral_quanta[t_previous] 

    #                 #quanta_decayed = decay_rate * room.viral_quanta[t] 
    #                 quanta_decayed = room.viral_quanta[t] * 
    #                                  (1.0 - exp(-1.0 * decay_rate * dt_minutes))
                    
    #                 # if room.viral_quanta[t] > 0.0
    #                 #     println(quanta_decayed)
    #                 # end


    #                 EntitiesTD.remove_quanta_from_room!(room, quanta_decayed, t) 

    #                #println("$(s.day_index), $(s.shift_index): $(room.viral_quanta[t])")
    #             end
    #             # add newly-shed quanta to room. 
    #             if !isempty(room.shedding_quanta)
    #                 if haskey(room.shedding_quanta, t)
    #                     #add shedding to viral quanta: 
    #                     shed_quanta = room.shedding_quanta[t]
    #                     EntitiesTD.add_quanta_to_room!(room, shed_quanta, t)
    #                 end
    #             end

    #         end

    #     end

    # end

#TODO: optimise this function - this is a major bottleneck. 
# the inner loop (over t in timesteps) is about 30% total run time. 
#function room_exposure_quanta(c::Config, room_id::Int64, timesteps::Array{Float64, 1})::Float64
function room_exposure_quanta(c::Config, room_id::Int64, timesteps::Array{Float64, 1})::Float64
    q = 0.0

    if haskey(c.rooms, room_id)
        room = c.rooms[room_id]
        vq = room.viral_quanta

        # if this for loop could be vectorised, it might speed up. 
        for t in timesteps
            if haskey(vq, t)
                quanta = vq[t]
                q += quanta
            end
        end

        q *= c.dt_days 
        q /= room.volume
      
   end
    return q
end


# function that produces mask effect - note that this parameter needs to be added to config. 


#TODO: consider updating floating-point time operations, low-precision integers will be faster
#TODO: implement mask wearing 
function exposure_quanta_for_shift(s::Shift_State, c::Config, traj::DataFrame, agent::EntitiesTD.AgentTD)::Float64
    q = 0.0
    n_events = size(traj, 1)
    for e_i in 1:n_events
        event = traj[e_i, :]
 
        room_id = event.location_id

        t_ind_1 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                    s.t_i,
                                                    c.minutes_per_day,
                                                    event.start_time)

        t_ind_2 = Utilities.shift_time_2_tstep_index(c.timesteps_g,
                                                    s.t_i,
                                                    c.minutes_per_day,
                                                    event.finish_time)

        timesteps = c.timesteps_g[t_ind_1:t_ind_2]

        m = evaluate_masking_effect(c, c.rooms[room_id], agent)

        q += room_exposure_quanta(c, room_id, timesteps) * (1.0 - m)


    
    end


    return q
end




#evaluate infection probability 
function eval_infection_prob!(ids_to_infect::Array{Int64, 1}, q::Float64, id::Int64, rng_infections::MersenneTwister)

    if q > 0.0
        #println("agent $id exposed to $exposure_quanta quanta")
        p_infect = 1.0 - exp(-1.0 * q)
        if rand(rng_infections) < p_infect
            push!(ids_to_infect, id)
            #println("going to infect $id (q = $q, p = $p_infect)")
        end
    end


end

# infect agents 
function infect_agents!(s::Shift_State, c::Config, a::SIR_Agents, ids_to_infect::Array{Int64, 1})

    for id in ids_to_infect
        #new infections 
        t_infect = s.t_f #initiate infection at end of shift (approximation)
        #TODO: consider timing of infections. 
        EntitiesTD.new_agent_infection!(a.all[id], c.pathogen, t_infect, c.rng_infections)
        a.infected[id] = a.all[id]
        delete!(a.susceptible, id)

        #add infection to linelist: 
        OutputStructures.fill_Infection_Data!(c.output_structure, #includes dict of infections 
                                              id, #agent_id,
                                              a.infected[id].infections[end],#infection,
                                              false,#is_index,
                                              t_infect,#t,
                                              s.day_index,#d,
                                              s.shift_index,#s,
                                              a.infected[id].label)#type)
        

    end

end

# compute exposure and new infections 
# TODO: add new infections to config's output structure. 
function update_exposure_hcw!(s::Shift_State, c::Config, a::SIR_Agents)
    #TODO: assert that a is a list of hcws, no patients. 
    ids_to_infect = Array{Int64, 1}()
    #println("-----")
    #println(a.susceptible)
    for (id, agent) in a.susceptible
        #NOTE: the two lines below differentiate hcw and patient functions 
        if in(id, s.ids_present_hcw) #_hcw
            agent_trajectory = s.trajectories_hcw.agent_trajectories[id] #_hcw
            exposure_quanta = exposure_quanta_for_shift(s, c, agent_trajectory, agent) #NOTE: agent needed for masking evaluation
            # if exposure_quanta > 0.0
                #     println("-----")
                #     println("on day: $(s.day_index), shift: $(s.shift_index)")
                #     println("hcw: $id exposed")
                #     println("quanta: $exposure_quanta")
                #     println("-----")
            # end
            eval_infection_prob!(ids_to_infect, exposure_quanta, id, c.rng_infections)
        end
    end
    infect_agents!(s, c, a, ids_to_infect)
end

# compute exposure and new infections
#TODO: patients get exposed the whole time, simplify 
function update_exposure_patient!(s::Shift_State, c::Config, a::SIR_Agents) 
    #TODO: assert that a is a list of patients, no hcws. 
    ids_to_infect = Array{Int64, 1}()
    for (id, agent) in a.susceptible
        #NOTE: the two lines below differentiate hcw and patient functions
        if in(id, s.ids_present_patient)
            #println("patient present")
            agent_trajectory = s.trajectories_patient.agent_trajectories[id]
            exposure_quanta = exposure_quanta_for_shift(s, c, agent_trajectory, agent)

            # if exposure_quanta > 0.0

            #     println("-----")
            #     println("on day: $(s.day_index), shift: $(s.shift_index)")
            #     println("patient: $id exposed")
            #     println("quanta: $exposure_quanta")
            #     println("-----")
            # end



            eval_infection_prob!(ids_to_infect, exposure_quanta, id, c.rng_infections)
        end
    end
    infect_agents!(s, c, a, ids_to_infect)
end


function simulate_transmission_forward_euler!(c::Config_T, hcws::SIR_Agents, patients::SIR_Agents) 

    for day in 1:c.n_days
        #println("***** day: $day *****")
        for shift in 1:c.shifts_per_day 
            #println("*** shift $shift ***")
            
            
            s = Shift_State()
            init_Shift_State!(s, c, day, shift)

            #println("shift: $shift")
            #println(s.ids_present_hcw)

            update_infections_hcw!(hcws, c, s)
            
            update_infections_patient!(patients, c, s)

            if c.q_vs_t_flag 
                println("flag: $(c.q_vs_t_flag)")
                update_quanta_forward_euler!(s, c)
            else
                update_quanta_forward_euler_no_q_vs_t!(s, c)
            end

            #println("**infecting hcws")
            update_exposure_hcw!(s, c, hcws)
            #println("**infecting patients")
            update_exposure_patient!(s, c, patients)

        end

        #check if anyone is infected: 
        if (isempty(hcws.infected) && isempty(patients.infected))
            break
        end

    end

end






#end module 
end