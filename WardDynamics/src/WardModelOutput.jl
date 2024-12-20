## structures and functions for output interpretation
# idea is to generate sequences of location trajectories
# for input into airflow and transmission models. 

module WardModelOutput

using DataFrames
using CSV 

import Main.Entities 
import Main.Events

#abstract types 
abstract type Event_Trajectory_T end 
abstract type Agent_Trajectory_T end 
abstract type Shift_Trajectory_T end 

abstract type Event_Summary_T end 

abstract type Timestamp_T end

#concrete types
Event_Sequence = Array{Event_Summary_T, 1}


mutable struct Timestamp<:Timestamp_T 
    day::Int64
    shift::Int64
    time::Int64 #in minutes from start of shift. 
    Timestamp() = new()
end

#for readable printing of output. 
function timestamp_2_vector(ts::Timestamp)::Array{Int, 1}
    outvec = [ts.day, ts.shift, ts.time]
    return outvec 
end

function default_init_timestamp!(t::Timestamp)
    t.day = -1
    t.shift = -1
    t.time = -1.0
end

function init_timestamp!(t::Timestamp, 
                        day::Int64, 
                        shift::Int64, 
                        time::Int64)
    t.day = day
    t.shift = shift
    t.time = time 
end


mutable struct Event_Summary<:Event_Summary_T 
    event_name::String
    location_name::String 
    location_id::Int64
    start_time::Timestamp
    finish_time::Timestamp
    Event_Summary() = new()
end

function default_init_event_summary!(e::Event_Summary)
    e.event_name = "not defined"
    e.location_name = "not defined"
    e.location_id = -1
    
    ts_1 = Timestamp()
    ts_2 = Timestamp()
    default_init_timestamp!(ts_1)
    default_init_timestamp!(ts_2)

    e.start_time = ts_1
    e.finish_time = ts_2
end

function init_event_summary(event::Events.Event_T, 
                             location_name::String, 
                             location_id::Int64, 
                             shift::Int64, 
                             day::Int64)::Event_Summary

    es = Event_Summary()
    default_init_event_summary!(es)
    es.event_name = event.name
    es.location_name = location_name
    es.location_id = location_id
    init_timestamp!(es.start_time, day, shift, event.start_time)
    init_timestamp!(es.finish_time, day, shift, event.finish_time)

    return es
end

# not sure I'll need this - keep it for now. 
mutable struct Shift_Trajectory<:Event_Trajectory_T 
    day_index::Int64
    shift_index::Int64
    event_sequence::Event_Sequence
    Shift_Trajectory() = new()
end

# contains a timeseries of events involving an agent. 
mutable struct Agent_Trajectory<:Agent_Trajectory_T
    agent_id::Int64
    event_sequence::Event_Sequence

    Agent_Trajectory() = new()
end

#maps agent id to trajectory structure. 
Agent_Trajectories = Dict{Int64, Agent_Trajectory}

function init_agent_trajectory(id::Int64)::Agent_Trajectory

    trajectory = Agent_Trajectory()
    trajectory.agent_id = id
    trajectory.event_sequence = Event_Sequence()

    return trajectory 

end

function init_agent_trajectories!(at::Agent_Trajectories,
                                  agents::Dict{Int64, <:Entities.Entity_T} )
    for (id, a) in agents
        at[id] = init_agent_trajectory(id)        
    end
end


# TODO: generalise this to a utility (this is super usful.)
function init_sequence_dataframe()::DataFrames.DataFrame
    #inititalise empty dataframe
    # with fields of event summary struct. 
    s = Event_Summary()
    default_init_event_summary!(s)
    fields = fieldnames(typeof(s))
    df = DataFrame()
    for f in fields 
        T = typeof(getfield(s, f))
        # TODO: decide whether or not this is needed, perhaps this is over-engineered 
        #special case for timestamps: 
        if T == Timestamp
            T = Array{Int64, 1}
        end
        df[:, f] = Array{T, 1}()
    end

    return df 
end

function add_event_summary_to_dataframe!(df::DataFrames.DataFrame, es::Event_Summary)
    fields = fieldnames(typeof(es))
    df_entry = []
    for f in fields 
        
        p = getproperty(es, f)
        #TODO: check that properties have correct indices in df. 
        if typeof(p) == Timestamp
            p = timestamp_2_vector(p)
        end
        push!(df_entry, p)
       
    end
    push!(df, df_entry)
end

function fill_sequence_dataframe!(df::DataFrames.DataFrame, sequence::Event_Sequence)
    for s in sequence
        add_event_summary_to_dataframe!(df, s)
    end
end

function trajectories_to_dataframe(trajectories::Agent_Trajectories)::DataFrames.DataFrame

    df_all = init_sequence_dataframe()
    # add the agent_id field: 
    df_all[:, :agent_id] = Array{Int64, 1}()

    #sort: 
    agent_ids = sort(collect(keys(trajectories)))

    id_post_burn_in = 0

    for id in agent_ids

        traj = trajectories[id]
        if !isempty(traj.event_sequence)
            id_post_burn_in += 1
            df_i = init_sequence_dataframe()
            sequence = traj.event_sequence
            fill_sequence_dataframe!(df_i, sequence)
            # add agent_id field 
        
            
            df_i[!, :agent_id] .= id_post_burn_in

            append!(df_all, df_i)
        end
    
        
    end

    return df_all 
end

function trajectories_to_CSV(agent_trajectories::Agent_Trajectories, filename::String)
    trajectories_df = trajectories_to_dataframe(agent_trajectories)
    CSV.write(filename, trajectories_df, delim = '\t')
end




## end module 
end