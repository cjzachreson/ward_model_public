module Utilities

using DataFrames
using CSV

# here, I'll put IO functions that I want to run in 
# compile mode during debugging (this can be changed 
# if there's problems that need to be addressed)


function parse_CSV_vec_string_2_Int64(str_in)::Array{Int64, 1}

    str = strip(str_in, [']','['])
    str = split(str, ',')
    
    vec = parse.(Int64, str)

    return vec 

end


function parse_trajectories_from_file(fname)::DataFrame

    #println("parsing trajectories frome input file")

    traj_raw = DataFrame(CSV.File(fname))

    t_start_str = traj_raw.start_time
    t_fin_str = traj_raw.finish_time 

    start_day = Array{Int64, 1}()
    finish_day = Array{Int64, 1}()

    start_shift = Array{Int64, 1}()
    finish_shift = Array{Int64, 1}()

    start_time = Array{Int64, 1}()
    finish_time = Array{Int64, 1}()



    for i in 1:size(traj_raw, 1)

        tstamp_start = parse_CSV_vec_string_2_Int64(t_start_str[i])
        tstamp_finish = parse_CSV_vec_string_2_Int64(t_fin_str[i])

        push!(start_day, tstamp_start[1])
        push!(finish_day, tstamp_finish[1])

        push!(start_shift, tstamp_start[2])
        push!(finish_shift, tstamp_finish[2])

        push!(start_time, tstamp_start[3])
        push!(finish_time, tstamp_finish[3])

    end

    df_out = select(traj_raw, Not([:start_time, :finish_time]))
    
    df_out[!, :start_day] = start_day
    df_out[!, :finish_day] = finish_day
    
    df_out[!, :start_shift] = start_shift
    df_out[!, :finish_shift] = finish_shift
    
    df_out[!, :start_time] = start_time
    df_out[!, :finish_time] = finish_time

    #println("done parsing trajectories frome input file")  

    return df_out 

end

function extract_first_location_id(s::DataFrame, id::Int64)::Int64
    loc_id = s[s.agent_id .== id, :].location_id[1]
    return loc_id 
end

function init_DF_from_struct(s::Any)::DataFrame
    #inititalise empty dataframe
    # with fields of event summary struct. 

    fields = fieldnames(typeof(s))
    df = DataFrame()
    for f in fields 
        T = typeof(getfield(s, f))
 
        df[:, f] = Array{T, 1}()
    end

    return df 
end

function DF_row_from_struct(s::Any)::DataFrame

    fields = fieldnames(typeof(s))
    df = init_DF_from_struct(s)

    for f in fields 
        push!(df[!, f],  getproperty(s, f))
    end

    return df

end

function write_DF_to_CSV(df::DataFrame, fname::String) 

    #TODO add any additional metadata etc. 
    CSV.write(fname, df)

end

# returns the closest absolute timestep to the input time 
# where the input time is the time in minutes from the 
# start of a shift. 
function shift_time_2_tstep(timesteps::Array{Float64, 1}, 
                   t_shift_start::Float64,
                   minutes_per_day::Int64,
                   time_in_shift::Int64)::Float64

    t_raw = t_shift_start + (Float64(time_in_shift) / Float64(minutes_per_day))

    t_step_index = argmin(broadcast(abs, timesteps .- t_raw))

    t_step = timesteps[t_step_index]

    return t_step

end

function shift_time_2_tstep_index(timesteps::Array{Float64, 1}, 
    t_shift_start::Float64,
    minutes_per_day::Int64,
    time_in_shift::Int64)::Int64

    t_raw = round(t_shift_start + (Float64(time_in_shift) / Float64(minutes_per_day)), digits = 10)

    # the below is slow - try something else. 
    #t_step_index = argmin(broadcast(abs, timesteps .- t_raw))

    t_step_index = searchsortedfirst(timesteps, t_raw)

    # if t_step_index > size(timesteps, 1)
    #     println("time index not found")
    #     println(last(timesteps))
    #     println(t_raw)
    # end

    return t_step_index

end



# end module 
end