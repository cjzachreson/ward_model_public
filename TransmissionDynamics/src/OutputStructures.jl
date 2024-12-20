module OutputStructures


import Main.EntitiesTD
import Main.Utilities

using DataFrames
using CSV


abstract type Infection_Data_T end

abstract type Output_Structure_T end 

abstract type Output_File_T end

mutable struct Infection_Data <: Infection_Data_T
    agent_id::Int64
    agent_type::String 
    time_infected::Float64
    day_infected::Int64
    shift_infected::Int64
    beta_max::Float64
    symptomatic::Bool 
    is_index_case::Bool

    Infection_Data() = new()
end

function default_init_Infection_Data!(idat::Infection_Data)
    idat.agent_id = -1
    idat.agent_type = "Not Specified"
    idat.time_infected = -1.0
    idat.day_infected = -1
    idat.shift_infected = -1
    idat.beta_max = -1.0
    idat.symptomatic = false
    idat.is_index_case = false 
end


mutable struct Output_Structure <: Output_Structure_T 
    # a viral load timeseries for each zone 
    viral_loads::Dict{Int64, Dict{Float64, Float64}}

    # a linelist of all infections on the ward
    # linelist row index -> infection data
    n_infections::Int64 
    infections::Dict{Int64, Infection_Data}

    Output_Structure() = new()
end

function default_init_Output_Structure!(out::Output_Structure)
    out.viral_loads = Dict{Int64, Dict{Float64, Float64}}()
    out.infections = Dict{Int64, Infection_Data}()
    out.n_infections = 0
end

function add_infection!(out::Output_Structure, idat::Infection_Data)
    out.n_infections += 1
    out.infections[out.n_infections] =  idat
end

function fill_Infection_Data!(out::Output_Structure, 
                              agent_id::Int64, 
                              infection::EntitiesTD.Infection_T, 
                              is_index::Bool,
                              t::Float64,
                              d::Int64,
                              s::Int64,
                              type::String)

    idat = Infection_Data()
    default_init_Infection_Data!(idat)
    
    idat.agent_id = agent_id 
    idat.agent_type = type 
    idat.time_infected = t
    idat.day_infected = d
    idat.shift_infected = s
    idat.beta_max = infection.beta_max
    idat.symptomatic = infection.symptomatic
    idat.is_index_case = is_index 

    add_infection!(out, idat)

end


function infection_linelist_to_DF(out::Output_Structure)::DataFrame

    infections = out.infections
    if haskey(infections, 1)
        df = Utilities.init_DF_from_struct(infections[1])
    else
        df = DataFrame() # return empty
        return df 
    end
    
    # add all 
    row_indices = sort(collect(keys(infections)))
    for row in row_indices
        df_i = Utilities.DF_row_from_struct(infections[row])
        append!(df, df_i)
    end


    return df

end


function linelist_to_timeseries(df::DataFrame,
    d_min::Int64,
    d_max::Int64)::DataFrame

    timeseries = DataFrame()
    timeseries[:, Symbol("new_infections")] = Array{Int64, 1}()
    timeseries[:, Symbol("day")] = Array{Int64, 1}()

    index = 0
    for d in d_min:d_max
        index += 1
        n_cases_d = sum(df.day_infected .== d)
        push!(timeseries[!, :new_infections], n_cases_d)
        push!(timeseries[!, :day], d)
    end

    return timeseries

end


function linelist_to_timeseries_nurse(df::DataFrame,
                                      d_min::Int64,
                                      d_max::Int64)::DataFrame

    nurse_filter(val::AbstractString) = val == "nurse"

    df_nurse = filter(:agent_type => nurse_filter, df)

    timeseries = DataFrame()
    timeseries[:, Symbol("new_infections")] = Array{Int64, 1}()
    timeseries[:, Symbol("day")] = Array{Int64, 1}()

    index = 0
    for d in d_min:d_max
        index += 1
        n_cases_d = sum(df_nurse.day_infected .== d)
        push!(timeseries[!, :new_infections], n_cases_d)
        push!(timeseries[!, :day], d)
    end

    return timeseries

end

function linelist_to_timeseries_patient(df::DataFrame,
                                        d_min::Int64,
                                        d_max::Int64)::DataFrame

    patient_filter(val::AbstractString) = val == "patient"

    df_patient = filter(:agent_type => patient_filter, df)

    timeseries = DataFrame()
    timeseries[:, Symbol("new_infections")] = Array{Int64, 1}()
    timeseries[:, Symbol("day")] = Array{Int64, 1}()

    index = 0
    for d in d_min:d_max
        index += 1
        n_cases_d = sum(df_patient.day_infected .== d)
        push!(timeseries[!, :new_infections], n_cases_d)
        push!(timeseries[!, :day], d)
    end

    return timeseries

end

function generate_infection_timeseries(op_struct::Output_Structure, d1::Int64, d2::Int64)::DataFrame

    linelist = infection_linelist_to_DF(op_struct)

    infection_timeseries = linelist_to_timeseries(linelist, d1, d2) 
    infection_timeseries_nurse = linelist_to_timeseries_nurse(linelist, d1, d2)
    infection_timeseries_patient = linelist_to_timeseries_patient(linelist, d1, d2);

    ts_out = DataFrame() 
    ts_out[:, :day] = infection_timeseries.day
    ts_out[:, :new_infections] = infection_timeseries.new_infections
    ts_out[:, :new_infections_nurse] = infection_timeseries_nurse.new_infections
    ts_out[:, :new_infections_patient] = infection_timeseries_patient.new_infections

    return ts_out 




end




function add_room_ids_to_q_vs_t!(out::Output_Structure, rooms::EntitiesTD.Rooms_T)

    for (room_id, room) in rooms
        out.viral_loads[room_id] = Dict{Float64, Float64}()
    end

end

function viral_load_timseries_to_df(out::Output_Structure)::DataFrame

    df_out = DataFrame() 


    for (room_id, ts) in out.viral_loads

        col_1 = "time_$room_id"
        col_2 = "q_$room_id"

        timesteps = sort(collect(keys(ts))) 
        q_vals = [ts[t] for t in timesteps]

        df_out[!, col_1] = timesteps 
        df_out[!, col_2] = q_vals

    end

    return df_out 


end


function update_zone_conc!(out::Output_Structure, zone_id::Int64, t::Float64, q::Float64)

    out.viral_loads[zone_id][t] = q

end


#end module 
end 