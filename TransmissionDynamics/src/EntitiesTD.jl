#entities module for the transmission dynamics layer

module EntitiesTD

using DataFrames
using CSV
using Distributions
using StatsBase
using Random

# abstract types
abstract type RoomTD_T end
abstract type AgentTD_T end 
abstract type Flow_Network_T end
abstract type Pathogen_T end # contains hyperparameters that specify infections 
abstract type Infection_T end # can pull this in from RACF model. 
abstract type Event_Trajectories_T end 

# parameters: (note, these will move to config system)
# seed_infections = rand(1:10000000)
# println("rng seed infections = $seed_infections")
# rng_infections = MersenneTwister(seed_infections)


# concrete types 
Rooms_T = Dict{Int64, RoomTD_T}
Agents_T = Dict{Int64, AgentTD_T}

## Agent (i.e., person)
mutable struct AgentTD<:AgentTD_T
    id::Int64
    label::String
    current_room_id::Int64
    current_activity::String
    infections::Array{<:Infection_T, 1}
    is_index_case::Bool 
    AgentTD() = new()
end
function default_init_AgentTD!(a::AgentTD)
    a.id = -1
    a.label = "not defined"
    a.current_room_id = -3
    a.current_activity = "not defined"
    a.infections = Array{Infection, 1}()
    a.is_index_case = false
end
function init_AgentTD!(a::AgentTD, id::Int64, label::String)
    a.id = id
    a.label = label
end
function new_agent_activity!(a::AgentTD, act::String)
    a.current_activity = act
end
function new_agent_infection!(a::AgentTD, covid::Pathogen_T, t::Float64, rng::MersenneTwister)
    infection = Infection()
    init_COVID_Infection!(infection, covid, t, rng)
    #print_FoI_vs_t(infection, Array{Float64, 1}(1:16))
    push!(a.infections, infection)
    
end
function move_agent!(a::AgentTD, new_room_id::Int64, rooms::Rooms_T)
    remove_agent_from_room!(rooms[a.current_room_id], a)
    add_agent_to_room!(rooms[new_room_id], a)
    a.current_room_id = new_room_id
end
# checks if the ids have been initialised and adds any new ones
# note: should only be applied to agents of the same type (same label)
function add_new_agents!(ids::Array{Int64, 1}, agents::Agents_T, label::String)
    for id in ids 
        if !haskey(agents, id)
            #initialise new agent 
            new_agent = AgentTD()
            default_init_AgentTD!(new_agent)
            init_AgentTD!(new_agent, id, label)
            agents[id] = new_agent 
        end
    end
end



## Flow network - keeps track of which rooms share airflows 
mutable struct Flow_Network<:Flow_Network_T
    # lookup table for source-target links to contaminant transfer coefficients 
    st_links_2_lambdas::Dict{Tuple{Int64, Int64}, Float64}
    Flow_Network() = new()
end

function default_init_Flow_Network!(F::Flow_Network)
    F.st_links_2_lambdas = Dict{Tuple{Int64, Int64}, Float64}()
end

function init_Flow_Network!(F::Flow_Network, df::DataFrame)
    for row in eachrow(df)
        s_vec = parse.(Int64, split(strip(row.adjacency_in, ['[',']']), ','))
        l_vec = parse.(Float64, split(strip(row.lambda, ['[',']']), ',')) 
        t = parse.(Int64, row.room_id)
        for i in 1:size(s_vec, 1)
            s = s_vec[i]
            lambda = l_vec[i]
            st = (s, t)
            F.st_links_2_lambdas[st] = lambda 
        end
    end
end

## Room 
mutable struct RoomTD<:RoomTD_T
    id::Int64
    name::String
    agents::Agents_T 
    volume::Float64
    #flows::Flow_Network 
    #a timeseries time (days) -> quanta
    viral_quanta::Dict{Float64, Float64} 
    #time->quanta shed by agents present 
    shedding_quanta::Dict{Float64, Float64} #set as Int time (in minutes) -> shedding rate (quanta per day)
    #TODO: check stability of these key values when computing viral quanta timeseries 
    
    RoomTD() = new()

end
function default_init_RoomTD!(r::RoomTD)
    r.id = -3 #note: -1 and -2 are special (supersource and supersink)
    r.name = "not defined"
    r.agents = Agents_T()
    r.volume = 0.0
    #flows = Flow_Network()
    #default_init_Flow_Network!(flows)
    #r.flows = flows
    r.viral_quanta = Dict{Float64, Float64}()#these will need to be reset at each shift 
    r.shedding_quanta = Dict{Float64, Float64}()
end
function init_RoomTD!(r::RoomTD, id::Int64, name::AbstractString, volume::Float64)
    r.id = id
    r.name = name
    r.volume = volume
end
function add_flows_to_rooms!(rooms::Rooms_T, flow_network::Flow_Network)
    for (st, lambda) in flow_network.st_links_2_lambdas
        s_id = st[1]
        t_id = st[2]
        rooms[s_id].flows.st_links_2_lambdas[st] = lambda
        rooms[t_id].flows.st_links_2_lambdas[st] = lambda
    end
end



function add_agent_to_room!(r::RoomTD, a::AgentTD)
    r.agents[a.id] = a
end
function remove_agent_from_room!(r::RoomTD, a::AgentTD)
    delete!(r.agents, a.id)
end
function add_quanta_to_room!(r::RoomTD, q::Float64, t::Float64)
    if !haskey(r.viral_quanta, t)
        r.viral_quanta[t] = q
    else
        r.viral_quanta[t] += q
    end
end
function remove_quanta_from_room!(r::RoomTD, q::Float64, t::Float64)
    r.viral_quanta[t] -= q
    if r.viral_quanta[t] < 0
        println("WARNING: Negative quanta detected in room $(r.id)")
    end
end
function add_shedding_to_room!(r::RoomTD, q::Float64, t::Float64)
    if !haskey(r.shedding_quanta, t)
        r.shedding_quanta[t] = q
    else
        r.shedding_quanta[t] += q
    end
end
function clear_room(r::RoomTD)::RoomTD

    r_clear = r

    r_clear.viral_quanta = Dict{Float64, Float64}()
    r_clear.shedding_quanta = Dict{Float64, Float64}()

    return r_clear
end
function clear_rooms(rooms::Rooms_T)::Rooms_T
    rooms_clear = Rooms_T()
    for (r_id, r) in rooms
        rooms_clear[r_id] = clear_room(r)
    end
    return rooms_clear 
end



## pathogen - an independent entity
mutable struct Pathogen<:Pathogen_T 
    ## beta, global transmission scaler
    beta::Float64 # 1
    
    # symptomatic fraction
    p_asymp::Float64 # 2

    #latent period 
    latent_period::Float64 #here initialised as dt (for asynchronous update) # 3

    ## incubation period, log-normal 
    inc_mu::Float64 # 4
    inc_sig::Float64 # 5
    inc_dist::LogNormal{Float64} # 6

    ## post-incubation (recovery) period 
    rec_min::Float64 # 7
    rec_max::Float64 # 8
    rec_dist::Uniform{Float64} # 9
    
    ## distribution of infectiousness 
    b_dispersion::Float64 # 10
    b_scale::Float64 # 11
    b_dist::Gamma{Float64} # 12

    Vmax::Float64 # scaling factor for growth and decline of viral load # 13
    inc_plat_fac::Float64 # proportion of incubation period in plateau viral load # 14
    rec_plat_fac::Float64 # proportion of recovery period in plateau viral load # 15

    #give the pathogen a name, for example
    name::String # 16

    #ranges of test sensitivity function parameters: 
    c_lims::Vector{Float64} # 17
    b1_lims::Vector{Float64} # 18
    b2_lims::Vector{Float64} # 19
    b3_lims::Vector{Float64} # 20

    #null constructor 
    Pathogen() = new()
end
function default_init_Pathogen!(P::Pathogen)
    # initialise disease (SARS-CoV-2)

    # new code:
    P.beta = 0.0
    # if multiple disease are used, this will need to be replaced with 
    # a vector in the Setup module. 

    P.p_asymp = 0.0
    P.latent_period = 0.0
    P.inc_mu = 0.0
    P.inc_sig = 0.0
    P.rec_min = -1.0#must be less than rec_max 
    P.rec_max = 0.0
    P.b_dispersion = 0.1
    P.Vmax = 0.0
    P.inc_plat_fac = 0.0
    P.rec_plat_fac = 0.0

    #test sensitivity function parameters: 
    P.c_lims = [0.0, 0.0]
    P.b1_lims = [0.0, 0.0]
    P.b2_lims = [0.0, 0.0]
    P.b3_lims = [0.0, 0.0]

    ## incubation period, log-normal 
    P.inc_dist = LogNormal(P.inc_mu, P.inc_sig) #inc_dist::LogNormal{Float64}

    ## post-incubation (recovery) period 
    P.rec_dist = Uniform(P.rec_min, P.rec_max) #rec_dist::Uniform{Float64}
        
    ## distribution of infectiousness 
    P.b_scale = P.beta / P.b_dispersion  #b_scale::Float64
    P.b_dist = Gamma(P.b_dispersion, P.beta/P.b_dispersion) #b_dist::Gamma{Float64}

    P.name = "Default"
end
function init_SARS_CoV_2!(P::Pathogen, beta::Float64)
    # initialise disease (SARS-CoV-2)
    # all parameters other than the transmission scalar
    # are taken from Zachreson et al. 2022 (Science advances 8 (14), eabm3624)

    # new code:
    #NOTE: in previous implementations, beta was a 'force of infection'
    # here, it's equivalent to expected maximum (replication competent) quanta emission per minute. 
    
    P.beta = beta #say this is <max> rate of infection per day 



    P.p_asymp = 0.33
    P.latent_period = 0.001 #in days #config.dt
    P.inc_mu = 1.62
    P.inc_sig = 0.418
    P.rec_min = 5.0
    P.rec_max = 10.0

    #for realism, dispersion between 0.1 and 0.15
    P.b_dispersion = 0.15
    #P.b_dispersion = 3.0


    P.Vmax = 7.0
    P.inc_plat_fac = 0.1
    P.rec_plat_fac = 0.0

    #test sensitivity function parameters: 
    P.c_lims = [1.0, 5.11]
    P.b1_lims = [0.8, 2.31]
    P.b2_lims = [1.26, 3.47]
    P.b3_lims = [1.05, 1.14]

    ## incubation period, log-normal 
    P.inc_dist = LogNormal(P.inc_mu, P.inc_sig) #inc_dist::LogNormal{Float64}

    ## post-incubation (recovery) period 
    P.rec_dist = Uniform(P.rec_min, P.rec_max) #rec_dist::Uniform{Float64}
        
    ## distribution of infectiousness 
    P.b_scale = P.beta / P.b_dispersion  #b_scale::Float64
    P.b_dist = Gamma(P.b_dispersion, P.beta/P.b_dispersion) #b_dist::Gamma{Float64}

    P.name = "SARS-CoV-2"
end


## infection - something an agent can have 
mutable struct Infection<:Infection_T 
    ## gets created during call to infect_agent!()
    ## stores properties and time-dependent characteristics of individual infections

    #pathogen object: 
    pathogen::Pathogen_T
    #pathogen name:
    pathogen_name::String
    # time since infection
    t_infected::Float64

    #latent period 
    t_latent::Float64

    # incubation period (time between exposure and symptom expression)
    t_inc::Float64
    # cdf of t_inc from incubation period dist, used to enforce correlations with test sensitivity and infectiousness 
    q_inc::Float64 
    # recovery period (time between peak viral load and recovery)
    t_rec::Float64
    # symptom expression (bool symptomatic or not) "Will they express symptoms?"
    symptomatic::Bool 
    # symptom expression (bool expressing symptoms or not) "are they currently expressing symptoms?"
    expressing_symptoms::Bool 

    # force of infection (i.e., baseline infection transmission rate as a function of time)
    beta_t::Float64 
    
    ## parameters determining force of infection
    #maximum infectiousness (peaks at symptom onset)
    beta_max::Float64 
    beta_min::Float64 #Def: beta_max/Vmax

    k_inc::Float64 #rate of exponential increase of infectiousness during incubation 
    k_rec::Float64 #rate of exponential decrease of infectiousness during recovery 

    # parameters determining detection probability
    # RAT detection probability, NOTE: decide whether this is an infection property, or something else... 
    # for now, I'm going to include these as infection properties, but they could potentially be implemented 
    # separately... 
    test_sensitivity::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    changepoint::Float64

    time_last_updated::Float64 # in days since the simulation started 

    recovered::Bool # has the infection been cleared?

    #null constructor: 
    Infection() = new()
end
function init_COVID_Infection!(I::Infection, P::Pathogen_T, t::Float64, rng_infections::MersenneTwister)
        #NOTE: rng_infections is used. 
        I.t_inc = rand(rng_infections, P.inc_dist)
        I.q_inc = cdf(P.inc_dist, I.t_inc)
        I.t_rec = rand(rng_infections, P.rec_dist)#t_rec::Float64
        I.beta_max = rand(rng_infections, P.b_dist)#beta_max::Float64

        
       #println("beta max = $(I.beta_max)")

        #
    
        ## pathogen reference
        I.pathogen = P
        ## name of pathogen 
        I.pathogen_name = P.name # pathogen::String 
        I.t_infected = 0.0 #t_infected::Float64
        ## latent period
        I.t_latent = P.latent_period # TODO: draw from distribution as with other intervals
        ## incubation period (time between exposure and symptom expression)
        #infection.t_inc = rand(config.rng_infections, pathogen.inc_dist)#t_inc::Float64
        #infection.q_inc = cdf(pathogen.inc_dist, infection.t_inc)#q_inc::Float64 NOTE: not sure if it will let me do this, need to check. *** Nope - need to pass as input to constructor. 
        ## recovery period (time between peak viral load and recovery)
        #infection.t_rec = rand(config.rng_infections, pathogen.rec_dist)#t_rec::Float64
        ## symptom expression (bool symptomatic or not) "Will they express symptoms?"
        I.symptomatic = rand(rng_infections) < (1.0 - P.p_asymp)#symptomatic::Bool 
        ## symptom expression (bool expressing symptoms or not) "are they currently expressing symptoms?"
        I.expressing_symptoms = false #expressing_symptoms::Bool 
        ## force of infection (i.e., baseline infection transmission rate as a function of time)
        I.beta_t = 0.0#beta_t::Float64 
        
        ### parameters determining force of infection
        ##maximum infectiousness (peaks at symptom onset)
        #infection.beta_max = rand(config.rng_infections, pathogen.b_dist)#beta_max::Float64 
        I.beta_min = I.beta_max / P.Vmax #beta_min::Float64 #Def: beta_max/Vmax #[NOTE: again, not sure if Julia will let me initialise in this way (does it know the value of beta_max?)] ***
    
        I.k_inc = compute_kinc(I.t_inc, P)#k_inc::Float64 #rate of exponential increase of infectiousness during incubation *** (may need to initialise these after construction)
        I.k_rec = compute_krec(I.t_rec, P)#k_rec::Float64 #rate of exponential decrease of infectiousness during recovery ***
    
        ## parameters determining detection probability
        ## RAT or PCR detection probability, NOTE: decide whether this is an infection property, or something else... 
        ## for now, I'm going to include these as infection properties, but they could potentially be implemented 
        ## separately... 
        
        I.test_sensitivity = 0.0#test_sensitivity::Float64 (initially set to 0)
        
        I.b1 = P.b1_lims[1] + 
                        rand(rng_infections) * (P.b1_lims[2] - P.b1_lims[1])#b1::Float64
        
        I.b2 = P.b2_lims[1] + 
                        (1.0 - I.q_inc) * (P.b2_lims[2] - P.b2_lims[1])#b2::Float64 ***quantile matching with incubation period. 
        
        I.b3 = P.b3_lims[1] + 
                        rand(rng_infections) * (P.b3_lims[2] - P.b3_lims[1])#b3::Float64
        
        I.changepoint = I.t_inc - 
                        (I.q_inc * (P.c_lims[2] - P.c_lims[1]))#changepoint::Float64

        I.time_last_updated = t

        I.recovered = false 
end


function step_infection!(I::Infection, t::Float64)
        dt = t - I.time_last_updated

        # increase time since infection
        I.t_infected += dt
    
        # check recovery status
        # what to do when an individual recovers? 
        # agent's immunity status should be updated, and the infection should 
        # be removed from memory. 
        if I.t_infected > (I.t_inc + I.t_rec) 
            I.recovered = true
        end
        
        # check symptom expression
        if (I.t_infected > I.t_inc) && I.symptomatic 
            I.expressing_symptoms = true
        end 
    
        # update force of infection 
        # TODO: allow asymptomatic cases to be less infectious?
        compute_FoI!(I)
    
        # update detection probability
        #TODO: don't need to do this every update, only if a test is about to be performed. 
        compute_test_sensitivity!(I)
    
        I.time_last_updated = t

end
#initialisation helper functions for computing some derived infection parameters: 
# maps t_plat_fac, t_inc (incubation period of agent i) 
#and Vmax to the growth rate of infectiousness during incubation
function compute_kinc(t_inc::Float64, pathogen::Pathogen_T) 
    t_plat_inc = pathogen.inc_plat_fac * t_inc
    k_inc = log(pathogen.Vmax)/(t_inc - t_plat_inc)
    return k_inc
end
# maps to the decay rate of infectiousness during 'recovery' i.e.,
# after symptom onset
function compute_krec(t_rec::Float64, pathogen::Pathogen_T) #t_rec 
    t_plat_rec = pathogen.rec_plat_fac * t_rec
    k_rec = log(1.0/pathogen.Vmax) / (t_rec - t_plat_rec)
    return k_rec
end
## update functions for computing time-dependent infection properties 
# compute FoI 
function compute_FoI!(I::Infection_T)

    if I.recovered
        I.beta_t = 0.0
    else 

        t_pt = 0.0 #local variable used for piecewise timepoint 
        #infection.beta_t = ?  #update force of infection based on t_infected 

        #check if latent: 
        if I.t_infected <= I.t_latent #latent 
            t_pt = I.t_infected
            I.beta_t = I.beta_min #note - this gets rescaled to 0 at the end 
        #check if agent is incubating 
        elseif I.t_infected <= I.t_inc # incubating
            t_pt = I.t_infected
            # fixed this - consider checking effects on previous results. 
            I.beta_t = I.beta_max * exp(I.k_inc * t_pt)/I.pathogen.Vmax

        # check if recovering 
        elseif I.t_infected <= (I.t_rec + I.t_inc) # recovering (i.e., symptomatic period if symptomatic)
            t_pt = I.t_infected - I.t_inc 

            I.beta_t = I.beta_max * exp(I.k_rec * t_pt)
        #if none of the above, then recovered. 
        else  
            I.beta_t = I.beta_min #note - this gets scaled to 0 in the rescaling step. 
        end

        #plateau implemented as cuttoff: 
        if I.beta_t > I.beta_max
            I.beta_t = I.beta_max
        end

        # rescale (missed this in the previous implementation)
        I.beta_t = (I.beta_t - I.beta_min)/(I.beta_max - I.beta_min) * I.beta_max


    end

end

function print_FoI_vs_t(I::Infection_T, timepoints::Array{Float64,})::DataFrame
    I_c = deepcopy(I)
    FoI_vs_t = DataFrame() 
    FoI_vals = Array{Float64, 1}()
    for t in timepoints

        I_c.t_infected = t
        compute_FoI!(I_c)
        push!(FoI_vals, I_c.beta_t)
    end
    FoI_vs_t.t = timepoints
    FoI_vs_t.beta_t = FoI_vals 
    println(FoI_vs_t)
    return FoI_vs_t
end


# compute test sensitivity 
function compute_test_sensitivity!(I::Infection_T)
    tau = I.t_infected - I.changepoint
    if tau < 0
        I.test_sensitivity = 
            1.0 / (1.0 + exp(-1.0 * 
            (I.b1 + I.b2 * tau )))
    else
        I.test_sensitivity = 
            1.0 / (1.0 + exp(-1.0 * 
            (I.b1 + I.b2 * tau + (-1.0 * I.b2*I.b3*tau) )))
    end
end




# link flow network data into rooms - this will speed up iteration
# during the airflow solve. 


function read_flow_network_from_file(fname::String)::Flow_Network
    input_table = DataFrame(CSV.File(fname, delim = "\t", type = String)) #string formatting for vector values
    flow_network = Flow_Network()
    default_init_Flow_Network!(flow_network)
    init_Flow_Network!(flow_network, input_table)
    return flow_network
        #parse.(Int64, split(workers_DF.roster[i], ','))
end

function read_rooms_from_file(fname::String)::Rooms_T
    input_table = DataFrame(CSV.File(fname, delim = "\t")) #string formatting for vector values
    rooms = Rooms_T()
    for r in 1:size(input_table, 1)
        room = RoomTD()
        default_init_RoomTD!(room)
        init_RoomTD!(room, 
                     input_table.room_id[r], 
                     input_table.name[r], 
                     Float64(input_table.volume[r]))
        rooms[input_table.room_id[r]] = room 
    end
    return rooms
end

# TODO: refactor so DataFrames are not used as dynamic objects. 
## trajectories - contains the list of events each agent is 
# involved in. 
mutable struct Shift_Trajectories<:Event_Trajectories_T  
    agent_trajectories::Dict{Int64, DataFrame}
    Shift_Trajectories() = new()
end

function extract_shift_trajectories(all::DataFrame, day::Int64, shift::Int64)::Shift_Trajectories
    day_df = all[all.start_day .== day, :]
    shift_df = day_df[day_df.start_shift .== shift, :]
    agent_ids = unique(shift_df.agent_id)
    shift_trajectories = Shift_Trajectories()
    shift_trajectories.agent_trajectories = Dict{Int64, DataFrame}()
    for id in agent_ids 
        # TODO: refactor dataframes (13% of runtime)
        df_i = shift_df[shift_df.agent_id .== id, :]
        shift_trajectories.agent_trajectories[id] = df_i
    end
    return shift_trajectories
end


## end module 
end