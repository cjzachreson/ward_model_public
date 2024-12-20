module AirflowDynamics

#using Plots
using Distributed

# for the deployment I'll use the structures defined in EntitiesTD (not needed for test case)
using Main.EntitiesTD
using Main.AirflowTransmissionInterface
#using Main.TransmissionDynamics

# notes: 
    # source of quanta from infected people: room.shedding_quanta            | Dict{Float64, Float64}              time -> quanta
    # viral quanta at each time t:           room.viral_quanta               | Dict{Float64, Float64}              time -> quanta 
    # flow network:                          Flow_Network.st_links_2_lambdas | Dict{Tuple{Int64, Int64}, Float64}  (source, target) -> flow rate (prop. of vol?)

    # need to find dq/dt for each room (q stands for quanta)
    # quanta concentration is q / room volume 

    # in any zone: 
    # dq/dt = [shedding rate + flow rate in] - [flow rate out + deactivation rate + exhaust rate]

    # for two zones A, B 

    # flow into A:

    # mass transport (constant density)
    # flow rate into A from B (B -> A, quanta/time ) = [q_B / vol_B] * lambda_BA
    # lambda is the volume of air in zone B at time t that moves into zone A by time t+dt 
    # units of lambda are volume/time (volumetric flow rate)

    # shedding rate into A = A.shedding_quanta / dt (NOTE: shedding_quanta is in units of quanta is computed over the interval dt from TransmissionDynamics)
    # see function TransmissionDynamics: compute_quanta (quanta = infection.beta_t * c.dt)

    # exhaust rate from A = q_A / vol_A * e_A 
    # e_A is the exhaust rate in volume/time 

    # deactivation rate in A = q_A * d_A 
    # d_A is a dimensionless deactivation rate (exponential)


    # initial states: 
    # v_A | volume of zone A  (m^3)
    # v_B | volume of zone B  (m^3)

    # q_A_i | initial quanta in zone A (quanta)
    # q_B_i | initial quanta in zone B (quanta) 

    # s_A | shedding rate in zone A (quanta/day)
    # s_B | shedding rate in zone B (quanta/day)

    # d_A | deactivation rate in zone A (dimensionless, proportion deactivated)
    # d_B | deactivation rate in zone B (dimensionless, proportion deactivated)

    # e_A | exhaust rate in zone A (volume / day) NOTE: for 6 ACH this is e_A = v_A (m^3 per change) * 6 (changes per hr) * 24 (hr per day) = 144 * v_A 
    # e_B | exhaust rate in zone B (volume / day)

    # lambda_AB | flow rate from zone A to zone B (volume / day)
    # lambda_BA | flow rate from zone B to zone A (volume / day)

    # dt | discrete time step (days)

    # updates: 

    #dq/dt = [shedding rate + flow rate in] - [flow rate out + deactivation rate + exhaust rate]

    # components: 
    # c_A = q_A / v_A | quanta concentration in zone A (quanta/vol)
    # c_B = q_B / v_B | quanta concentration in zone B (quanta/vol)

    # f_AB = lambda_AB * c_A | flow rate of quanta from A to B (quanta/time)
    # f_BA = lambda_BA * c_B | flow rate of quanta from B to A (quanta/time)
    # f_AE = e_A * c_A | flow rate of quanta from A to exhaust (quanta/time)
    # f_BE = e_B * c_B | flow rate of quanta from B to exhaust (quanta/time)

    # q_A(t + dt) = q_A(t) + [s_A + f_BA]*dt - [f_AB + f_AE + q_A * d_A]*dt
    # q_B(t + dt) = q_B(t) + [s_B + f_AB]*dt - [f_BA + f_BE + q_B * d_B]*dt

    #NOTE: no input flows are modelled, we assume ventilation input balances exhaust. 

function test_case(dt::Float64 #time resolution of test case 
                   )::Dict{AbstractString, Array{Float64, 1}} #output is a timeseries for each zone. 
    
    t_i = 0.0
    #t_f = 0.0035
    t_f = 0.1
    timepoints = t_i + dt:dt:t_f

    timeseries_out = Dict{AbstractString, Array{Float64, 1}}()

    # set up initial conditions for test case: 
    v_A = 10.0 # v_A | volume of zone A  (m^3)
    v_B = 10.0 # v_B | volume of zone B  (m^3)

    q_A_i = 0.0 # q_A_i | initial quanta in zone A (quanta)
    q_B_i = 0.0 # q_B_i | initial quanta in zone B (quanta) 

    s_A = 10.0 # s_A | shedding rate in zone A (quanta/day)
    s_B = 0.0 # s_B | shedding rate in zone B (quanta/day)

    d_A = 0.0 # d_A | deactivation rate in zone A (dimensionless, proportion deactivated)
    d_B = 0.0 # d_B | deactivation rate in zone B (dimensionless, proportion deactivated)

    e_A = 144.0 * v_A # e_A | exhaust rate in zone A (volume / day) NOTE: for 6 ACH this is e_A = v_A (m^3 per change) * 6 (changes per hr) * 24 (hr per day) = 144 * v_A 
    e_B = 144.0 * v_B # e_B | exhaust rate in zone B (volume / day)

    lambda_AB = 0.0#10.0 * v_B # lambda_AB | flow rate from zone A to zone B (volume / day)
    lambda_BA = 0.0#10.0 * v_B # lambda_BA | flow rate from zone B to zone A (volume / day)

    #NOTE: no input flows are modelled, we assume ventilation input balances exhaust. 

    #initialise: 
    timeseries_out["A"] = Array{Float64, 1}()
    timeseries_out["B"] = Array{Float64, 1}()

    push!(timeseries_out["A"], q_A_i)
    push!(timeseries_out["B"], q_B_i)

    for t in timepoints

        q_A = last(timeseries_out["A"])
        q_B = last(timeseries_out["B"])

        # components: 
        c_A = q_A / v_A #| quanta concentration in zone A (quanta/vol)
        c_B = q_B / v_B #| quanta concentration in zone B (quanta/vol)

        #NOTE: here, lambda is expressed as an absolute volume
        f_AB = lambda_AB * c_A #| flow rate of quanta from A to B (quanta/time)
        f_BA = lambda_BA * c_B #| flow rate of quanta from B to A (quanta/time)
        #NOTE: here, exhaust rate is expressed as absolute volume. 
        f_AE = e_A * c_A #| flow rate of quanta from A to exhaust (quanta/time)
        f_BE = e_B * c_B #| flow rate of quanta from B to exhaust (quanta/time)


        q_A_next = q_A + (s_A + f_BA)*dt - (f_AB + f_AE + q_A * d_A)*dt
        q_B_next = q_B + (s_B + f_AB)*dt - (f_BA + f_BE + q_B * d_B)*dt

        push!(timeseries_out["A"], q_A_next)
        push!(timeseries_out["B"], q_B_next)

    end

    return timeseries_out 

end

#TODO: conder if lambdas should be expressed as relative to source room volume
# to have some equivalence to ACH (i.e., 144 -> [super sink] equates to 6 ACH)
# input: config contains Rooms, Flow network, and timestep from scheduler
# Here, lambda is expressed in the input files as a multiple of source zone volume. 
# in the example above, it's expressed as an absolute volume. 
function quanta_flow_forward_euler!(c::Config_T,
                                    t_i::Float64, 
                                    t_f::Float64)::Dict{Int64, Array{Float64, 1}}#modifies 


    # set fine timescale 
    dt_fine = c.dt_days/10 #might not be small enough. 
    timepoints = t_i + dt_fine : dt_fine : t_f

    # set up fine-scale timeseries Dict{Int64, Array{Float64}}()
    # add a dict entry for each room (each gets its own timeseries)
    timeseries = Dict{Int64, Array{Float64}}()
    for (id, room) in c.rooms
        timeseries[id] = Array{Float64, 1}()
        #set initial values: 
        if haskey(room.viral_quanta, t_i)
            push!(timeseries[id], room.viral_quanta[t_i])
        else
            push!(timeseries[id], 0.0)
        end
    end
    # TODO: once stability is tested and assured
    # the timeseries arrays can be adjusted to store only the
    # most recent update (less memory intensive)


    q_next = Dict{Int64, Float64}()
    # iterate through time: 
    for t in timepoints 

        for (id, room) in c.rooms            # iterate through rooms r_i: 
            # solve q_r_i(t + dt)
            q_i = last(timeseries[id])
            v_i = room.volume
            c_i = q_i / v_i #| quanta concentration in zone A (quanta/vol)

            if haskey(room.shedding_quanta, t_i)
                s_i = room.shedding_quanta[t_i] / c.dt_days # TODO: remove the dt term here and from the initial calcualtion of shedding in TransmissionDynamics
            else
                s_i = 0.0
            end
            
            f_in = 0.0
            f_out = 0.0 

            for (st, lambda) in c.flow_network.st_links_2_lambdas # tuples of integers source id -> target id 
                # compute in-flows 
                if st[2] == id 
                    r_j = st[1]
                    q_j = last(timeseries[r_j])
                    v_j = c.rooms[r_j].volume 
                    c_j = q_j / v_j #| quanta concentration in zone B (quanta/vol)
                    lambda_ji = lambda
                    f_ji = lambda_ji * c_j * v_j #| flow rate of quanta from B to A (quanta/time)
                    #NOTE: here, lambda is expressed a a fraction of source volume per day. 
                    f_in += f_ji 
                
                # compute out-flow 
                elseif st[1] == id
                    lambda_ij = lambda
                    f_ij = lambda_ij * c_i * v_i#| flow rate of quanta from B to A (quanta/time)
                    # NOTE: here, lambda is expressed as a multiple of source volume per day. 
                    f_out += f_ij 
                end

            end

            # exhaust will be accounted for by 'super sink' target above 
            #f_AE = e_A * c_A #| flow rate of quanta from A to exhaust (quanta/time)
            # TODO: include deactivation terms too if necessary. 
            
            # update quanta for next fine-scale iteration: 
            q_next[id] = q_i + (s_i + f_in)*dt_fine - f_out*dt_fine
            

        end

        # update fine-scale timeseries: 
        for (id, rm) in c.rooms
            push!(timeseries[id], q_next[id])
            q_next[id] = 0.0 # zero next value. 
        end

    end

    # update rooms 
    for (id, room) in c.rooms
        room.viral_quanta[t_f] = last(timeseries[id])
    end

    return timeseries

end

function quanta_flow_forward_euler_optimised!(c::Config_T,
                                              t_i::Float64, 
                                              t_f::Float64,
                                              timescale_fac::Float64)::Dict{Int64, Float64}#modifies 


    # set fine timescale 
    dt_fine = Float64(c.dt_days/timescale_fac) #might not be small enough. 

    # add a dict entry for each room 
    q = Dict{Int64, Float64}()
    q_next = Dict{Int64, Float64}()
    s_i = Dict{Int64, Float64}()

    for (id, room) in c.rooms
        q[id] = 0.0
        q_next[id] = 0.0
        #set initial values: 
        if haskey(room.viral_quanta, t_i)
            q[id] = room.viral_quanta[t_i]
        else
            q[id] = 0.0
            #println("no time key for room $id at time $t_i")
        end
    end


    # zero shedding quanta for rooms without any:
    for (id, room) in c.rooms
        if !haskey(room.shedding_quanta, t_i)
            room.shedding_quanta[t_i] = 0.0 
        end
        s_i[id] = room.shedding_quanta[t_i]
        
    end

    # iterate through time: 

    t = t_i 

    while t < t_f 
        
        t += dt_fine 

        for (id, room) in c.rooms            # iterate through rooms r_i: 
            q_next[id] += s_i[id] * dt_fine  # add shedding
            q_next[id] += q[id] # add previous value 
        end

        for (st, lambda) in c.flow_network.st_links_2_lambdas # tuples of integers source id -> target id 
            
            s_id = st[1]
            
            if q[s_id] != 0.0
            
                t_id = st[2]
                
                flow = lambda * q[s_id] * dt_fine 

                q_next[s_id] -= flow 
                q_next[t_id] += flow
            end 

        end

        # update fine-scale timeseries: 
        for (id, room) in c.rooms
            q[id] = q_next[id]
            q_next[id] = 0.0 # zero next value. 
        end

    end


    # update rooms 
    for (id, room) in c.rooms
        room.viral_quanta[t_f] = q[id]
        c.output_structure.viral_loads[id][t_f] = q[id]
    end

    return q

end

function quanta_flow_forward_euler_optimised_DEBUG!(c::Config_T,
                                              t_i::Float64, 
                                              t_f::Float64,
                                              timescale_fac::Float64)::Dict{Int64, Float64}#modifies 


    # set fine timescale 
    dt_fine = Float64(c.dt_days/timescale_fac) #might not be small enough. 

    # add a dict entry for each room 
    q = Dict{Int64, Float64}()
    q_next = Dict{Int64, Float64}()
    s_i = Dict{Int64, Float64}()

    for (id, room) in c.rooms
        q[id] = 0.0
        q_next[id] = 0.0
        #set initial values: 
        if haskey(room.viral_quanta, t_i)
            q[id] = room.viral_quanta[t_i]
        else
            q[id] = 0.0
            println("no time key for room $id at time $t_i")
        end
    end


    # zero shedding quanta for rooms without any:
    for (id, room) in c.rooms
        if !haskey(room.shedding_quanta, t_i)
            room.shedding_quanta[t_i] = 0.0 
        end
        s_i[id] = room.shedding_quanta[t_i]
        
    end

    # iterate through time: 

    t = t_i 

    while t < t_f 
        
        t += dt_fine 

        for (id, room) in c.rooms            # iterate through rooms r_i: 
            q_next[id] += s_i[id] * dt_fine  # add shedding
            q_next[id] += q[id] # add previous value 
        end

        for (st, lambda) in c.flow_network.st_links_2_lambdas # tuples of integers source id -> target id 
            
            s_id = st[1]
            
            if q[s_id] != 0.0
            
                t_id = st[2]
                
                flow = lambda * q[s_id] * dt_fine 

                q_next[s_id] -= flow 
                q_next[t_id] += flow
            end 

        end

        # update fine-scale timeseries: 
        for (id, room) in c.rooms
            q[id] = q_next[id]
            q_next[id] = 0.0 # zero next value. 
        end

    end


    # update rooms 
    for (id, room) in c.rooms
        room.viral_quanta[t_f] = q[id]
        c.output_structure.viral_loads[id][t_f] = q[id]
    end

    return q

end



function quanta_flow_forward_euler_optimised_no_q_vs_t!(c::Config_T,
                                              t_i::Float64, 
                                              t_f::Float64,
                                              timescale_fac::Float64)::Dict{Int64, Float64}#modifies 


    # set fine timescale 
    dt_fine = Float64(c.dt_days/timescale_fac) #might not be small enough. 

    # add a dict entry for each room 
    q = Dict{Int64, Float64}()
    q_next = Dict{Int64, Float64}()
    s_i = Dict{Int64, Float64}()

    for (id, room) in c.rooms
        q[id] = 0.0
        q_next[id] = 0.0
        #set initial values: 
        if haskey(room.viral_quanta, t_i)
            q[id] = room.viral_quanta[t_i]
        else
            q[id] = 0.0
            #println("no time key for room $id at time $t_i")
        end
    end


    # zero shedding quanta for rooms without any:
    for (id, room) in c.rooms
        if !haskey(room.shedding_quanta, t_i)
            room.shedding_quanta[t_i] = 0.0 
        end
        s_i[id] = room.shedding_quanta[t_i]
        
    end

    # iterate through time: 

    t = t_i 

    while t < t_f 
        
        t += dt_fine 

        for (id, room) in c.rooms            # iterate through rooms r_i: 
            q_next[id] += s_i[id] * dt_fine  # add shedding
            q_next[id] += q[id] # add previous value 
        end

        for (st, lambda) in c.flow_network.st_links_2_lambdas # tuples of integers source id -> target id 
            
            s_id = st[1]
            
            if q[s_id] != 0.0
            
                t_id = st[2]
                
                flow = lambda * q[s_id] * dt_fine 

                q_next[s_id] -= flow 
                q_next[t_id] += flow
            end 

        end

        # update fine-scale timeseries: 
        for (id, room) in c.rooms
            q[id] = q_next[id]
            q_next[id] = 0.0 # zero next value. 
        end

    end


    # update rooms 
    for (id, room) in c.rooms
        room.viral_quanta[t_f] = q[id]
        #c.output_structure.viral_loads[id][t_f] = q[id]
    end

    return q

end




# function quanta_flow_forward_euler_parallel!(c::Config_T,
#                                               t_i::Float64, 
#                                               t_f::Float64,
#                                               timescale_fac::Float64)::Dict{Int64, Float64}#modifies 


    # set fine timescale 
    #dt_fine = Float64(c.dt_days/timescale_fac) #might not be small enough. 

#     # add a dict entry for each room 
#     q = Dict{Int64, Float64}()
#     q_next = Dict{Int64, Float64}()
#     s_i = Dict{Int64, Float64}()

#     for (id, room) in c.rooms
#         q[id] = 0.0
#         q_next[id] = 0.0
#         #set initial values: 
#         if haskey(room.viral_quanta, t_i)
#             q[id] = room.viral_quanta[t_i]
#         else
#             q[id] = 0.0
#         end
#     end


#     # zero shedding quanta for rooms without any:
#     for (id, room) in c.rooms
#         if !haskey(room.shedding_quanta, t_i)
#             room.shedding_quanta[t_i] = 0.0 
#         end
#         s_i[id] = room.shedding_quanta[t_i]
        
#     end

#     # iterate through time: 

#     t = t_i 

#     #edges = collect(keys(c.flow_network.st_links_2_lambdas))

#     while t < t_f 
        
#         t += dt_fine 

#         for (id, room) in c.rooms            # iterate through rooms r_i: 
#             q_next[id] += s_i[id] * dt_fine  # add shedding
#             q_next[id] += q[id] # add previous value 
#         end

#         Threads.@threads for (st, lambda) in c.flow_network.st_links_2_lambdas # tuples of integers source id -> target id 
            
#             #println(Threads.threadid())

#             #lambda = c.flow_network.st_links_2_lambdas[st]

#             s_id = st[1]
            
#             if q[s_id] != 0.0
            
#                 t_id = st[2]
                
#                 flow = lambda * q[s_id] * dt_fine 

#                 q_next[s_id] -= flow 
#                 q_next[t_id] += flow
#             end 

#         end

#         # update fine-scale timeseries: 
#         for (id, room) in c.rooms
#             q[id] = q_next[id]
#             q_next[id] = 0.0 # zero next value. 
#         end

#     end

#     # update rooms 
#     for (id, room) in c.rooms
#         room.viral_quanta[t_f] = q[id]
#     end

#     return q

# end



## end module 
end






# test module: 
# test_series = AirflowDynamics.test_case(0.0001)

# plot(test_series["A"])
# plot!(test_series["B"])