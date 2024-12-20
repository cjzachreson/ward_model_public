using DataFrames 
using CSV 


function add_zone_entry!(df::DataFrame,
                         zone_id::Int64,
                         zone_name::AbstractString,
                         volume::Float64,
                         adjacency_in::Array{Int64, 1},
                         lambdas_in::Array{Float64, 1} )

    #add room id
    push!(df[!, :room_id], zone_id)
    # add room name 
    push!(df[!, :name], zone_name)
    #volume (the only one that currently matters)
    push!(df[!, :volume], volume) # for airflow dynamics, all nodes must have volume (or nan values may be produced)
    #adjacency in (specified as a single edge from the super source)
    push!(df[!, :adjacency_in], adjacency_in) 

    # add back-flow from external zone to rooms: 
    #lambda - flow parameter from supersource
    push!(df[!, :lambda], lambdas_in)

    # placeholders: 
    #add width 
    push!(df[!, :width], NaN)
    #length 
    push!(df[!, :length], NaN)
    #height 
    push!(df[!, :height], NaN)

end

function main()


    #set constant parameters 
    Vfac_BR = 4.0 
    Vfac_EZ = 10.0

    # set parameter ranges for varied properties: 
    #proportion shared rooms: 
    p_vals = [0.0, 1.0]
    # air mixing rate (scaled to input ACH)
    f_vals = [0.0]
    # air changes for patient rooms 
    ACH_PR_vals = [6.0]
    # air changes for break room/lunch room
    ACH_BR_vals = [2.0]
    # factor relating volume of double room to vol of single room 
    # note v_double = Vfac * v_single 
    Vfac_double_vals = [1.5]#[1.0, 2.0]


    # set input-output directories 
    base_directory = dirname(@__FILE__)

    for layer in 1:3
        base_directory = dirname(base_directory)
    end

    input_dir = joinpath(base_directory, 
                        "WardDynamics",
                        "input",
                        "shared_room_ratio_experiments",
                        "room_tables",
                        "40_bed_ward")

    output_dir = joinpath(base_directory, 
                        "TransmissionDynamics", 
                        "input",
                        "shared_room_ratio_experiments",
                        "ward_struct",
                        "40_bed_ward")

    if !isdir(output_dir)
        mkpath(output_dir)
    end
    

    # add the break room id to the room id list 
    break_room_id = Int64(-10)
    # supersource id 
    source_id = Int64(-2) 
    #sink id 
    sink_id = Int64(-1) 
    # external zone id
    external_zone_id = Int64(-100)

    for p in p_vals 

        p_str = string(round(p, digits=2))
        p_label = "pShared_" * replace(p_str, '.' => 'p')

        for f in f_vals 

            f_str = string(round(f, digits=2))
            f_label = "f_" * replace(f_str, '.' => 'p')

            for ACH_patient_rooms in ACH_PR_vals

                ACH_pr_str = string(round(ACH_patient_rooms, digits=2))
                ACH_pr_label = "ACHpr_" * replace(ACH_pr_str, '.' => 'p')

                for ACH_break_room in ACH_BR_vals

                    ACH_br_str = string(round(ACH_break_room, digits = 2))
                    ACH_br_label = "ACHbr_" * replace(ACH_br_str, '.' => 'p')

                    for Vfac_double in Vfac_double_vals

                        Vfac_str = string(round(Vfac_double, digits = 2))
                        Vfac_label = "VfacDbl_" * replace(Vfac_str, '.' => 'p')


                        #TODO: include ACH and Vfacs in label (or nest if too extensive)
                        label = "ward_struct_" * 
                                 p_label * "_" * 
                                 f_label * "_" * 
                                 ACH_pr_label * "_" *
                                 ACH_br_label * "_" *
                                 Vfac_label 

                        output_fname = joinpath(output_dir, label * ".csv")

                        input_fname = joinpath(input_dir, "room_table_" * p_label * ".csv")

                        df_input = DataFrame(CSV.File(input_fname))

                        room_ids = unique(df_input.room_id) 

                        push!(room_ids, break_room_id) # add id to list of room ids 

                        df_output = DataFrame(
                                    room_id = Array{Int64, 1}(),
                                    name = Array{String, 1}(),
                                    width = Array{Float64, 1}(),
                                    length = Array{Float64, 1}(),
                                    height = Array{Float64, 1}(),
                                    volume = Array{Float64, 1}(),
                                    adjacency_in = Array{Array{Int64, 1}, 1}(),
                                    lambda = Array{Array{Float64, 1}, 1}() 
                        )



                        # set parameters 
                        # ACH_patient_rooms = 6.0
                        # ACH_break_room = 6.0
                        # Vfac_double = 1.0

                        Vfac_BR = 4.0 
                        Vfac_EZ = 10.0
                        
                        lambda_source_to_patient_rooms = ACH_patient_rooms * 24.0
                        lambda_source_to_break_room = ACH_break_room * 24.0     

                        #lambda_from_source_to_rooms = 144.0 # 6 ACH to rooms from source (relative to room volume)

                        room_volume_single = 25.0

                        room_volume_double = room_volume_single * Vfac_double
                        external_zone_volume = Vfac_EZ * room_volume_single
                        break_room_volume = Vfac_BR * room_volume_single
                        source_volume = room_volume_single

                        # set mixing parameter f 
                        lambda_backflow_external_zone_to_patient_rooms = lambda_source_to_patient_rooms * f # relative to room volume 
                        lambda_backflow_external_zone_to_break_room = lambda_source_to_break_room * f 
                        #NOTE: the above are relative to the source volume, which is adjusted below
                        # so that x ACH of room volume is the final input. 

                        # add rooms - adjust source input to room volume 
                        # also build lambda vectors for intput and output to/from external zone
                        # so flows can be balanced from external zone to sink (exhaust)
                        lambdas_rooms_to_external = Array{Float64, 1}()
                        flow_to_external = Array{Float64, 1}()
                        lambdas_external_to_rooms = Array{Float64, 1}()
                        flow_from_external_to_rooms = Array{Float64, 1}()
                        net_input = 0.0


                        #for good bookkeeping, start with source input to rooms (mechanical ventilation input)
                        # this maps to air changes for the ward rooms, wich is meant to exceed 6 ACH 

                        # add source 

                        add_zone_entry!(df_output,
                                    source_id,
                                    "source",
                                    source_volume,
                                    [source_id],
                                    [0.0] )



                        for r in room_ids

                            room_name = "patient room $r"
                            
                            n_beds = size(df_input.room_id[df_input.room_id .== r], 1) #check this code. 


                            # compute lambdas: 
                            # l_source is the ACH relative to source volume 
                            # l_ext is ACH relative to external zone volume 

                            if n_beds == 2 
                                v = room_volume_double
                                l_source = (v / source_volume) * lambda_source_to_patient_rooms
                                l_ext = (v / external_zone_volume) * lambda_backflow_external_zone_to_patient_rooms  
                                l_room_to_ext = lambda_source_to_patient_rooms * (1.0 + f)
                            elseif n_beds == 1 
                                v = room_volume_single
                                l_source = (v / source_volume) * lambda_source_to_patient_rooms 
                                l_ext = (v / external_zone_volume) * lambda_backflow_external_zone_to_patient_rooms  
                                l_room_to_ext = lambda_source_to_patient_rooms * (1.0 + f)
                            end

                            if r == break_room_id
                                room_name = "break room"
                                v = break_room_volume 
                                #input
                                l_source = (v / source_volume) * lambda_source_to_break_room
                                #backflow
                                l_ext = (v / external_zone_volume) * lambda_backflow_external_zone_to_break_room
                                # input + backflow (balanced, to external zone)
                                l_room_to_ext = lambda_source_to_break_room * (1.0 + f)
                            end

                            net_input += l_source * source_volume 

                            push!(lambdas_rooms_to_external, l_room_to_ext)
                            push!(lambdas_external_to_rooms, l_ext)

                            # keep track of net input and output for balancing into sink. 
                            push!(flow_to_external, l_room_to_ext * v)
                            push!(flow_from_external_to_rooms, l_ext * external_zone_volume)
                        
                            add_zone_entry!(df_output, # output df
                                            r, # zone id
                                            room_name, # zone name 
                                            v, # volume 
                                            [source_id, external_zone_id], # adjacency
                                            [l_source, l_ext] ) # lambda in 


                        end

                        # add external mixing zone (i.e. halls etc.)
                        # include input from each room, which is a function of the backflow rate
                        # and the input rate from the source to the rooms. 

                        add_zone_entry!(df_output, # output df
                                        external_zone_id, # zone id
                                        "external zone", # zone name 
                                        external_zone_volume, # volume 
                                        room_ids, # adjacency
                                        lambdas_rooms_to_external ) # lambda in 


                        # add sink (from external zone): 

                        #lambda - flow parameter from supersource
                        # has to balance the input volumes 
                        l_external_zone_to_sink = (sum(flow_to_external) - sum(flow_from_external_to_rooms)) / external_zone_volume

                        add_zone_entry!(df_output, # output df
                                        sink_id, # zone id
                                        "sink", # zone name 
                                        25.0, # volume 
                                        [external_zone_id], # adjacency
                                        [l_external_zone_to_sink] ) # lambda in 


                        net_output = l_external_zone_to_sink * external_zone_volume 

                        # final check, allowing for floating point errors: 
                        if abs((net_input - net_output))/net_input > 0.000001
                            println("WARNING: flows are not balanced")
                            break
                        end

                        CSV.write(output_fname, df_output, delim = '\t')
                    end
                end
            end
        end
    end

# end function
end 

main()
