
include("init_test.jl")


# output label 
output_label = "full_run_pShared_N4_same_2024_11_25"


# need to input input and output directories as arguments 

function configure_constants(model_type::String, 
                             index_case_type::String,
                             experiment_type::String,
                             fname_ward_structure::String,
                             dirname_trajectories::String,
                             trajectory_label::String,
                             output_dirname::String,
                             detailed_output_flag::Bool)::TransmissionDynamics.Experiment_Constants

    # altering for cross-platform compatibility: 
    #input and output directories: 

    if !isdir(output_dirname)
        mkpath(output_dirname)
    end

    #event lists: 
    fname_trajectories_hcw = joinpath(dirname_trajectories, 
                                      "trajectory_hcw_" * trajectory_label * ".csv")
    fname_trajectories_patient = joinpath(dirname_trajectories, 
                                          "trajectory_patient_" * trajectory_label * ".csv")

    ex = TransmissionDynamics.Experiment_Constants()
    TransmissionDynamics.init_Experiment_Constants!(ex,
                                                   fname_trajectories_hcw, 
                                                   fname_trajectories_patient, 
                                                   fname_ward_structure,
                                                   output_dirname,
                                                   model_type,
                                                   index_case_type,
                                                   experiment_type,
                                                   detailed_output_flag)

    return ex 

end

# the configure funcion can be customised for different experiments,
# all other functions and modules should be general (including flags
# and multiple dispatch etc., to handle different special cases.) 
function configure_run(seed_infections::Int64, 
                       seed_index_cases::Int64,
                       beta::Float64, 
                       ex::TransmissionDynamics.Experiment_Constants,
                       q_vs_t_flag::Bool)::TransmissionDynamics.Config
    
    output_dirname = ex.output_dirname 

    config = TransmissionDynamics.Config()



    #TODO: consider splitting this into separate functions, 
    # a minimal initialiser with default args, 
    # followed by custom modifiers for different experiments 
    TransmissionDynamics.init_config!(config, 
                                      ex,
                                      output_dirname,
                                      seed_infections,
                                      seed_index_cases,
                                      beta,
                                      q_vs_t_flag)
    
    if ex.experiment_type == "R0"
        config.R0_calib_flag = true
    else
        config.R0_calib_flag = false
    end 

    ## TODO: undo this after detailed test. 
    #config.n_days = 1
    
    return config 
end

function configure_run_index_IDs(seed_infections::Int64, 
                                seed_index_cases::Int64,
                                beta::Float64, 
                                ex::TransmissionDynamics.Experiment_Constants,
                                q_vs_t_flag::Bool,
                                index_case_ids_hcws::Array{Int64, 1},
                                index_case_ids_patients::Array{Int64, 1})::TransmissionDynamics.Config

    output_dirname = ex.output_dirname 

    config = TransmissionDynamics.Config()



    #TODO: consider splitting this into separate functions, 
    # a minimal initialiser with default args, 
    # followed by custom modifiers for different experiments 
    TransmissionDynamics.init_config!(config, 
                                      ex,
                                      output_dirname,
                                      seed_infections,
                                      seed_index_cases,
                                      beta,
                                      q_vs_t_flag)

    if ex.experiment_type == "R0"
        config.R0_calib_flag = true
    else
        config.R0_calib_flag = false
    end 

    # set specific index cases to infect: 
    config.index_case_ids_hcws = index_case_ids_hcws
    config.index_case_ids_patients = index_case_ids_patients

    return config 
end




function run!(config)

    #println("generating agents")
    hcws = TransmissionDynamics.SIR_Agents()
    patients = TransmissionDynamics.SIR_Agents()

    #println("initialising agents")
    TransmissionDynamics.default_init_SIR_Agents!(hcws)
    TransmissionDynamics.default_init_SIR_Agents!(patients)

    TransmissionDynamics.init_SIR_Agents!(hcws, config.hcw_ids, "nurse")
    TransmissionDynamics.init_SIR_Agents!(patients, config.patient_ids, "patient")
    # needed to ensure patients are always registered to a location even if they're not involved in events
    # on a given shift. 

    TransmissionDynamics.init_Patient_Locations!(patients, config.trajectories_patient)

    #TODO: make the functions below more flexible and transparent (setting initial conditions will be essential)
    #println("setting initial infections") - use the index case type string as a function argument
    # rather than as a function selector. 
    if config.index_case_type == "nurse"
        TransmissionDynamics.initialise_infections_hcw!(config, hcws, 1)
        TransmissionDynamics.initialise_infections_patient!(config, patients, 0)
    elseif config.index_case_type == "patient"
        TransmissionDynamics.initialise_infections_patient!(config, patients, 1)
        TransmissionDynamics.initialise_infections_hcw!(config, hcws, 0)
    else
        prinln("Error: index case not specified!")
    end

    TransmissionDynamics.simulate_transmission_forward_euler!(config, hcws, patients)
    #TransmissionDynamics.simulate_transmission!(config, hcws, patients)

end

function run_index_IDs!(config)

    #println("generating agents")
    hcws = TransmissionDynamics.SIR_Agents()
    patients = TransmissionDynamics.SIR_Agents()

    #println("initialising agents")
    TransmissionDynamics.default_init_SIR_Agents!(hcws)
    TransmissionDynamics.default_init_SIR_Agents!(patients)

    TransmissionDynamics.init_SIR_Agents!(hcws, config.hcw_ids, "nurse")
    TransmissionDynamics.init_SIR_Agents!(patients, config.patient_ids, "patient")
    # needed to ensure patients are always registered to a location even if they're not involved in events
    # on a given shift. 

    TransmissionDynamics.init_Patient_Locations!(patients, config.trajectories_patient)

    #TODO: make the functions below more flexible and transparent (setting initial conditions will be essential)
    #println("setting initial infections") - use the index case type string as a function argument
    # rather than as a function selector. 
    if !isempty(config.index_case_ids_hcws)
        TransmissionDynamics.initialise_infections_by_ID!(config, 
                                                          hcws, 
                                                          config.index_case_ids_hcws)
        TransmissionDynamics.initialise_infections_patient!(config, patients, 0)                                                   
    elseif !isempty(config.index_case_ids_patients)
        TransmissionDynamics.initialise_infections_by_ID!(config, 
                                                          patients, 
                                                          config.index_case_ids_patients)
    
        TransmissionDynamics.initialise_infections_hcw!(config, hcws, 0)
    else
        prinln("Error: index case not specified!")
    end

    TransmissionDynamics.simulate_transmission_forward_euler!(config, hcws, patients)
    #TransmissionDynamics.simulate_transmission!(config, hcws, patients)

end


function main(fname_tag::AbstractString) 


    base_dirname = dirname(@__FILE__)
    for i in 1:2
        base_dirname = dirname(base_dirname) #should be 'ward_model'
    end

    dirname_ward_struct = joinpath(base_dirname,
                                    "TransmissionDynamics",
                                    "input",
                                    "shared_room_ratio_experiments",
                                    "ward_struct",
                                    "4_bed_ward")

    dirname_trajectories = joinpath(base_dirname, 
                                "WardDynamics", 
                                "output", 
                                "shared_room_ratio", 
                                "4_bed_ward")


    n_runs = 5000#2500
    # make sure this is enforced. 
    # the ids of index cases should be specified 
    # in the configuration file because these will be 
    # crucial for the initial conditions. 
    n_index_cases = 2

    Random.seed!(4)

    # seeds for rng_infections 
    #seeds = rand(1:1000000, n_runs)
    betas = [10000.0]#[0.15]


    index_case_types = ["patient"]

    index_case_ids_patients = [1, 2] # single rooms, or double rooms same room
    #index_case_ids_patients = [1, 3] # double rooms, different rooms
    index_case_ids_hcws = Array{Int64, 1}() 

    #println("seed:  $seeds")

    # specify experiment type: 
    #experiment_type = "R0"
    #experiment_type = "R0"
    # iterate through experiment scenarios: 

    model_type = "heterogeneous"

    q_vs_t_flag = false # flag saying to write quanta timeseries to output (high storage requirement)
    detailed_output_flag = false # flag saying to write detailed timeseries output for each run. 

    for experiment_type in ["R0"]#, "OB"]

        for index_case_type in index_case_types 

            p_vals = [0.0, 1.0]
            f_vals = [0.0]
            ACHpr_vals = [6.0]
            ACHbr_vals = [2.0]
            VfacDbl_vals = [1.5]

            df_out =  DataFrame(
                p = Array{Float64, 1}(),
                f = Array{Float64, 1}(),
                beta = Array{Float64, 1}(),
                attack_rate = Array{Float64, 1}(),
                attack_rate_patient = Array{Float64, 1}(),
                attack_rate_nurse = Array{Float64, 1}()
            )

            tag = fname_tag

            output_label = experiment_type * "_p_vs_f_" * index_case_type * "_N$n_runs" *"_"*tag

            output_dirname = joinpath(base_dirname, 
                                    "TransmissionDynamics", 
                                    "output", 
                                    "shared_room_experiments",
                                    output_label)

            output_dirname_detailed = joinpath(output_dirname, "detailed_output")
            if !isdir(output_dirname_detailed)
                mkpath(output_dirname_detailed)
            end


            output_fname = joinpath(output_dirname, 
                                    "AR_summary.csv")

            for p in p_vals 

                p_str = replace(string(round(p, digits=2)), '.' => 'p')

                trajectory_label = "pShared_" * p_str

                for f in f_vals 

                    f_str = replace(string(round(f, digits=2)), '.' => 'p')

                    for ACHpr in ACHpr_vals

                        ACHpr_str = replace(string(round(ACHpr, digits=2)), '.' => 'p')

                        for ACHbr in ACHbr_vals

                            ACHbr_str = replace(string(round(ACHbr, digits=2)), '.' => 'p')

                            for VfacDbl in VfacDbl_vals

                                VfacDbl_str = replace(string(round(VfacDbl, digits = 2)), '.' => 'p')

                                struct_label = "pShared_" * p_str * 
                                               "_f_" * f_str * 
                                               "_ACHpr_" * ACHpr_str *
                                               "_ACHbr_" * ACHbr_str *
                                               "_VfacDbl_" * VfacDbl_str 

                                fname_ward_structure_pf = joinpath(dirname_ward_struct,
                                                                "ward_struct_" * struct_label * ".csv")
                                                                
                                
                                dirname_trajectories_pf = joinpath(dirname_trajectories, 
                                                                    "pShared_" * p_str)

                                ex = configure_constants(model_type, 
                                                        index_case_type, 
                                                        experiment_type,
                                                        fname_ward_structure_pf,
                                                        dirname_trajectories_pf,
                                                        trajectory_label,
                                                        output_dirname,
                                                        detailed_output_flag)



                                #need to spit out a comprehensive, human-readable config file too. 

                                for b in betas 

                                    
                                    beta_label = "beta_$(replace("$b", "." => "p" ))"
                                    
                                    println("**** beta = $b, p_shared = $p, flowback = $f ****")

                                    ar = ones(n_runs)
                                    ar_patient = ones(n_runs)
                                    ar_nurse = ones(n_runs)

                                    batch_label = experiment_type * 
                                                "_" * index_case_type *
                                                "_p_" * p_str *
                                                "_f_" * f_str *
                                                "_" * beta_label *
                                                "_N$n_runs" 

                                    output_dirname_batch = joinpath(output_dirname_detailed, "runs_" * batch_label)

                                    if ex.detailed_output_flag
                                        if !isdir(output_dirname_batch)
                                            mkdir(output_dirname_batch)
                                        end
                                    end

                                    Threads.@threads :static for i in 1:n_runs 
                                    #for i in 1:n_runs

                                        seed_i = Random.rand(1:1000000)

                                        println("run: $i, seed: $(seed_i), thread id: $(Threads.threadid())")
                                        
                                        c = configure_run_index_IDs(seed_i, 
                                                                    seed_i, 
                                                                    b, 
                                                                    ex, 
                                                                    q_vs_t_flag,
                                                                    index_case_ids_hcws,
                                                                    index_case_ids_patients)

                                        if ex.detailed_output_flag
                                            output_dirname_run = joinpath(output_dirname_batch, "run_$i")
                                            if !isdir(output_dirname_run)
                                                mkdir(output_dirname_run)
                                            end
                                        end

                                        run_index_IDs!(c)

                                        infection_linelist = 
                                        OutputStructures.infection_linelist_to_DF(c.output_structure)

                                        ts_out = OutputStructures.generate_infection_timeseries(c.output_structure, 1, c.n_days)

                                        attack_rate = sum(ts_out.new_infections) - n_index_cases
                                        attack_rate_nurse = sum(ts_out.new_infections_nurse)
                                        attack_rate_patient = sum(ts_out.new_infections_patient)

                                        #TODO: generalise to mixed initialisation 
                                        if index_case_type == "nurse"
                                            attack_rate_nurse -= n_index_cases
                                        end

                                        if index_case_type == "patient"
                                            attack_rate_patient -= n_index_cases
                                        end

                                        println("attack rate: $attack_rate")
                                        #println("attack rate (nurses): $attack_rate_nurse")
                                        #println("attack rate (patients): $attack_rate_patient")

                                        #display(plot(case_timeseries.day, case_timeseries.new_infections))

                                        ar[i] = attack_rate
                                        ar_nurse[i] = attack_rate_nurse
                                        ar_patient[i] = attack_rate_patient 

                                        if ex.detailed_output_flag
                                            if q_vs_t_flag 
                                                df_q_vs_t = OutputStructures.viral_load_timseries_to_df(c.output_structure)
                                                output_fname_q_vs_t = joinpath(output_dirname_run, "q_vs_t.csv")
                                                CSV.write(output_fname_q_vs_t, df_q_vs_t)
                                            end

                                        
                                            output_fname_linelist = joinpath(output_dirname_run, "infection_linelist.csv")
                                            output_fname_infection_ts = joinpath(output_dirname_run, "infection_timeseries.csv")

                                            CSV.write(output_fname_linelist, infection_linelist)
                                            CSV.write(output_fname_infection_ts, ts_out)

                                        end

                                        


                                    end

                                    output_fname_pf = joinpath(output_dirname_detailed, 
                                                                struct_label *
                                                                "_" * beta_label * ".csv")

                                    df_out_pf = DataFrame(AR = ar, AR_nurse = ar_nurse, AR_patient = ar_patient)

                                    CSV.write(output_fname_pf, df_out_pf)

                                    println("**** beta = $b, mean AR = $(mean(ar)) ****")

                                    push!(df_out[!, :attack_rate], mean(ar))
                                    push!(df_out[!, :attack_rate_nurse], mean(ar_nurse))
                                    push!(df_out[!, :attack_rate_patient], mean(ar_patient))
                                    push!(df_out[!, :beta], b)
                                    push!(df_out[!, :p], p)
                                    push!(df_out[!, :f], f)


                                end
                            end
                        end
                    end

                end

            end

            CSV.write(output_fname, df_out)

        end

    end

end

main(output_label)
