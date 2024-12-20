

include("./setup.jl")


# #####

input_dirname = joinpath(base_directory,
                        "input",
                        "shared_room_ratio_experiments")

p_vals = [0.0, 1.0]


for p in p_vals

    p_str = string(round(p, digits=2))

    p_str = replace(p_str, '.' => 'p' )

    run_label = "pShared_" * p_str

    output_dirname = joinpath(base_directory,
                    "output",
                    "shared_room_ratio",
                    "4_bed_ward",
                    run_label)

    if !isdir(output_dirname)
        mkpath(output_dirname)
    end
    
    input_fname = joinpath(input_dirname,
                "room_tables",
                "4_bed_ward",
                "room_table_" * run_label * ".csv")

    #### run parameters ####
    
    homogeneous_flag = false 
    admissions_per_day = 4
    
    # mean and variance
    length_of_stay_mean = 42 # 42 shifts => 14 days 
    #length_of_stay_var = 81 # std is 9 shifts, 3 days. 
    
    ## if we want a constant: 
    # Uniform (will be rounded up to Dirac delta dist. in UpdateWard) 
    length_of_stay = Uniform(length_of_stay_mean-0.5, length_of_stay_mean)
    
    ## Gamma dist. 
    #Gamma (will be rounded up to integer values after sampling)
    #shp = (length_of_stay_mean^2) / length_of_stay_var 
    #scl = length_of_stay_var / length_of_stay_mean
    
    #length_of_stay = Gamma(shp, scl)
    ## 


    n_roster_periods_burn_in = 0
    n_roster_periods = 1 # should be 1 set of patients
    rng_seed = 1
    
    ####
    
    c = Run.configure(input_fname, 
                output_dirname,
                run_label, 
                homogeneous_flag,
                admissions_per_day,
                length_of_stay,
                n_roster_periods,
                n_roster_periods_burn_in,
                rng_seed)
    
    Run.run(c)


end


#  #####




