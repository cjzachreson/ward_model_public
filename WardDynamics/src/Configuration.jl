module Configuration 

import Random 
import Distributions

#export

abstract type Config_T end

mutable struct Config <: Config_T
    rng_seed::Int64
    dt::Int64 #timestep in minutes
    output_dirname::String 
    output_fname_hcw_traj::String
    output_fname_patient_traj::String

    room_table_fname::String

    n_roster_periods::Int64

    n_roster_periods_burn_in::Int64

    days_per_roster::Int64
    shifts_per_day::Int64
    shifts_per_nurse_per_roster::Int64

    staff_deficit::Int64

    max_consecutive_shifts::Int64

    allow_back_to_back_shifts::Bool 
    allow_late_earlies::Bool

    admissions_per_day::Float64
    
    homogeneous::Bool

    length_of_stay::Distributions.Sampleable #how many shifts do patients stay for 

    rng_LoS::Random.MersenneTwister

    rng_roster::Random.MersenneTwister

    rng_trajectories::Random.MersenneTwister

    Config() = new()
    
end

function default_init_Config!(c::Config)

    c.rng_seed = 1

    c.dt = 5 #timestep in minutes


    c.output_dirname = joinpath(dirname(pwd()),
                                "output",
                                "test")

    c.output_fname_hcw_traj = joinpath(c.output_dirname,"trajectory_hcw.csv")
    c.output_fname_patient_traj = joinpath(c.output_dirname, "trajectory_patient.csv")

    c.room_table_fname = joinpath(dirname(pwd()),
                            "input",
                            "room_table_example.csv")


    c.n_roster_periods = 10

    c.n_roster_periods_burn_in = 5

    c.days_per_roster = 14
    c.shifts_per_day = 3
    c.shifts_per_nurse_per_roster = 8

    c.staff_deficit = 0

    c.max_consecutive_shifts = 6

    c.allow_back_to_back_shifts = false

    c.allow_late_earlies = false 

    c.admissions_per_day = 6

    c.homogeneous = false 

    c.length_of_stay = Distributions.Binomial(18, 0.5) #default mean is 3 days or 9 shifts

    c.rng_LoS = Random.MersenneTwister(1)

    c.rng_roster = Random.MersenneTwister(1)

    c.rng_trajectories = Random.MersenneTwister(1)
    

end

function default_homogeneous_Config!(c::Config)
    c.rng_seed = 1
    c.dt = 5
    c.output_dirname = joinpath(dirname(pwd()),
                                "output",
                                "test")

    c.output_fname_hcw_traj = joinpath(c.output_dirname,"trajectory_hom_N100_hcw.csv")
    c.output_fname_patient_traj = joinpath(c.output_dirname, "trajectory_hom_N100_patient.csv")

    c.room_table_fname = joinpath(dirname(pwd()),
                                  "input",
                                  "room_table_hom_N100.csv")

    c.n_roster_periods = 10

    c.n_roster_periods_burn_in = 5

    c.days_per_roster = 14
    c.shifts_per_day = 3
    c.shifts_per_nurse_per_roster = c.days_per_roster * c.shifts_per_day

    c.staff_deficit = 0

    c.max_consecutive_shifts = c.n_roster_periods * c.days_per_roster * c.shifts_per_day

    c.allow_back_to_back_shifts = true

    c.allow_late_earlies = true
    #late-early restriction, make sure to turn it off. 

    c.admissions_per_day = 300

    c.homogeneous = true 

    # all patients stay the whole time (fixed population)
    c.length_of_stay = Distributions.Dirac(c.n_roster_periods * c.days_per_roster * c.shifts_per_day)

    c.rng_LoS = Random.MersenneTwister(1)

    c.rng_roster = Random.MersenneTwister(1)

    c.rng_trajectories  = Random.MersenneTwister(1)

end

function set_input_fname!(c::Config, input_fname::String)
    c.room_table_fname = input_fname
end

function set_output_fnames!(c::Config, output_dirname::String, label::String)
    c.output_dirname = output_dirname
    c.output_fname_hcw_traj = joinpath(c.output_dirname, "trajectory_hcw_" * label * ".csv")
    c.output_fname_patient_traj = joinpath(c.output_dirname, "trajectory_patient_" * label * ".csv")
end

function set_n_roster_periods!(c::Config, n::Int64)
    c.n_roster_periods = n
end

function set_burn_in!(c::Config, n::Int64)
    c.n_roster_periods_burn_in = n
end

function set_admissions_per_day!(c::Config, n::Int64)
    c.admissions_per_day = n
end

#TODO: decouple seed control for roseter from seed control for patient bed assignments 
function set_rng_seed!(c::Config, n::Int64)
    c.rng_seed = n
    c.rng_LoS = Random.MersenneTwister(n)
    c.rng_roster = Random.MersenneTwister(n)
    c.rng_trajectories = Random.MersenneTwister(n)
end

function set_length_of_stay!(c::Config, LoS::Distributions.ContinuousUnivariateDistribution)
    c.length_of_stay = LoS
end

# default ward config flag 

# ward data input file

# n roster periods 

# days per roster period 

# shifts per day 

# shifts per nurse per roster 

# max number of consecutive shifts (i.e., shifts on consecutive days)

# flag for allowing back-to-back shifts (doubles etc.)




end