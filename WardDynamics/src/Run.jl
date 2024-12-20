module Run

import Main.Configuration
import Main.GenerateTrajectories
import Main.InitialiseWard 
import Main.WardModelOutput

import Random 

import Distributions 

# TODO: make sure important control parameters are specified
# in a place that's not confusing - the parameter settings below are 
# acting like defaults. 
function configure(input_fname::String, 
                   output_dirname::String, 
                   output_label::String,
                   h_flag::Bool,
                   admissions_per_day::Int64,
                   length_of_stay::Distributions.ContinuousUnivariateDistribution,
                   n_roster_periods::Int64,
                   n_roster_periods_burn_in::Int64,
                   rng_seed::Int64)::Configuration.Config

    c = Configuration.Config()

    homogeneous = h_flag  

    if homogeneous
        Configuration.default_homogeneous_Config!(c)
    else
        Configuration.default_init_Config!(c)
    end 

    Configuration.set_rng_seed!(c, rng_seed)
    Configuration.set_input_fname!(c, input_fname)
    Configuration.set_output_fnames!(c, output_dirname, output_label)
    Configuration.set_n_roster_periods!(c, n_roster_periods)#NOTE: ensure this is updated.
    Configuration.set_burn_in!(c, n_roster_periods_burn_in) 
    Configuration.set_admissions_per_day!(c, admissions_per_day)
    Configuration.set_length_of_stay!(c, length_of_stay)

    
    return c

end

function run(c::Configuration.Config)

    #TODO: try removing this and see if anything changes. 
    #Random.seed!(c.rng_seed)

    ward = InitialiseWard.Ward()
    InitialiseWard.default_init_Ward!(ward)
    InitialiseWard.configure_Ward!(ward, c)

    # discrete-time dynamics should run for each shift. 

    #output structures: 
    nurse_trajectories = WardModelOutput.Agent_Trajectories()
    patient_trajectories = WardModelOutput.Agent_Trajectories()

    GenerateTrajectories.generate_trajectories!(ward, 
                                                c, 
                                                nurse_trajectories, 
                                                patient_trajectories)


    GenerateTrajectories.insert_patient_wait_times!(patient_trajectories, c)

    #record trajectory output: 
    WardModelOutput.trajectories_to_CSV(nurse_trajectories, c.output_fname_hcw_traj)
    WardModelOutput.trajectories_to_CSV(patient_trajectories, c.output_fname_patient_traj)

end

#end module
end