#using Random 
using Plots 
using DataFrames
using CSV
using Statistics
using Distributions


base_directory = dirname(dirname((@__FILE__)))

println(base_directory)

#TODO: go through and check all rng calls 
#Random.seed!(11)

include(joinpath(base_directory, "src","Configuration.jl"))
include(joinpath(base_directory, "src","Utilities.jl"))
include(joinpath(base_directory, "src","EventEntityInterface.jl"))
include(joinpath(base_directory, "src","Events.jl"))
include(joinpath(base_directory, "src","Entities.jl"))
include(joinpath(base_directory, "src","WardModelOutput.jl"))
include(joinpath(base_directory, "src","NurseRosterGenerator.jl"))
include(joinpath(base_directory, "src","InitialiseWard.jl"))
include(joinpath(base_directory, "src","PatientFunctions.jl"))
include(joinpath(base_directory, "src","NurseFunctions.jl"))
include(joinpath(base_directory, "src","UpdateWard.jl"))
include(joinpath(base_directory, "src","GenerateTrajectories.jl"))
include(joinpath(base_directory, "src","Run.jl"))

import Main.run

