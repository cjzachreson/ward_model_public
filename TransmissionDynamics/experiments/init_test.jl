using Random 
#using Plots 
using DataFrames
using CSV
using Statistics


base_directory = dirname(dirname((@__FILE__)))

println(base_directory)

#TODO: go through and check all rng calls 
Random.seed!(11)

#TODO: work on path navigation 
include(joinpath(base_directory, "src", "Utilities.jl"))
include(joinpath(base_directory, "src","EntitiesTD.jl"))
include(joinpath(base_directory, "src","OutputStructures.jl"))
include(joinpath(base_directory, "src", "AirflowTransmissionInterface.jl"))
include(joinpath(base_directory, "src", "AirflowDynamics.jl"))
include(joinpath(base_directory, "src", "TransmissionDynamics.jl"))


#cd(joinpath(base_directory, "experiments"))

import Main.EntitiesTD 
import Main.TransmissionDynamics
import Main.OutputStructures