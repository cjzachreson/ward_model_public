# defines some abstract types that are used 
# in both TransmissionDynamics.jl and AirflowDynamics.jl 
# note- in downstream implementations, 


# all calls to TransmissionDynamics.Config_T 
# will refer to the same definition contained here. (e.g., AirflowTransmissionInterface.Config_T)
module AirflowTransmissionInterface


export Config_T

abstract type Config_T end


##### end module 
end