# defines some abstract types that are used 
# in both Events.jl and Entities.jl 
# note- in downstream implementations, 
# all calls to Events.Entity_T, Events.HCW_T, Entities.HCW_T or Entities.Entity_T etc. 
# will refer to the same definition contained here. (e.g., EventEntityInterface.Entity_T)
module EventEntityInterface 


export Entity_T, HCW_T, Patient_T

abstract type Entity_T end

abstract type HCW_T <: Entity_T end
abstract type Patient_T <: Entity_T end



##### end module 
end