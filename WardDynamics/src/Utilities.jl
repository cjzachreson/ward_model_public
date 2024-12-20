module Utilities

function push_to_dict_array!(dict, key, newval)

    if haskey(dict, key)
        push!(dict[key], newval)
    else
        dict[key] = [newval]
    end
end

function push_vector_to_set!(set, vector)
    for i in vector
        push!(set, i)
    end
end


### end module
end