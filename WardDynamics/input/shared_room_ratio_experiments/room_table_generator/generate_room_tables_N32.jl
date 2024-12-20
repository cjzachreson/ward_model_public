using DataFrames 
using CSV 


function add_bed!(df::DataFrame, 
                  bed_id::Int64, 
                  room_id::Int64, 
                  block_id::Int64)

    push!(df[!, :bed_id], bed_id)
    push!(df[!, :room_id], room_id)
    push!(df[!, :block_id], block_id)

end

function add_double_room!(df::DataFrame, room_id::Int64, n_beds_per_block::Int64)
    # adds two beds to the dataframe 
    for i in 1:2
        if isempty(df.bed_id)
            bed_id = 1
        else
        bed_id = last(df.bed_id) + 1
        end
        block_id = Int64(ceil(bed_id/n_beds_per_block))
        add_bed!(df, bed_id, room_id, block_id)
    end
end

function add_single_room!(df::DataFrame, room_id::Int64, n_beds_per_block::Int64)
    if isempty(df.bed_id)
        bed_id = 1
    else
        bed_id = last(df.bed_id) + 1
    end
    block_id = Int64(ceil(bed_id/n_beds_per_block))
    add_bed!(df, bed_id, room_id, block_id)
end

n_beds = 32

label = "$(n_beds)_bed_ward"

base_directory = dirname(dirname((@__FILE__)))

output_dirname = joinpath(base_directory, "room_tables", label)

if !isdir(output_dirname)
    mkpath(output_dirname)
end


n_beds_per_block = 4 

n_beds_per_shared_room = 2 

# proportion of double-occupancy rooms 
p_vec = 0:0.1:1

#p_vec = [0.1]

for p in p_vec 


    # make empty dataframe for output: 
    df_output = DataFrame(
            bed_id = Array{Int64, 1}(),
            room_id = Array{Int64, 1}(),
            block_id = Array{Int64, 1}()
            
    )


    p_str = string(round(p, digits=2))

    label_p = "p_shared_" * replace(p_str, '.' => 'p') 

    output_fname = joinpath(output_dirname, "room_table_" * label_p * ".csv")

    # number of beds in shared rooms
    n_double = Int64(floor(p * n_beds)) 

    # to ensure we get the correct number of beds. 
    if isodd(n_double)
        n_double -= 1
    end


    # number of beds in single rooms 
    n_single = Int64(n_beds - n_double)

    n_double_rooms = Int64(ceil(n_double/n_beds_per_shared_room))

    current_room_id = 0

    n_beds_created = 0



    for i in 1:n_double_rooms

        current_room_id += 1

        # generate double-occupancy rooms 
        add_double_room!(df_output, current_room_id, n_beds_per_block )

        n_beds_created += 2

    end

    for i in 1:n_single

        current_room_id += 1 
        # generate remaining single-occupancy rooms 
        add_single_room!(df_output, current_room_id, n_beds_per_block )

        n_beds_created += 1

    end

    if n_beds_created > n_beds
        println("WARNING: created too many beds")
    elseif n_beds_created < n_beds
        println("WARNING: created too few beds")
    end



    CSV.write(output_fname, df_output, delim = '\t')

end