#TODO: this is the NurseRosterGenerator - it's not generalisable to other hcws. 


# TODO: limit shift scheduling with only one shift in between. (i.e., late-early)

module NurseRosterGenerator


import Random

using DataFrames
using CSV


import Main.Entities
import Main.Utilities


function init_BedTable_from_CSV!(bedtable::Entities.BedTable, input_fname::String )
    
    input_table = DataFrame(CSV.File(input_fname, delim = "\t"))
    
    append!(bedtable.data.room_id, input_table.room_id)
    append!(bedtable.data.bed_id, input_table.bed_id)
    append!(bedtable.data.block_id, input_table.block_id)

end

function compute_total_nursing_staff(n_bedblocks::Int64, 
                                     days_per_roster::Int64,
                                     shifts_per_day::Int64,
                                     shifts_per_nurse_per_roster::Int64 )::Int64
    n_tot = ceil(n_bedblocks * shifts_per_day  * days_per_roster / shifts_per_nurse_per_roster)
    return n_tot 
end

#TODO: add late-early restriction.
function update_nurse_eligibility!(nurse::Entities.Nurse, 
                                   max_consecutive_shifts::Int64, 
                                   max_shifts_per_roster::Int64,
                                   allow_back_to_back_shifts::Bool,
                                   allow_late_earlies::Bool)

    if (nurse.worked_last_shift && !allow_back_to_back_shifts)
        nurse.is_eligible = false
        #println("nurse not eligible: worked last shift")

    elseif (nurse.worked_shift_before_last && !allow_late_earlies)
        nurse.is_eligible = false
        #println("nurse not eligible: worked shift before last (late early)") 

    elseif (nurse.n_consecutive_shifts >= max_consecutive_shifts)
        nurse.is_eligible = false
        #println("nurse not eligible: too many consecutive days: $(nurse.n_consecutive_shifts)")

    elseif (nurse.n_shifts_this_roster >= max_shifts_per_roster)
        nurse.is_eligible = false
    else
        nurse.is_eligible = true
    end
end

#TODO: add late-early restriction.
function update_eligibility_status!(hcws::Dict{Int64, Entities.Nurse},
                                    max_consecutive_shifts::Int64, 
                                    max_shifts_per_roster::Int64,
                                    allow_back_to_back_shifts::Bool,
                                    allow_late_earlies::Bool)
    for i in keys(hcws)
        update_nurse_eligibility!(hcws[i], 
                                  max_consecutive_shifts, 
                                  max_shifts_per_roster,
                                  allow_back_to_back_shifts,
                                  allow_late_earlies)
    end
end

function nurse_is_eligible(nurse::Entities.Nurse)::Bool
    
    if nurse.is_eligible
        return true
    else
        return false
    end
end

function find_eligible_ids(hcws::Dict{Int64, Entities.Nurse})::Array{Int64, 1}
    eligible_ids = Array{Int64, 1}()
    for i in keys(hcws) # NOTE: fix random ordering. 
        if nurse_is_eligible(hcws[i])
            push!(eligible_ids, hcws[i].id)
        end
    end
    return eligible_ids
end

function find_ineligible_ids(hcws::Dict{Int64, Entities.Nurse})::Array{Int64, 1}
    ineligible_ids = Array{Int64, 1}()
    for i in keys(hcws)
        if !nurse_is_eligible(hcws[i])
            push!(ineligible_ids, hcws[i].id)
        end
    end
    return ineligible_ids
end

# add a function to detect late-early criteria 


function reset_consecutive_shifts!(hcws::Dict{Int64, Entities.Nurse}, 
                                   ids::Set{Int64})
    for i in ids
        hcws[i].n_consecutive_shifts = 0
    end
end

function add_consecutive_shifts!(hcws::Dict{Int64, Entities.Nurse}, 
                                 ids::Set{Int64})
    for i in ids
        hcws[i].n_consecutive_shifts += 1
    end
end

# implementing late-early check
# TODO: test this works 
function reset_last_shift_flag!(hcws::Dict{Int64, Entities.Nurse})
    for i in keys(hcws)
        if hcws[i].worked_last_shift
            hcws[i].worked_shift_before_last = true
            hcws[i].worked_last_shift = false
        else
            hcws[i].worked_shift_before_last = false
        end
    end
end

#TODO: use different terminology for bedblock. 
function reset_bedblock_ids!(hcws::Dict{Int64, Entities.Nurse})
    for i in keys(hcws)
        hcws[i].bedblock_id = -1
    end
end

function reset_shifts_this_roster!(hcws::Dict{Int64, Entities.Nurse})
    for i in keys(hcws)
        hcws[i].n_shifts_this_roster = 0
    end
end

function get_workload_array(hcws::Dict{Int64, Entities.Nurse},
                            ids::Array{Int64, 1})::Array{Int64, 1}
    workload = Array{Int64, 1}()
    for i in ids
        push!(workload, hcws[i].n_shifts_this_roster)
    end
    return workload
end

function get_consecutive_shift_array(hcws::Dict{Int64, Entities.Nurse},
                                     ids::Array{Int64, 1})::Array{Int64, 1}
    workload = Array{Int64, 1}()
    for i in ids
        push!(workload, hcws[i].n_consecutive_shifts)
    end
    return workload
end

function sort_ids_by_workload(ids::Array{Int64, 1}, 
                              workload::Array{Int64, 1},
                              rng_roster::Random.MersenneTwister)::Array{Int64, 1}

    #for non-deterministic sorting of equal elements. 
    #NOTE: call to random
    workload_perturbed = workload + rand(rng_roster, size(workload, 1)) .* 0.0001

    perm = sortperm(workload_perturbed)
    ids_out = ids[perm]
    return ids_out
end

function update_time_of_day(time_of_day::String)::String
    new_time_of_day = ""
    if time_of_day == "morning"
        new_time_of_day = "daytime"
    elseif time_of_day == "daytime"
        new_time_of_day = "nighttime"
    else
        new_time_of_day = "morning"
    end

    return new_time_of_day
end



#TODO: add late-early restriction. 

function generate_roster!(nurses::Dict{Int64, Entities.Nurse},
                          shifts::Dict{Int64, Entities.Shift},
                          beds::Dict{Int64, Entities.Bed},
                          bedblocks::Dict{Int64, Entities.BedBlock}, 
                          max_consecutive_shifts::Int64, 
                          shifts_per_nurse_per_roster::Int64,
                          n_roster_periods::Int64,
                          days_per_roster::Int64,
                          shifts_per_day::Int64,
                          allow_back_to_back_shifts::Bool,
                          allow_late_earlies::Bool,
                          rng_roster::Random.MersenneTwister)::Entities.Roster

    nurse_ids = collect(keys(nurses))

    n_bedblocks = size(collect(keys(bedblocks)), 1)


    #initialise roster 
    roster = Entities.Roster()
    Entities.default_init_Roster!(roster)

    #iterate through roster periods

    shift_index = 0
    day_index = 0
    time_of_day = "nighttime" #will update so first shift is morning

    for i in 1:n_roster_periods

        #println("roster period $i")
        reset_shifts_this_roster!(nurses)

        for j in 1:days_per_roster

            day_index += 1

            #println("roster day $j")


            # make a list of hcws assigned today. 
            # reset consecutive shift counter for all those not assigned. 

            IDs_assigned_today = Set{Int64}()
            IDs_not_assigned_today = Set(nurse_ids)

            for k in 1:shifts_per_day
                #println("shift $k")

                shift_index += 1

                #create a new shift: 
                time_of_day = update_time_of_day(time_of_day)
                shift_k = Entities.Shift()
                Entities.default_init_Shift!(shift_k)
                shift_k.id = shift_index
                shift_k.time_of_day = time_of_day
                shift_k.day_index = day_index
                shift_k.roster_day = j
                shifts[shift_index] = shift_k

                roster.assignments[shift_index] = Array{Int64, 1}()

                update_eligibility_status!(nurses, 
                                            max_consecutive_shifts, 
                                            shifts_per_nurse_per_roster,
                                            allow_back_to_back_shifts,
                                            allow_late_earlies)

                eligible_IDs = find_eligible_ids(nurses)
                #ineligible_IDs = find_ineligible_ids(nurses)
                #check if the number of eligible nurses is 
                #sufficient to cover all bedblocks:
                if size(eligible_IDs, 1) < n_bedblocks
                    println(" ROSTER $(i) NOT VALID: there are " * 
                             "$(size(eligible_IDs, 1)) nurses and " *
                             "$n_bedblocks bed blocks \n")
                    return 1 
                end

                # reset worked_last_shift for all nurses (this will be updated for those assigned)
                # NOTE: now includes update to worked_shift_before_last flag. 
                reset_last_shift_flag!(nurses)
                #reset bedblock IDs: 
                reset_bedblock_ids!(nurses)

                #shuffle!(eligible_IDs)

                #instead of random assignment, assign first to those with the lowest workload. 
                workload = get_workload_array(nurses, eligible_IDs)
                n_consecutive = get_consecutive_shift_array(nurses, eligible_IDs)
                sorted_eligible_IDs = sort_ids_by_workload(eligible_IDs, workload, rng_roster)
                #pop should pull the lowest one first, so reverse the order
                #println("$(workload)")
                #println("$n_consecutive")
                #println("$(eligible_IDs)")            
                sorted_eligible_IDs = reverse(sorted_eligible_IDs) # uses pop!() below, so reverse order
                assigned_ids = Array{Int64, 1}() #record which nurses were assigned 


                for bb in collect(keys(bedblocks))
                    #draw a nurse at random and assign to the next bedblock
                    nurse_id = pop!(sorted_eligible_IDs) #pops from the top, so lowest workload first. 
                    nurses[nurse_id].bedblock_id = bb
                    nurses[nurse_id].worked_last_shift = true
                    nurses[nurse_id].n_shifts_this_roster += 1
                    push!(assigned_ids, nurse_id)
                    #record that nurse was assigned today:
                    push!(IDs_assigned_today, nurse_id)
                    delete!(IDs_not_assigned_today, nurse_id)
                    #add assignment to roster
                    push!(roster.assignments[shift_index], nurse_id)


                    Utilities.push_to_dict_array!(shifts[shift_index].nurse_to_bbs, 
                                                  nurse_id, 
                                                  bb)

                    Utilities.push_to_dict_array!(shifts[shift_index].bb_to_nurses,
                                                  bb,
                                                  nurse_id)

                    #beds in bedblock are assigned:
                    for b in bedblocks[bb].bed_ids
                        beds[b].is_assigned = true
                    end

                end

            end

            # all eligible hcws who were not assigned a bedblock today get 
            # their consecutive shift counter reset. 
            reset_consecutive_shifts!(nurses, IDs_not_assigned_today)
            add_consecutive_shifts!(nurses, IDs_assigned_today)

        end

    end

    println("completed roster")
    return roster


end

function check_roster(roster::Entities.Roster, # Int -> Array(Int), links shift index to array of HCW IDs 
                      days_per_roster::Int64, # also called roster_period
                      shifts_per_day)::String

    results_string = "" 

    # to through roster and check for issues: 

    # check for double shifts: 

    # check for late-earlies:

    # check distribution of days worked for each roster period: 

    # TODO: add any other summary statistics of interest for validation 



    return results_string

end


### end of module

end