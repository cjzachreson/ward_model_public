module InitialiseWard

import Main.Configuration 

import Main.Entities

import Main.NurseRosterGenerator



abstract type Ward_T end

mutable struct Ward <: Ward_T
    nurses::Dict{Int64, Entities.Nurse}
    roster::Entities.Roster
    shifts::Dict{Int64, Entities.Shift}
    bedblocks::Dict{Int64, Entities.BedBlock}
    bedblock_ids::Array{Int64, 1}
    bedtable::Entities.BedTable 
    beds::Dict{Int64, Entities.Bed}
    #TODO: make sure the term bedblock is replaced with something unambiguous. 
    n_bedblocks::Int64 
    staff_deficit::Int64 
    admitted_patients::Entities.PatientList_T

    shift_indices::Array{Int64, 1}
    n_shifts::Int64
    n_days::Int64 

    admissions_per_shift::Float64 

    shift_index::Int64 

    Ward() = new()
end

function default_init_Ward!(w::Ward)
    w.nurses = Dict{Int64, Entities.Nurse}()
    w.roster = Entities.Roster()
    w.shifts = Dict{Int64, Entities.Shift}()
    w.bedblocks = Dict{Int64, Entities.BedBlock}()
    w.bedblock_ids = Array{Int64, 1}()
    w.bedtable = Entities.BedTable()
    w.beds = Dict{Int64, Entities.Bed}()
    w.n_bedblocks = 0 
    w.admitted_patients = Entities.PatientList_T()

    w.shift_indices = Array{Int64, 1}()
    w.n_shifts = 0
    w.n_days = 0 

    w.admissions_per_shift = 0.0

    w.shift_index = 0


end

function configure_Ward!(w::Ward, c::Configuration.Config)

    NurseRosterGenerator.init_BedTable_from_CSV!(w.bedtable, c.room_table_fname)
    w.n_bedblocks = length(unique(w.bedtable.data.block_id))

    #determines the required number of nursing staff (rounding up)
    # TODO: does this allow late-early restrictions? 
    n_nurses_tot = NurseRosterGenerator.compute_total_nursing_staff(w.n_bedblocks,
                                                            c.days_per_roster,
                                                            c.shifts_per_day,
                                                            c.shifts_per_nurse_per_roster)
    
    n_nurses_tot -= c.staff_deficit 
    generate_nurses!(w.nurses, n_nurses_tot)


    generate_beds!(w)

    generate_bedblocks!(w)

    n_roster_periods_tot = c.n_roster_periods + c.n_roster_periods_burn_in

    w.roster = NurseRosterGenerator.generate_roster!(w.nurses, 
                                                    w.shifts,
                                                    w.beds,
                                                    w.bedblocks,
                                                    c.max_consecutive_shifts, 
                                                    c.shifts_per_nurse_per_roster, 
                                                    n_roster_periods_tot,
                                                    c.days_per_roster,
                                                    c.shifts_per_day,
                                                    c.allow_back_to_back_shifts,
                                                    c.allow_late_earlies,
                                                    c.rng_roster)


    w.shift_indices = sort(collect(keys(w.shifts)))
    w.n_shifts = size(w.shift_indices, 1)
    w.n_days = Int64(ceil(w.n_shifts / c.shifts_per_day)) 

    w.admissions_per_shift = c.admissions_per_day / c.shifts_per_day

    w.shift_index = 0

end

#initialise beds: 
function generate_beds!(w::Ward)
    n_beds = size(w.bedtable.data, 1)
    for i in 1:n_beds
        index = w.bedtable.data.bed_id[i] #should just be i 
        w.beds[index] = Entities.Bed()
        Entities.default_init_Bed!(w.beds[index])
        w.beds[index].id = index
        w.beds[index].room_id = w.bedtable.data.room_id[i]
        w.beds[index].block_id = w.bedtable.data.block_id[i]
        # not touching other params, because they will be updated dynamically. 
    end
end

# initialise bed blocks 
function generate_bedblocks!(w::Ward)
    w.bedblock_ids = unique(w.bedtable.data.block_id)
    for i in w.bedblock_ids
        w.bedblocks[i] = Entities.BedBlock()
        Entities.default_init_BedBlock!(w.bedblocks[i])
        w.bedblocks[i].id = i
        w.bedblocks[i].bed_ids =  w.bedtable.data[w.bedtable.data.block_id .== i, :].bed_id
    end
end

# initialise HCWs 
function generate_nurses!(nurses::Dict{Int64, Entities.Nurse}, n_nurses_tot::Int64)
    for i in 1:n_nurses_tot
        nurses[i] = Entities.Nurse()
        Entities.default_init_Nurse!(nurses[i])
        nurses[i].id = i
    end
end

### end module
end