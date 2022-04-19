using Setfield
using ModelParameters
import ModelParameters: update!, Model


abstract type EcoModel end


function set!(p::Param, val)
    if hasproperty(p, :ptype) && p.ptype == "integer" && !isinteger(val)
        # For integer/categorical parameters, take floor of v + 1, capping to the upper bound
        val = floor(val + 1)
        val = min(val, p.bounds[2])
    end

    # @set! p.val = val
    return val
end


function update!(m::Model, vals::Union{Vector,Tuple,Array})::Nothing
    m[:val] = map((x) -> set!(x...), zip(params(m), vals))

    return
end


Base.@kwdef struct Intervention{N,P} <: EcoModel
    # Intervention Parameters
    guided::N = Param(0, ptype="integer", bounds=(0, 5)) # Guided, choice of MCDA approach
    seed_TA::N = Param(0, ptype="integer", bounds=(0, 500000)) # Seed1, integer, number of Enhanced TA to seed
    seed_TC::N = Param(0, ptype="integer", bounds=(0, 500000)) # Seed2, integer, number of Enhanced TC to seed
    fogging::P = Param(0.2, ptype="real", bounds=(0.0, 0.3)) # fogging, float, assumed percent reduction in bleaching mortality
    SRM::P = Param(0.0, ptype="real", bounds=(0.0, 12.0)) # SRM, float, reduction in DHWs due to shading
    a_adapt::P = Param(0.0, ptype="real", bounds=(0.0, 12.0)) # Aadpt, float, float, increased adaptation rate
    n_adapt::P = Param(0.025, ptype="real", bounds=(0.0, 0.1)) # Natad, float, natural adaptation rate
    seed_years::N = Param(10, ptype="integer", bounds=(5, 16)) # Seedyrs, integer, years into simulation during which seeding is considered
    shade_years::N = Param(10, ptype="integer", bounds=(5, 74)) # Shadeyrs, integer, years into simulation during which shading is considered
    seed_freq::N = Param(5, ptype="integer", bounds=(0, 6)) # Seedfreq, integer, yearly intervals to adjust seeding site selection (0 is set and forget)
    shade_freq::N = Param(1, ptype="integer", bounds=(0, 6)) # Shadefreq, integer, yearly intervals to adjust shading (fogging) site selection (0 is set and forget)
    seed_year_start::N = Param(2, ptype="integer", bounds=(2, 26)) # Seedyr_start, integer, seed intervention start offset from simulation start
    shade_year_start::N = Param(2, ptype="integer", bounds=(2, 26)) # Shadeyr_start, integer, shade intervention start offset from simulation start
end


Base.@kwdef struct Criteria{P} <: EcoModel
    wave_stress::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    heat_stress::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    shade_connectivity::P = Param(0.0, ptype="real", bounds=(0.0, 1.0))
    seed_connectivity::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    coral_cover_high::P = Param(0.0, ptype="real", bounds=(0.0, 1.0))
    coral_cover_low::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    seed_priority::P = Param(1.0, ptype="real", bounds=(0.0, 1.0))
    shade_priority::P = Param(0.0, ptype="real", bounds=(0.0, 1.0))
    deployed_coral_risk_tol::P = Param(5.0, ptype="real", bounds=(0.0, 1.0))
    depth_min::P = Param(5.0, ptype="real", bounds=(3.0, 5.0))     # minimum depth
    depth_offset::P = Param(5.0, ptype="real", bounds=(5.0, 6.0))  # offset from minimum depth to indicate maximum depth**
end
# **This is simply to avoid parameterization/implementation
#   that requires one parameter to be greater than another.


"""
    GenerateCoralStruct(field_defs::Dict)::Nothing

Helper function to dynamically create structs

https://stackoverflow.com/a/27084705/2694952
https://stackoverflow.com/questions/27083816/is-it-possible-to-create-types-in-julia-at-runtime


Example
-------
julia> GenerateCoralStruct(Dict("a"=>200))

julia> y = Corals()
Corals{Int64}(200)

julia> y.a
    200
"""
function GenerateCoralStruct(field_defs::Dict)::Nothing
    s::String = "Base.@kwdef struct Corals{P} <: EcoModel\n"

    for (f, v) in field_defs
        s = s * "$(f)::P = $(v)\n"
    end

    eval(Meta.parse(s * "end"))

    return
end


function CoralParams(bounds=(0.9, 1.1))
    _, base_coral_params, x = coral_spec()

    coral_ids = x.coral_id

    struct_fields = Dict()
    for c_id in coral_ids
        for p in base_coral_params
            f_name = c_id * "_" * p
            f_val = x[x.coral_id.==c_id, p][1]
            struct_fields[f_name] = Param(f_val, ptype="real", bounds=(f_val * bounds[1], f_val * bounds[2]))
        end
    end

    GenerateCoralStruct(struct_fields)

    return
end


function to_spec(m::Model)::DataFrame
    _, pnames, spec = coral_spec()

    val_df = DataFrame(m)
    res = copy(spec)
    fnames = m[:fieldname]
    for p in pnames
        # tmp = String
        target = [occursin(p, String(fn)) for fn in fnames]
        target_names = map(String, fnames[target])
        for tn in target_names
            c_id = rsplit(tn, "_$p", keepempty=false)
            res[res.coral_id.==c_id, p] = val_df[val_df.fieldname.==Symbol(tn), :val]
        end
    end

    return res
end


struct Domain
    TP_data     # site connectivity data
    site_ranks  # site rank
    strongpred  # strongest predecessor
    site_data   # table of site data (depth, carrying capacity, etc)
    init_coral_cover  # initial coral cover dataset
    connectivity_site_ids  # Site IDs as specified by the connectivity dataset (indicates order of `TP_data`)
    removed_sites  # indices of sites that were removed. Used to align site_data, DHW, connectivity, etc.
    dhw_scens  # DHW scenarios
    wave_scens # wave scenarios
end
