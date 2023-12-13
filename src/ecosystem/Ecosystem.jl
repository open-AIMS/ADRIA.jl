using Setfield

using DataFrames
using ModelParameters
import ModelParameters: update!, Model
using Distributions

abstract type EcoModel end

const DISCRETE_FACTOR_TYPES = ["integer", "categorical"]

"""
    _check_discrete(p_type::String)::Bool
    _check_discrete(dom, fieldname::Symbol)::Bool

Check ptype for discrete variable types. Returns true if discrete, false otherwise.

# Arguments
- `ptype` : String representing variable type.
"""
function _check_discrete(p_type::String)::Bool
    return p_type ∈ DISCRETE_FACTOR_TYPES
end
function _check_discrete(dom, fieldname::Symbol)::Bool
    model::Model = dom.model
    param_filter::BitVector = collect(model[:fieldname]) .== fieldname
    ptype::String = model[:ptype][param_filter][1]
    return _check_discrete(ptype)
end

"""Set a model parameter value directly."""
function set(p::Param, val::Union{Int64,Float64})
    if hasproperty(p, :ptype)
        if _check_discrete.(p.ptype) && !isinteger(val)
            val = map_to_discrete(val, p.bounds[2])
        end
    end

    return val
end

"""
    map_to_discrete(v::Union{Int64,Float64}, u::Union{Int64,Float64})::Int64

For integer/categorical parameters, take floor of `v`, capping to `u - 1`
"""
function map_to_discrete(v::Union{Int64,Float64}, u::Union{Int64,Float64})::Int64
    return Int64(min(floor(v), u - 1))
end

"""
    map_to_discrete!(df::Union{DataFrame,SubDataFrame}, u::Union{AbstractVector{Union{Int64,Float64}},Tuple})::Nothing

Update a dataframe of parameters.
Length of `u` (the upper bounds) is expected to match number of columns in `df`.
"""
function map_to_discrete!(
    df::Union{DataFrame,SubDataFrame},
    u::Union{AbstractVector{<:Union{Int64,Float64}},Tuple},
)::Nothing
    for (idx, b) in enumerate(u)
        df[!, idx] .= map_to_discrete.(df[!, idx], b)
    end
    return nothing
end

struct EnvironmentalLayer{P<:Param} <: EcoModel
    dhw_scenario::P
    wave_scenario::P
    cyclone_mortality_scenario::P
end

function EnvironmentalLayer(
    dhw::AbstractArray{T}, wave::AbstractArray{T2}, cyclone_mortality::AbstractArray{<:Real}
)::EnvironmentalLayer where {
    T<:Union{Missing,Float32,Float64},T2<:Union{Missing,Float32,Float64}
}
    return EnvironmentalLayer(
        Param(
            1;
            bounds=(1.0, Float64(size(dhw, 3)) + 1.0),
            default_bounds=(1.0, Float64(size(dhw, 3)) + 1.0),
            ptype="integer",
            dists="unif",
            criteria_keywords=(""),
            name="DHW Scenario",
            description="DHW scenario member identifier.",
        ),
        Param(
            1;
            bounds=(1.0, Float64(size(wave, 3)) + 1.0),
            default_bounds=(1.0, Float64(size(wave, 3)) + 1.0),
            ptype="integer",
            dists="unif",
            criteria_keywords=(""),
            name="Wave Scenario",
            description="Wave scenario member identifier.",
        ),
        Param(
            1;
            bounds=(1.0, Float64(size(cyclone_mortality, 4)) + 1.0),
            default_bounds=(1.0, Float64(size(cyclone_mortality, 4)) + 1.0),
            ptype="integer",
            dists="unif",
            criteria_keywords=(""),
            name="Cyclone Mortality",
            description="Cyclone mortality scenario identifier.",
        ),
    )
end
