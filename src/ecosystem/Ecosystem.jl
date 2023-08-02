using Setfield

using DataFrames
using ModelParameters
import ModelParameters: update!, Model
using Distributions

abstract type EcoModel end


"""Set a model parameter value directly."""
function set(p::Param, val::Union{Int64,Float64})
    if hasproperty(p, :ptype) && p.ptype == "integer" && !isinteger(val)
        val = map_to_discrete(val, p.bounds[2])
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
    map_to_discrete!(df::DataFrame, u::AbstractArray)::Nothing

Update a dataframe of parameters.
Length of `u` is expected to match number of columns in `df`.
"""
function map_to_discrete!(df::DataFrame, u::Union{AbstractVector,Tuple})::Nothing
    for (idx, b) in enumerate(u)
        df[!, idx] .= map_to_discrete.(df[!, idx], b)
    end
end

struct EnvironmentalLayer{P} <: EcoModel
    dhw_scenario::P
    wave_scenario::P
end

function EnvironmentalLayer(dhw::AbstractArray, wave::AbstractArray)
    return EnvironmentalLayer(
        Param(1, bounds=(1.0, Float64(size(dhw, 3)) + 1.0), ptype="integer", dists="unif", name="DHW Scenario", description="DHW scenario member identifier."),
        Param(1, bounds=(1.0, Float64(size(wave, 3)) + 1.0), ptype="integer", dists="unif", name="Wave Scenario", description="Wave scenario member identifier.")
    )
end