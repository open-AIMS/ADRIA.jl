using Setfield

using DataFrames
using ModelParameters
import ModelParameters: update!, Model
using Distributions

abstract type EcoModel end

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
