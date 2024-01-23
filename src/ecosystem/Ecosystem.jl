using Setfield

using DataFrames
using ModelParameters
import ModelParameters: update!, Model
using Distributions
using StatsFuns

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
        Factor(
            1;
            ptype="unordered categorical",
            dist=DiscreteUniform,
            dist_params=(1.0, Float64(size(dhw, 3))),
            name="DHW Scenario",
            description="DHW scenario member identifier.",
        ),
        Factor(
            1;
            ptype="unordered categorical",
            dist=DiscreteUniform,
            dist_params=(1.0, Float64(size(wave, 3))),
            name="Wave Scenario",
            description="Wave scenario member identifier.",
        ),
        Factor(
            1;
            ptype="unordered categorical",
            dist=DiscreteUniform,
            dist_params=(1.0, Float64(size(cyclone_mortality, 4))),
            name="Cyclone Mortality",
            description="Cyclone mortality scenario identifier.",
        ),
    )
end

"""
    truncated_standard_normal_mean(a::Float64, b::Float64)::Float64

Calculates the mean of the truncated standard normal distribution. Implementation taken
from Distribution.jl excluding unused error checks

# Arguements
- 'a' : lower bound of the truncated distribution
- 'b' : upper bound of the truncated distribution

# Returns
- Float64, mean of the truncated standard normal distribution

"""
function truncated_standard_normal_mean(a::Float64, b::Float64)::Float64
    if abs(a) > abs(b)
        return - truncated_standard_normal_mean(-b, -a)
    elseif (a == b) 
        return a
    end
    mid = (a + b) / 2 
    Δ = (b - a) * mid
    a′ = a * invsqrt2
    b′ = b * invsqrt2
    if a ≤ 0 ≤ b
        m = expm1(-Δ) * exp(-a^2 / 2) / erf(b′, a′)
    elseif 0 < a < b
        z = exp(-Δ) * erfcx(b′) - erfcx(a′)
        iszero(z) && return mid
        m = expm1(-Δ) / z
    end
    return clamp(m / sqrthalfπ, a, b)
end

"""
    truncated_normal_mean(
        normal_mean::Float64, 
        normal_stdev::Float64, 
        lower_bound::Float64, 
        upper_bound::Float64
    )::Float64   

Calculates the mean of the truncated normal distribution. Implementation is taken from
Distributions.jl excluding unused error checks

# Arguements
- 'normal_mean' : mean of the underlying (untruncated) normal distribution
- 'normal_stdev' : standard deviation of the underlying (untruncated) normal distribution
- 'lower_bound' : lower bound of the truncated normal distribution
- 'upper_bound' : upper bound of the truncated normal distributionk

# Returns
- Float64, mean of the truncated normal distribution

"""
function truncated_normal_mean(normal_mean::Float64, normal_stdev::Float64, lower_bound::Float64, upper_bound::Float64)::Float64

    alpha::Float64 = (lower_bound - normal_mean) / normal_stdev
    beta::Float64 = (upper_bound - normal_mean) / normal_stdev

    return normal_mean + truncated_standard_normal_mean(alpha, beta) * normal_stdev
end

# --------------------- need to fix numerical stability below -----------------------------

"""
    standard_normal_pdf(x::Float64)

Evaluate the probability density function (pdf) of the stanrdard normal distribution at
the given point.

"""
function standard_normal_pdf(x::Float64)::Float64
    return 1/√(2π) * ℯ^(-0.5 * x^2)
end

"""
    truncared_normal_cdf(x::Float64, normal_mean::Float64, normal_stdev::Float64,
        lower_bound::Float64, upper_bound::Float64)::Float64

Calculates the cdf of a truncated normal distribution given the mean and standard deviation
of the normal distribution, and the lower and upper bounds of the truncated distribution.

Note: this function uses the aproximation of the error function (erf) provided by
Julia SpecialFunctions.

return μ + (ϕ(c) - ϕ(b))/(Φ(b) - Φ(a)) * σ

where,

a = (lower_bound - μ)/σ
b = (upper_bound - μ)/σ
c = (x - μ)/σ

# Arguments
- x : value to evaluate cdf at
- mean : mean of underlying normal distribution
- stdev : standard deviation of underlying normal Distribution
- lower_bound : lower bound of truncated distribution
- upper_bound : upper bound of truncated distribution

# Returns
- mean of truncated normal distribution
"""
function truncated_normal_cdf(
    x::Float64,
    normal_mean::Float64,
    normal_stdev::Float64,
    lower_bound::Float64,
    upper_bound::Float64
)::Float64

    alpha::Float64 = (lower_bound - normal_mean) / normal_stdev
    beta::Float64 = (upper_bound - normal_mean) / normal_stdev
    zeta::Float64 = (x - normal_mean) / normal_stdev

    sqrt2::Float64 = √2

    Z::Float64 = 0.5 * (erf(beta / sqrt2) - erf(alpha / sqrt2))

    return 0.5 * (erf(zeta / sqrt2) - erf(alpha / sqrt2)) / Z
end
