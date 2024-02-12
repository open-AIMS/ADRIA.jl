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
    rational_erf(x::Float64)::Float64

Rational approximation of the error function only using elementary functions [1].
Maximum error of 1.5 × 10^{-7}.

# References
1. Abramowitz, Milton; Stegun, Irene Ann, eds. (1983) [June 1964]. "Chapter 7". 
   Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables. 
   Applied Mathematics Series. Vol. 55 (Ninth reprint with additional corrections of 
   tenth original printing with corrections (December 1972); first ed.). Washington D.C.; 
   New York: United States Department of Commerce, National Bureau of Standards; 
   Dover Publications. p. 297
"""
function rational_erf(x::Float64)::Float64
    coef::Float64 = 1.0
    if (x < 0)
        x *= -1
        coef = -1
    end

    # Use multiplication to avoid slow power function (x^n)
    x2::Float64 = x * x
    x3::Float64 = x2 * x
    x4::Float64 = x3 * x
    x5::Float64 = x4 * x
    x6::Float64 = x5 * x

    a1::Float64 = 0.0705230784
    a2::Float64 = 0.0422820123
    a3::Float64 = 0.0092705272
    a4::Float64 = 0.0001520143
    a5::Float64 = 0.0002765672
    a6::Float64 = 0.0000430638

    denom = 1.0 + a1 * x + a2 * x2 + a3 * x3 + a4 * x4 + a5 * x5 + a6 * x6

    denom = denom * denom # power 2
    denom = denom * denom # power 4
    denom = denom * denom # power 8
    denom = denom * denom # power 16

    return coef * (1 - 1.0 / denom)
end

"""
    rational_erfcx(x::Float64)::Float64

Approximation of erfcx using a rational approximation of the error function. 

erfcx(x) = e^{x^2} ⋅ (1 - erf(x))
"""
function rational_erfcx(x::Float64)::Float64
    return exp(x * x) * (1 - erf(x))
end

"""
    truncated_standard_normal_mean(lb::Float64, ub::Float64)::Float64

Calculates the mean of the truncated standard normal distribution. Implementation taken
from Distributions.jl [1] excluding unused error checks.

Note: This implementation attempts to avoid NaNs where possible. In cases where `lb` > `ub`,
      `ub` is returned.

# Arguments
- `lb` : lower bound of the truncated distribution
- `ub` : upper bound of the truncated distribution

# Returns
The mean of the truncated standard normal distribution or the upper bound if `lb` > `ub`.

# References
1. Documentation of (Distributions.jl)[https://juliastats.org/Distributions.jl/stable/],
   Specific code taken from (truncated normal directory)[https://github.com/JuliaStats/Distributions.jl/blob/c1705a3015d438f7e841e82ef5148224813831e8/src/truncated/normal.jl#L24-L46] on 24/01/2024
"""
function truncated_standard_normal_mean(lb::Float64, ub::Float64)::Float64
    if abs(lb) > abs(ub)
        return -truncated_standard_normal_mean(-ub, -lb)
    elseif (lb == ub)
        return lb
    end

    mid = (lb + ub) / 2
    Δ = (ub - lb) * mid
    lb′ = lb * StatsFuns.invsqrt2
    ub′ = ub * StatsFuns.invsqrt2

    m = ub
    if lb ≤ 0 ≤ ub
        m = expm1(-Δ) * exp(-lb^2 / 2) / (rational_erf(lb′) - rational_erf(ub′))
    elseif 0 < lb < ub
        z = exp(-Δ) * rational_erfcx(ub′) - rational_erfcx(lb′)
        iszero(z) && return mid
        m = expm1(-Δ) / z
    end

    return clamp(m / StatsFuns.sqrthalfπ, lb, ub)
end

"""
    truncated_normal_mean(
        normal_mean::Float64,
        normal_stdev::Float64,
        lower_bound::Float64,
        upper_bound::Float64
    )::Float64

Calculates the mean of the truncated normal distribution.

# Arguments
- `normal_mean` : mean of the underlying (untruncated) normal distribution
- `normal_stdev` : standard deviation of the underlying (untruncated) normal distribution
- `lower_bound` : lower bound of the truncated normal distribution
- `upper_bound` : upper bound of the truncated normal distributionk

# Returns
The mean of the truncated normal distribution
"""
function truncated_normal_mean(normal_mean::Float64, normal_stdev::Float64, lower_bound::Float64, upper_bound::Float64)::Float64
    alpha::Float64 = (lower_bound - normal_mean) / normal_stdev
    beta::Float64 = (upper_bound - normal_mean) / normal_stdev

    return normal_mean + truncated_standard_normal_mean(alpha, beta) * normal_stdev
end

"""
    truncated_normal_cdf(x::Float64, normal_mean::Float64, normal_stdev::Float64,
        lower_bound::Float64, upper_bound::Float64)::Float64

Calculates the cdf of a truncated normal distribution given the mean and standard deviation
of the normal distribution, and the lower and upper bounds of the truncated distribution.

# Arguments
- `x` : value to evaluate cdf at
- `mean` : mean of underlying normal distribution
- `stdev` : standard deviation of underlying normal Distribution
- `lower_bound` : lower bound of the truncated distribution
- `upper_bound` : upper bound of the truncated distribution

# Returns
The cdf of truncated normal distribution evaluated at `x`
"""
function truncated_normal_cdf(
    x::Float64,
    normal_mean::Float64,
    normal_stdev::Float64,
    lower_bound::Float64,
    upper_bound::Float64
)::Float64
    if x <= lower_bound
        return 0.0
    elseif x >= upper_bound
        return 1.0
    end

    alpha::Float64 = (lower_bound - normal_mean) / normal_stdev
    beta::Float64 = (upper_bound - normal_mean) / normal_stdev
    zeta::Float64 = (x - normal_mean) / normal_stdev

    # Large errors occurs when bounds deviate from the mean significantly and
    # are close together relative to the standard deviation.
    threshold = 3
    if (alpha > threshold && beta > threshold) || (alpha < -threshold && beta < -threshold)
        @debug "Possible loss of accuracy: the given truncated normal distribution bounds \
            are more than 5 standard deviations from the normal mean. \
            \nLower and upper bounds of the truncated normal distribution are \
            $(alpha) and $(beta) standard deviations from the normal mean respectively. \
            Falling back to more accurate calculation."
        
        return erf(alpha * StatsFuns.invsqrt2, zeta * StatsFuns.invsqrt2) /
               erf(alpha * StatsFuns.invsqrt2, beta * StatsFuns.invsqrt2)
    end

    # Store error function of alpha to avoid duplicate calculations
    erf_alpha = rational_erf(alpha * StatsFuns.invsqrt2)


    return (rational_erf(zeta * StatsFuns.invsqrt2) - erf_alpha) /
           (rational_erf(beta * StatsFuns.invsqrt2) - erf_alpha)
end
