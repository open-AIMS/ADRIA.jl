module analysis

using Statistics, DataFrames
using ADRIA
import ADRIA: ResultSet


"""
    normalize(data::AbstractArray)::AbstractArray
    normalize(data::Vector)::Vector

Normalize a matrix (∈ [0, 1]) on a per-column basis.
"""
function normalize(data::AbstractMatrix{T})::AbstractMatrix{T} where {T<:Real}
    return hcat(normalize.(eachcol(data))...)
end
function normalize(data::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    limits = extrema(data)
    scaled = let (mi, ma) = limits
        (data .- mi) ./ (ma - mi)
    end

    replace!(scaled, NaN => 0.0)
    return scaled
end

"""
    discretize_outcomes(y; S=20)

Normalize outcomes and discretize them into \$S\$ bins.

Classify as 1 to S where:
S := 1.0 - 0.9
S-1 := 0.9 - 0.8
etc
"""
function discretize_outcomes(y; S=20)
    steps = 0.0:(1/S):1.0

    y_s_hat = normalize(y)
    y_disc = zeros(size(y)...)
    for i in axes(steps, 1)[2:end]
        for j in size(y_s_hat, 2)
            y_disc[steps[i-1].<y_s_hat[:, j].<=steps[i], j] .= steps[i-j]
        end
    end

    return y_disc
end


include("pareto.jl")
include("sensitivity.jl")


end