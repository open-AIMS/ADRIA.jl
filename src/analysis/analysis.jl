module analysis

using Statistics, DataFrames
using ADRIA
import ADRIA: ResultSet


"""
    col_normalize(data::AbstractArray)::AbstractArray
    col_normalize(data::Vector)::Vector

Normalize a matrix or vector on a per-column basis (∈ [0, 1]).
"""
function col_normalize(data::AbstractMatrix{T})::AbstractMatrix{T} where {T<:Real}
    return hcat(col_normalize.(eachcol(data))...)
end
function col_normalize(data::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    limits = extrema(data)
    scaled = let (mi, ma) = limits
        (data .- mi) ./ (ma - mi)
    end

    replace!(scaled, NaN => 0.0)
    return scaled
end

"""
    discretize_outcomes(y; S=20)

Normalize outcomes (column wise) and discretize them into \$S\$ bins.

Classify as 1 to S where:
S := 1.0 - 0.9
S-1 := 0.9 - 0.8
etc
"""
function discretize_outcomes(y; S=20)
    steps = 0.0:(1/S):1.0

    y_s_hat = col_normalize(y)
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