module analysis

using Statistics, DataFrames
using NamedDims, AxisKeys
import ADRIA: ResultSet
using ADRIA.metrics: nds


"""
    col_normalize(data::AbstractArray)::AbstractArray

Normalize a matrix on a per-column basis (∈ [0, 1]).
"""
function col_normalize(data::AbstractMatrix{T})::AbstractMatrix{T} where {T<:Real}
    d = copy(data)
    Threads.@threads for ax in axes(d, 2)
        d[:, ax] .= normalize!(d[:, ax])
    end

    return d
end
function col_normalize(data::AbstractVector{T})::AbstractVector{T} where {T<:Real}
    return normalize!(copy(data))
end

"""
    normalize(data::AbstractArray{T})::AbstractArray{T} where {T<:Real}

Normalize a matrix or vector (∈ [0, 1]).
"""
function normalize(data::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    return normalize!(copy(data))
end

"""
    normalize!(data::AbstractArray{T})::AbstractArray{T} where {T<:Real}

Normalize a matrix or vector (∈ [0, 1]) in place.
"""
function normalize!(data::AbstractArray{T})::AbstractArray{T} where {T<:Real}
    (mi, ma) = extrema(data)
    data .= (data .- mi) ./ (ma - mi)

    replace!(data, NaN => 0.0)
    return data
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
        Threads.@threads for j in size(y_s_hat, 2)
            y_disc[steps[i-1].<y_s_hat[:, j].<=steps[i], j] .= steps[i-j]
        end
    end

    return y_disc
end

include("pareto.jl")
include("sensitivity.jl")

end