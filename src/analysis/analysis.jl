module analysis

using ADRIA: ZeroDataCube, ResultSet, n_locations
using ADRIA.metrics: nds

using JuliennedArrays: Slices
using Statistics, DataFrames
using YAXArrays

"""
    col_normalize(data::AbstractArray)::AbstractArray

Normalize a matrix on a per-column basis (∈ [0, 1]).
"""
function col_normalize(
    data::AbstractMatrix{T}
)::AbstractMatrix{T} where {T<:Union{Missing,Real}}
    local d
    try
        d = copy(data)
    catch
        d = copy(data.data)
    end

    Threads.@threads for ax in axes(d, 2)
        @inbounds d[:, ax] .= normalize!(d[:, ax])
    end

    return d
end
function col_normalize(
    data::AbstractVector{T}
)::AbstractVector{T} where {T<:Union{Missing,Real}}
    return normalize(data)
end

"""
    normalize(data::AbstractArray{T})::AbstractArray{T} where {T<:Real}

Normalize a matrix or vector (∈ [0, 1]).
"""
function normalize(data::AbstractArray{T})::AbstractArray{T} where {T<:Union{Missing,Real}}
    d = copy(data)
    normalize!(d)

    return d
end

"""
    normalize!(data::AbstractArray{T})::AbstractArray{T} where {T<:Real}

Normalize a matrix or vector (∈ [0, 1]) in place.
"""
function normalize!(data::AbstractArray{T})::AbstractArray{T} where {T<:Union{Missing,Real}}
    if count(!ismissing, data) == 0
        data .= zeros(size(data)...)
        return data
    end

    (mi, ma) = extrema(skipmissing(data))
    if mi == ma
        pos = findall(!ismissing, data)
        data[pos] .= zeros(typeof(data[pos][1]), size(data[pos])...)
        return data
    end

    data .= (data .- mi) ./ (ma - mi)
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
    steps = 0.0:(1 / S):1.0

    y_s_hat = col_normalize(y)
    y_disc = zeros(size(y)...)
    for i in axes(steps, 1)[2:end]
        Threads.@threads for j in size(y_s_hat, 2)
            y_disc[steps[i - 1] .< y_s_hat[:, j] .<= steps[i], j] .= steps[i - j]
        end
    end

    return y_disc
end

"""
    series_confint(data::AbstractMatrix; agg_dim::Symbol=:scenarios)::Matrix{Float64}

Computes confidence interval for series of data.

# Arguments
- `data` : Matrix with series of data
- `agg_dim` : Dimension used to aggregate data if a YAXArray is passed

# Returns
Confidence interval (lower bound, median and higher bound) for each series step
"""
function series_confint(data::AbstractMatrix; agg_dim::Symbol=:scenarios)::Matrix{Float64}
    slice_dim = data isa YAXArray ? YAXArrays.findAxis(agg_dim, data) : 2
    return quantile.(Slices(data, slice_dim), [0.025 0.5 0.975])
end

include("pareto.jl")
include("screening.jl")
include("intervention.jl")
include("clustering.jl")
include("rule_extraction.jl")
include("scenario.jl")
include("feature_set.jl")
include("data_envelopment.jl")

end
