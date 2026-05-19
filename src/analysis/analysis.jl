module analysis

using ADRIA: ZeroDataCube, ResultSet, n_locations
using ADRIA.metrics: nds

using Clustering, Distances
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

# Clustering

function _complexity(x::AbstractMatrix{<:Real})::Vector{Float64}
    return vec(sqrt.(sum(diff(Matrix(x); dims=1) .^ 2; dims=1)) .+ 1)
end

function correction_factor(ce_i::T, ce_j::T)::Float64 where {T<:Real}
    return max(ce_i, ce_j) / min(ce_i, ce_j)
end

function _complexity_invariance(
    data_i, data_j, complexity_i, complexity_j, dist_fn
)::Float64
    ed = dist_fn(data_i, data_j)
    cf = correction_factor(complexity_i, complexity_j)
    return ed * cf
end

"""
    complexity_invariance_distance(data::AbstractMatrix{<:Real}; distance=:euclidean)::AbstractMatrix{Float64}

Compute Complexity Invariance Distance (CID) between every pair of columns in `data`
(shape T × S, where T is timesteps and S is scenarios). Returns an S × S distance matrix.

# Arguments
- `data` : Matrix of shape T × S
- `distance` : `:euclidean` or `:weuclidean` (weighted Euclidean). Defaults to `:euclidean`
"""
function complexity_invariance_distance(
    data::AbstractMatrix{<:Real};
    distance=:euclidean
)::AbstractMatrix{Float64}
    complexity::Vector{Float64} = _complexity(data)

    n_timesteps, n_scenarios = size(data)
    cid_matrix::Matrix{Float64} = zeros(n_scenarios, n_scenarios)

    local weights::Vector{Float64}
    if distance == :weuclidean
        weights = sqrt.(1 ./ (1:n_timesteps))
    end
    dist_fn(x, y) = (distance == :euclidean) ? euclidean(x, y) : weuclidean(x, y, weights)

    for i in 1:n_scenarios
        Threads.@threads for j in (i + 1):n_scenarios
            cid_matrix[i, j] =
                cid_matrix[j, i] = _complexity_invariance(
                    data[:, i], data[:, j], complexity[i], complexity[j], dist_fn
                )
        end
    end

    return cid_matrix
end

"""
    cluster_series(data::AbstractMatrix{<:Real}, n_clusters::Int64; method::Symbol=:kmedoids, distance::Symbol=:euclidean)::Vector{Int64}

Cluster S scenarios with T timesteps each using Complexity Invariance Distance.

# Arguments
- `data` : Matrix of shape T × S
- `n_clusters` : Number of clusters
- `method` : `:kmedoids` (default) or `:hclust`
- `distance` : `:euclidean` (default) or `:weuclidean`

# Returns
Vector of cluster assignments (length S).

# References
1. Steinmann et al. (2020). Behavior-based scenario discovery using time series clustering.
   https://doi.org/10.1016/j.techfore.2020.120052
2. Batista et al. (2014). CID: an efficient complexity-invariant distance for time series.
   https://doi.org/10.1007/s10618-013-0312-3
"""
function cluster_series(
    data::AbstractMatrix{<:Real},
    n_clusters::Int64;
    method::Symbol=:kmedoids,
    distance::Symbol=:euclidean
)::Vector{Int64}
    distances = complexity_invariance_distance(data; distance=distance)

    if method == :kmedoids
        return kmedoids(distances, n_clusters).assignments
    end

    dendogram = hclust(distances; linkage=:average)
    return cutree(dendogram; k=n_clusters)
end

"""
    cluster_scenarios(data::AbstractArray{<:Real}, n_clusters::Int64; method::Symbol=:kmedoids, distance::Symbol=:euclidean)::Array{Int64}

Cluster scenarios across one or more metrics. For a 2D matrix, equivalent to `cluster_series`.
For a 3D array (timesteps × scenarios × metrics), clusters each metric slice independently
and returns an S × M matrix of assignments.

# Arguments
- `data` : Matrix (T × S) or 3D array (T × S × M)
- `n_clusters` : Number of clusters
- `method` : `:kmedoids` (default) or `:hclust`
- `distance` : `:euclidean` (default) or `:weuclidean`
"""
function cluster_scenarios(
    data::AbstractArray{<:Real},
    n_clusters::Int64;
    method::Symbol=:kmedoids,
    distance::Symbol=:euclidean
)::Array{Int64}
    ndims(data) == 2 && return cluster_series(data, n_clusters; method=method, distance=distance)

    _, n_scenarios, n_metrics = size(data)

    clusters = zeros(Int64, n_scenarios, n_metrics)
    for m in 1:n_metrics
        clusters[:, m] = cluster_series(
            data[:, :, m], n_clusters; method=method, distance=distance
        )
    end

    return clusters
end

end
