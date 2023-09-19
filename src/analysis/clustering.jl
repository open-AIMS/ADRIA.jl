using Distances
using Clustering
using FLoops
using ADRIA
using JuliennedArrays

"""
    _complexity(x::AbstractMatrix{T})::AbstractMatrix{Float64} where {T <: Real}

Compute Complexity (CE) of an Matrix `x` of shape \$T ⋅ S\$, where \$T\$ is total number of
time steps and \$S\$ is number of scenarios.

# Arguments
- `x` : series matrix of shape \$T ⋅ S\$

# Returns
Vector of \$N\$ elements

# Examples
```julia-repl
julia> CE([[1; 2; 3] [1; 3; 4]))
Vector{Float64}:
 2
 5
"""
function _complexity(x::AbstractMatrix{T})::Vector{Float64} where {T<:Real}
    return vec(sqrt.(sum(diff(Matrix(x); dims=1) .^ 2; dims=1)) .+ 1)
end

"""
    correction_factor(ce_i::T, ce_j::T)::Float64 where {T<:Real}

Compute Correction Factor (CF) between two time series complexities `ce_i` and `ce_j`.

# Arguments
- `ce_i` : Time series `i`
- `ce_j` : Time series `j`

# Returns
Float64

# Examples
```julia-repl
julia> ce = CE([[1; 2; 3] [1; 3; 4]))
julia> CF(ce[1], ce[2])
Float64:
 2.5
"""
function correction_factor(ce_i::T, ce_j::T)::Float64 where {T<:Real}
    return max(ce_i, ce_j) / min(ce_i, ce_j)
end

"""
    complexity_invariance_distance(data::AbstractMatrix{T})::AbstractMatrix{Float64} where {T<:Real}

Compute Complexity Invariance Distance (CID) between every matrix column pairs (`data`) of
shape \$T ⋅ S\$. The distance between every two series is the weighted euclidian distanced
multiplied by the correction factor, which takes into account the ration between the two
series complexities. Returns a matrix of distances (\$S ⋅ S\$).

# Arguments
- `data` : Matrix of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is number of scenarios

# Returns
Matrix of complexity invariance distances
"""
function complexity_invariance_distance(
    data::AbstractMatrix{T}
)::AbstractMatrix{Float64} where {T<:Real}
    # Compute complexity vector
    complexity = _complexity(data)

    # Create empty Matrix
    data_size = size(data, 2)
    cid_matrix::AbstractMatrix{Float64} = zeros(data_size, data_size)

    # [1, 1/2, 1/3, ..., 1/n]
    weights = sqrt.(1 ./ (1:size(data, 1)))

    #? do we want to normalize the amplitudes of all series?

    # Iterate over data matrix to compute CID (Complexity Invariance Distance)
    for i in axes(data, 2)
        @floop for j in axes(data, 2)
            # Weight the WEuclidian Distance (ed) using the Correction Factor (cf)
            ed = weuclidean(data[:, i], data[:, j], weights)
            cf = correction_factor(complexity[i], complexity[j])
            cid_matrix[i, j] = cid_matrix[j, i] = ed * cf
        end
    end

    return cid_matrix
end

"""
    cluster_series(data::AbstractMatrix{T}, n_clusters::Int64)::Vector{Int64} where {T<:Real}

Hierarchically cluster \$S\$ scenarios with \$T\$ time steps each.


# Arguments
- `data` : Matrix of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is
  number of scenarios
- `n_clusters` : Number of clusters determined _a priori_.

# Returns
- `Vector` : Cluster ids indicating which cluster each scenario belongs to.

# References
1. Steinmann, P., Auping, W.L., Kwakkel, J.H., 2020.
   Behavior-based scenario discovery using time series clustering.
   Technological Forecasting and Social Change 156, 120052.
   https://doi.org/10.1016/j.techfore.2020.120052

2. Batista, G.E.A.P.A., Keogh, E.J., Tataw, O.M., de Souza, V.M.A., 2014.
   CID: an efficient complexity-invariant distance for time series.
   Data Min Knowl Disc 28, 634-669.
   https://doi.org/10.1007/s10618-013-0312-3
"""
function cluster_series(
    data::AbstractMatrix{T}, n_clusters::Int64
)::Vector{Int64} where {T<:Real}
    # Create dendogram using distantes matrix
    distances = complexity_invariance_distance(data)
    dendogram = hclust(distances; linkage=:average)

    # Return hierarchical clustering with n_clusters
    return cutree(dendogram; k=n_clusters)
end

"""
    cluster_scenarios(data::AbstractMatrix{T}, n_clusters::Int64)::Vector{Int64} where {T<:Real}

Alias to cluster_series.

# Arguments
- `data` : Matrix of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is
  number of scenarios
- `n_clusters` : Number of clusters determined _a priori_.

# Returns
- `Vector` : Cluster ids indicating which cluster each scenario belongs to.
"""
function cluster_scenarios(
    data::AbstractMatrix{T}, n_clusters::Int64
)::Vector{Int64} where {T<:Real}
    return cluster_series(data, n_clusters)
end

"""
    target_clusters(clusters::Vector{T}, outcomes::AbstractMatrix{F}; metric=temporal_variability, size_limit=0.01) where {T<:Int64,F<:Real}

Cluster scenarios into target and non target based on median outcome temporal variability of
previous time series cluster.

# Arguments
- `clusters` : Vector with outcome cluster indexes
- `outcomes` : AbstractMatrix of scenario outcomes
- `metric` : Metric used to aggregate outcomes for each cluster
- `size_limit` : This function will iteratively merge the best cluster with the second best
    if the fraction of scenarios inside it is below `size_limit`

# Returns
Vector containing 1's for target and 0's for non-target clusters
"""
function target_clusters(
    clusters::Vector{T},
    outcomes::AbstractMatrix{F};
    metric=temporal_variability,
    size_limit=0.01,
)::Vector{T} where {T<:Int64,F<:Real}

    # Compute statistic for each cluster
    clusters_statistics::Vector{Float64} = []
    for cluster in unique(clusters)
        normalized_outcomed = outcomes[:, clusters .== cluster] ./ maximum(outcomes)
        statistic = median(metric(normalized_outcomed))
        push!(clusters_statistics, statistic)
    end

    target_index = argmax(clusters_statistics)
    target_indexes = [target_index]

    # Merge target cluster if it is below 1% of size
    sizes = [size(outcomes[:, clusters .== c], 2) for c in unique(clusters)]
    target_size = sizes[target_index] / sum(sizes)
    while target_size < size_limit
        # Nullify target_index to find the next argmax
        clusters_statistics[target_index] = 0

        # Find next best cluster and add to target_indexes
        target_index = argmax(clusters_statistics)
        push!(target_indexes, target_index)

        # Update target_size with next best cluster size
        target_size += sizes[target_index] / sum(sizes)
    end

    # Return new clusters vector with only 1 and 0 for target and non-target clusters
    return [c ∈ target_indexes ? 1 : 0 for c in clusters]
end
