using Distances
using Clustering
using FLoops
using ADRIA
using JuliennedArrays

"""
    _complexity(x::AbstractMatrix{<:Real})::AbstractMatrix{Float64}

Compute Complexity (CE) of an Matrix `x` of shape \$T ⋅ S\$, where \$T\$ is total number of
time steps and \$S\$ is number of scenarios.

# Arguments
- `x` : series matrix of shape \$T ⋅ S\$

# Returns
Vector of \$N\$ elements
"""
function _complexity(x::AbstractMatrix{<:Real})::Vector{Float64}
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
```julia
julia> ce = complexity([[1, 2, 3] [1, 3, 4]])
julia> correction_factor(ce[1], ce[2])
Float64:
 2.5
 ```
"""
function correction_factor(ce_i::T, ce_j::T)::Float64 where {T<:Real}
    return max(ce_i, ce_j) / min(ce_i, ce_j)
end

"""
    complexity_invariance_distance(data::AbstractMatrix{<:Real})::AbstractMatrix{Float64}

Compute Complexity Invariance Distance (CID) between every matrix column pairs (`data`) of
shape \$T ⋅ S\$. The distance between every two series is the weighted euclidian distanced
multiplied by the correction factor, which takes into account the ration between the two
series complexities. Returns a matrix of distances (\$S ⋅ S\$).

# Arguments
- `data` : Matrix of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is
number of scenarios

# Returns
Matrix of complexity invariance distances.
"""
function complexity_invariance_distance(
    data::AbstractMatrix{<:Real}
)::AbstractMatrix{Float64}
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
            # ed = weuclidean(data[:, i], data[:, j], weights)
            ed = euclidean(data[:, i], data[:, j])
            cf = correction_factor(complexity[i], complexity[j])
            cid_matrix[i, j] = cid_matrix[j, i] = ed * cf
        end
    end

    return cid_matrix
end

"""
    cluster_series(data::AbstractMatrix{<:Real}, n_clusters::Int64)::Vector{Int64}

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
function cluster_series(data::AbstractMatrix{<:Real}, n_clusters::Int64;
    method=:kmedoids)::Vector{Int64}
    # Calculate distantes matrix
    distances = complexity_invariance_distance(data)

    if method == :kmedoids
        return kmedoids(distances, n_clusters).assignments
    end

    # Return hierarchical clustering with n_clusters
    dendogram = hclust(distances; linkage=:average)
    return cutree(dendogram; k=n_clusters)
end

"""
    cluster_scenarios(data::AbstractArray{<:Real}, n_clusters::Int64)::Array{Int64}

Alias to cluster_series.

# Arguments
- `data` : Matrix of \$T ⋅ S\$, where \$T\$ is total number of time steps and \$S\$ is
  number of scenarios
- `n_clusters` : Number of clusters determined _a priori_.

# Returns
- `Vector` : Cluster ids indicating which cluster each scenario belongs to.

# Examples
One can cluster scenarios based on a single Metric, passing a Matrix of outcomes for each
timestep and scenario:

```julia
# Matrix of outcomes
s_tac = ADRIA.metrics.scenario_total_cover(rs)

# Cluster scenarios
n_cluster = 6
clusters = ADRIA.analysis.cluster_series(s_tac, n_clusters)
```

And perform multiple clusterings, based on multiple Metrics, passing a 3-dimensional Array
(or YAXArray) of outcomes for each timestep, scenario and Metric.

```julia
metrics::Vector{ADRIA.metrics.Metric} = [
    ADRIA.metrics.scenario_total_cover,
    ADRIA.metrics.scenario_asv,
    ADRIA.metrics.scenario_absolute_juveniles,
]

# 3-dimensional array of outcomes
outcomes = ADRIA.metrics.scenario_outcomes(rs, metrics)

# Cluster scenarios
num_clusters = 6
outcomes_clusters = ADRIA.analysis.cluster_scenarios(outcomes, num_clusters)
```
"""
function cluster_scenarios(data::AbstractArray{<:Real}, n_clusters::Int64)::Array{Int64}
    ndims(data) == 2 && return cluster_series(data, n_clusters)

    _, n_scenarios, n_metrics = size(data)

    clusters = zeros(Int64, n_scenarios, n_metrics)
    for m in 1:n_metrics
        clusters[:, m] = cluster_series(data[:, :, m], n_clusters)
    end

    return clusters
end

"""
    target_clusters(clusters::Vector{T}, outcomes::AbstractMatrix{<:Real}; metric=temporal_variability, size_limit=0.01) where {T<:Int64}

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
    outcomes::AbstractMatrix{<:Real};
    metric=temporal_variability,
    size_limit=0.01,
)::Vector{T} where {T<:Int64}

    # Compute statistic for each cluster
    clusters_statistics::Vector{Float64} = []
    for cluster in unique(clusters)
        normalized_outcomed = outcomes[:, clusters.==cluster] ./ maximum(outcomes)
        statistic = median(metric(normalized_outcomed))
        push!(clusters_statistics, statistic)
    end

    target_index = argmax(clusters_statistics)
    target_indexes = [target_index]

    # Merge target cluster if it is below 1% of size
    sizes = [size(outcomes[:, clusters.==c], 2) for c in unique(clusters)]
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

"""
    find_scenarios(outcomes::AbstractMatrix{<:Real}, clusters::Vector{Int64}, filter_func::Function, aggregation_func::Function=temporal_variability)::BitVector
    find_scenarios(outcomes::AbstractArray{<:Real,3}, clusters::AbstractMatrix{Int64}, filter_funcs::Vector{Function}; aggregation_function::Function=temporal_variability)::BitVector
    find_scenarios(outcomes::AbstractArray{<:Real,3}, clusters::AbstractMatrix{Int64}, filter_func::Function; aggregation_function::Function=temporal_variability)::BitVector

If outcomes is Matrix of scenario outcomes and clusters is a Vector of clusters:
- Computes a median series for each cluster
- Use aggregation_func to compute a summary statistics for each median series
- Select scenarios for which `filter_func` returns true

If outcomes is a 3-dimensional array of scenario outcomes:
- Computes a median series for each outcome cluster
- Use aggregation_func to compute a summary statistics for each median series
- Select scenarios for which `filter_func` returns true for each matrix of outcomes
- Select scenarios that were selected for all outcomes

# Arguments
- `outcomes` : Outcomes for one or more scenario metrics
- `clusters` : Clusters for one or more scenario metric outcomes
- `filter_funcs` : Function used to filter/target clusters
- `aggregation_function` : Function used to aggregate each median temporal series into a
    single number (default is temporal_variability)

# Returns
BitVector with true/false for selected/not selected scenarios

# Examples
```julia
metrics = [
    ADRIA.metrics.scenario_total_cover,
    ADRIA.metrics.scenario_asv
]

# Get outcomes
outcomes = ADRIA.metrics.scenario_outcomes(rs, metrics)
num_clusters = 6

# Cluster scenarios based on outcomes
outcomes_clusters = ADRIA.analysis.cluster_scenarios(outcomes, num_clusters)

# Find scenarios above 0.25-quantile for all metrics
robustness_func(x) = x .>= quantile(x, 0.25)
robust_scens = ADRIA.analysis.find_scenarios(outcomes, outcomes_clusters, robustness_func)

# Find scenarios in the three highest clusters for all metrics
highest_clusters(x) = x .>= x .∈ [sort(x; rev=true)[1:3]]
high_scens = ADRIA.analysis.find_scenarios(outcomes, outcomes_clusters, highest_clusters)
```
"""
function find_scenarios(
    outcomes::AbstractMatrix{<:Real},
    clusters::AbstractVector{Int64},
    filter_func::Function;
    aggregation_func::Function=temporal_variability,
)::BitVector
    clusters_summary::Vector{Float64} = zeros(length(unique(clusters)))

    for (idx_c, c) in enumerate(unique(clusters))
        cluster_metric = outcomes[:, clusters.==c]

        # Median series for current cluster
        tf = axes(cluster_metric, :timesteps)
        timesteps_slices = JuliennedArrays.Slices(cluster_metric[timesteps=tf], 2)
        median_series = median.(timesteps_slices)

        # Summary statistics for that cluster metric
        clusters_summary[idx_c] = aggregation_func(median_series)
    end

    return clusters .∈ [unique(clusters)[filter_func(clusters_summary)]]
end
function find_scenarios(
    outcomes::AbstractArray{<:Real,3},
    clusters::AbstractMatrix{Int64},
    filter_funcs::Vector{Function};
    aggregation_func::Function=temporal_variability,
)::BitVector
    scenarios = trues(size(clusters, 1))

    # Find scenarios for each clustered outcomes
    for (idx, clust) in enumerate(eachcol(clusters))
        scenarios =
            scenarios .& find_scenarios(
                outcomes[:, :, idx],
                collect(clust),
                filter_funcs[idx];
                aggregation_func=aggregation_func,
            )
    end

    return scenarios
end
function find_scenarios(
    outcomes::AbstractArray{<:Real,3},
    clusters::AbstractMatrix{Int64},
    filter_func::Function;
    aggregation_func::Function=temporal_variability,
)::BitVector
    filter_funcs::Vector{Function} = fill(filter_func, size(clusters, 2))
    return find_scenarios(
        outcomes, clusters, filter_funcs; aggregation_func=aggregation_func
    )
end
