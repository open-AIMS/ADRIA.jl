using ADRIA.performance: temporal_variability

"""
    target_clusters(clusters::Vector{<:Integer}, outcomes::AbstractMatrix{<:Real}; metric=temporal_variability, size_limit=0.01)

Classify scenario clusters into target and non-target groups based on median outcome
temporal variability of each cluster.

# Arguments
- `clusters` : Vector of cluster index assignments for each scenario
- `outcomes` : Matrix of scenario outcomes (timesteps × scenarios)
- `metric` : Function used to aggregate outcomes within each cluster. Defaults to
  `temporal_variability`
- `size_limit` : Minimum fraction of scenarios required in the target cluster. If the best
  cluster is smaller, it is iteratively merged with the next-best cluster until the
  size threshold is met

# Returns
`Vector{T}` with 1 for target scenarios and 0 for non-target scenarios
"""
function target_clusters(
    clusters::Vector{<:Integer},
    outcomes::AbstractMatrix{<:Real};
    metric=temporal_variability,
    size_limit=0.01
)::Vector{Int64}
    clusters_statistics::Vector{Float64} = []
    for cluster in unique(clusters)
        normalized_outcomed = outcomes[:, clusters .== cluster] ./ maximum(outcomes)
        statistic = median(metric(normalized_outcomed))
        push!(clusters_statistics, statistic)
    end

    target_index = argmax(clusters_statistics)
    target_indexes = [target_index]

    sizes = [size(outcomes[:, clusters .== c], 2) for c in unique(clusters)]
    target_size = sizes[target_index] / sum(sizes)
    while target_size < size_limit
        clusters_statistics[target_index] = 0
        target_index = argmax(clusters_statistics)
        push!(target_indexes, target_index)
        target_size += sizes[target_index] / sum(sizes)
    end

    return [c ∈ target_indexes ? 1 : 0 for c in clusters]
end

"""
    find_scenarios(outcomes::AbstractArray{<:Real}, clusters::AbstractVector{Int64}, filter_func::Function; aggregation_func::Function=temporal_variability)::BitVector
    find_scenarios(outcomes::AbstractArray{<:Real,3}, clusters::AbstractMatrix{Int64}, filter_funcs::Vector{Function}; aggregation_func::Function=temporal_variability)::BitVector
    find_scenarios(outcomes::AbstractArray{<:Real,3}, clusters::AbstractMatrix{Int64}, filter_func::Function; aggregation_func::Function=temporal_variability)::BitVector

Select scenarios based on cluster-level summary statistics.

For a 2D `outcomes` matrix (timesteps × scenarios):
- Computes a median series for each cluster
- Applies `aggregation_func` to each median series
- Returns scenarios for which `filter_func(summary_statistics)` is true

For a 3D `outcomes` array (timesteps × scenarios × metrics):
- Applies the 2D logic per metric slice
- Returns scenarios selected across all metric slices (intersection)

# Arguments
- `outcomes` : Outcomes for one or more scenario metrics
- `clusters` : Cluster assignments (Vector for 2D, Matrix for 3D outcomes)
- `filter_func` / `filter_funcs` : Function(s) to select target clusters from summary stats
- `aggregation_func` : Aggregates each cluster's median series to a scalar.
  Defaults to `temporal_variability`

# Returns
`BitVector` — `true` for selected scenarios, `false` otherwise

# Examples
```julia
metrics = [ADRIA.metrics.scenario_total_cover, ADRIA.metrics.scenario_asv]
outcomes = ADRIA.metrics.scenario_outcomes(rs, metrics)
clusters = ADRIA.analysis.cluster_scenarios(outcomes, 6)
robust_scens = ADRIAanalysis.find_scenarios(outcomes, clusters, x -> x .>= quantile(x, 0.25))
```
"""
function find_scenarios(
    outcomes::AbstractArray{<:Real},
    clusters::AbstractVector{Int64},
    filter_func::Function;
    aggregation_func::Function=temporal_variability
)::BitVector
    clusters_summary::Vector{Float64} = zeros(length(unique(clusters)))

    for (idx_c, c) in enumerate(unique(clusters))
        cluster_metric = outcomes[:, clusters .== c]
        median_series = vec(
            median(cluster_metric[timesteps = axes(cluster_metric, :timesteps)]; dims=2)
        )
        clusters_summary[idx_c] = aggregation_func(median_series)
    end

    return clusters .∈ [unique(clusters)[filter_func(clusters_summary)]]
end
function find_scenarios(
    outcomes::AbstractArray{<:Real,3},
    clusters::AbstractMatrix{Int64},
    filter_funcs::Vector{Function};
    aggregation_func::Function=temporal_variability
)::BitVector
    scenarios = trues(size(clusters, 1))

    for (idx, clust) in enumerate(eachcol(clusters))
        scenarios =
            scenarios .& find_scenarios(
                outcomes[:, :, idx],
                collect(clust),
                filter_funcs[idx];
                aggregation_func=aggregation_func
            )
    end

    return scenarios
end
function find_scenarios(
    outcomes::AbstractArray{<:Real,3},
    clusters::AbstractMatrix{Int64},
    filter_func::Function;
    aggregation_func::Function=temporal_variability
)::BitVector
    filter_funcs::Vector{Function} = fill(filter_func, size(clusters, 2))
    return find_scenarios(
        outcomes, clusters, filter_funcs; aggregation_func=aggregation_func
    )
end
