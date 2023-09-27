include("pareto.jl")

"""
    function find_scenarios(outcomes::AbstractArray, clusters::AbstractMatrix, filter_functions::Vector{Function}; aggregation_function::Function=temporal_variability)::BitVector
    function find_scenarios(outcomes::AbstractArray, clusters::AbstractMatrix, filter_function::Function; aggregation_function::Function=temporal_variability)::BitVector

Find scenarios across a list of given metrics outcomes using clustering strategy:
For each metric outcome:
    - A median series is computed for each previously computed cluster (passed as argument)
    - A summary statistics is computed for each median series (default is temporal variance)
    - Selected clusters are those whose summary statistics is true for given function
Selected scenarios are, then, those that belong to the selected cluster for all outcomes

# Arguments
- `outcomes_clusters` : Matrix where each col is a cluster vector for a different outcome,
outcomes::AbstractArray,
robustness_funcs::Vector{Function};
temporal_aggregation_func::Function=temporal_variability,

# Returns
BitVector with true for robust scenarios and false for non-robust

# Example
```julia
metrics::Vector{ADRIA.metrics.Metric} = [
    ADRIA.metrics.scenario_total_cover,
    ADRIA.metrics.scenario_asv
]

# Get outcomes
outcomes = ADRIA.metrics.scenario_outcomes(rs, metrics)
num_clusters = 6

# Cluster scenarios based on outcomes
outcomes_clusters::AbstractMatrix{Int64} = ADRIA.analysis.cluster_scenarios(
    outcomes, num_clusters
)

robustness_limit = 0.25
robustness_func(x) = x .>= quantile(x, robustness_limit)
robust_scenarios = ADRIA.analysis.find_scenarios(outcomes, outcomes_clusters, robustness_func)
```
"""
function find_scenarios(
    outcomes::AbstractArray,
    clusters::AbstractMatrix,
    filter_functions::Vector{Function};
    aggregation_function::Function=temporal_variability,
)::BitVector
    robust_scenarios = trues(size(clusters, 1))

    # Compute summary statistics for each metric
    for (idx, cluster) in enumerate(eachcol(clusters))
        robust_clusters::BitVector = _find_clusters(
            outcomes[:, :, idx],
            collect(cluster),
            filter_functions[idx],
            aggregation_function,
        )

        # Select scenarios that are robust across all metrics
        aux_rc::BitVector = cluster .∈ Vector{Int64}[unique(cluster)[robust_clusters]]

        robust_scenarios = robust_scenarios .& aux_rc
    end

    return robust_scenarios
end
function find_scenarios(
    outcomes::AbstractArray,
    clusters::AbstractMatrix,
    filter_function::Function;
    aggregation_function::Function=temporal_variability,
)::BitVector
    filter_functions::Vector{Function} = fill(filter_function, size(clusters, 2))
    return find_scenarios(
        outcomes, clusters, filter_functions; aggregation_function=aggregation_function
    )
end

function _find_clusters(
    outcomes::AbstractArray,
    clusters::Vector{Int64},
    filter_function::Function,
    aggregation_function::Function,
)::BitVector
    clusters_summary::Vector{Float64} = zeros(length(unique(clusters)))

    for (idx_c, c) in enumerate(unique(clusters))
        cluster_metric = outcomes[:, clusters .== c]

        # Median series for current cluster
        tf = axes(cluster_metric, :timesteps)
        timesteps_slices = JuliennedArrays.Slices(
            cluster_metric[timesteps=tf], NamedDims.dim(cluster_metric, :scenarios)
        )
        median_series = median.(timesteps_slices)

        # Summary statistics for that cluster metric
        clusters_summary[idx_c] = aggregation_function(median_series)
    end

    # Robust clusters are the ones with summary statistics above the limit
    return filter_function(clusters_summary)
end
