include("pareto.jl")

"""
    find_robust_clustering(result_set::ResultSet, metrics::Vector{ADRIA.metrics.Metric}; num_clusters::Int64=6)

Find scenarios that perform well across a list of given metrics using clustering strategy:
For each metric outcome:
    - The scenarios are clusterized for each metric outcome.
    - A median series is computed for each cluster.
    - A summary statistics is computed for each median series (temporal variance).
    -  Robust Clusters are those whose summary statistics is above the median.
Robust Scenarios are, then, those that belong to the Robust Cluster for all metrics.

# Arguments
- `result_set` : ResultSet
- `metrics` : Vector of metrics to consider when looking for robust clusters
- `num_clusters` : Number of clusters to consider when clustering scenarios for each metric

# Returns
BitVector with true for robust scenarios and false for non-robust

# Example
```julia
rs::ResultSet = ADRIA.run_scenarios(samples, domain, "45")

_metrics::Vector{ADRIA.metrics.Metric} = [
    ADRIA.metrics.scenario_total_cover,
    ADRIA.metrics.scenario_asv,
    ADRIA.metrics.scenario_absolute_juveniles,
]

robustness_limit::Float64 = 0.25

robust_scens = ADRIA.analysis.find_robust_clustering(
    rs, _metrics; robustness_limit=robustness_limit
)
```
"""
function find_robust_clustering(
    result_set::ResultSet,
    metrics::Vector{ADRIA.metrics.Metric};
    num_clusters::Int64=6,
    robustness_limit=0.2,
)::BitVector
    robust_scenarios = trues(size(result_set.inputs, 1))

    # Compute summary statistics for each metric
    for metric in metrics
        metric_outcome::NamedDims.NamedDimsArray = metric(result_set)

        clusters::Vector{Int64} = ADRIA.analysis.cluster_scenarios(
            metric_outcome, num_clusters
        )

        robust_clusters::BitVector = _find_robust(
            clusters, metric_outcome, robustness_limit
        )

        # Select scenarios that are robust across all metrics
        aux_rc::BitVector = clusters .∈ Vector{Int64}[unique(clusters)[robust_clusters]]

        robust_scenarios = robust_scenarios .& aux_rc
    end

    return robust_scenarios
end

function _find_robust(
    clusters::Vector{Int64},
    metric_outcome::NamedDims.NamedDimsArray,
    robustness_limit::Float64,
)::BitVector
    clusters_summary = zeros(length(unique(clusters)))

    for (idx_c, c) in enumerate(unique(clusters))
        cluster_metric = metric_outcome[:, clusters .== c]

        # Median series for current cluster
        tf = axes(cluster_metric, :timesteps)
        timesteps_slices = JuliennedArrays.Slices(
            cluster_metric[timesteps=tf], NamedDims.dim(cluster_metric, :scenarios)
        )
        median_series = median.(timesteps_slices)

        # Summary statistics for that cluster metric
        clusters_summary[idx_c] = temporal_variability(median_series)
    end

    # Robust clusters are the ones with summary statistics above the limit
    return clusters_summary .>= quantile(clusters_summary, robustness_limit)
end
