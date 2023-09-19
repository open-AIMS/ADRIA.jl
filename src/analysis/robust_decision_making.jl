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
"""
function find_robust_clustering(
    result_set::ResultSet,
    metrics::Vector{ADRIA.metrics.Metric};
    num_clusters::Int64=6,
    robustness_limit=0.2,
)::BitVector
    _robust_scenarios = trues(size(result_set.inputs, 1))

    # Compute summary statistics for each metric
    for metric in metrics
        rs_metric = metric(result_set)

        clusters = ADRIA.analysis.cluster_scenarios(rs_metric, num_clusters)

        # Where statistics for each cluster will be saved
        clusters_summary = zeros(num_clusters)

        for (idx_c, c) in enumerate(unique(clusters))
            cluster_metric = rs_metric[:, clusters .== c]

            # Compute median series for current cluster
            # This could be extracted to a separate function
            tf = axes(cluster_metric, :timesteps)
            timesteps_slices = JuliennedArrays.Slices(
                cluster_metric[timesteps=tf], NamedDims.dim(cluster_metric, :scenarios)
            )
            median_series = median.(timesteps_slices)

            # Summary statistics for that cluster metric
            # clusters_summary[idx_c] = sum(cumsum(median_series))
            clusters_summary[idx_c] = temporal_variability(median_series)
        end

        # Robust clusters are the ones with summary statistics above the limit
        robust_clusters = clusters_summary .>= quantile(clusters_summary, robustness_limit)

        # Select scenarios that are robust across all metrics
        aux_rc::BitVector = clusters .∈ Vector{Int64}[unique(clusters)[robust_clusters]]

        _robust_scenarios = _robust_scenarios .& aux_rc
    end

    return _robust_scenarios
end
