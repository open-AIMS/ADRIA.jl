"""
    ts_cluster(data::AbstractMatrix, clusters::Vector{Int64}; fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    ts_cluster!(g::Union{GridLayout,GridPosition}, data::AbstractMatrix, clusters::Vector{Int64}; axis_opts::Dict=Dict())

Visualize clustered time series of scenarios.

# Arguments
- `data` : Matrix of scenario data
- `clusters` : Vector of numbers corresponding to clusters

# Returns
Figure
"""
function ADRIA.viz.ts_cluster(data::AbstractMatrix, clusters::Vector{Int64};
    fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.ts_cluster!(g, data, clusters; axis_opts=axis_opts)

    return f
end
function ADRIA.viz.ts_cluster!(g::Union{GridLayout,GridPosition}, data::AbstractMatrix,
    clusters::Vector{Int64}; axis_opts::Dict=Dict())
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(data)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(
        g[1, 1],
        xticks=xtick_vals,
        xticklabelrotation=xtick_rot;
        axis_opts...
    )

    # Filter clusters and data for non-zero clusters
    clusters_filtered = filter(c -> c != 0, clusters)
    data_filtered = data[:, clusters.>0]

    # Compute cluster colors
    clusters_colors = _clusters_colors(clusters_filtered)
    unique_cluster_colors = unique(clusters_colors)

    leg_entry = Any[]
    for cluster in unique(clusters_filtered)
        cluster_color = _cluster_color(unique_cluster_colors, clusters_filtered, cluster)
        push!(
            leg_entry,
            series!(
                ax,
                data_filtered[:, clusters_filtered.==cluster]',
                solid_color=cluster_color
            )
        )
    end

    # Plot Legend
    n_clusters = length(unique(clusters_filtered))
    Legend(g[1, 2], leg_entry, "Cluster " .* string.(1:n_clusters), framevisible=false)

    return g
end

"""
    map(rs::Union{Domain,ResultSet}, data::AbstractMatrix, clusters::Vector{Int64}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    map(g, rs, data, clusters; opts::Dict=Dict(), axis_opts::Dict=Dict())

Visualize clustered time series for each site and map.

# Arguments
- `rs` : ResultSet
- `data` : Matrix of scenario data for each location
- `clusters` : Vector of numbers corresponding to clusters
- `opts` : Options specific to this plotting method
    - `highlight` : Vector of colors indicating cluster membership for each location.
    - `summary` : function (which must support the `dims` keyword) to summarize data with.
                  Default: `mean`

# Returns
Figure
"""
function ADRIA.viz.map(rs::Union{Domain,ResultSet}, data::AbstractMatrix,
    clusters::Vector{Int64}; opts::Dict=Dict(), fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict())

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.map!(g, rs, data, clusters; opts=opts, axis_opts=axis_opts)

    return f
end
function ADRIA.viz.map!(g::Union{GridLayout,GridPosition},
    rs::Union{Domain,ResultSet}, data::AbstractMatrix, clusters::Vector{Int64};
    opts::Dict=Dict(), axis_opts::Dict=Dict())

    # Vector of summary statistics (default is mean) computed over timesteps
    opts[:summary] = get(opts, :summary, mean)
    data_stats = ADRIA.metrics.per_loc(opts[:summary], data)

    cluster_colors = _clusters_colors(clusters)
    legend_params = _cluster_legend_params(clusters, cluster_colors, data_stats)

    opts[:highlight] = get(opts, :highlight, cluster_colors)
    opts[:legend_params] = get(opts, :legend_params, legend_params)

    ADRIA.viz.map!(g, rs, data_stats; opts=opts, axis_opts=axis_opts)

    return g
end

"""
    _cluster_colors(clusters::Vector{Int64})

Vector of cluster colors.

# Arguments
- `clusters` : Vector of numbers corresponding to clusters

# Returns
Vector{RGBA{Float32}}
"""
function _clusters_colors(clusters::Vector{Int64})::Vector{RGBA{Float32}}
    # Number of non-zero clusters
    n_clusters = length(unique(filter(cluster -> cluster != 0, clusters)))

    # Vector of clusters colors for non-zero clusters
    clusters_colors = categorical_colors(:seaborn_bright, n_clusters)

    # Assign color "black" to cluster == 0
    rgba_black = parse(RGBA{Float32}, "transparent")
    return [cluster == 0 ? rgba_black : clusters_colors[cluster] for cluster in clusters]
end

"""
    _cluster_color(unique_cluster_colors::Vector{RGBA{FLoat32}}, clusters::Vector{Int64}, cluster::Int64)

Color parameter for current cluster weighted by number of scenarios

# Arguments
- `unique_cluster_colors` : Vector of all cluster options that are being used
- `clusters` : Vector of numbers corresponding to clusters
- `cluster` : Current cluster

# Returns
Tuple{RGBA{Float32}, Float64}
"""
function _cluster_color(unique_cluster_colors::Vector{RGBA{Float32}},
    clusters::Vector{Int64}, cluster::Int64)::Tuple{RGBA{Float32},Float64}
    # Number of scenarios for a cluster
    n_scens = count(clusters .== cluster)

    # Compute line weight
    color_weight = max(min((1.0 / (n_scens * 0.05)), 0.6), 0.1)

    return (unique_cluster_colors[cluster], color_weight)
end

"""
    _cluster_legend_params(clusters::Vector{Int64}, clusters_colors::Vector{RGBA{Float32}}, data_statistics::Vector{Float32})

Color parameter for current cluster weighted by number of scenarios

# Arguments
- `clusters` : Vector of numbers corresponding to clusters
- `cluster_colors` : Vector of all cluster options that are being used
- `data_statistics` : Vector of statistics for each site (default is mean)

# Returns
Tuple{RGBA{Float32}, Float64}
"""
function _cluster_legend_params(clusters::Vector{Int64},
    clusters_colors::Vector{RGBA{Float32}}, data_statistics::Vector{Float32})::Tuple
    # Filter non-zero clusters from clusters, colors and data
    non_zero_clusters = clusters .!= 0
    clusters_filtered = clusters[non_zero_clusters]
    colors_filtered = clusters_colors[non_zero_clusters]
    statistics_filtered = data_statistics[non_zero_clusters]

    # Fill legend entries colors
    legend_entries = [PolyElement(color=color, strokecolor=:transparent) for
                      color in unique(colors_filtered)]

    # Fill legend labels
    legend_labels = String[]
    clusters_numbers = unique(clusters_filtered)
    for cluster in clusters_numbers
        stat_mean = mean(statistics_filtered[cluster.==clusters_filtered])
        stat_mean_formatted = @sprintf "%.1e" stat_mean
        push!(legend_labels, "Cluster $(cluster): $stat_mean_formatted")
    end

    legend_title = "Clusters mean value"
    return (legend_entries, legend_labels, legend_title)
end
