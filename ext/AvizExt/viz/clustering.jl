using JuliennedArrays: Slices
using Statistics

"""
    clustered_scenarios(data::AbstractMatrix, clusters::Vector{Int64}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    clustered_scenarios!(g::Union{GridLayout,GridPosition}, data::AbstractMatrix, clusters::Vector{Int64}; opts::Dict=Dict(), axis_opts::Dict=Dict())

Visualize clustered time series of scenarios.

# Arguments
- `data` : Matrix of time series data for several scenarios or sites
- `clusters` : Vector of numbers corresponding to clusters
- `opts` : Aviz options
    - `summarize` : plot confidence interval. Defaults to true

# Returns
Figure
"""
function ADRIA.viz.clustered_scenarios(
    data::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.clustered_scenarios!(g, data, clusters; axis_opts=axis_opts, opts=opts)

    return f
end
function ADRIA.viz.clustered_scenarios!(
    g::Union{GridLayout,GridPosition},
    data::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)::Union{GridLayout,GridPosition}
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(data)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    if get(opts, :summarize, true)
        _plot_clusters_confint!(g, ax, data, clusters; opts=opts)
    else
        _plot_clusters_series!(g, ax, data, clusters; opts=opts)
    end

    return g
end

function _plot_clusters_series!(
    g::Union{GridLayout,GridPosition},
    ax::Axis,
    data::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
)::Nothing
    sorted_clusters = sort(clusters)
    alphas = cluster_alphas(sorted_clusters)
    colors = unique(cluster_colors(sorted_clusters))

    for (idx, cluster) in enumerate(unique(sorted_clusters))
        series!(ax, data[:, clusters .== cluster]'; solid_color=(colors[idx], alphas[idx]))
    end

    _render_clustered_scenarios_legend(g, cluster_labels(clusters), colors)

    return nothing
end

function _plot_clusters_confint!(
    g::Union{GridLayout,GridPosition},
    ax::Axis,
    data::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
)::Nothing
    if :timesteps ∉ dimnames(data)
        throw(ArgumentError("data does not have :timesteps axis"))
    elseif ndims(data) != 2
        throw(ArgumentError("data does not have two dimensions"))
    end

    sorted_clusters = sort(clusters)
    colors = unique(cluster_colors(sorted_clusters))

    x_timesteps::UnitRange{Int64} = 1:length(timesteps(data))
    data_dims = dimnames(data)
    slice_dimension = data_dims[findfirst(data_dims .!= :timesteps)]

    for (idx_c, cluster) in enumerate(unique(sorted_clusters))
        y_lower, y_median, y_upper = eachcol(
            ADRIA.analysis.series_confint(
                data[:, clusters .== cluster]; agg_dim=slice_dimension
            ),
        )

        band!(ax, x_timesteps, y_lower, y_upper; color=(colors[idx_c], 0.5))
        scatterlines!(ax, y_median; color=colors[idx_c], markersize=5)
    end

    _render_clustered_scenarios_legend(g, cluster_labels(clusters), colors)
    return nothing
end

function _render_clustered_scenarios_legend(
    g::Union{GridLayout,GridPosition}, labels::Vector{String}, colors::Vector{RGBA{Float32}}
)::Nothing
    line_elems = [LineElement(; color=c, linestyle=nothing) for c in colors]
    Legend(g[1, 2], line_elems, labels; framevisible=false)

    return nothing
end

"""
    map(rs::Union{Domain,ResultSet}, data::AbstractMatrix, clusters::Vector{Int64}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    map(g, rs, data, clusters; opts::Dict=Dict(), axis_opts::Dict=Dict())

Visualize clustered time series for each site and map.

# Arguments
- `rs` : ResultSet
- `data` : Vector of summary statistics data for each location
- `clusters` : Vector of numbers corresponding to clusters
- `opts` : Options specific to this plotting method
    - `highlight` : Vector of colors indicating cluster membership for each location.
    - `summary` : function (which must support the `dims` keyword) to summarize data with.
                  Default: `mean`

# Returns
Figure
"""
function ADRIA.viz.map(
    rs::Union{Domain,ResultSet},
    data::AbstractArray{<:Real},
    clusters::Vector{Int64};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.map!(g, rs, data, clusters; opts=opts, axis_opts=axis_opts)

    return f
end
function ADRIA.viz.map!(
    g::Union{GridLayout,GridPosition},
    rs::Union{Domain,ResultSet},
    data::AbstractVector{<:Real},
    clusters::Vector{Int64};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)::Union{GridLayout,GridPosition}
    colors = cluster_colors(clusters)
    legend_params = _cluster_legend_params(clusters, unique(colors), data)

    opts[:highlight] = get(opts, :highlight, colors)
    opts[:legend_params] = get(opts, :legend_params, legend_params)

    ADRIA.viz.map!(g, rs, data; opts=opts, axis_opts=axis_opts)

    return g
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
function _cluster_legend_params(
    clusters::Vector{Int64}, clusters_colors::Vector, data::AbstractArray
)::Tuple
    # Filter non-zero clusters from clusters, colors and data
    non_zero_clusters = clusters .!= 0
    clusters_filtered = clusters[non_zero_clusters]
    statistics_filtered = data[non_zero_clusters]

    # Fill legend entries colors
    legend_entries = [
        PolyElement(; color=color, strokecolor=:transparent) for color in clusters_colors
    ]

    # Fill legend labels
    legend_labels = String[]
    clusters_numbers = unique(clusters_filtered)
    for cluster in clusters_numbers
        stat_mean = mean(statistics_filtered[cluster .== clusters_filtered])
        stat_mean_formatted = @sprintf "%.1e" stat_mean
        push!(legend_labels, "Cluster $(cluster): $stat_mean_formatted")
    end

    legend_title = "Clusters mean value"
    return (legend_entries, legend_labels, legend_title)
end
