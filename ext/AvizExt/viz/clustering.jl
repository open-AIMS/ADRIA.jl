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

    _render_clustered_scenarios_legend(g, cluster_labels(sorted_clusters), colors)

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
    cluster_ids = unique(sorted_clusters)

    colors = unique(cluster_colors(sorted_clusters))
    alpha = get(opts, :alpha, 0.4)

    n_timesteps = length(timesteps(data))
    x_timesteps::UnitRange{Int64} = 1:n_timesteps
    data_dims = dimnames(data)
    slice_dimension = data_dims[findfirst(data_dims .!= :timesteps)]

    confints = zeros(n_timesteps, length(cluster_ids), 3)
    for (idx_c, cluster) in enumerate(cluster_ids)
        confints[:, idx_c, :] = ADRIA.analysis.series_confint(
            data[:, clusters .== cluster]; agg_dim=slice_dimension
        )
    end

    for idx in eachindex(cluster_ids)
        y_lower, y_upper = confints[:, idx, 1], confints[:, idx, 3]
        band!(ax, x_timesteps, y_lower, y_upper; color=(colors[idx], alpha))
    end

    series!(ax, confints[:, :, 2]'; solid_color=colors)

    _render_clustered_scenarios_legend(g, cluster_labels(sorted_clusters), colors)
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
    clusters::Union{BitVector,Vector{Int64}};
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
    clusters::Union{BitVector,Vector{Int64}};
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

Color parameter for current cluster weighted by number of scenarios.

# Arguments
- `clusters` : Vector of numbers corresponding to clusters
- `colors` : Vector of all cluster options that are being used
- `data` : Vector of some metric outcome for each site

# Returns
Tuple of legend params to be passed to map! containing legend_entries, legend_labels and
legend_title (in that order).
"""
function _cluster_legend_params(
    clusters::Union{BitVector,Vector{Int64}},
    colors::Vector{RGBA{Float32}},
    data::AbstractVector{<:Real},
)::Tuple
    legend_entries = [PolyElement(; color=c, strokecolor=:transparent) for c in colors]
    legend_labels = cluster_labels(clusters, data)
    legend_title = "Clusters mean value"

    return (legend_entries, legend_labels, legend_title)
end
