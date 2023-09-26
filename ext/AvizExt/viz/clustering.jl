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
    clusters::Vector{Int64};
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
    clusters::Vector{Int64};
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
    clusters::Vector{Int64};
    opts::Dict=Dict(),
)::Nothing
    unique_clusters = sort(unique(clusters))
    alphas = _get_alphas(clusters)
    clusters_colors = unique(_get_colors(clusters, alphas))

    leg_entry::Vector{Any} = [
        series!(ax, data[:, clusters .== cluster]'; solid_color=clusters_colors[idx_c]) for
        (idx_c, cluster) in enumerate(unique_clusters)
    ]

    Legend(g[1, 2], leg_entry, "Cluster " .* string.(unique_clusters); framevisible=false)

    return nothing
end

function _plot_clusters_confint!(
    g::Union{GridLayout,GridPosition},
    ax::Axis,
    data::NamedDimsArray,
    clusters::Vector{Int64};
    opts::Dict=Dict(),
)::Nothing
    if :timesteps ∉ dimnames(data)
        throw(ArgumentError("data does not have :timesteps axis"))
    elseif ndims(data) != 2
        throw(ArgumentError("data does not have two dimensions"))
    end

    unique_clusters = sort(unique(clusters))

    band_colors = unique(_get_colors(clusters, 0.5))
    line_colors = unique(_get_colors(clusters))
    legend_entry = Vector{Any}(undef, length(unique_clusters))

    x_timesteps::UnitRange{Int64} = 1:length(timesteps(data))
    slice_dimension = filter(x -> x != :timesteps, dimnames(data))[1]

    for (idx_c, cluster) in enumerate(unique_clusters)
        band_color = band_colors[idx_c]
        line_color = line_colors[idx_c]

        y_lower, y_median, y_upper = confint(data[:, clusters .== cluster], slice_dimension)

        band!(ax, x_timesteps, y_lower, y_upper; color=band_color)
        legend_entry[idx_c] = scatterlines!(ax, y_median; color=line_color, markersize=5)
    end

    Legend(
        g[1, 2], legend_entry, "Cluster " .* string.(unique_clusters); framevisible=false
    )

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
    cluster_colors = _get_colors(clusters)
    legend_params = _cluster_legend_params(clusters, cluster_colors, data)

    opts[:highlight] = get(opts, :highlight, cluster_colors)
    opts[:legend_params] = get(opts, :legend_params, legend_params)

    ADRIA.viz.map!(g, rs, data; opts=opts, axis_opts=axis_opts)

    return g
end

"""
    _get_colors(clusters::Vector{Int64}, alphas::Vector{Float64})::Vector{RGBA{Float32}}
    _get_colors(clusters::Vector{Int64}, alpha::Float64)::Vector{RGBA{Float32}}
    _get_colors(clusters::Vector{Int64})::Vector{RGBA{Float32}}
Vector of cluster colors.

# Arguments
- `clusters` : Vector of numbers corresponding to clusters
- `alphas` : Vector of alphas to be used for cluster color

# Returns
Vector{RGBA{Float32}}
"""
function _get_colors(
    clusters::Vector{Int64}, alphas::Vector{Float64}
)::Vector{RGBA{Float32}}
    # Number of non-zero clusters
    n_clusters = length(unique(filter(cluster -> cluster != 0, clusters)))

    # Vector of clusters colors for non-zero clusters
    colors = categorical_colors(:seaborn_bright, n_clusters)

    # Apply alpha to cluster colors
    if !isempty(alphas)
        colors = [RGBA(c.r, c.g, c.b, a) for (c, a) in zip(colors, alphas)]
    end

    # Assign color "black" to cluster == 0
    rgba_black = parse(RGBA{Float32}, "transparent")
    return [cluster == 0 ? rgba_black : colors[cluster] for cluster in clusters]
end
function _get_colors(clusters::Vector{Int64}, alpha::Float64)::Vector{RGBA{Float32}}
    alphas::Vector{Float64} = fill(alpha, length(clusters))
    return _get_colors(clusters, alphas)
end
function _get_colors(clusters::Vector{Int64})::Vector{RGBA{Float32}}
    return _get_colors(clusters, 1.0)
end

"""
    _get_alphas(clusters::Vector{Int64})::Vector{Float64}

Vector of color alphas for each clusters weighted by number of scenarios

# Arguments
- `clusters` : Vector with scenario cluster numbers

# Returns
Vector with one color alpha for each cluster
"""
function _get_alphas(clusters::Vector{Int64})::Vector{Float64}
    alphas::Vector{Float64} = zeros(Float64, length(unique(clusters)))

    for (i, cluster) in enumerate(unique(clusters))
        n_scens = count(clusters .== cluster)
        base_alpha = 1.0 / (n_scens * 0.05)
        alphas[i] = max(min(base_alpha, 0.6), 0.1)
    end

    return alphas
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
    clusters::Vector{Int64},
    clusters_colors::Vector{RGBA{Float32}},
    data_statistics::AbstractArray,
)::Tuple
    # Filter non-zero clusters from clusters, colors and data
    non_zero_clusters = clusters .!= 0
    clusters_filtered = clusters[non_zero_clusters]
    colors_filtered = clusters_colors[non_zero_clusters]
    statistics_filtered = data_statistics[non_zero_clusters]

    # Fill legend entries colors
    legend_entries = [
        PolyElement(; color=color, strokecolor=:transparent) for
        color in unique(colors_filtered)
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
