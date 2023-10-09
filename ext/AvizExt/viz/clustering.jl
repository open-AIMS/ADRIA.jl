using JuliennedArrays: Slices
using Statistics

function ADRIA.viz.scenarios(
    outcomes::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.scenarios!(
        g, outcomes, clusters; opts=opts, axis_opts=axis_opts, series_opts=series_opts
    )

    return f
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    outcomes::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Union{GridLayout,GridPosition}
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(outcomes)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    scen_groups = ADRIA.analysis.scenario_clusters(clusters)
    opts::Dict{Symbol,Any} = Dict(:histogram => false)
    return ADRIA.viz.scenarios!(
        g,
        ax,
        outcomes,
        scen_groups;
        opts=opts,
        axis_opts=axis_opts,
        series_opts=series_opts,
    )
end

"""
    clustered_scenarios(outcomes::AbstractMatrix, clusters::Vector{Int64}; opts::Dict=Dict(), fig_opts::Dict=Dict(), axis_opts::Dict=Dict())
    clustered_scenarios!(g::Union{GridLayout,GridPosition}, outcomes::AbstractMatrix, clusters::Vector{Int64}; opts::Dict=Dict(), axis_opts::Dict=Dict())

Visualize clustered time series of scenarios.

# Arguments
- `outcomes` : Matrix of outcomes for several scenarios or sites
- `clusters` : Vector of numbers corresponding to clusters
- `opts` : Aviz options
    - `summarize` : plot confidence interval. Defaults to true

# Returns
Figure
"""
function ADRIA.viz.clustered_scenarios(
    outcomes::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.scenarios!(g, outcomes, clusters; axis_opts=axis_opts, opts=opts)

    return f
end
function ADRIA.viz.clustered_scenarios!(
    g::Union{GridLayout,GridPosition},
    outcomes::NamedDimsArray,
    clusters::Union{BitVector,Vector{Int64}};
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
)::Union{GridLayout,GridPosition}
    return ADRIA.viz.scenarios!(g, outcomes, clusters; axis_opts=axis_opts, opts=opts)
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
