using Statistics

"""
    scenarios(outcomes::AbstractMatrix, clusters::Vector{Int64}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    scenarios!(g::Union{GridLayout,GridPosition}, outcomes::AbstractMatrix, clusters::Vector{Int64}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Visualize clustered time series of scenarios.

# Arguments
- `outcomes` : AbstractMatrix of outcomes for several scenarios or locations
- `clusters` : Vector of numbers corresponding to clusters
- `opts` : Aviz options
    - `summarize` : plot confidence interval. Defaults to true

# Returns
Figure
"""
function ADRIA.viz.scenarios(
    outcomes::AbstractMatrix{<:Real},
    clusters::Union{BitVector,AbstractVector{Int64}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.scenarios!(g, outcomes, clusters; opts=opts, axis_opts=axis_opts)

    return f
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    outcomes::AbstractMatrix{<:Real},
    clusters::Union{BitVector,Vector{Int64}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(outcomes)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    scen_groups = ADRIA.analysis.scenario_clusters(clusters)
    opts[:histogram] = false
    opts[:legend_labels] = sort(collect(keys(scen_groups)))

    return ADRIA.viz.scenarios!(
        g,
        ax,
        outcomes,
        scen_groups;
        opts=opts,
        axis_opts=axis_opts,
        series_opts=series_opts
    )
end

"""
    clustered_scenarios(outcomes::AbstractMatrix{<:Real}, clusters::Vector{Int64}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    clustered_scenarios!(g::Union{GridLayout,GridPosition}, outcomes::AbstractMatrix{<:Real}, clusters::Vector{Int64}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Visualize clustered time series of scenarios.

# Arguments
- `outcomes` : Matrix of outcomes for several scenarios or locations
- `clusters` : Vector of numbers corresponding to clusters
- `opts` : Aviz options
    - `summarize` : plot confidence interval. Defaults to true

# Returns
Figure
"""
function ADRIA.viz.clustered_scenarios(
    outcomes::AbstractMatrix{<:Real},
    clusters::Union{BitVector,Vector{Int64}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    if !haskey(axis_opts, :title)
        axis_opts[:title] = "$(outcome_title(outcomes)) Clusters"
    end

    return ADRIA.viz.scenarios(
        outcomes, clusters; opts=opts, fig_opts=fig_opts, axis_opts=axis_opts
    )
end
function ADRIA.viz.clustered_scenarios!(
    g::Union{GridLayout,GridPosition},
    outcomes::AbstractMatrix{<:Real},
    clusters::Union{BitVector,Vector{Int64}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    return ADRIA.viz.scenarios!(g, outcomes, clusters; axis_opts=axis_opts, opts=opts)
end

"""
    map(rs::Union{Domain,ResultSet}, data::AbstractMatrix, clusters::AbstractVector{Int64}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    map(g, rs, data, clusters; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())

Visualize clustered time series for each location and map.

# Arguments
- `rs` : ResultSet
- `loc_outcomes` : Vector of summary statistics data for each location
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
    loc_outcomes::AbstractVector{<:Real},
    clusters::Union{BitVector,AbstractVector{Int64}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=set_axis_defaults(DEFAULT_OPT_TYPE())
)::Figure
    # Setup fig_opts size before the other default fig_opts
    fig_opts[:size] = get(fig_opts, :size, (800, 700))
    set_figure_defaults(fig_opts)

    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    set_plot_opts!(loc_outcomes, opts, :colorbar_label)

    ADRIA.viz.map!(g, rs, loc_outcomes, clusters; opts=opts, axis_opts=axis_opts)

    return f
end
function ADRIA.viz.map!(
    g::Union{GridLayout,GridPosition},
    rs::Union{Domain,ResultSet},
    loc_outcomes::AbstractVector{<:Real},
    clusters::Union{BitVector,Vector{Int64}};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    # Although this function is called scenario_clusters, here we have locations clusters
    loc_groups::Dict{Symbol,BitVector} = ADRIA.analysis.scenario_clusters(clusters)
    group_colors::Dict{Symbol,Union{Symbol,RGBA{Float32}}} = colors(loc_groups)

    legend_params::Tuple = _cluster_legend_params(loc_outcomes, loc_groups, group_colors)

    _colors::Vector{Union{Symbol,RGBA{Float32}}} = Vector{Union{Symbol,RGBA{Float32}}}(
        undef, length(clusters)
    )
    for (idx, filt) in loc_groups
        _colors[filt] .= group_colors[idx]
    end

    # Highlight is a vector of stroke colors for each location
    opts[:highlight] = get(opts, :highlight, _colors)
    opts[:legend_params] = get(opts, :legend_params, legend_params)

    ADRIA.viz.map!(g, rs, loc_outcomes; opts=opts, axis_opts=axis_opts)

    return g
end

"""
    _cluster_legend_params(data::AbstractVector{<:Real}, scen_groups::Dict{Symbol,BitVector}, group_colors::Dict{Symbol,Union{Symbol,RGBA{Float32}}})::Tuple

Color parameter for current cluster weighted by number of scenarios.

# Arguments
- `data` : Vector of some metric outcome for each location
- `loc_groups` : Dictionary of (group_names => filter), where filter is a BitVector to
select locations that belong to each group
- `group_colors` : Dictionary of (group_names => colors), where colors can be Symbols or
RGBA{Float32}

# Returns
Tuple of legend params to be passed to map! containing legend_entries, legend_labels and
legend_title (in that order).
"""
function _cluster_legend_params(
    data::AbstractVector{<:Real},
    loc_groups::Dict{Symbol,BitVector},
    group_colors::Dict{Symbol,Union{Symbol,RGBA{Float32}}}
)::Tuple
    group_keys = sort(collect(keys(group_colors)))
    colors = [group_colors[key] for key in group_keys]
    legend_entries = [PolyElement(; color=c, strokecolor=:transparent) for c in colors]

    label_means::Vector{Float64} = zeros(length(group_keys))
    for (idx_key, key) in enumerate(group_keys)
        label_means[idx_key] = mean(data[loc_groups[key]])
    end

    legend_labels =
        labels(group_keys) .* ": " .* ADRIA.to_scientific.(label_means, digits=2)
    legend_title = "Cluster mean $(outcome_label(data; label_case=lowercase))"

    return (legend_entries, legend_labels, legend_title)
end
