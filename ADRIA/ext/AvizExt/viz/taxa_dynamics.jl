using ADRIA.analysis: series_confint
using ADRIA: functional_group_names, human_readable_name

"""
    ADRIA.viz.taxonomy(rs::ResultSet; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Figure
    ADRIA.viz.taxonomy(scenarios::DataFrame, relative_taxa_cover::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Figure
    ADRIA.viz.taxonomy!(g::Union{GridLayout,GridPosition}, relative_taxa_cover::YAXArray, scen_groups::Dict{Symbol, BitVector}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Union{GridLayout,GridPosition}

Plot relative cover divided by taxonomy over time.

# Examples
```julia
# Plot relative cover divided by taxonomy or functional group
ADRIA.viz.taxonomy(rs)

# Plot relative cover divided by taxonomy with custom scenario and relative cover inputs
ADRIA.viz.taxonomy(scenarios, relative_taxa_cover)
```

# Arguments
- `rs` : ADRIA result set
- `scenarios` : Scenario specification
- `relative_taxa_cover` : YAXArray of dimensions [timesteps ⋅ groups ⋅ scenarios]
- `opts` : Aviz options
    - `by_RCP` : Split plots by RCP otherwise split by scenario type. Defaults to false.
    - `by_functional_groups` : If true, split plots by scenario types, otherwise split by taxonomy. Defaults to true.
    - `show_confints` : Show confidence intervals around series. Defaults to true.
    - `colors` : Colormap for each taxonomy or scenario type. Defaults to Set1_5 for groups and ADRIA defaults for scenario type.
- `axis_opts` : Additional options to pass to adjust Axis attributes.
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
Figure or GridPosition
"""
function ADRIA.viz.taxonomy(
    rs::ResultSet;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    if !haskey(rs.outcomes, :relative_taxa_cover)
        throw(
            ArgumentError(
                "Unable to found relative_taxa_cover in outcomes. This variable may be passed manually."
            )
        )
    end
    return ADRIA.viz.taxonomy(
        rs.inputs,
        rs.outcomes[:relative_taxa_cover];
        opts=opts,
        fig_opts=fig_opts,
        axis_opts=axis_opts,
        series_opts=series_opts
    )
end
function ADRIA.viz.taxonomy(
    scenarios::DataFrame,
    relative_taxa_cover::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    _scenarios = copy(scenarios[1:end .∈ [relative_taxa_cover.scenarios], :])
    scen_groups = if get(opts, :by_RCP, false)
        ADRIA.analysis.scenario_rcps(_scenarios)
    else
        ADRIA.analysis.scenario_types(_scenarios)
    end

    by_functional_groups::Bool = get(opts, :by_functional_groups, true)

    fig_height::Int64 = _fig_height_taxonomy(
        by_functional_groups, scen_groups, relative_taxa_cover
    )

    fig_opts[:size] = get(fig_opts, :size, (1200, fig_height))
    f = Figure(; fig_opts...)

    g = f[1, 1] = GridLayout()

    ADRIA.viz.taxonomy!(
        g,
        relative_taxa_cover,
        scen_groups;
        opts=opts,
        axis_opts=axis_opts,
        series_opts=series_opts
    )

    return f
end
function ADRIA.viz.taxonomy!(
    g::Union{GridLayout,GridPosition},
    relative_taxa_cover::YAXArray,
    scen_groups::Dict{Symbol,BitVector};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Year")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Relative Cover")

    show_confints::Bool = get(opts, :show_confints, true)
    by_functional_groups::Bool = get(opts, :by_functional_groups, true)
    if by_functional_groups
        # Create colors
        n_functional_groups::Int64 = length(relative_taxa_cover.groups)
        default_color = Symbol("Set1_" * string(n_functional_groups))
        color = get(opts, :colors, default_color)
        _colors = categorical_colors(color, n_functional_groups)

        # Plot results
        taxonomy_by_intervention!(
            g,
            relative_taxa_cover,
            scen_groups,
            _colors;
            show_confints=show_confints,
            axis_opts,
            series_opts
        )
    else
        # Use default ADRIA colors for scenario type if use not specified
        n_scenario_groups::Int64 = length(keys(scen_groups))
        color = get(opts, :colors, nothing)
        _colors =
            if isnothing(color)
                [
                    COLORS[scen_name] for scen_name in keys(scen_groups)
                ]
            else
                categorical_colors(color, n_scenario_groups)
            end

        # Plot results
        intervention_by_taxonomy!(
            g,
            relative_taxa_cover,
            scen_groups,
            _colors;
            show_confints=show_confints,
            axis_opts,
            series_opts
        )
    end

    Label(
        g[1, :, Top()], "Taxa Dynamics"; padding=(0, 0, 30, 0), font=:bold, valign=:bottom
    )

    return g
end

"""
    taxonomy_by_intervention!(g::Union{GridLayout, GridPosition}, relative_taxa_cover::YAXArray, scen_groups::Dict{Symbol}, colors::Union{Vector{Symbol},Vector{RGBA{T}}}; show_confints=true, axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}
    taxonomy_by_intervention!(ax::Axis, relative_taxa_cover::YAXArray, colors::Union{Vector{Symbol},Vector{RGBA{T}}}; show_confints=true, show_legend=true, series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}

Create a plot for each scenario group, displaying the relative coral cover split between
functional groups.
"""
function taxonomy_by_intervention!(
    g::Union{GridLayout,GridPosition},
    relative_taxa_cover::YAXArray,
    scen_groups::Dict{Symbol},
    colors::Union{Vector{Symbol},Vector{RGBA{T}}};
    show_confints=true,
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Nothing where {T<:Float32}
    # Get taxonomy names for legend
    taxa_names = human_readable_name(functional_group_names(); title_case=true)
    series_opts[:labels] = get(series_opts, :labels, taxa_names)

    # Get default axis options
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(relative_taxa_cover)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    for (idx, scen_name) in enumerate(keys(scen_groups))
        ax = Axis(
            g[idx, 1];
            title=titlecase(String(scen_name)),
            xticks=xtick_vals,
            xticklabelrotation=xtick_rot,
            axis_opts...
        )

        taxonomy_by_intervention!(
            ax,
            relative_taxa_cover[scenarios=scen_groups[scen_name]],
            colors;
            show_confints=show_confints,
            series_opts=series_opts
        )
    end
    return nothing
end
function taxonomy_by_intervention!(
    ax::Axis,
    relative_taxa_cover::YAXArray,
    colors::Union{Vector{Symbol},Array{RGBA{T},1}};
    show_confints=true,
    show_legend=true,
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Nothing where {T<:Float32}
    n_timesteps::Int64 = length(relative_taxa_cover.timesteps)
    functional_groups::Vector{Int64} = ADRIA.axis_labels(relative_taxa_cover, :groups)
    n_functional_groups::Int64 = length(functional_groups)

    # Plot and calculate confidence intervals
    confints = zeros(n_timesteps, n_functional_groups, 3)
    for (idx, group) in enumerate(functional_groups)
        confints[:, idx, :] = series_confint(relative_taxa_cover[groups=At(group)])
        if show_confints
            band!(
                ax, 1:n_timesteps, confints[:, idx, 1], confints[:, idx, 3];
                color=(colors[idx], 0.4)
            )
        else
            nothing
        end
    end

    # Plot series
    series!(ax, 1:n_timesteps, confints[:, :, 2]'; solid_color=colors, series_opts...)
    show_legend ? axislegend(ax; position=:lt) : nothing

    return nothing
end

"""
    intervention_by_taxonomy!( g::Union{GridLayout, GridPosition}, relative_taxa_cover::YAXArray, scen_groups::Dict{Symbol,BitVector}, colors::Union{Vector{Symbol},Vector{RGBA{T}}}; show_confints::Bool=true, axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}
    intervention_by_taxonomy!( ax::Axis, relative_taxa_cover::YAXArray, colors::Union{Vector{Symbol},Vector{RGBA{T}}}, scen_groups::Dict{Symbol,BitVector}; show_confints::Bool=true, show_legend::Bool=true, series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}

Plot relative cover, comparing taxa cover between different scenario groups. Each plot is a
different coral taxonomy or functional group.
"""
function intervention_by_taxonomy!(
    g::Union{GridLayout,GridPosition},
    relative_taxa_cover::YAXArray,
    scen_groups::Dict{Symbol,BitVector},
    colors::Union{Vector{Symbol},Vector{RGBA{T}}};
    show_confints::Bool=true,
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Nothing where {T<:Float32}
    taxa_names = human_readable_name(functional_group_names(); title_case=true)

    scenario_group_names::Vector{Symbol} = collect(keys(scen_groups))
    series_opts[:labels] = get(
        series_opts, :labels, titlecase.(String.(scenario_group_names))
    )

    for (idx, taxa_name) in enumerate(taxa_names)
        xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(relative_taxa_cover)))
        xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
        ax = Axis(
            g[idx, 1];
            title=taxa_name,
            xticks=xtick_vals,
            xticklabelrotation=xtick_rot,
            axis_opts...
        )

        intervention_by_taxonomy!(
            ax,
            relative_taxa_cover[groups=idx],
            colors,
            scen_groups;
            show_confints=show_confints,
            series_opts=series_opts
        )
    end
    return nothing
end
function intervention_by_taxonomy!(
    ax::Axis,
    relative_taxa_cover::YAXArray,
    colors::Union{Vector{Symbol},Vector{RGBA{T}}},
    scen_groups::Dict{Symbol,BitVector};
    show_confints::Bool=true,
    show_legend::Bool=true,
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Nothing where {T<:Float32}
    scenario_group_names::Vector{Symbol} = collect(keys(scen_groups))

    n_timesteps::Int64 = length(relative_taxa_cover.timesteps)
    n_scenario_groups::Int64 = length(scenario_group_names)

    confints = zeros(n_timesteps, n_scenario_groups, 3)
    for (idx, scen) in enumerate(scenario_group_names)
        confints[:, idx, :] = series_confint(
            relative_taxa_cover[scenarios=scen_groups[scen]]
        )
        if show_confints
            band!(
                ax, 1:n_timesteps, confints[:, idx, 1], confints[:, idx, 3];
                color=(colors[idx], 0.4)
            )
        else
            nothing
        end
    end

    # Plot series
    series!(ax, 1:n_timesteps, confints[:, :, 2]'; solid_color=colors, series_opts...)
    show_legend ? axislegend(ax; position=:lt) : nothing

    return nothing
end

"""
    _fig_height_taxonomy(
        by_functional_groups::Bool,
        scen_groups::Dict{Symbol,BitVector},
        relative_taxa_cover::YAXArray
    )::Int64

Auto-adjust taxonomy figure height based on number of plots.
"""
function _fig_height_taxonomy(
    by_functional_groups::Bool,
    scen_groups::Dict{Symbol,BitVector},
    relative_taxa_cover::YAXArray
)::Int64
    fig_base_height::Int64 = 1200
    fig_height = if by_functional_groups
        n_scen_groups = length(scen_groups)
        base_n_scen_groups = 3
        fig_base_height * (n_scen_groups / base_n_scen_groups)
    else
        n_functional_groups = length(relative_taxa_cover.groups)
        base_n_functional_groups = 5
        fig_base_height * (n_functional_groups / base_n_functional_groups)
    end
    return Int64(ceil(fig_height))
end
