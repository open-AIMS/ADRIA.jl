using ADRIA.analysis: series_confint
using ADRIA: functional_group_names, human_readable_name

"""
    ADRIA.viz.taxonomy(rs::ResultSet; opts::Dict=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Figure
    ADRIA.viz.taxonomy(scenarios::DataFrame, relative_taxa_cover::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Figure
    ADRIA.viz.taxonomy!(g::Union{GridLayout,GridPosition}, relative_taxa_cover::YAXArray, scen_groups::Dict{Symbol, BitVector}; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Union{GridLayout,GridPosition}

Plot relative cover divided by taxonomy over time.

# Examples
```julia
# Plot relative cover divided by taxonomy or functonal group
ADRIA.viz.taxonomy(rs)

# Plot relative cover divided by taxonomy with custom scenario and relative cover inputs
ADRIA.viz.taxonomy(scenarios, relative_taxa_cover)
```

# Arguments
- `rs` : ADRIA result set
- `scenarios` : Scenario specification
- `relative_taxa_cover` : YAXArray of dimensions [timesteps ⋅ taxa ⋅ scenarios]

# Keyword Arguments
- `opts` : Aviz options
    - `by_RCP` : Split plots by RCP otherwise split by scenario type. Defaults to false.
    - `compare_taxa` : If true, split plots by scenario types, otherwise split by taxonomy. Defaults to true.
    - `show_confints` : Show confidence intervals around series. Defaults to true.
    - `colors` : Colormap for each taxonomy or scenario type. Defaults to Set1_5 for taxa and ADRIA defaults for scenario type.
- `axis_opts` : Additional options to pass to adjust Axis attributes.
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
Figure or GridPosition
"""
function ADRIA.viz.taxonomy(
    rs::ResultSet;
    opts::Dict=Dict{Symbol,Any}(),
    fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    if !haskey(rs.outcomes, :relative_taxa_cover)
        throw(ArgumentError("Unable to found relative_taxa_cover in outcomes. This variable may be passed manually."))
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
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    fig_opts[:size] = get(fig_opts, :size, (1200, 1200))
    f = Figure(; fig_opts...)

    g = f[1, 1] = GridLayout()

    _scenarios = copy(scenarios[1:end .∈ [relative_taxa_cover.scenarios], :])
    scen_groups = if get(opts, :by_RCP, false)
        ADRIA.analysis.scenario_rcps(_scenarios)
    else
        ADRIA.analysis.scenario_types(_scenarios)
    end

    ADRIA.viz.taxonomy!(
        g,
        relative_taxa_cover,
        scen_groups,
        opts=opts,
        axis_opts=axis_opts,
        series_opts=series_opts
    )

    return f
end
function ADRIA.viz.taxonomy!(
    g::Union{GridLayout,GridPosition},
    relative_taxa_cover::YAXArray,
    scen_groups::Dict{Symbol, BitVector};
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Union{GridLayout,GridPosition}
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Year")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Relative Cover")

    show_confints::Bool = get(opts, :show_confints, true)
    compare_taxa::Bool = get(opts, :compare_taxa, true)
    if compare_taxa
        n_groups::Int64 = length(relative_taxa_cover.taxa)
        default_color = Symbol("Set1_" * string(n_groups))
        color = get(opts, :colors, default_color)
        _colors = categorical_colors(color, n_groups)

        taxonomy_series_confint!(
            g,
            relative_taxa_cover,
            scen_groups,
            _colors;
            show_confints=show_confints,
            axis_opts,
            series_opts
        )
    else
        _colors = [COLORS[scen_name] for scen_name in keys(scen_groups)]

        taxonomy_intervened_series_confint!(
            g,
            relative_taxa_cover,
            scen_groups,
            _colors;
            show_confints=show_confints,
            axis_opts,
            series_opts
        )
    end

    return g
end

"""
    taxonomy_series_confint!(g::Union{GridLayout, GridPosition}, relative_taxa_cover::YAXArray, scen_groups::Dict{Symbol}, colors::Union{Vector{Symbol},Vector{RGBA{T}}}; show_confints=true, axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}
    taxonomy_series_confint!( ax::Axis, relative_taxa_cover::YAXArray, colors::Union{Vector{Symbol},Vector{RGBA{T}}}; show_confints=true, show_legend=true, series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}

Plot relative cover divided between coral taxonomy with each grid cell a different type of
scenario.
"""
function taxonomy_series_confint!(
    g::Union{GridLayout, GridPosition},
    relative_taxa_cover::YAXArray,
    scen_groups::Dict{Symbol},
    colors::Union{Vector{Symbol},Vector{RGBA{T}}};
    show_confints=true,
    axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Nothing where {T<:Float32}
    taxa_names = human_readable_name(String.(functional_group_names()), title_case=true)
    series_opts[:labels] = get(series_opts, :labels, taxa_names)

    for (idx, scen_name) in enumerate(keys(scen_groups))
        xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(relative_taxa_cover)))
        xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
        ax = Axis(g[idx, 1]; title=String(scen_name), xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

        taxonomy_series_confint!(
            ax,
            relative_taxa_cover[scenarios=scen_groups[scen_name]],
            colors;
            show_confints=show_confints,
            series_opts=series_opts
        )

    end
    return nothing
end
function taxonomy_series_confint!(
    ax::Axis,
    relative_taxa_cover::YAXArray,
    colors::Union{Vector{Symbol},Array{RGBA{T},1}};
    show_confints=true,
    show_legend=true,
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Nothing where {T<:Float32}
    n_timesteps::Int64 = length(relative_taxa_cover.timesteps)
    n_groups::Int64 = length(relative_taxa_cover.taxa)

    # Plot and calculate confidence intervals
    confints = zeros(n_timesteps, n_groups, 3)
    for (idx, taxa) in enumerate(relative_taxa_cover.taxa)
        confints[:, idx, :] = series_confint(relative_taxa_cover[taxa=At(taxa)])
        show_confints ? band!(
            ax, 1:n_timesteps, confints[:, idx, 1], confints[:, idx, 3];
            color=(colors[idx], 0.4)
        ) : nothing
    end

    # Plot series
    series!(ax, 1:n_timesteps, confints[:, :, 2]'; solid_color=colors, series_opts...)
    show_legend ? axislegend(ax) : nothing

    return nothing
end

"""
    taxonomy_intervened_series_confint!( g::Union{GridLayout, GridPosition}, relative_taxa_cover::YAXArray, scen_groups::Dict{Symbol,BitVector}, colors::Union{Vector{Symbol},Vector{RGBA{T}}}; show_confints::Bool=true, axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}
    taxonomy_intervened_series_confint!( ax::Axis, relative_taxa_cover::YAXArray, colors::Union{Vector{Symbol},Vector{RGBA{T}}}, scen_groups::Dict{Symbol,BitVector}; show_confints::Bool=true, show_legend::Bool=true, series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Nothing where {T<:Float32}

Plot relative cover, comparing taxa cover between different scenario groups. Each plot is a
different coral taxonomy or functional group.
"""
function taxonomy_intervened_series_confint!(
    g::Union{GridLayout, GridPosition},
    relative_taxa_cover::YAXArray,
    scen_groups::Dict{Symbol,BitVector},
    colors::Union{Vector{Symbol},Vector{RGBA{T}}};
    show_confints::Bool=true,
    axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Nothing where {T<:Float32}
    taxa_names = human_readable_name(String.(functional_group_names()), title_case=true)

    scenario_names::Vector{Symbol} = collect(keys(scen_groups))
    series_opts[:labels] = get(series_opts, :labels, String.(scenario_names))

    for (idx, taxa_name) in enumerate(taxa_names)
        xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(relative_taxa_cover)))
        xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
        ax = Axis(g[idx, 1]; title=taxa_name, xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

        taxonomy_intervened_series_confint!(
            ax,
            relative_taxa_cover[taxa=idx],
            colors,
            scen_groups;
            show_confints=show_confints,
            series_opts=series_opts
        )

    end
    return nothing
end
function taxonomy_intervened_series_confint!(
    ax::Axis,
    relative_taxa_cover::YAXArray,
    colors::Union{Vector{Symbol},Vector{RGBA{T}}},
    scen_groups::Dict{Symbol,BitVector};
    show_confints::Bool=true,
    show_legend::Bool=true,
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Nothing where {T<:Float32}
    scenario_names::Vector{Symbol} = collect(keys(scen_groups))

    n_timesteps::Int64 = length(relative_taxa_cover.timesteps)
    n_scen_types::Int64 = length(scenario_names)

    confints = zeros(n_timesteps, n_scen_types, 3)
    for (idx, scen) in enumerate(scenario_names)
        confints[:, idx, :] = series_confint(
            relative_taxa_cover[scenarios=scen_groups[scen]]
        )
        show_confints ? band!(
            ax, 1:n_timesteps, confints[:, idx, 1], confints[:, idx, 3];
            color=(colors[idx], 0.4)
        ) : nothing
    end

    # Plot series
    series!(ax, 1:n_timesteps, confints[:, :, 2]'; solid_color=colors, series_opts...)
    show_legend ? axislegend(ax) : nothing

    return nothing
end
