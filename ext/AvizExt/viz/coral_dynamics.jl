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
    - `show_confints` : Show confidence intervals around series. Defaults to true.
    - `colors` : Color for each taxonomy. Defaults to categorical colors.
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
    for (idx, scen_name) in enumerate(keys(scen_groups))
        xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(relative_taxa_cover)))
        xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
        ax = Axis(g[idx, 1]; title=String(scen_name), xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

        ADRIA.viz.taxonomy!(
            g,
            ax,
            relative_taxa_cover[scenarios=scen_groups[scen_name]],
            opts=opts,
            series_opts=series_opts
        )

    end

    return g
end
function ADRIA.viz.taxonomy!(
    g::Union{GridLayout,GridPosition},
    ax::Axis,
    relative_taxa_cover::YAXArray;
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(),
    series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Union{GridLayout,GridPosition}
    n_groups::Int64 = length(relative_taxa_cover.taxa)
    n_timesteps::Int64 = length(relative_taxa_cover.timesteps)

    # Allow user to specify taxonomy colors
    default_color = Symbol("Set1_" * string(n_groups))
    color = get(opts, :colors, default_color)
    _colors = categorical_colors(color, n_groups)

    # Plot confidence intervals
    show_confints = get(opts, :show_confints, true)
    if show_confints
        confints = zeros(n_timesteps, n_groups, 3)
        for (idx, taxa) in enumerate(relative_taxa_cover.taxa)
            confints[:, idx, :] = series_confint(relative_taxa_cover[taxa=At(taxa)])
            band!(
                ax, 1:n_timesteps, confints[:, idx, 1], confints[:, idx, 3];
                color=(_colors[idx], 0.4)
            )
        end
    end

    # Plot lines
    taxa_names = human_readable_name(String.(functional_group_names()), title_case=true)
    series_opts[:labels] = get(series_opts, :labels, taxa_names)

    series!(ax, 1:n_timesteps, confints[:, :, 2]'; solid_color=_colors, series_opts...)
    axislegend(ax)
    return g
end
