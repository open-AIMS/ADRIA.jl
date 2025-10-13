using JuliennedArrays: Slices
using ADRIA.analysis: series_confint
using ADRIA: axes_names, RMEResultSet

"""
    ADRIA.viz.scenarios(rs::ADRIA.ResultSet, outcomes::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(by_RCP => false), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())
    ADRIA.viz.scenarios(scenarios::DataFrame, outcomes::YAXArray; opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(:by_RCP => false), fig_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), axis_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}(), series_opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}())::Figure
    ADRIA.viz.scenarios!(g::Union{GridLayout,GridPosition}, scenarios::DataFrame, outcomes::YAXArray; opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

Plot scenario outcomes over time.

# Examples
```julia
scens = ADRIA.sample(dom, 64)

# ... run scenarios ...

s_tac = ADRIA.metrics.scenario_total_cover(rs)

# Plot scenario outcomes
ADRIA.viz.scenarios(rs.inputs, s_tac)

# Plot outcomes of scenarios where SRM < 1.0
ADRIA.viz.scenarios(rs.inputs, s_tac[:, scens.SRM .< 1.0])
```

# Arguments
- `scenarios` : Scenario specification
- `outcomes` : Results of scenario metric
- `opts` : Aviz options
    - `by_RCP` : color by RCP otherwise color by scenario type. Defaults to false.
    - `legend` : show legend. Defaults to true.
    - `summarize` : plot confidence interval. Defaults to true
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
Figure or GridPosition
"""
function ADRIA.viz.scenarios(
    rs::ResultSet,
    outcomes::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(:by_RCP => false),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(:size => (800, 300)),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    return ADRIA.viz.scenarios(
        rs.inputs,
        outcomes;
        opts=opts,
        fig_opts=fig_opts,
        axis_opts=axis_opts,
        series_opts=series_opts
    )
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    outcomes::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(:by_RCP => false),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    opts[:histogram] = get(opts, :histogram, false)

    return ADRIA.viz.scenarios!(
        g,
        rs.inputs,
        outcomes;
        opts=opts,
        axis_opts=axis_opts,
        series_opts=series_opts
    )
end
function ADRIA.viz.scenarios(
    rs::RMEResultSet,
    outcomes::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(:by_RCP => false),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(outcomes)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)
    scen_groups = rs.scenario_groups

    ADRIA.viz.scenarios!(
        g,
        ax,
        outcomes,
        scen_groups;
        opts=opts,
        axis_opts=axis_opts,
        series_opts=series_opts
    )
    return f
end
function ADRIA.viz.scenarios(
    scenarios::DataFrame,
    outcomes::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(:by_RCP => false),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(:size => (800, 300)),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.scenarios!(
        g, scenarios, outcomes; opts=opts, axis_opts=axis_opts, series_opts=series_opts
    )

    return f
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    scenarios::DataFrame,
    outcomes::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(outcomes)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)

    if !haskey(axis_opts, :title)
        axis_opts[:title] = outcome_title(outcomes)
    end

    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    _scenarios = copy(scenarios[1:end .∈ [outcomes.scenarios], :])
    scen_groups = if get(opts, :by_RCP, false)
        ADRIA.analysis.scenario_rcps(_scenarios)
    else
        ADRIA.analysis.scenario_types(_scenarios)
    end

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
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    ax::Axis,
    outcomes::YAXArray,
    scen_groups::Dict{Symbol,BitVector};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    by_RCP::Bool = get(opts, :by_RCP, false)
    sort_by::Symbol = get(opts, :sort_by, :default)
    default_names::Vector{Symbol} = get(opts, :legend_labels, [])
    group_names::Vector{Symbol} = _sort_keys(
        scen_groups;
        by_RCP=by_RCP, by=sort_by, outcomes=outcomes, default_names=default_names
    )

    if get(opts, :summarize, true)
        scenarios_confint!(ax, outcomes, scen_groups, group_names)
    else
        scenarios_series!(ax, outcomes, scen_groups, group_names; series_opts=series_opts)
    end

    get(opts, :histogram, true) && scenarios_hist(g, outcomes, scen_groups, group_names)

    if get(opts, :legend, true)
        legend_position = get(opts, :histogram, true) ? (1, 3) : (1, 2)
        legend_labels = get(opts, :legend_labels, group_names)
        _render_legend(g, scen_groups, legend_position, legend_labels)
    end

    ax.xlabel = "Year"
    ax.ylabel = outcome_label(outcomes)
    return g
end

# TODO Add support for YAXArrays (`Makie.lines!` expects an Array). For now this works.
"""
    scenario_stratified(data::AbstractArray{Float64, 4})

Plot a single scenario, averaged across all locations, stratified by size class and
functional group.

# Arguments
- `data` : Either a 4D array comprising (timesteps, size_classes, functional_groups and
locations) or a 3D array comprising (timesteps, size_groups, locations).
- `axis_opts` : Opts to be passed to Makie.Axis.

# Example
```
# Run a single scenario directly using `run_model` directly
m_rs = ADRIA.run_model(dom, samples[1, :])

n_timesteps, n_size_groups, n_locs = size(m_rs.raw)
n_groups = dom.coral_growth.n_groups
n_sizes = dom.coral_growth.n_sizes

# Plot scenario cover stratified by size class and functional group
ADRIA.viz.scenario_stratified(reshape(m_rs.raw, n_timesteps, n_sizes, n_groups, n_locations))

# Plot scenario bleaching mortality stratified by size class and functional group
ADRIA.viz.scenario_stratified(
    permutedims(m_rs.bleaching_mortality, (1,3,2,4));
    axis_opts=Dict{Symbol,Any}(:limits=(nothing, (0.0, 1.0)))
)

# Plot scenario cover
```
"""
function ADRIA.viz.scenario_stratified(
    data::Array{Float64,4};
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    n_timesteps, n_sizes, n_groups, n_locations = size(data)
    f = Figure(; size=(1000, 1200), fig_opts...)
    xdata = 1:n_timesteps

    # TODO plot confints
    mean_data = dropdims(mean(data; dims=4); dims=4)

    limits = pop!(
        axis_opts,
        :limits,
        (nothing, (nothing, ceil(maximum(mean_data); digits=2)))
    )
    ylabel = pop!(axis_opts, :ylabel, "")

    # This isn't technically an axis parameter because this title is going to be applied on
    # as a Label to the whole figure, but left it like this for simplicity
    title = pop!(axis_opts, :title, "")

    for fg in 1:n_groups
        for sc in 1:n_sizes
            ax = Axis(
                f[sc, fg];
                title="Functional group = $fg\nSize class = $sc",
                ylabel=ylabel,
                xlabel="Year",
                limits=limits,
                axis_opts...
            )
            ydata = mean_data[:, sc, fg]
            lines!(ax, xdata, ydata)
        end
    end
    Label(f[0, :], title; fontsize=24)
    return f
end

function _confints(
    outcomes::YAXArray,
    scen_groups::Dict{Symbol,BitVector},
    group_names::Vector{Symbol}
)::Array{Float64}
    n_timesteps::Int64 = size(outcomes, 1)
    n_groups::Int64 = length(group_names)

    # Compute confints
    confints::Array{Float64} = zeros(n_timesteps, n_groups, 3)
    agg_dim = symdiff(axes_names(outcomes), [:timesteps])[1]
    for (idx, group) in enumerate(group_names)
        confints[:, idx, :] = series_confint(
            outcomes.data[:, scen_groups[group]]; agg_dim=agg_dim
        )
    end

    return confints
end

function scenarios_confint!(
    ax::Axis,
    confints::AbstractArray,
    ordered_groups::Vector{Symbol},
    _colors::Dict{Symbol,T};
    x_vals::Union{Vector{Int64},Vector{Float64}}=collect(1:size(confints, 1))
)::Nothing where {T<:Union{RGBA{Float32},String,Symbol}}
    for idx in eachindex(ordered_groups)
        band_color = (_colors[ordered_groups[idx]], 0.4)
        y_lower, y_upper = confints[:, idx, 1], confints[:, idx, 3]
        band!(ax, x_vals, y_lower, y_upper; color=band_color)
    end

    series_colors = [_colors[group] for group in ordered_groups]
    series!(ax, x_vals, confints[:, :, 2]'; solid_color=series_colors)

    return nothing
end
function scenarios_confint!(
    ax::Axis,
    outcomes::YAXArray,
    scen_groups::Dict{Symbol,BitVector},
    group_names::Vector{Symbol}
)::Nothing
    _colors::Dict{Symbol,Union{Symbol,RGBA{Float32},String}} = colors(scen_groups)
    confints = _confints(outcomes, scen_groups, group_names)
    return scenarios_confint!(
        ax,
        confints,
        group_names,
        _colors;
        x_vals=collect(1:size(confints, 1))
    )
end

function scenarios_series!(
    ax::Axis,
    outcomes::YAXArray,
    scen_groups::Dict{Symbol,BitVector},
    group_names::Vector{Symbol};
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    x_vals::T=collect(1:size(outcomes, 1))
)::Nothing where {T<:Union{Vector{Int64},Vector{Float64}}}
    _colors::Dict{Symbol,Union{Symbol,RGBA{Float32}}} = colors(scen_groups)
    _alphas::Dict{Symbol,Float64} = alphas(scen_groups)

    for group in group_names
        color = (_colors[group], _alphas[group])
        scens = outcomes[:, scen_groups[group]]'
        series!(ax, x_vals, scens.data; solid_color=color, series_opts...)
    end

    return nothing
end

function scenarios_hist(
    g::Union{GridLayout,GridPosition},
    outcomes::YAXArray,
    scen_groups::Dict{<:Any,BitVector},
    group_names::Vector{Symbol}
)::Nothing
    scen_dist = dropdims(mean(outcomes; dims=:timesteps); dims=:timesteps)
    ax_hist = Axis(g[1, 2]; width=100)
    _colors = colors(scen_groups)

    for group in group_names
        color = (_colors[group], 0.7)
        dist = scen_dist[scen_groups[group]]
        hist!(ax_hist, dist.data; direction=:x, color=color, bins=30, normalization=:pdf)
    end

    hidedecorations!(ax_hist)
    hidespines!(ax_hist)
    ylims!(
        ax_hist,
        minimum(scen_dist) - quantile(scen_dist, 0.05),
        maximum(scen_dist) + quantile(scen_dist, 0.05)
    )

    return nothing
end

function _render_legend(
    g::Union{GridLayout,GridPosition},
    scen_groups::Dict{Symbol,BitVector},
    legend_position::Tuple{Int64,Int64},
    legend_labels::Vector{Symbol}
)::Nothing
    _colors = colors(scen_groups)
    line_els::Vector{LineElement} = [LineElement(; color=_colors[n]) for n in legend_labels]

    Legend(g[legend_position...], line_els, labels(legend_labels); framevisible=false)

    return nothing
end

"""
    _sort_keys(scenario_types::Dict{Symbol, BitVector}, outcomes::YAXArray)::Vector{Symbol}

Sort scenario types (counterfactual, unguided and guided).

# Arguments
- `outcomes` : Results of scenario metric
- `scenario_types` : NamedTuple of BitVectors to filter scenarios for each scenario type of:
    - :guided
    - :unguided
    - :counterfactual
- `by` : Sort criteria. Can be either `:default`, `:variance` or `:size`
"""
function _sort_keys(
    scenario_types::Dict{Symbol,BitVector};
    by_RCP::Bool,
    by::Symbol=:default,
    outcomes::AbstractArray=[],
    default_names::Vector{Symbol}=[]
)::Vector{Symbol}
    scen_types::Vector{Symbol} = collect(keys(scenario_types))
    if by == :default
        by_RCP && return sort(collect(keys(scenario_types)))
        !isempty(default_names) && return default_names

        default_keys = [:counterfactual, :unguided, :guided]
        return default_keys[default_keys .∈ [scen_types]]
    elseif by == :variance
        msg = "When sorting by variance, optional parameter `outcomes` must be provided"
        isempty(outcomes) && throw(ArgumentError(msg))
        return sort(
            scen_types;
            by=type -> sum(var(outcomes[:, scenario_types[type]]; dims=2)),
            rev=true
        )
    elseif by == :size
        msg = "When sorting by size, optional parameter `outcomes` must be provided"
        isempty(outcomes) && throw(ArgumentError(msg))
        return sort(
            scen_types; by=type -> size(outcomes[:, scenario_types[type]], 2), rev=true
        )
    else
        throw(
            ArgumentError(
                "Invalid 'by' option. Must be one of: [:default, :variance, :size]"
            )
        )
    end
end
