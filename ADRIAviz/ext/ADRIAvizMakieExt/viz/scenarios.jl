using YAXArrays
using ADRIA.analysis: series_confint
using ADRIA: axes_names, RMEResultSet, AnnotatedOutcomes
using OrderedCollections

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
    outcomes::YAXArray,
    scen_groups::Dict{Symbol,BitVector};
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(:size => (800, 300)),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(outcomes)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    ADRIA.viz.scenarios!(
        g, ax, outcomes, scen_groups;
        opts=opts, axis_opts=axis_opts, series_opts=series_opts
    )
    return f
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
        _scenario_rcps(_scenarios)
    else
        _scenario_types(_scenarios)
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
        _render_legend(g[legend_position...], scen_groups, legend_labels)
    end

    ax.xlabel = get(axis_opts, :xlabel, "Year")
    ax.ylabel = get(axis_opts, :ylabel, outcome_label(outcomes))
    return g
end

function ADRIA.viz.scenarios_legend!(
    g::GridPosition, rs::ResultSet, outcomes::YAXArray; opts::OPT_TYPE=Dict{Symbol,Any}(),
    legend_opts::OPT_TYPE=Dict{Symbol,Any}()
)
    _scenarios = rs.inputs #copy(@view(scenarios[1:end .∈ [outcomes.scenarios], :]))
    by_RCP::Bool = get(opts, :by_RCP, false)
    scen_groups = if by_RCP
        _scenario_rcps(_scenarios)
    else
        _scenario_types(_scenarios)
    end
    return ADRIA.viz.scenarios_legend!(
        g, scen_groups, outcomes; opts=opts, legend_opts=legend_opts
    )
end
function ADRIA.viz.scenarios_legend!(
    g::GridPosition, scen_groups::Dict, outcomes::YAXArray;
    opts::OPT_TYPE=Dict{Symbol,Any}(),
    legend_opts::OPT_TYPE=Dict{Symbol,Any}()
)
    by_RCP::Bool = get(opts, :by_RCP, false)
    sort_by::Symbol = get(opts, :sort_by, :default)
    default_names::Vector{Symbol} = get(opts, :legend_labels, [])

    group_names::Vector{Symbol} = _sort_keys(
        scen_groups;
        by_RCP=by_RCP, by=sort_by, outcomes=outcomes, default_names=default_names
    )

    legend_labels = get(opts, :legend_labels, group_names)
    return _render_legend(
        g, scen_groups, legend_labels; legend_opts=legend_opts
    )
end

# TODO Add support for YAXArrays (`Makie.lines!` expects an Array). For now this works.
"""
    scenario_by_group_and_size(
        data::Array{Float64,4};
        agg_fn::Function=median,
        fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
        axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
    )
    scenario_by_group_and_size(
        data::Array{Float64,3};
        fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
        axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
    )

Plot a single scenario, averaged across all locations, stratified by size class and
functional group. When a 4 dimensional Array is passed, the last dimension is assumed to be
`locations` and the median data is aggregated across this dimension.

# Arguments
- `data` : A 4D array comprising (`timesteps`, `functional_groups`, `size_classes` and
`locations`) or a 3D array comprising (`timesteps`, `functional_groups`, `size_classes`).
- `agg_fn` : Function used to aggregate `location` dimension data. Defaults to `median`.
- `axis_opts` : Options to be passed to `Makie.Axis()`.

# Example
```
# Run a single scenario using `run_model` directly
m_rs = ADRIA.run_model(dom, scens[1, :])

# Plot scenario cover by size class and functional group
ADRIA.viz.scenario_by_group_and_size(m_rs.raw)

# Plot scenario bleaching mortality by size class and functional group
ADRIA.viz.scenario_by_group_and_size(
    m_rs.bleaching_mortality;
    axis_opts=Dict{Symbol,Any}(:limits=>(nothing, (0.0, 1.0)))
)
```
"""
function ADRIA.viz.scenario_by_group_and_size(
    data::Array{Float64,4};
    agg_fn::Function=median,
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    # TODO plot confints
    agg_data = dropdims(agg_fn(data; dims=4); dims=4)

    return ADRIA.viz.scenario_by_group_and_size(
        agg_data; fig_opts=fig_opts, axis_opts=axis_opts
    )
end
function ADRIA.viz.scenario_by_group_and_size(
    data::Array{Float64,3};
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    n_timesteps, n_groups, n_sizes = size(data)
    fig_size = pop!(fig_opts, :size, (1000, 1200))
    fig = Figure(; size=fig_size, fig_opts...)
    xdata = 1:n_timesteps

    limits = pop!(
        axis_opts,
        :limits,
        (nothing, (nothing, ceil(maximum(data); digits=2)))
    )
    ylabel = pop!(axis_opts, :ylabel, "")

    # This isn't technically an axis parameter because this title is going to be applied on
    # as a Label to the whole figure, but left it like this for simplicity
    title = pop!(fig_opts, :title, "")

    for fg in 1:n_groups
        for sc in 1:n_sizes
            ax = Axis(
                fig[sc, fg];
                title="Functional group = $fg\nSize class = $sc",
                ylabel=ylabel,
                xlabel="Year",
                limits=limits,
                axis_opts...
            )
            ydata = data[:, fg, sc]
            lines!(ax, xdata, ydata)
        end
    end

    Label(fig[0, :], title; fontsize=24)

    return fig
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
    outcomes_arr = Array(outcomes)
    for (idx, group) in enumerate(group_names)
        confints[:, idx, :] = series_confint(
            outcomes_arr[:, scen_groups[group]]; agg_dim=agg_dim
        )
    end

    return confints
end

function scenarios_confint!(
    ax::Axis,
    confints::AbstractArray,
    ordered_groups::Vector{Symbol},
    _colors::Dict{Symbol,COLOR_TYPE};
    x_vals::Union{Vector{Int64},Vector{Float64}}=collect(1:size(confints, 1))
)::Nothing
    for idx in eachindex(ordered_groups)
        band_color = (_colors[ordered_groups[idx]], 0.4)
        y_lower, y_upper = confints[:, idx, 1], confints[:, idx, 3]
        band!(ax, x_vals, y_lower, y_upper; color=band_color)
    end

    series_colors = [_colors[group] for group in ordered_groups]
    series!(ax, x_vals, Matrix(confints[:, :, 2]'); solid_color=series_colors)

    return nothing
end
function scenarios_confint!(
    ax::Axis,
    outcomes::YAXArray,
    scen_groups::Dict{Symbol,BitVector},
    group_names::Vector{Symbol}
)::Nothing
    _colors::Dict{Symbol,COLOR_TYPE} = colors(scen_groups)
    x_vals = collect(1:size(outcomes, 1))
    outcomes_arr = Array(outcomes)

    single_groups = filter(g -> sum(scen_groups[g]) == 1, group_names)
    for group in single_groups
        idx = findfirst(scen_groups[group])
        lines!(ax, x_vals, outcomes_arr[:, idx]; color=_colors[group])
    end

    multi_groups = filter(g -> sum(scen_groups[g]) > 1, group_names)
    isempty(multi_groups) && return nothing

    confints = _confints(outcomes, scen_groups, multi_groups)
    return scenarios_confint!(ax, confints, multi_groups, _colors; x_vals=x_vals)
end

function scenarios_series!(
    ax::Axis,
    outcomes::YAXArray,
    scen_groups::Dict{Symbol,BitVector},
    group_names::Vector{Symbol};
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    x_vals::T=collect(1:size(outcomes, 1))
)::Nothing where {T<:Union{Vector{Int64},Vector{Float64}}}
    _colors::Dict{Symbol,COLOR_TYPE} = colors(scen_groups)
    _alphas::Dict{Symbol,Float64} = alphas(scen_groups)

    outcomes_arr = Array(outcomes)
    for group in group_names
        color = (_colors[group], _alphas[group])
        scens = Matrix(outcomes_arr[:, scen_groups[group]]')
        series!(ax, x_vals, scens; solid_color=color, series_opts...)
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
    g::Union{GridLayout,GridPosition,GridSubposition},
    scen_groups::Dict{Symbol,BitVector},
    legend_labels::Vector{Symbol};
    legend_opts::Dict{Symbol,Any}=Dict{Symbol,Any}()
)::Nothing
    _colors = colors(scen_groups)
    line_els::Vector{LineElement} = [LineElement(; color=_colors[n]) for n in legend_labels]

    title = pop!(legend_opts, :title, "Scenarios")
    Legend(g, line_els, labels(legend_labels), title; framevisible=false, legend_opts...)

    return nothing
end

"""
    ADRIA.viz.scenarios(outcomes::YAXArray; opts=Dict(), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())::Figure
    ADRIA.viz.scenarios!(g::Union{GridLayout,GridPosition}, outcomes::YAXArray; opts=Dict(), axis_opts=Dict(), series_opts=Dict())

Plot scenario outcomes over time directly from a `YAXArray`, without requiring a result set
or scenario `DataFrame`. All scenarios are treated as a single group.

# Examples
```julia
# outcomes is a (timesteps × scenarios) YAXArray, e.g. from an external source
ADRIA.viz.scenarios(outcomes)

# Plot individual lines instead of the confidence interval band
ADRIA.viz.scenarios(outcomes; opts=Dict(:summarize => false))
```

# Arguments
- `outcomes` : Results of scenario metric, with dimensions `(timesteps, scenarios)`
- `opts` : Aviz options
    - `summarize` : plot confidence interval band. Defaults to true.
    - `legend` : show legend. Defaults to false.
    - `histogram` : show marginal histogram. Defaults to false.
- `fig_opts` : Additional options to pass to `Figure`
- `axis_opts` : Additional options to pass to `Axis`
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to `series!`
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
Figure or GridPosition
"""
function ADRIA.viz.scenarios(
    outcomes::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(:size => (800, 300)),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.scenarios!(
        g, outcomes; opts=opts, axis_opts=axis_opts, series_opts=series_opts
    )
    return f
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    outcomes::YAXArray;
    opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    series_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)::Union{GridLayout,GridPosition}
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(outcomes)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)

    if !haskey(axis_opts, :title)
        axis_opts[:title] = outcome_title(outcomes)
    end

    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    n_scenarios = size(outcomes, 2)
    scen_groups = Dict{Symbol,BitVector}(:scenarios => trues(n_scenarios))
    opts[:legend] = get(opts, :legend, false)
    opts[:histogram] = get(opts, :histogram, false)

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

        default_keys = [:counterfactual, :interventions, :unguided, :guided]
        filtered = default_keys[default_keys .∈ [scen_types]]
        return isempty(filtered) ? scen_types : filtered
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

"""
    scenario_bands!(ax, data, scen_groups; kwargs...)

Render grouped time-series data as median lines with CI bands, one band per scenario
group. When `summarize=false`, draws individual member lines instead.

# Arguments
- `ax`: a Makie `Axis`
- `data`: `YAXArray` with a scenarios dimension
- `scen_groups`: `OrderedDict{Symbol,BitVector}` mapping group name to member mask

# Keyword arguments
- `summarize::Bool=true` — aggregate to median + CI band; `false` draws all member lines
- `colormap=:tableau_10`
- `alpha::Float64=0.3` — opacity of CI band fill
- `legend_labels::Union{AbstractDict{Symbol,String},Nothing}=nothing` — optional override labels keyed by group name
"""
function scenario_bands!(
    ax,
    data::AbstractArray,
    scen_groups::AbstractDict{Symbol,BitVector};
    summarize::Bool=true,
    colormap=:tableau_10,
    alpha::Float64=0.3,
    legend_labels::Union{AbstractDict{Symbol,String},Nothing}=nothing
)
    _colors = categorical_colors(colormap, length(scen_groups))
    x_vals = collect(1:size(data, 1))
    _resolve_label =
        isnothing(legend_labels) ?
        (label -> string(label)) :
        (label -> get(legend_labels, label, string(label)))
    for (i, (label, mask)) in enumerate(scen_groups)
        display_label = _resolve_label(label)
        members = data isa YAXArray ? Array(data[scenarios = mask]) : data[:, mask]
        if summarize
            probs = [0.025, 0.5, 0.975]
            n_t = size(members, 1)
            confints = Matrix{Float64}(undef, n_t, 3)
            scratch = Vector{Float64}(undef, size(members, 2))
            for (t, sl) in enumerate(eachslice(members; dims=1))
                copyto!(scratch, sl)
                confints[t, :] .= quantile!(scratch, probs)
            end

            band!(ax, x_vals, confints[:, 1], confints[:, 3]; color=(_colors[i], alpha))
            lines!(ax, x_vals, confints[:, 2]; color=_colors[i], label=display_label)
        else
            for s in axes(members, 2)
                lines!(ax, x_vals, members[:, s];
                    color=(_colors[i], alpha), label=s == 1 ? display_label : nothing)
            end
        end
    end
end

"""
    ADRIA.viz.scenarios(ao::AnnotatedOutcomes; kwargs...) -> Figure

Plot scenario outcomes over time from an `AnnotatedOutcomes`, grouping lines by scenario
type or RCP and optionally rendering a legend panel.
"""
function ADRIA.viz.scenarios(
    ao::AnnotatedOutcomes;
    by_RCP::Bool=false,
    summarize::Bool=true,
    legend::Bool=true,
    sort_by::Symbol=:default,
    legend_labels::Union{AbstractDict{Symbol,String},Nothing}=nothing,
    size::Tuple{Int,Int}=(800, 300),
    title::AbstractString="",
    xlabel::AbstractString="Year",
    ylabel::AbstractString="",
    xticks=nothing,
    xticklabelrotation::Float64=π / 2
)::Figure
    scen_groups = _get_scenario_groups(ao; by_RCP)
    f = Figure(; size)
    g = f[1, 1] = GridLayout()
    xtick_vals = isnothing(xticks) ? _time_labels(timesteps(ao.data)) : xticks
    ax = Axis(g[1, 1];
        title, xlabel, ylabel,
        xticks=xtick_vals, xticklabelrotation
    )
    scenario_bands!(ax, ao.data, scen_groups;
        summarize, legend_labels
    )
    legend && _scenarios_legend_from_groups!(g[1, 2], scen_groups, ao.data;
        by_RCP, sort_by, legend_labels
    )
    return f
end

"""
    ADRIA.viz.scenarios!(g, ao::AnnotatedOutcomes; kwargs...) -> Union{GridLayout,GridPosition}

Render scenario outcome bands from `ao` into an existing grid position `g`.
"""
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    ao::AnnotatedOutcomes;
    by_RCP::Bool=false,
    summarize::Bool=true,
    legend::Bool=true,
    sort_by::Symbol=:default,
    legend_labels::Union{AbstractDict{Symbol,String},Nothing}=nothing,
    title::AbstractString="",
    xlabel::AbstractString="Year",
    ylabel::AbstractString="",
    xticks=nothing,
    xticklabelrotation::Float64=π / 2
)::Union{GridLayout,GridPosition}
    scen_groups = _get_scenario_groups(ao; by_RCP)
    xtick_vals = isnothing(xticks) ? _time_labels(timesteps(ao.data)) : xticks
    ax = Axis(g[1, 1];
        title, xlabel, ylabel,
        xticks=xtick_vals, xticklabelrotation
    )
    scenario_bands!(ax, ao.data, scen_groups;
        summarize, legend_labels
    )
    legend && _scenarios_legend_from_groups!(g[1, 2], scen_groups, ao.data;
        by_RCP, sort_by, legend_labels
    )
    return g
end

"""
    _scenarios_legend_from_groups!(g, scen_groups, data; kwargs...)

Build a scenario legend panel from pre-computed `scen_groups`, avoiding a redundant
`_get_scenario_groups` call when groups are already available at the call site.
"""
function _scenarios_legend_from_groups!(
    g::GridPosition,
    scen_groups::OrderedDict{Symbol,BitVector},
    data::YAXArray;
    by_RCP::Bool=false,
    sort_by::Symbol=:default,
    legend_labels::Union{AbstractDict{Symbol,String},Nothing}=nothing,
    legend_title::AbstractString="Scenarios"
)
    opts = Dict{Symbol,Any}(:by_RCP => by_RCP, :sort_by => sort_by)
    isnothing(legend_labels) || (opts[:legend_labels] = Symbol.(keys(legend_labels)))
    legend_opts = Dict{Symbol,Any}(:title => legend_title)
    return ADRIA.viz.scenarios_legend!(g, Dict{Symbol,BitVector}(scen_groups), data;
        opts=opts, legend_opts=legend_opts
    )
end

"""
    ADRIA.viz.scenarios_legend!(g, ao::AnnotatedOutcomes; kwargs...)

Add a scenario group legend panel to grid position `g` using metadata from `ao`.
"""
function ADRIA.viz.scenarios_legend!(
    g::GridPosition,
    ao::AnnotatedOutcomes;
    by_RCP::Bool=false,
    sort_by::Symbol=:default,
    legend_labels::Union{AbstractDict{Symbol,String},Nothing}=nothing,
    legend_title::AbstractString="Scenarios"
)
    scen_groups = _get_scenario_groups(ao; by_RCP)
    return _scenarios_legend_from_groups!(g, scen_groups, ao.data;
        by_RCP, sort_by, legend_labels, legend_title
    )
end
