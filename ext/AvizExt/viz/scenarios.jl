using JuliennedArrays: Slices
using ADRIA.analysis: series_confint

"""
    ADRIA.viz.scenarios(rs::ADRIA.ResultSet, outcomes::NamedDimsArray; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.scenarios(rs_inputs::DataFrame, outcomes::NamedDimsArray; opts::Dict=Dict(:by_RCP => false), fig_opts::Dict=Dict(), axis_opts::Dict=Dict(), series_opts::Dict=Dict())::Figure
    ADRIA.viz.scenarios!(g::Union{GridLayout,GridPosition}, rs_inputs::DataFrame, outcomes::NamedDimsArray; opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

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
- `rs_input` : DataFrame with ResultSet inputs
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
    outcomes::NamedDimsArray;
    opts::Dict=Dict(:by_RCP => false),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Figure
    return ADRIA.viz.scenarios(
        rs.inputs,
        outcomes;
        opts=opts,
        fig_opts=fig_opts,
        axis_opts=axis_opts,
        series_opts=series_opts,
    )
end
function ADRIA.viz.scenarios(
    rs_inputs::DataFrame,
    outcomes::NamedDimsArray;
    opts::Dict=Dict(:by_RCP => false),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.scenarios!(
        g, rs_inputs, outcomes; opts=opts, axis_opts=axis_opts, series_opts=series_opts
    )

    return f
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    rs_inputs::DataFrame,
    outcomes::NamedDimsArray;
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Union{GridLayout,GridPosition}
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(outcomes)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    _rs_inputs = copy(rs_inputs[1:end .== outcomes.scenarios, :])
    scen_groups = if get(opts, :by_RCP, false)
        ADRIA.analysis.scenario_rcps(_rs_inputs)
    else
        ADRIA.analysis.scenario_types(_rs_inputs)
    end

    if get(opts, :summarize, true)
        _plot_scenarios_confint!(ax, outcomes, scen_groups)
    else
        _plot_scenarios_series!(ax, outcomes, scen_groups, series_opts)
    end

    get(opts, :histogram, true) ? _plot_scenarios_hist(g, outcomes, scen_groups) : nothing

    if get(opts, :legend, true)
        legend_position = get(opts, :histogram, true) ? (1, 3) : (1, 2)
        _render_legend(g, scen_groups, legend_position)
    end

    ax.xlabel = "Year"
    # ax.ylabel = metric_label(metric)

    return g
end

function _plot_scenarios_confint!(
    ax::Axis, outcomes::NamedDimsArray, scen_groups::Dict{Symbol,BitVector}
)::Nothing
    ordered_groups::Vector{Symbol} = _sort_keys(scen_groups, outcomes)
    n_timesteps::Int64 = size(outcomes, 1)
    n_scens::Int64 = length(ordered_groups)

    # Compute confints
    confints::Array{Float64} = zeros(n_timesteps, n_scens, 3)
    for (idx, group) in enumerate(ordered_groups)
        confints[:, idx, :] = series_confint(outcomes[:, scen_groups[group]])
    end

    _colors::Dict{Symbol,Symbol} = colors(scen_groups)

    for idx in eachindex(ordered_groups)
        band_color = (_colors[ordered_groups[idx]], 0.4)
        y_lower, y_upper = confints[:, idx, 1], confints[:, idx, 3]
        band!(ax, 1:n_timesteps, y_lower, y_upper; color=band_color)
    end

    series_colors = [_colors[group] for group in keys(scen_groups)]
    series!(ax, confints[:, :, 2]'; solid_color=series_colors)

    return nothing
end

function _plot_scenarios_series!(
    ax::Axis,
    outcomes::NamedDimsArray,
    scen_groups::Dict{Symbol,BitVector},
    series_opts::Dict,
)::Nothing
    _colors = colors(scen_groups)
    _alphas = alphas(scen_groups)

    for group in _sort_keys(scen_groups, outcomes)
        color = (_colors[group], _alphas[group])
        scens = outcomes[:, scen_groups[group]]'
        series!(ax, scens; solid_color=color, series_opts...)
    end

    return nothing
end

function _plot_scenarios_hist(
    g::Union{GridLayout,GridPosition},
    outcomes::NamedDimsArray,
    scen_groups::Dict{<:Any,BitVector},
)::Nothing
    scen_dist = dropdims(mean(outcomes; dims=:timesteps); dims=:timesteps)
    ax_hist = Axis(g[1, 2]; width=100)
    _colors = colors(scen_groups)

    for group in _sort_keys(scen_groups, outcomes)
        color = (_colors[group], 0.7)
        dist = scen_dist[scen_groups[group]]
        hist!(ax_hist, dist; direction=:x, color=color, bins=30, normalization=:pdf)
    end

    hidedecorations!(ax_hist)
    hidespines!(ax_hist)
    ylims!(
        ax_hist,
        minimum(scen_dist) - quantile(scen_dist, 0.05),
        maximum(scen_dist) + quantile(scen_dist, 0.05),
    )

    return nothing
end

function _render_legend(
    g::Union{GridLayout,GridPosition},
    scen_groups::Dict{<:Any,BitVector},
    legend_position::Tuple{Int64,Int64},
)::Nothing
    _colors = colors(scen_groups)
    line_els::Vector{LineElement} = [
        LineElement(; color=_colors[group]) for group in keys(scen_groups)
    ]
    labels::Vector{String} = string.(keys(scen_groups))

    Legend(
        g[legend_position...],
        line_els,
        labels;
        halign=:left,
        valign=:top,
        margin=(5, 5, 5, 5),
    )

    return nothing
end

"""
    _sort_keys(scenario_types::Dict{Symbol, BitVector}, data::NamedDimsArray)::Vector{Symbol}

Sort types by variance in reverse order.

# Arguments
- `data` : Results of scenario metric
- `scenario_types` : NamedTuple of BitVectors to filter scenarios for each scenario type of:
    - :guided
    - :unguided
    - :counterfactual
"""
function _sort_keys(
    scenario_types::Dict{Symbol,BitVector}, data::NamedDimsArray
)::Vector{Symbol}
    scen_types::Vector{Symbol} = collect(keys(scenario_types))
    return sort(
        scen_types; by=type -> sum(var(data[:, scenario_types[type]]; dims=2)), rev=true
    )
end
