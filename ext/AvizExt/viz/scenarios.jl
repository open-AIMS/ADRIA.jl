using JuliennedArrays: Slices

"""
    ADRIA.viz.scenarios(rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), fig_opts=Dict(), axis_opts=Dict(), series_opts=Dict())
    ADRIA.viz.scenarios!(g::Union{GridLayout,GridPosition}, rs::ADRIA.ResultSet, y::NamedDimsArray; opts=Dict(by_RCP => false), axis_opts=Dict(), series_opts=Dict())

Plot scenario outcomes over time.

# Examples
```julia
scens = ADRIA.sample(dom, 64)

# ... run scenarios ...

s_tac = ADRIA.metrics.scenario_total_cover(rs)

# Plot scenario outcomes
ADRIA.viz.scenarios(rs, s_tac)

# Plot outcomes of scenarios where SRM < 1.0
ADRIA.viz.scenarios(rs, s_tac[:, scens.SRM .< 1.0])
```

# Arguments
- `rs` : ResultSet
- `data` : results of scenario metric
- `opts` : Aviz options
    - `by_RCP` : color by RCP otherwise color by scenario type. Defaults to false.
    - `legend` : show legend. Defaults to true.
    - `summarize` : plot confidence interval. Defaults to true
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
GridPosition
"""
function ADRIA.viz.scenarios(
    rs::ResultSet,
    data::NamedDimsArray;
    opts::Dict=Dict(:by_RCP => false),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.scenarios!(g, rs, data; opts, axis_opts, series_opts)

    return f
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    data::NamedDimsArray;
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Union{GridLayout,GridPosition}
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(data)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)

    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    # Set series colors
    merge!(series_opts, _get_series_opt_colors(rs, data, opts, series_opts))

    if get(opts, :summarize, true)
        _plot_scenarios_confint!(ax, rs, data)
    else
        _plot_scenarios_series!(ax, rs, data, series_opts)
    end

    # Plot Histograms when opts[:histogram] is true
    get(opts, :histogram, true) && _plot_scenarios_hist(g, rs, data)

    # Render legend
    legend_position = get(opts, :histogram, true) ? (1, 3) : (1, 2)
    _render_scenarios_legend(g, rs, legend_position, opts)

    ax.xlabel = "Year"
    # ax.ylabel = metric_label(metric)

    return g
end

function _plot_scenarios_confint!(ax::Axis, rs::ResultSet, data::NamedDimsArray)::Nothing
    n_timesteps = size(data, 1)
    x_timesteps::UnitRange{Int64} = 1:n_timesteps
    scenario_types = scenario_type(rs)
    ordered_types = _order_by_variance(data, scenario_types)

    selected_scenarios = [scenario_types[type] for type in ordered_types]
    colors = [scenario_colors(rs)[scenario][1][1] for scenario in selected_scenarios]

    confints = zeros(n_timesteps, length(scenario_types), 3)
    for (idx_s, scenario) in enumerate(selected_scenarios)
        confints[:, idx_s, :] = series_confint(data[:, scenario]; agg_dim=:scenarios)
    end

    for idx in eachindex(ordered_types)
        band_alpha = max(0.7 - idx * 0.1, 0.4)
        band_color = (colors[idx], band_alpha)
        y_lower, y_upper = confints[:, idx, 1], confints[:, idx, 3]
        band!(ax, x_timesteps, y_lower, y_upper; color=band_color)
    end

    series!(ax, confints[:, :, 2]'; solid_color=colors)

    return nothing
end

function _plot_scenarios_series!(
    ax::Axis, rs::ResultSet, data::NamedDimsArray, series_opts::Dict
)::Nothing
    series_colors = pop!(series_opts, :color)
    scenario_types = scenario_type(rs)
    for type in _order_by_variance(data, scenario_types)
        selected_scenarios = scenario_types[type]
        _color = series_colors[selected_scenarios]

        series!(ax, data[:, selected_scenarios]'; solid_color=_color, series_opts...)
    end
    return nothing
end

function _plot_scenarios_hist(
    g::Union{GridLayout,GridPosition}, rs::ResultSet, data::NamedDimsArray
)::Nothing
    scen_match = 1:nrow(rs.inputs) .∈ [_dimkeys(data).scenarios]
    scen_types = scenario_type(rs; scenarios=scen_match)
    scen_dist = dropdims(mean(data; dims=:timesteps); dims=:timesteps)

    hist_color_weights = (counterfactual=0.8, unguided=0.7, guided=0.6)

    ax_hist = Axis(g[1, 2]; width=100)
    for type in keys(scen_types)
        if !isempty(scen_types[type])
            hist!(
                ax_hist,
                scen_dist[scen_types[type]];
                direction=:x,
                color=(COLORS[type], hist_color_weights[type]),
                bins=30,
                normalization=:pdf,
            )
        end
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

function _get_series_opt_colors(
    rs::ResultSet, data::NamedDimsArray, opts::Dict, series_opts::Dict
)::Dict{Symbol,Vector{Tuple{Symbol,Float64}}}
    if get(opts, :by_RCP, false)
        rcp::Vector{Symbol} = Symbol.(:RCP, Int64.(rs.inputs[:, :RCP]))
        return Dict(:color => map(x -> (COLORS[x], _color_weight(data)), rcp))
    else
        hide_idx = get(series_opts, :hide_series, BitVector())
        return Dict(:color => scenario_colors(rs, _color_weight(data), hide_idx))
    end
end

function _color_weight(data::NamedDimsArray)::Float64
    min_step::Float64 = (1.0 / 0.05)
    return max(min((1.0 / (size(data, 2) / min_step)), 0.6), 0.05)
end

function _render_scenarios_legend(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    legend_position::Tuple{Int64,Int64},
    opts::Dict,
)::Nothing
    labels::Vector{String} = Vector{String}(undef, 0)
    line_elements::Vector{LineElement} = Vector{LineElement}(undef, 0)

    if get(opts, :by_RCP, false)
        rcp::Vector{Symbol} = sort((unique(Symbol.(:RCP, Int64.(rs.inputs[:, :RCP])))))

        rcp_colors = [COLORS[r] for r in rcp]
        line_elements = [LineElement(; color=c, linestyle=nothing) for c in rcp_colors]
        labels = String.(rcp)
    else
        cf::LineElement = LineElement(; color=COLORS[:counterfactual], linestyle=nothing)
        ug::LineElement = LineElement(; color=COLORS[:unguided], linestyle=nothing)
        gu::LineElement = LineElement(; color=COLORS[:guided], linestyle=nothing)

        line_elements = [cf, ug, gu]
        labels = ["No Intervention", "Unguided", "Guided"]
    end

    # Add legend
    if get(opts, :legend, true)
        Legend(
            g[legend_position...],
            line_elements,
            labels;
            halign=:left,
            valign=:top,
            margin=(5, 5, 5, 5),
        )
    end

    return nothing
end

"""
    function _order_by_variance(data::NamedDimsArray, scenario_types::NamedTuple)::Tuple{Symbol,Symbol,Symbol}

Sort types by variance in reverse order to plot highest variances first

# Arguments
- `data` : Results of scenario metric
- `scenario_types` : Named tuple with BitVectors to filter scenarios of each type (
    :guided, :unguided or :counterfactual)
"""
function _order_by_variance(
    data::NamedDimsArray, scenario_types::NamedTuple
)::Tuple{Symbol,Symbol,Symbol}
    return sort(
        keys(scenario_types);
        by=type -> sum(var(data[:, scenario_types[type]]; dims=2)),
        rev=true,
    )
end
