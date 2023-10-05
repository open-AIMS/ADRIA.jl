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

    rs_inputs = copy(rs.inputs)
    keepat!(rs_inputs, data.scenarios)

    # Set series colors
    by_rcp = get(opts, :by_RCP, false)
    colors = by_rcp ? _RCP_colors(rs_inputs, data) : _type_colors(rs_inputs)
    merge!(series_opts, colors)

    if get(opts, :summarize, true)
        _plot_scenarios_confint!(ax, rs_inputs, data)
    else
        _plot_scenarios_series!(ax, rs_inputs, data, series_opts)
    end

    # Plot Histograms when opts[:histogram] is true
    get(opts, :histogram, true) ? _plot_scenarios_hist(g, rs, data) : nothing

    # Render legend
    legend_position = get(opts, :histogram, true) ? (1, 3) : (1, 2)
    _render_scenarios_legend(g, rs, legend_position, opts)

    ax.xlabel = "Year"
    # ax.ylabel = metric_label(metric)

    return g
end

function _plot_scenarios_confint!(
    ax::Axis, rs_inputs::DataFrame, data::NamedDimsArray
)::Nothing
    n_timesteps = size(data, 1)
    x_timesteps::UnitRange{Int64} = 1:n_timesteps
    scenario_types = ADRIA.analysis.scenario_types(rs_inputs)
    ordered_types = _sort(scenario_types, data)

    selected_scenarios = [scenario_types[type] for type in ordered_types]
    colors = [COLORS[scenario] for scenario in keys(scenario_types)]

    confints = zeros(n_timesteps, length(scenario_types), 3)
    for (idx_s, scenario) in enumerate(selected_scenarios)
        confints[:, idx_s, :] = ADRIA.analysis.series_confint(
            data[:, scenario]; agg_dim=:scenarios
        )
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
    ax::Axis, rs_inputs::DataFrame, data::NamedDimsArray, series_opts::Dict
)::Nothing
    scenario_types = ADRIA.analysis.scenario_types(rs_inputs)

    colors = pop!(series_opts, :color)
    _alphas = alphas(scenario_types)

    for type in _sort(scenario_types, data)
        selected_scenarios = scenario_types[type]
        color = (colors[selected_scenarios][1], _alphas[type])

        series!(ax, data[:, selected_scenarios]'; solid_color=color, series_opts...)
    end
    return nothing
end

function _plot_scenarios_hist(
    g::Union{GridLayout,GridPosition}, rs::ResultSet, data::NamedDimsArray
)::Nothing
    scen_match = 1:nrow(rs.inputs) .∈ [_dimkeys(data).scenarios]
    scen_types = ADRIA.analysis.scenario_types(rs; scenarios=scen_match)
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

# TODO Move these two to theme.jl
function _RCP_colors(rs_input::DataFrame, data::NamedDimsArray)
    rcp::Vector{Symbol} = Symbol.(:RCP, Int64.(rs_input[:, :RCP]))
    return Dict(:color => map(x -> (COLORS[x], _color_weight(data)), rcp))
end

function _type_colors(rs_input::DataFrame)
    return Dict(:color => scenario_colors(rs_input))
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
    _sort(scenario_types::Dict{Symbol, BitVector}, data::NamedDimsArray)::Vector{Symbol}

Sort types by variance in reverse order to plot highest variances first

# Arguments
- `data` : Results of scenario metric
- `scenario_types` : NamedTuple of BitVectors to filter scenarios for each scenario type of:
    - :guided
    - :unguided
    - :counterfactual
"""
function _sort(scenario_types::Dict{Symbol,BitVector}, data::NamedDimsArray)::Vector{Symbol}
    scen_types::Vector{Symbol} = collect(keys(scenario_types))
    return sort(
        scen_types; by=type -> sum(var(data[:, scenario_types[type]]; dims=2)), rev=true
    )
end
