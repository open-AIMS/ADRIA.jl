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
- `y` : results of scenario metric
- `opts` : Aviz options
    - `by_RCP` : color by RCP otherwise color by scenario type. Defaults to false.
    - `legend` : show legend. Defaults to true.
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!

# Returns
GridPosition
"""
function ADRIA.viz.scenarios(
    rs::ResultSet,
    y::NamedDimsArray;
    opts::Dict=Dict(:by_RCP => false),
    fig_opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Figure
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()
    ADRIA.viz.scenarios!(g, rs, y; opts, axis_opts, series_opts)

    return f
end
function ADRIA.viz.scenarios!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    y::NamedDimsArray;
    opts::Dict=Dict(),
    axis_opts::Dict=Dict(),
    series_opts::Dict=Dict(),
)::Union{GridLayout,GridPosition}
    # Ensure last year is always shown in x-axis
    xtick_vals = get(axis_opts, :xticks, _time_labels(timesteps(y)))
    xtick_rot = get(axis_opts, :xticklabelrotation, 2 / π)

    ax = Axis(g[1, 1]; xticks=xtick_vals, xticklabelrotation=xtick_rot, axis_opts...)

    min_step = (1.0 / 0.05)
    color_weight = max(min((1.0 / (size(y, 2) / min_step)), 0.6), 0.05)

    # Create legend for each guided scenario type
    if :color ∉ keys(series_opts)
        hide_idx = get(series_opts, :hide_series, BitVector())

        if get(opts, :by_RCP, false)
            rcp::Vector{Symbol} = Symbol.(:RCP, Int64.(rs.inputs[:, :RCP]))
            unique_rcp = sort((unique(rcp)))

            line_elements::Vector{LineElement} = LineElement[
                LineElement(; color=_c, linestyle=nothing) for
                _c in [COLORS[r] for r in unique_rcp]
            ]
            LineElement.(; color=[COLORS[r] for r in unique_rcp], linestyle=nothing)

            series_opts = merge(
                series_opts, Dict(:color => map(x -> (COLORS[x], color_weight), rcp))
            )
            labels = String.(unique_rcp)
        else
            series_opts = merge(
                series_opts, Dict(:color => scenario_colors(rs, color_weight, hide_idx))
            )

            cf = LineElement(; color=COLORS[:counterfactual], linestyle=nothing)
            ug = LineElement(; color=COLORS[:unguided], linestyle=nothing)
            gu = LineElement(; color=COLORS[:guided], linestyle=nothing)

            line_elements = [cf, ug, gu]
            labels = ["No Intervention", "Unguided", "Guided"]
        end

        # Add legend
        if get(opts, :legend, true)
            Legend(
                g[1, 3],
                line_elements,
                labels;
                halign=:left,
                valign=:top,
                margin=(5, 5, 5, 5),
            )
        end
    end

    _plot_scenarios_series!(ax, y; series_opts)
    _plot_scenarios_hist(g, rs, y)

    return g
end

function _plot_scenarios_series!(ax::Axis, data::NamedDimsArray; series_opts)::Nothing
    series!(ax, data'; series_opts...)

    # ax.ylabel = metric_label(metric)
    ax.xlabel = "Year"

    return nothing
end

function _plot_scenarios_hist(
    g::Union{GridLayout,GridPosition}, rs::ResultSet, data::NamedDimsArray
)::Nothing
    scen_match = 1:nrow(rs.inputs) .∈ [_dimkeys(data).scenarios]
    scen_types = scenario_type(rs; scenarios=scen_match)
    scen_dist = dropdims(mean(data; dims=:timesteps); dims=:timesteps)

    color_weights = (counterfactual=0.8, unguided=0.7, guided=0.6)

    ax_hist = Axis(g[1, 2]; width=100)
    for type in keys(scen_types)
        if !isempty(scen_types[type])
            hist!(
                ax_hist,
                scen_dist[scen_types[type]];
                direction=:x,
                color=(COLORS[type], color_weights[type]),
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
