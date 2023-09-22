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
)
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
)
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

    series!(ax, y'; series_opts...)

    # Density (TODO: Separate into own function)
    scen_match = 1:nrow(rs.inputs) .∈ [_dimkeys(y).scenarios]
    scen_types = scenario_type(rs; scenarios=scen_match)
    scen_dist = dropdims(mean(y; dims=:timesteps); dims=:timesteps)
    ax2 = Axis(g[1, 2]; width=100)
    if count(scen_types.counterfactual) > 0
        hist!(
            ax2,
            scen_dist[scen_types.counterfactual];
            direction=:x,
            color=(COLORS[:counterfactual], 0.8),
            bins=30,
            normalization=:pdf,
        )
    end

    if count(scen_types.unguided) > 0
        hist!(
            ax2,
            scen_dist[scen_types.unguided];
            direction=:x,
            color=(COLORS[:unguided], 0.7),
            bins=30,
            normalization=:pdf,
        )
    end

    if count(scen_types.guided) > 0
        hist!(
            ax2,
            scen_dist[scen_types.guided];
            direction=:x,
            color=(COLORS[:guided], 0.6),
            bins=30,
            normalization=:pdf,
        )
    end

    hidedecorations!(ax2)
    hidespines!(ax2)
    ylims!(
        ax2,
        minimum(scen_dist) - quantile(scen_dist, 0.05),
        maximum(scen_dist) + quantile(scen_dist, 0.05),
    )

    # ax.ylabel = metric_label(metric)
    ax.xlabel = "Year"

    return g
end
