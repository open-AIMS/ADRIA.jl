using ADRIA.analysis: DEAResult

"""
    ADRIA.viz.data_envelopment_analysis(rs::ResultSet, DEA_output::DEAResult; axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(), fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(), opts::OPT_TYPE=DEFAULT_OPT_TYPE())
    ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition}, rs::ResultSet, DEA_output::DEAResult; axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(), opts::OPT_TYPE=DEFAULT_OPT_TYPE())
    ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition}, DEA_output::DEAResult; axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(), opts::OPT_TYPE=DEFAULT_OPT_TYPE())

Plot results from a DEA analysis. Plots the first 2 dimensions of the effificency frontier,
along side the technical and scale efficiencies.

# Examples
```
# Run scenarios
scens = ADRIA.sample(dom, 2^12)
rs = ADRIA.run_scenarios(dom, scens, ["45"])

# Compute cost an mean metrics for each scenario
cost = cost_function(scens)
s_tac::Vector{Float64} = Array(
    dropdims(
        mean(ADRIA.metrics.scenario_total_cover(rs); dims=:timesteps); dims=:timesteps
    )
)
s_sv::Vector{Float64} = Array(
    dropdims(
        mean(mean(ADRIA.metrics.absolute_shelter_volume(rs); dims=:timesteps); dims=:sites);
        dims=(:timesteps, :sites)
    )
)

# Apply DEA analysis
DEA_output = ADRIA.data_envelopment_analysis(cost, s_tac, s_sv)

# Plot frontier, scale and technical efficiencies
ADRIA.viz.data_envelopment_analysis(rs, DEA_output)
```

# Arguments
- `rs` : ResultSet
- `DEA_output` : output structure from a DEA analysis, carried out using `data_envelopment_analysis`
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure
- `opts` : Aviz options
    - `frontier_type`, type of frontier to plot: "CRS","VRS" or "FDH".
    - `line_color`, color to use for best practice frontier.
    - `data_color`, color to use for data cloud when plotting efficiency frontier.
    - `frontier_name`, text to label efficiency frontier.
    - `data_name`, text to label data cloud used to plot efficiency frontier.
    - `metrics_x_lab` : Name of metric on x axis.
    - `metrics_y_lab` : Name of metric on y axis.
    - `scale_eff_y_lab` : Label for scale efficiency.
    - `scale_eff_y_lab` : Label for technical efficiency.
"""
function ADRIA.viz.data_envelopment_analysis(
    rs::ResultSet,
    DEA_output::DEAResult;
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.data_envelopment_analysis!(
        g, rs, DEA_output; opts=opts, axis_opts=axis_opts
    )

    return f
end
function ADRIA.viz.data_envelopment_analysis!(
    g::Union{GridLayout,GridPosition},
    rs::ResultSet,
    DEA_output::DEAResult;
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    return ADRIA.viz.data_envelopment_analysis!(
        g, DEA_output; opts=opts, axis_opts=axis_opts
    )
end
function ADRIA.viz.data_envelopment_analysis!(
    g::Union{GridLayout,GridPosition},
    DEA_output::DEAResult;
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    opts::OPT_TYPE=DEFAULT_OPT_TYPE())
    scatter_colors = get(opts, :scatter_colors, [:red, :black])
    legend_names = get(
        opts, :legend_names, ["Best practice frontier", "Scenario data cloud"]
    )
    scale_eff_y_lab = get(opts, :scale_eff_y_lab, L"$\frac{eff_{vrs}}{eff_{crs}}$")
    tech_eff_y_lab = get(opts, :tech_eff_y_lab, L"$\frac{1}{eff_{vrs}}$")
    metrics_x_lab = get(opts, :metrics_x_lab, L"$metric 1$")
    metrics_y_lab = get(opts, :metrics_y_lab, L"$metric 2$")

    # Determines which returns to scale approach is used to select scenario peers
    # (most efficient scenarios)
    frontier_type = get(opts, :frontier_type, :vrs_peers)

    Y = DEA_output.Y # Output values

    # Find points on best practice frontier
    best_practice_scens = getfield(DEA_output, frontier_type).J

    scale_efficiency = DEA_output.crs_vals ./ DEA_output.vrs_vals

    # Plot efficiency frontier and data cloud
    ax_a = Axis(g[1, 1]; xlabel=metrics_x_lab, ylabel=metrics_y_lab, axis_opts...)
    scatter!(ax_a, Y[:, 1], Y[:, 2]; color=scatter_colors[2])
    scatter!(
        ax_a, Y[best_practice_scens, 1], Y[best_practice_scens, 2]; color=scatter_colors[1]
    )
    _render_marker_legend(g, (1, 2), legend_names, scatter_colors)

    # Plot the scale efficiency (ratio of efficiencies assuming CRS vs. assuming VRS)
    ax_b = Axis(g[2, 1]; title="Scale efficiency", ylabel=scale_eff_y_lab, axis_opts...)
    scatter!(ax_b, scale_efficiency; color=scatter_colors[2])
    scatter!(
        ax_b,
        best_practice_scens,
        scale_efficiency[best_practice_scens];
        color=scatter_colors[1]
    )

    # Plot the technical efficiency (inverse VRS efficiencies)
    ax_c = Axis(g[3, 1]; title="Technical efficiency", ylabel=tech_eff_y_lab, axis_opts...)
    scatter!(ax_c, DEA_output.vrs_vals; color=scatter_colors[2])
    scatter!(
        ax_c,
        best_practice_scens,
        DEA_output.vrs_vals[best_practice_scens];
        color=scatter_colors[1]
    )

    return g
end

function _render_marker_legend(g::Union{GridLayout,GridPosition},
    legend_position::Tuple{Int64,Int64},
    legend_labels::Union{Vector{Symbol},Vector{String}},
    colors::Union{Vector{Symbol},RGBA{Float32}}
)::Nothing
    marker_els::Vector{MarkerElement} = [
        MarkerElement(; color=color, marker=:circle) for color in colors
    ]
    Legend(g[legend_position...], marker_els, legend_labels; framevisible=false)

    return nothing
end
