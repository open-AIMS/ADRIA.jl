using ADRIA: DEAResult

"""
    ADRIA.viz.data_envelopment_analysis(rs::ResultSet, DEA_output::DEAResult;axis_opts=Dict(),
        fig_opts=Dict(), opts=Dict())
    ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition},rs::ResultSet,
        DEA_output::DEAResult; axis_opts=Dict(),opts=Dict())
    ADRIA.viz.data_envelopment_analysis(g::Union{GridLayout,GridPosition}, DEA_output::DEAResult;
        axis_opts=Dict(), opts=Dict())

Plot results from a DEA analysis. Plots the first 2 dimensions of the effificency frontier,
along side the technical and scale efficiencies.

# Examples
```
# Run scenarios
scens = ADRIA.sample(dom, 2^12)
rs = ADRIA.run_scenarios(dom, scens, ["45"])

# Compute cost an mean metrics for each scenario
CAD_cost = ADRIA.CAD_cost(scens)
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
DEA_output = ADRIA.data_envelopment_analysis(CAD_cost, s_tac, s_sv)

# Plot frontier, sclae and technical efficiencies
ADRIA.viz.data_envelopment_analysis(rs, DEA_output)
```

# Arguments
- `rs` : ResultSet
- `DEA_output` : output structure from a DEA analysis, carried out using `data_envelopment_analysis`
- `opts` : Aviz options
    - `frontier_type`, type of frontier to plot: "CRS","VRS" or "FDH".
    - `line_color`, color to use for best practice frontier.
    - `data_color`, color to use for data cloud when plotting efficiency frontier.
    - `frontier_name`, text to label efficiency frontier.
    - `data_name`, text to label data cloud used to plot efficiency frontier.
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `fig_opts` : Additional options to pass to adjust Figure creation
  See: https://docs.makie.org/v0.19/api/index.html#Figure
"""
function ADRIA.viz.data_envelopment_analysis(
    rs::ResultSet, DEA_output::DEAResult;
    axis_opts=Dict(), fig_opts=Dict(),
    opts=Dict()
)
    f = Figure(; fig_opts...)
    g = f[1, 1] = GridLayout()

    ADRIA.viz.data_envelopment_analysis!(
        g, rs, DEA_output; opts=opts, axis_opts=axis_opts
    )

    return f
end
function ADRIA.viz.data_envelopment_analysis!(g::Union{GridLayout,GridPosition},
    rs::ResultSet, DEA_output::DEAResult; axis_opts=Dict(),
    opts=Dict()
)
    return ADRIA.viz.data_envelopment_analysis!(
        g, DEA_output; opts=opts, axis_opts=axis_opts
    )
end
function ADRIA.viz.data_envelopment_analysis!(g::Union{GridLayout,GridPosition},
    DEA_output::DEAResult; axis_opts=Dict(), opts=Dict())
    line_color = get(opts, :line_color, :red)
    data_color = get(opts, :data_color, :black)
    frontier_name = get(opts, :frontier_name, "Best practice frontier")
    data_name = get(opts, :data_name, "Scenario data cloud")

    # Determines which returns to scale approach is used to select scenario peers
    # (most efficient scenarios)
    frontier_type = get(opts, :frontier_type, :vrs_peers)

    Y = DEA_output.Y # Output values

    # Find points on best practice frontier
    best_practice_scens = DEA_output.$frontier_type.J

    scale_efficiency = DEA_output.crs_vals ./ DEA_output.vrs_vals

    # Plot efficiency frontier and data cloud
    axa = Axis(g[1, 1]; axis_opts...)
    scatter!(axa, Y[:, 1], Y[:, 2]; color=data_color)
    scatter!(
        axa, Y[best_practice_scens, 1], Y[best_practice_scens, 2]; color=line_color
    )
    #Legend(g[1, 2], [frontier, data], [frontier_name, data_name])

    # Plot the scale efficiency (ratio of efficiencies assuming CRS vs. assuming VRS)
    axb = Axis(g[2, 1]; axis_opts...)
    scatter!(axb, scale_efficiency; color=data_color)#, title="Scale efficiency")

    # Plot the technical efficiency (inverse VRS efficiencies)
    axc = Axis(g[3, 1]; axis_opts...)
    scatter!(axc, DEA_output.vrs_vals; color=data_color)#, title="Technical efficiency")
    return g
end
