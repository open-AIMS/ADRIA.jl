using Statistics

"""
    ADRIA.viz.dhw_scenario(
        dom::Domain,
        scen_id::Int64;
        fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
        axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
    )

Display the trajectory of a single DHW scenario.

# Arguments
- `dom` : ADRIA Domain
- `scen_id` : ID of DHW scenario to display
- `fig_opts` : Additional options to pass to adjust Figure creation
- `axis_opts` : Additional options to pass to adjust Axis attributes

# Returns
Figure
"""
function ADRIA.viz.dhw_scenario(
    dom::Domain,
    scen_id::Int64;
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    loc_scens = dom.dhw_scens[:, :, scen_id]
    mean_dhw_scen = dropdims(mean(loc_scens; dims=2); dims=2)

    fig_opts[:size] = get(axis_opts, :size, (800, 400))
    f = Figure(; fig_opts...)

    axis_opts[:title] = get(axis_opts, :title, "DHW Scenario $(scen_id)")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Year")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "DHW")
    axis_opts[:xticks] = get(axis_opts, :xticks, _time_labels(dom.env_layer_md.timeframe))
    axis_opts[:xticklabelrotation] = get(axis_opts, :xticklabelrotation, 2 / π)
    ax = Axis(f[1, 1]; axis_opts...)

    each_loc = series!(ax, loc_scens'; solid_color=(:red, 0.03))
    scen_mean = lines!(ax, mean_dhw_scen; color=(:black, 0.5))

    Legend(
        f[1, 2],
        [scen_mean, each_loc],
        ["Scenario Mean", "Locations"]
    )

    return f
end

"""
    dhw_scenarios(
        dom::Domain;
        ci_level=0.95,
        fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
        axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
    )

Display all DHW scenarios for the currently set RCP/SSP.

The set RCP/SSP can be changed with `switch_RCPs!()`.

# Arguments
- `dom` : ADRIA Domain
- `ci_level` : Level of Confidence Interval
- `fig_opts` : Additional options to pass to adjust Figure creation
- `axis_opts` : Additional options to pass to adjust Axis attributes

# Returns
Figure
"""
function ADRIA.viz.dhw_scenarios(
    dom::Domain;
    scens=(:),
    ci_level=0.95,
    fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(),
    axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
)
    data = dom.dhw_scens[:, :, scens]

    # Extract dimension information
    years = dom.env_layer_md.timeframe

    fig_size = get(fig_opts, :size, (1200, 600))

    axis_opts[:title] = get(axis_opts, :title, "DHW Scenarios")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Year")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "DHW")
    axis_opts[:xticklabelrotation] = get(axis_opts, :xticklabelrotation, 2 / π)
    fig = Figure(; size=fig_size)
    ax = Axis(fig[1, 1]; axis_opts...)

    # Determine CI bounds
    α = (1 - ci_level) / 2

    means = zeros(length(years))
    lower_ci = zeros(length(years))
    upper_ci = zeros(length(years))

    for t in collect(axes(data, 1))
        # Calculate statistics across scenarios for each timestep
        means[t] = mean(data[t, :, :])
        lower_ci[t] = quantile(data[t, :, :][:], α)
        upper_ci[t] = quantile(data[t, :, :][:], 1 - α)
    end

    ci_band = band!(ax, years, lower_ci, upper_ci; color=(:red, 0.3))
    each_scen = series!(
        ax, years, dropdims(mean(data; dims=:sites); dims=:sites).data';
        solid_color=(:red, 0.2)
    )
    mean_line = lines!(ax, years, means; color=(:black, 0.7))

    # Add legend
    Legend(
        fig[1, 2],
        [ci_band, each_scen, mean_line],
        ["$(floor(Int64, ci_level*100))% CI", "Individual trajectories", "Mean"]
    )

    return fig
end
