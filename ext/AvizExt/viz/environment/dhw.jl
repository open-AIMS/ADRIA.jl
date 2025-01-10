using Statistics

function ADRIA.viz.dhw_scenario(
    dom::Domain, scen_id::Int64; fig_opts::OPT_TYPE=DEFAULT_OPT_TYPE(), axis_opts::OPT_TYPE=DEFAULT_OPT_TYPE()
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
