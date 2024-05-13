using Statistics



function ADRIA.viz.dhw_scenario(dom::Domain, scen_id::Int64; fig_opts=Dict(), axis_opts=Dict())
    loc_scens = dom.dhw_scens[:, :, scen_id]
    mean_dhw_scen = dropdims(mean(loc_scens, dims=2), dims=2)

    ts = dom.env_layer_md.timeframe[1]:dom.env_layer_md.timeframe[end]

    fig_opts[:size] = get(axis_opts, :size, (800, 400))
    f = Figure(; fig_opts...)

    axis_opts[:xticks] = get(axis_opts, :xticks, _time_labels(ts))
    axis_opts[:xticklabelrotation] = get(axis_opts, :xticklabelrotation, 2 / π)
    axis_opts[:title] = get(axis_opts, :title, "DHW Scenario $(scen_id)")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Year")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "DHW")
    ax = Axis(f[1, 1]; axis_opts...)

    each_loc = series!(ax, loc_scens'; solid_color=(:red, 0.05))
    scen_mean = lines!(ax, mean_dhw_scen; color=(:black, 0.8))

    Legend(
        f[1, 2],
        [scen_mean, each_loc],
        ["Scenario Mean", "Locations"]
    )

    return f
end