using Statistics



function ADRIA.viz.cyclone_scenario(dom::Domain, scen_id::Int64; fig_opts=Dict(), axis_opts=Dict())
    taxa_mean = dropdims(
        mean(dom.cyclone_mortality_scens[:, :, :, scen_id], dims=:species), dims=:species
    )
    mean_cyc_scen = dropdims(mean(taxa_mean, dims=:locations), dims=:locations)

    fig_opts[:size] = get(axis_opts, :size, (800, 400))
    f = Figure(; fig_opts...)

    axis_opts[:xticks] = get(axis_opts, :xticks, collect(taxa_mean.timesteps[1:5:end]))
    axis_opts[:xticklabelrotation] = get(axis_opts, :xticklabelrotation, 2 / π)
    axis_opts[:title] = get(axis_opts, :title, "Cyclone Mortality Scenario $(scen_id)")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Timesteps [years]")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Cyclone Mortality")
    ax = Axis(f[1, 1]; axis_opts...)

    each_loc = series!(ax, taxa_mean'; solid_color=(:red, 0.05))
    scen_mean = lines!(ax, mean_cyc_scen; color=(:black, 0.8))

    Legend(
        f[1, 2],
        [scen_mean, each_loc],
        ["Scenario Mean", "Locations"]
    )

    return f
end