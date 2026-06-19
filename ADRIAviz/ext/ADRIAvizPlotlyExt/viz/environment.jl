using Statistics

# ──────────────────────────────────────────────────────────────────────────────
# cyclone_scenario
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.cyclone_scenario(cyclone_scens::YAXArray, scen_id::Int; kwargs...)
    ADRIA.viz.cyclone_scenario(dom::Domain, scen_id::Int; kwargs...)

Plot per-location cyclone mortality trajectory for a single scenario.
`cyclone_scens` is expected to have dimensions `(timesteps × locations × species × scenarios)`.
Each location is drawn as a semi-transparent line; the cross-location mean is overlaid.

# Keyword arguments
- `xlabel` : x-axis label (default `"Timestep"`)
- `ylabel` : y-axis label (default `"Cyclone Mortality [%]"`)
- `title`  : plot title (default `"Cyclone Mortality Scenario <scen_id>"`)
"""
function ADRIA.viz.cyclone_scenario(
    cyclone_scens::YAXArray,
    scen_id::Int;
    xlabel::String="Timestep",
    ylabel::String="Cyclone Mortality [%]",
    title::String="Cyclone Mortality Scenario $(scen_id)",
    kwargs...
)::PlotlyBase.Plot
    # Materialise and slice: (timesteps × locations × species × scenarios) → (T × L × Sp)
    mat_full = collect(cyclone_scens) .* 100.0
    scen_mat = mat_full[:, :, :, scen_id]            # (T × L × Sp)
    taxa_mean = dropdims(mean(scen_mat; dims=3); dims=3)  # (T × L) — avg over species
    mean_vals = vec(mean(taxa_mean; dims=2))          # (T,) — avg over locations

    x_vals = collect(cyclone_scens.timesteps)         # Dim → plain Vector
    n_locs = size(taxa_mean, 2)

    traces = PlotlyBase.AbstractTrace[]

    # Per-location traces (semi-transparent)
    for i = 1:n_locs
        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_vals, y=taxa_mean[:, i],
                mode="lines",
                line=PlotlyBase.attr(; color="rgba(220,50,50,0.12)", width=1),
                showlegend=false, hoverinfo="skip",
                name="loc_$(i)", type="scatter"
            )
        )
    end

    # Cross-location mean
    push!(
        traces,
        PlotlyBase.scatter(;
            x=x_vals, y=mean_vals,
            mode="lines",
            line=PlotlyBase.attr(; color="rgba(30,30,30,0.85)", width=2),
            name="Scenario Mean", type="scatter"
        )
    )

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text=xlabel),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end

"""
    ADRIA.viz.cyclone_scenario(dom::Domain, scen_id::Int; kwargs...)

Convenience dispatch: extracts `dom.cyclone_mortality_scens` and delegates.
"""
function ADRIA.viz.cyclone_scenario(
    dom::Domain, scen_id::Int; kwargs...
)::PlotlyBase.Plot
    return ADRIA.viz.cyclone_scenario(
        dom.cyclone_mortality_scens, scen_id; kwargs...
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# dhw_scenario (single scenario)
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.dhw_scenario(dhw_scens::YAXArray, scen_id::Int; timeframe=nothing, kwargs...)
    ADRIA.viz.dhw_scenario(dom::Domain, scen_id::Int; kwargs...)

Plot per-location DHW trajectory for a single scenario.
`dhw_scens` is expected to have dimensions `(timesteps × sites × scenarios)`.
Each location is drawn as a semi-transparent line; the cross-location mean is overlaid.

# Keyword arguments
- `timeframe`  : optional vector of calendar years for the x-axis; falls back to timestep integers
- `xlabel`     : x-axis label (default `"Year"`)
- `ylabel`     : y-axis label (default `"DHW"`)
- `title`      : plot title (default `"DHW Scenario <scen_id>"`)
"""
function ADRIA.viz.dhw_scenario(
    dhw_scens::YAXArray,
    scen_id::Int;
    timeframe=nothing,
    xlabel::String="Year",
    ylabel::String="DHW",
    title::String="DHW Scenario $(scen_id)",
    kwargs...
)::PlotlyBase.Plot
    mat_all = collect(dhw_scens)             # (T × S × Scens) plain Array
    loc_scens = mat_all[:, :, scen_id]       # (T × S) for selected scenario
    mean_dhw = vec(mean(loc_scens; dims=2))  # (T,) avg over sites

    x_vals = isnothing(timeframe) ? collect(dhw_scens.timesteps) : collect(timeframe)
    n_sites = size(loc_scens, 2)

    traces = PlotlyBase.AbstractTrace[]

    # Per-site traces (semi-transparent)
    for i = 1:n_sites
        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_vals, y=loc_scens[:, i],
                mode="lines",
                line=PlotlyBase.attr(; color="rgba(220,50,50,0.12)", width=1),
                showlegend=false, hoverinfo="skip",
                name="site_$(i)", type="scatter"
            )
        )
    end

    # Cross-site mean
    push!(
        traces,
        PlotlyBase.scatter(;
            x=x_vals, y=mean_dhw,
            mode="lines",
            line=PlotlyBase.attr(; color="rgba(30,30,30,0.85)", width=2),
            name="Scenario Mean", type="scatter"
        )
    )

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text=xlabel),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end

"""
    ADRIA.viz.dhw_scenario(dom::Domain, scen_id::Int; kwargs...)

Convenience dispatch: extracts `dom.dhw_scens` and `dom.env_layer_md.timeframe`.
"""
function ADRIA.viz.dhw_scenario(
    dom::Domain, scen_id::Int; kwargs...
)::PlotlyBase.Plot
    return ADRIA.viz.dhw_scenario(
        dom.dhw_scens, scen_id; timeframe=dom.env_layer_md.timeframe, kwargs...
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# dhw_scenarios (all scenarios)
# ──────────────────────────────────────────────────────────────────────────────

"""
    ADRIA.viz.dhw_scenarios(dhw_scens::YAXArray; ci_level=0.95, timeframe=nothing, kwargs...)
    ADRIA.viz.dhw_scenarios(dom::Domain; ci_level=0.95, kwargs...)

Plot all DHW scenarios with a confidence-interval band and mean trajectory.
`dhw_scens` is expected to have dimensions `(timesteps × sites × scenarios)`.

# Keyword arguments
- `ci_level`   : confidence interval level (default `0.95`)
- `timeframe`  : optional vector of calendar years; falls back to timestep integers
- `xlabel`     : x-axis label (default `"Year"`)
- `ylabel`     : y-axis label (default `"DHW"`)
- `title`      : plot title (default `"DHW Scenarios"`)
"""
function ADRIA.viz.dhw_scenarios(
    dhw_scens::YAXArray;
    ci_level::Real=0.95,
    timeframe=nothing,
    xlabel::String="Year",
    ylabel::String="DHW",
    title::String="DHW Scenarios",
    kwargs...
)::PlotlyBase.Plot
    mat_all = collect(dhw_scens)                               # (T × S × Scens)
    n_timesteps = size(mat_all, 1)
    n_scens = size(mat_all, 3)

    # Mean over sites for each (timestep, scenario)
    site_means = dropdims(mean(mat_all; dims=2); dims=2)       # (T × Scens)

    α = (1 - ci_level) / 2
    means = vec(mean(site_means; dims=2))                      # (T,)
    lower_ci = [quantile(site_means[t, :], α) for t = 1:n_timesteps]
    upper_ci = [quantile(site_means[t, :], 1 - α) for t = 1:n_timesteps]

    x_vals = isnothing(timeframe) ? collect(dhw_scens.timesteps) : collect(timeframe)
    ci_pct = floor(Int, ci_level * 100)

    traces = PlotlyBase.AbstractTrace[]

    # Individual scenario trajectories (semi-transparent)
    for i = 1:n_scens
        push!(
            traces,
            PlotlyBase.scatter(;
                x=x_vals, y=site_means[:, i],
                mode="lines",
                line=PlotlyBase.attr(; color="rgba(220,50,50,0.2)", width=1),
                showlegend=false, hoverinfo="skip",
                name="scen_$(i)", type="scatter"
            )
        )
    end

    # CI band
    fill_x = vcat(x_vals, reverse(x_vals))
    fill_y = vcat(upper_ci, reverse(lower_ci))
    push!(
        traces,
        PlotlyBase.scatter(;
            x=fill_x, y=fill_y,
            fill="toself", fillcolor="rgba(220,50,50,0.3)",
            line_color="rgba(0,0,0,0)",
            showlegend=true, hoverinfo="skip",
            name="$(ci_pct)% CI", type="scatter"
        )
    )

    # Overall mean
    push!(
        traces,
        PlotlyBase.scatter(;
            x=x_vals, y=means,
            mode="lines",
            line=PlotlyBase.attr(; color="rgba(30,30,30,0.85)", width=2),
            name="Mean", type="scatter"
        )
    )

    layout = PlotlyBase.Layout(;
        ADRIA_LAYOUT_DEFAULTS...,
        title_text=title,
        xaxis=PlotlyBase.attr(; title_text=xlabel),
        yaxis=PlotlyBase.attr(; title_text=ylabel)
    )
    return PlotlyBase.Plot(traces, layout)
end

"""
    ADRIA.viz.dhw_scenarios(dom::Domain; kwargs...)

Convenience dispatch: extracts `dom.dhw_scens` and `dom.env_layer_md.timeframe`.
"""
function ADRIA.viz.dhw_scenarios(
    dom::Domain; kwargs...
)::PlotlyBase.Plot
    return ADRIA.viz.dhw_scenarios(
        dom.dhw_scens; timeframe=dom.env_layer_md.timeframe, kwargs...
    )
end
