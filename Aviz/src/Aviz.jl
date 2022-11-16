module Aviz

using RelocatableFolders, FileIO
using GLMakie, GeoMakie
using GLMakie.GeometryBasics
using Statistics, Distributions, ThreadsX, Random

using DataFrames, Bootstrap, DecisionTree

using GeoInterface
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import GeoMakie.GeoJSON.FeatureCollection as FC

using ADRIA


Random.seed!(101)

const ASSETS = @path joinpath(@__DIR__, "../assets")
const LOGO = @path joinpath(ASSETS, "imgs", "ADRIA_logo.png")


include("./plotting.jl")
include("./layout.jl")
include("./theme.jl")
include("./spatial.jl")
include("./rf_analysis.jl")
include("./analysis.jl")


"""Main entry point for app."""
function julia_main()::Cint
    if "explore" in ARGS
        rs_pkg = ARGS[2]
        explore(rs_pkg)
        return 0
    end

    if "run" in ARGS
        domain_loc = ARGS[2]
        rcp_id = ARGS[3]
        input_set = ARGS[4]

        dom = ADRIA.load_domain(domain_loc, rcp_id)
        p_df = ADRIA.load_scenarios(dom, input_set)

        dom = ADRIA.run_scenarios(p_df, dom)
        rs = ADRIA.load_results(dom)
        explore(rs)
        return 0
    end

    main_menu()

    return 0
end

function main_menu()
    f = Figure()

    logo = image(
        f[1, 1],
        rotr90(load(convert(String, LOGO))),
        axis=(aspect=DataAspect(),)
    )

    hidedecorations!(f.content[1])
    hidespines!(f.content[1])

    Label(f[2, 1], "Enter ADRIA Result Set to analyze")
    rs_path_tb = Textbox(f[3, 1], placeholder="./Moore_RS")
    rs_path_tb.stored_string[] = "./Moore_RS"
    status_label = Label(f[4, 1], "")

    launch_button = Button(f[5, 1], label="Explore")

    on(launch_button.clicks) do c
        rs_path = rs_path_tb.stored_string[]
        if !isnothing(rs_path) && ispath(rs_path)
            status_label.text[] = "Loading Result Set..."
            rs = nothing
            try
                rs = ADRIA.load_results(rs_path)
            catch
                rs_path_tb.bordercolor = :red
                status_label.text[] = "Invalid ADRIA Result Set"
            else
                # Clear current figure and launch new display
                empty!(f)
                explore(rs)
            end
        else
            rs_path_tb.bordercolor = :red
            status_label.text[] = "Invalid path"
        end
    end

    gl_screen = display(f)
    wait(gl_screen)
end


function _get_seeded_sites(seed_log, ts, scens; N=5)
    t = dropdims(sum(seed_log[timesteps=ts, scenarios=scens], dims=:timesteps), dims=:timesteps)
    site_scores = dropdims(sum(t, dims=:scenarios), dims=:scenarios)

    # @info "Scores", site_scores
    if length(unique(site_scores)) == 1
        return zeros(Int64, N)
    end

    return sortperm(site_scores)[1:N]
end


function comms(rs::ADRIA.ResultSet)
    layout = comms_layout(resolution=(1920, 1080))

    f = layout.figure
    # controls = layout.controls
    traj_display = layout.trajectory.temporal
    traj_outcome_sld = layout.trajectory.outcome_slider
    traj_time_sld = layout.trajectory.time_slider

    # Generate trajectory
    tac_scens = ADRIA.metrics.scenario_total_cover(rs)
    tac_data = Matrix(tac_scens')
    tac_min_max = (minimum(tac_scens), maximum(tac_scens))

    # Generate trajectory controls
    num_steps = Int(ceil((tac_min_max[2] - tac_min_max[1]) + 1))
    tac_slider = IntervalSlider(traj_outcome_sld[2, 1],
        range=LinRange(floor(Int64, tac_min_max[1]) - 1, ceil(Int64, tac_min_max[2]) + 1, num_steps),
        startvalues=tac_min_max,
        horizontal=false
    )

    # Dynamic label text for TAC slider
    tac_bot_val = Observable(floor(tac_min_max[1]) - 1)
    tac_top_val = Observable(ceil(tac_min_max[2]) + 1)
    Label(traj_outcome_sld[1, 1], @lift("$(round($tac_top_val / 1e6, digits=2)) M (m²)"))
    Label(traj_outcome_sld[3, 1], @lift("$(round($tac_bot_val / 1e6, digits=2)) M (m²)"))

    # Time slider
    years = timesteps(rs)
    year_range = first(years), last(years)
    time_slider = IntervalSlider(
        traj_time_sld[1, 2:3],
        range=LinRange(year_range[1], year_range[2], (year_range[2] - year_range[1]) + 1),
        startvalues=year_range
    )

    # Dynamic label text for TAC slider
    left_year_val = Observable("$(year_range[1])")
    right_year_val = Observable("$(year_range[2])")
    Label(traj_time_sld[1, 1], left_year_val)
    Label(traj_time_sld[1, 4], right_year_val)

    # Generate map
    map_display = layout.map

    # Get bounds to display
    centroids = rs.site_centroids
    mean_rc_sites = ADRIA.metrics.relative_cover(rs)
    obs_mean_rc_sites = Observable(vec(mean(mean_rc_sites, dims=(:scenarios, :timesteps))))

    # Placeholder store to control which trajectories are visible
    color_map = scenario_colors(rs)
    obs_color = Observable(color_map)
    scen_types = scenario_type(rs)

    seed_log = rs.seed_log[:, 1, :, :]

    # Trajectories
    series!(traj_display, years, tac_data, color=obs_color)

    # Density (TODO: Separate into own function)
    tac_scen_dist = dropdims(mean(tac_scens, dims=:timesteps), dims=:timesteps)
    obs_cf_scen_dist = Observable(tac_scen_dist[scen_types.counterfactual])
    obs_ug_scen_dist = Observable(tac_scen_dist[scen_types.unguided])
    obs_g_scen_dist = Observable(tac_scen_dist[scen_types.guided])

    # Color transparency for density plots
    # Note: Density plots currently cannot handle empty datasets
    #       as what might happen if user selects a region with no results.
    #       so instead we set alpha to 0.0 to hide it.
    cf_hist_alpha = Observable(0.5)
    ug_hist_alpha = Observable(0.5)
    g_hist_alpha = Observable(0.5)

    # Legend(traj_display)  legend=["Counterfactual", "Unguided", "Guided"]
    scen_hist = layout.scen_hist
    density!(scen_hist, obs_cf_scen_dist, direction=:y, color=(:red, cf_hist_alpha))
    density!(scen_hist, obs_ug_scen_dist, direction=:y, color=(:green, ug_hist_alpha))
    density!(scen_hist, obs_g_scen_dist, direction=:y, color=(:blue, g_hist_alpha))
    hidedecorations!(scen_hist)
    hidespines!(scen_hist)

    # Prep seed log
    seed_rank_log = rs.ranks[:, :, 1, :]
    n_sites = size(rs.site_data, 1)

    # Random forest stuff
    # Feature importance
    # layout.outcomes
    ft_import = layout.importance

    # https://github.dev/JuliaAI/DecisionTree.jl
    # X = Matrix(rs.inputs)
    X = rs.inputs
    p = outcome_probability(tac_scen_dist)
    # model = build_forest(p, X, ceil(Int, sqrt(size(X, 1))), 30, 0.7, -1; rng=101)
    # p_tbl = probability_table(model, X, p)
    # @time ft_tbl = ft_importance(model, rs.inputs, p; rng=101)

    asv_scens = ADRIA.metrics.scenario_asv(rs)
    asv_scen_dist = dropdims(mean(asv_scens, dims=:timesteps), dims=:timesteps)

    interv_names = ADRIA.component_params(rs.model_spec, Intervention).fieldname
    interv_idx = findall(x -> x in interv_names, names(rs.inputs))
    @time begin
        mean_tac_med = relative_sensitivities(X, Array(tac_scen_dist))
        mean_tac_med = mean_tac_med[interv_idx]

        mean_asv_med = relative_sensitivities(X, Array(asv_scen_dist))
        mean_asv_med = mean_asv_med[interv_idx]
    end

    sample_cv = std(tac_scen_dist) ./ mean(tac_scen_dist)
    cf_cv = std(tac_scen_dist[scen_types.counterfactual]) ./ mean(tac_scen_dist[scen_types.counterfactual])
    ug_cv = std(tac_scen_dist[scen_types.unguided]) ./ mean(tac_scen_dist[scen_types.unguided])
    g_cv = std(tac_scen_dist[scen_types.guided]) ./ mean(tac_scen_dist[scen_types.guided])

    ft_import = Axis(
        layout.importance[1, 1],
        xticks=([1, 2], ["Mean TAC", "Mean ASV"]),
        yticks=(1:length(interv_names), interv_names)
    )
    ft_import.yreversed = true

    sensitivities = Observable(hcat(mean_tac_med, mean_asv_med)')
    heatmap!(ft_import, sensitivities)
    Colorbar(layout.importance[1, 2]; colorrange=(0.0, 1.0))

    outcomes_ax = layout.outcomes
    # barplot!(
    #     outcomes_ax,
    #     eachindex(p_tbl.mean),
    #     p_tbl.mean,
    #     bar_labels=p_tbl.Outcome,
    #     flip_labels_at=maximum(p_tbl.mean) - minimum(p_tbl.mean),
    #     direction=:x
    # )
    hideydecorations!(outcomes_ax)

    # xticks!(outcomes, 1:size(p_tbl,1), p_tbl[:, :Outcome])
    # ytick=(1:size(p_tbl,1), p_tbl[:, :Outcome])
    # xaxis="Probability", yaxis="TAC Outcome", legend=false
    # xerror!(a1, p_tbl[:, :mean], 1:nrow(p_tbl), xerror=(p_tbl[:, :lower], p_tbl[:, :upper]))

    # tac_ts = vec(mean(tac_outcomes[timesteps=50:60], dims=:timesteps)')
    # tac_dist = fit(Normal, tac_ts)
    # p_tac_outcomes = cdf.(tac_dist, tac_ts)


    # TODO: Separate this out into own function
    # Make temporary copy of GeoPackage as GeoJSON
    tmpdir = mktempdir()

    local geo_fn = joinpath(tmpdir, "Aviz_$(rs.name).geojson")
    try
        GDF.write(geo_fn, rs.site_data; driver="geojson")
    catch
        GDF.write(geo_fn, rs.site_data; geom_columns=(:geom,), driver="geojson")
    end
    geodata = GeoMakie.GeoJSON.read(read(geo_fn))

    map_disp = create_map!(map_display, geodata, obs_mean_rc_sites, (:black, 0.05), centroids)
    curr_highlighted_sites = _get_seeded_sites(seed_log, (:), (:))

    obs_site_sel = FC(geodata[curr_highlighted_sites, :][:])
    obs_site_sel = Observable(obs_site_sel)
    obs_site_highlight = Observable((:green, 1.0))
    overlay_site = poly!(map_disp, obs_site_sel, color=(:white, 0.0), strokecolor=obs_site_highlight, strokewidth=0.75, overdraw=true)

    onany(time_slider.interval, tac_slider.interval) do time_val, tac_val
        # Update slider labels
        left_year_val[] = "$(Int(floor(time_val[1])))"
        right_year_val[] = "$(Int(ceil(time_val[2])))"
        tac_bot_val[] = tac_val[1]
        tac_top_val[] = tac_val[2]

        # Convert time ranges to index values
        timespan = floor(Int, time_val[1] - (year_range[1]) + 1):ceil(Int, time_val[2] - (year_range[1]) + 1)
        hide_idx = vec(all((tac_val[1] .<= tac_scens[timespan, :] .<= tac_val[2]) .== 0, dims=1))
        show_idx = Bool.(ones(Int64, length(hide_idx)) .⊻ hide_idx)  # inverse of hide

        scen_dist = dropdims(mean(tac_scens[timespan, :], dims=:timesteps), dims=:timesteps)

        # Update map
        obs_mean_rc_sites[] = vec(mean(mean_rc_sites[timesteps=timespan][scenarios=findall(show_idx)], dims=(:scenarios, :timesteps)))

        seeded_sites = _get_seeded_sites(seed_log, (:), show_idx)
        if seeded_sites != curr_highlighted_sites
            # Highlight seeded sites
            if !all(seeded_sites .== 0.0) && !all(show_idx .== 0)
                obs_site_sel[] = FC(geodata[seeded_sites, :][:])
                obs_site_highlight[] = (:green, 1.0)
                curr_highlighted_sites .= seeded_sites
            else
                obs_site_highlight[] = (:green, 0.0)
            end
        else
            if all(seeded_sites .== 0.0) || all(show_idx .== 0)
                obs_site_highlight[] = (:green, 0.0)
            end
        end

        # Boolean index of scenarios to hide (inverse of tac_idx)
        if !all(hide_idx .== 0)
            # Hide scenarios that were filtered out
            cf_dist = scen_dist[show_idx.&scen_types.counterfactual]
            ug_dist = scen_dist[show_idx.&scen_types.unguided]
            g_dist = scen_dist[show_idx.&scen_types.guided]
        else
            cf_dist = scen_dist[scen_types.counterfactual]
            ug_dist = scen_dist[scen_types.unguided]
            g_dist = scen_dist[scen_types.guided]
        end

        # Update scenario density plot
        if !isempty(cf_dist)
            obs_cf_scen_dist[] = cf_dist
            cf_hist_alpha[] = 0.5
        else
            cf_hist_alpha[] = 0.0
        end

        if !isempty(ug_dist)
            obs_ug_scen_dist[] = ug_dist
            ug_hist_alpha[] = 0.5
        else
            ug_hist_alpha[] = 0.0
        end

        if !isempty(g_dist)
            obs_g_scen_dist[] = g_dist
            g_hist_alpha[] = 0.5
        else
            g_hist_alpha[] = 0.0
        end

        # Determine level of transparency for each line (maximum of 0.6)
        min_step = (1.0 / 0.05)
        color_weight = min((1.0 / (count(show_idx .> 0) / min_step)), 0.6)

        # Update visible trajectories
        scenario_colors!(obs_color, scen_types, color_weight, hide_idx)

        # Update sensitivities
        if !all(show_idx .== 0)
            sel_tac_scens = dropdims(mean(tac_scens[timesteps=timespan, scenarios=show_idx], dims=:timesteps), dims=:timesteps)
            mean_tac_med = relative_sensitivities(X[show_idx, :], Array(sel_tac_scens))
            mean_tac_med = mean_tac_med[interv_idx]

            sel_asv_scens = dropdims(mean(asv_scens[timesteps=timespan, scenarios=show_idx], dims=:timesteps), dims=:timesteps)
            mean_asv_med = relative_sensitivities(X[show_idx, :], Array(sel_asv_scens))
            mean_asv_med = mean_asv_med[interv_idx]
        else
            # Display nothing if no data is available
            mean_tac_med = fill(NaN, length(interv_idx))
            mean_asv_med = fill(NaN, length(interv_idx))
        end

        sensitivities[] = hcat(mean_tac_med, mean_asv_med)'
    end


    gl_screen = display(f)
    DataInspector()

    wait(gl_screen)
end
function comms(rs_path::String)
    comms(ADRIA.load_results(rs_path))
end


function explore(rs::ADRIA.ResultSet)
    layout = modeler_layout(resolution=(1920, 1080))

    f = layout.figure
    # controls = layout.controls
    traj_display = layout.trajectory.temporal
    traj_outcome_sld = layout.trajectory.outcome_slider
    traj_time_sld = layout.trajectory.time_slider

    scen_hist = layout.scen_hist
    map_display = layout.map

    interv_pcp_display = layout.interv_pcp
    pair_display = layout.pairplot
    outcome_pcp_display = layout.outcome_pcp

    colsize!(f.layout, 1, Fixed(400))
    colsize!(f.layout, 2, Fixed(1100))

    color_map = scenario_colors(rs)
    obs_color = Observable(color_map)

    n_visible_scenarios = Observable(size(rs.inputs, 1))

    seed_log = rs.seed_log[:, 1, :, :]

    # Controls
    ms = rs.model_spec
    # intervention_components = ms[ms.component .== "Intervention", [:fieldname, :full_bounds]]

    tac_scens = ADRIA.metrics.scenario_total_cover(rs)
    # rc_scens = ADRIA.metrics.scenario_relative_cover(rs)
    mean_tac_outcomes = vec(mean(tac_scens, dims=1))
    # mean_rc_outcomes = vec(mean(rc_scens, dims=1))
    tac_min_max = (minimum(tac_scens), maximum(tac_scens))

    # tac_label = Label(traj_outcome_sld[1,1], "Mean TAC (m²)")  # , rotation = pi/2
    num_steps = Int(ceil((tac_min_max[2] - tac_min_max[1]) + 1))
    tac_slider = IntervalSlider(traj_outcome_sld[2, 1],
        range=LinRange(floor(Int64, tac_min_max[1]) - 1, ceil(Int64, tac_min_max[2]) + 1, num_steps),
        startvalues=tac_min_max,
        horizontal=false
        # width=350
    )

    # Dynamic label text for TAC slider
    tac_bot_val = Observable(floor(tac_min_max[1]) - 1)
    tac_top_val = Observable(ceil(tac_min_max[2]) + 1)

    tac_bot = @lift("$(round($tac_bot_val / 1e6, digits=2))")
    tac_top = @lift("$(round($tac_top_val / 1e6, digits=2))")
    Label(traj_outcome_sld[1, 1], tac_top)
    Label(traj_outcome_sld[3, 1], tac_bot)

    # Time slider
    years = timesteps(rs)
    year_range = first(years), last(years)
    time_slider = IntervalSlider(
        traj_time_sld[1, 2],
        range=LinRange(year_range[1], year_range[2], (year_range[2] - year_range[1]) + 1),
        startvalues=year_range
    )

    # Dynamic label text for TAC slider
    left_year_val = Observable("$(year_range[1])")
    right_year_val = Observable("$(year_range[2])")
    Label(traj_time_sld[1, 1], left_year_val)
    Label(traj_time_sld[1, 3], right_year_val)

    tac = ADRIA.metrics.scenario_total_cover(rs)
    tac_data = Matrix(tac')

    asv = ADRIA.metrics.scenario_asv(rs)

    # Histogram/Density plot
    scen_types = scenario_type(rs)

    scen_dist = dropdims(mean(tac, dims=:timesteps), dims=:timesteps)
    # scen_dist = vec(mean(scen_dist, dims=:sites))

    cf_scen_dist = scen_dist[scen_types.counterfactual]
    ug_scen_dist = scen_dist[scen_types.unguided]
    g_scen_dist = scen_dist[scen_types.guided]

    obs_cf_scen_dist = Observable(cf_scen_dist)
    obs_ug_scen_dist = Observable(ug_scen_dist)
    obs_g_scen_dist = Observable(g_scen_dist)

    # Color transparency for density plots
    # Note: Density plots currently cannot handle empty datasets
    #       as what might happen if user selects a region with no results.
    #       so instead we set alpha to 0.0 to hide it.
    cf_hist_alpha = Observable(0.5)
    ug_hist_alpha = Observable(0.5)
    g_hist_alpha = Observable(0.5)

    # Get intervention/criteria inputs for each scenario
    interv_criteria = ms[(ms.component.=="EnvironmentalLayer").|(ms.component.=="Intervention").|(ms.component.=="Criteria"), [:fieldname, :full_bounds]]
    input_names = vcat(["RCP", interv_criteria.fieldname...])
    in_pcp_data = normalize(Matrix(rs.inputs[:, input_names]))
    # in_pcp_lines = Observable(in_pcp_data)


    # Get mean outcomes for each scenario
    outcome_pcp_data = hcat([
        mean_tac_outcomes,
        vec(mean(ADRIA.metrics.scenario_asv(rs), dims=1)),
        vec(mean(ADRIA.metrics.scenario_rsv(rs), dims=1))
    ]...)
    disp_names = ["TAC", "ASV", "RSV"]

    out_pcp_data = normalize(outcome_pcp_data)

    # Specify interactive elements and behavior
    # TODO: Lift on data to be plotted
    #       Controls simply update transparency settings etc, and update the dataset to be plotted
    #       All other elements simply update when the underlying dataset updates
    # https://discourse.julialang.org/t/interactive-plot-with-makielayout/48843
    onany(time_slider.interval, tac_slider.interval) do time_val, tac_val
        # Update slider labels
        tac_bot_val[] = tac_val[1]
        tac_top_val[] = tac_val[2]

        left_year_val[] = "$(Int(floor(time_val[1])))"
        right_year_val[] = "$(Int(ceil(time_val[2])))"

        # Trajectories
        # tac_idx = (mean_tac_outcomes .>= tac_val[1]-0.5) .& (mean_tac_outcomes .<= tac_val[2]+0.5)

        # Convert time ranges to index values
        t_idx = Int(time_val[1] - (year_range[1]) + 1), Int(time_val[2] - (year_range[1]) + 1)

        hide_idx = vec(all((tac_val[1] .<= tac_scens[t_idx[1]:t_idx[2], :] .<= tac_val[2]) .== 0, dims=1))
        show_idx = Bool.(zeros(Int64, length(hide_idx)) .⊻ hide_idx)  # inverse of hide

        scen_dist = dropdims(mean(tac[timesteps=t_idx[1]:t_idx[2]], dims=:timesteps), dims=:timesteps)

        # Boolean index of scenarios to hide (inverse of tac_idx)
        # hide_idx = Bool.(ones(Int64, length(tac_idx)) .⊻ tac_idx)
        if !all(hide_idx .== 0)
            # Hide scenarios that were filtered out
            cf_dist = scen_dist[show_idx.&scen_types.counterfactual]
            ug_dist = scen_dist[show_idx.&scen_types.unguided]
            g_dist = scen_dist[show_idx.&scen_types.guided]
        else
            cf_dist = scen_dist[scen_types.counterfactual]
            ug_dist = scen_dist[scen_types.unguided]
            g_dist = scen_dist[scen_types.guided]
        end

        # Update scenario density plot
        if !isempty(cf_dist)
            obs_cf_scen_dist[] = cf_dist
            cf_hist_alpha[] = 0.5
        else
            cf_hist_alpha[] = 0.0
        end

        if !isempty(ug_dist)
            obs_ug_scen_dist[] = ug_dist
            ug_hist_alpha[] = 0.5
        else
            ug_hist_alpha[] = 0.0
        end

        if !isempty(g_dist)
            obs_g_scen_dist[] = g_dist
            g_hist_alpha[] = 0.5
        else
            g_hist_alpha[] = 0.0
        end

        # Determine level of transparency for each line (maximum of 0.6)
        min_step = (1 / 0.05)
        color_weight = min((1.0 / (count(show_idx .> 0) / min_step)), 0.6)

        obs_color[] = scenario_colors(rs, color_weight, hide_idx)
    end

    # Trajectories
    series!(traj_display, timesteps(rs), tac_data, color=@lift($obs_color[:]))

    # Legend(traj_display)  legend=["Counterfactual", "Unguided", "Guided"]
    density!(scen_hist, @lift($obs_cf_scen_dist[:]), direction=:y, color=(:red, @lift($cf_hist_alpha[])))
    density!(scen_hist, @lift($obs_ug_scen_dist[:]), direction=:y, color=(:green, @lift($ug_hist_alpha[])))
    density!(scen_hist, @lift($obs_g_scen_dist[:]), direction=:y, color=(:blue, @lift($g_hist_alpha[])))

    hidedecorations!(scen_hist)
    hidespines!(scen_hist)

    # TODO: Separate this out into own function
    # Make temporary copy of GeoPackage as GeoJSON
    tmpdir = mktempdir()
    geo_fn = GDF.write(joinpath(tmpdir, "Aviz_$(rs.name).geojson"), rs.site_data; driver="geojson")
    geodata = GeoMakie.GeoJSON.read(read(geo_fn))

    # Get bounds to display
    centroids = rs.site_centroids
    lon = [c[1] for c in centroids]
    lat = [c[2] for c in centroids]

    # Display map
    mean_rc_sites = mean(ADRIA.metrics.relative_cover(rs), dims=(:scenarios, :timesteps))
    map_buffer = 0.005
    spatial = GeoAxis(
        map_display;  # any cell of the figure's layout
        lonlims=(minimum(lon) - map_buffer, maximum(lon) + map_buffer),
        latlims=(minimum(lat) - map_buffer, maximum(lat) + map_buffer),
        xlabel="Long",
        ylabel="Lat"
    )
    # datalims!(spatial)  # auto-adjust limits (doesn't work if there are Infs...)
    poly!(spatial, geodata, color=vec(mean_rc_sites), colormap=:plasma)

    # Fill pairplot
    # Get mean outcomes for each scenario
    pairplot!(pair_display, outcome_pcp_data, disp_names)

    # Parallel Coordinate Plot
    pcp!(interv_pcp_display, in_pcp_data, input_names; color=@lift($obs_color[:]))
    pcp!(outcome_pcp_display, out_pcp_data, disp_names; color=@lift($obs_color[:]))

    gl_screen = display(f)
    DataInspector()

    wait(gl_screen)
end
function explore(rs_path::String)
    explore(ADRIA.load_results(rs_path))
end

end


# Allow use from terminal if this file is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    if "explore" in ARGS
        rs_pkg = ARGS[2]
        Aviz.explore(rs_pkg)
    end
end
