module AvizExt

using Base.Iterators
using Reexport

using RelocatableFolders
@reexport using GeoMakie

using Statistics, Distributions, FLoops, Random
using DataFrames, Bootstrap

using ImageIO, GeoInterface

import GeoMakie.GeoJSON.FeatureCollection as FC

import ADRIA.FileIO, ADRIA.GFT
using ADRIA
using ADRIA:
    load_results, load_domain, load_scenarios,
    ResultSet, run_scenarios, metrics, viz.explore, GDF

import ADRIA: timesteps as AD_timesteps
import ADRIA.viz: explore

Random.seed!(101)

const ASSETS = @path joinpath(pkgdir(ADRIA), "assets")
const LOGO = @path joinpath(ASSETS, "imgs", "ADRIA_logo.png")
const LOADER = @path joinpath(ASSETS, "imgs", "ADRIA_loader.gif")

include("./plotting.jl")
include("./layout.jl")
include("./theme.jl")
include("./analysis.jl")
include("./viz/viz.jl")

"""Main entry point for app."""
function julia_main()::Cint
    if "explore" in ARGS
        rs_pkg = ARGS[2]
        ADRIA.viz.explore(rs_pkg)
        return 0
    end

    if "run" in ARGS
        domain_loc = ARGS[2]
        rcp_id = ARGS[3]
        input_set = ARGS[4]

        dom = load_domain(domain_loc, rcp_id)
        p_df = load_scenarios(dom, input_set)

        rs = run_scenarios(dom, p_df)
        ADRIA.viz.explore(rs)
        return 0
    end

    main_menu()

    return 0
end

function main_menu()
    f = Figure()

    logo = image(
        f[1, 1],
        rotr90(load(convert(String, LOGO)));
        axis=(aspect=DataAspect(),)
    )

    hidedecorations!(f.content[1])
    hidespines!(f.content[1])

    Label(f[2, 1], "Enter ADRIA Result Set to analyze")
    rs_path_tb = Textbox(f[3, 1]; placeholder="./Moore_RS")
    rs_path_tb.stored_string[] = "./Moore_RS"
    status_label = Label(f[4, 1], "")

    launch_button = Button(f[5, 1]; label="Explore")

    on(launch_button.clicks) do c
        rs_path = rs_path_tb.stored_string[]
        if !isnothing(rs_path) && ispath(rs_path)
            status_label.text[] = "Loading Result Set..."
            rs = nothing
            try
                rs = load_results(rs_path)
            catch
                rs_path_tb.bordercolor = :red
                status_label.text[] = "Invalid ADRIA Result Set"
            else
                # Clear current figure and launch new display
                empty!(f)
                ADRIA.explore(rs)
            end
        else
            rs_path_tb.bordercolor = :red
            status_label.text[] = "Invalid path"
        end
    end

    gl_screen = display(f)
    return wait(gl_screen)
end

function _get_seeded_sites(seed_log, ts, scens; N=10)
    t = dropdims(
        sum(seed_log[timesteps=ts, scenarios=scens]; dims=:timesteps); dims=:timesteps
    )
    site_scores = dropdims(sum(t; dims=:scenarios); dims=:scenarios)

    # @info "Scores", site_scores
    if length(unique(site_scores)) == 1
        return zeros(Int64, N)
    end

    return sortperm(site_scores)[1:N]
end

function display_loader(fig, anim)
    a = image(fig[1, 1], anim[:, :, 1])
    hidedecorations!(a.axis)
    hidespines!(a.axis)

    for i in cycle(axes(anim, 3))
        image!(a.axis, anim[:, :, i])

        sleep(0.1)
    end
end
function remove_loader(fig, task)
    Base.throwto(task, InterruptException())
    return empty!(fig)
end

"""
    ADRIA.viz.explore(rs::String)
    ADRIA.viz.explore(rs::ResultSet)

Display GUI for quick visualization and analysis of results.
"""
function ADRIA.viz.explore(rs::ResultSet)
    @info "Creating display"
    layout = comms_layout(; size=(1920, 1080))

    f = layout.figure
    traj_display = layout.trajectory.temporal
    traj_outcome_sld = layout.trajectory.outcome_slider
    traj_time_sld = layout.trajectory.time_slider

    # Generate trajectory
    @info "Extracting core scenario outcomes"
    tac_scens = metrics.scenario_total_cover(rs)
    tac_data = Matrix(tac_scens')
    tac_min_max = (minimum(tac_scens), maximum(tac_scens))

    mean_rc_sites = metrics.relative_cover(rs)
    obs_rc = vec(mean(mean_rc_sites; dims=(:scenarios, :timesteps)))
    obs_mean_rc_sites = Observable(obs_rc)

    asv_scens = metrics.scenario_asv(rs)
    asv_scen_dist = dropdims(mean(asv_scens; dims=:timesteps); dims=:timesteps)

    juves_scens = metrics.scenario_relative_juveniles(rs)
    juves_scen_dist = dropdims(mean(juves_scens; dims=:timesteps); dims=:timesteps)

    # Generate trajectory controls
    @info "Creating controls"
    num_steps = Int(ceil((tac_min_max[2] - tac_min_max[1]) + 1))
    tac_slider = IntervalSlider(traj_outcome_sld[2, 1];
        range=LinRange(
            floor(Int64, tac_min_max[1]) - 1, ceil(Int64, tac_min_max[2]) + 1, num_steps
        ),
        startvalues=tac_min_max,
        horizontal=false
    )

    # Dynamic label text for TAC slider
    tac_bot_val = Observable(floor(tac_min_max[1]) - 1)
    tac_top_val = Observable(ceil(tac_min_max[2]) + 1)
    Label(traj_outcome_sld[1, 1], @lift("$(round($tac_top_val / 1e6, digits=2)) M (m²)"))
    Label(traj_outcome_sld[3, 1], @lift("$(round($tac_bot_val / 1e6, digits=2)) M (m²)"))

    # Time slider
    years = AD_timesteps(rs)
    year_range = first(years), last(years)
    time_slider = IntervalSlider(
        traj_time_sld[1, 2:3];
        range=LinRange(year_range[1], year_range[2], (year_range[2] - year_range[1]) + 1),
        startvalues=year_range
    )

    # Dynamic label text for TAC slider
    left_year_val = Observable("$(year_range[1])")
    right_year_val = Observable("$(year_range[2])")
    Label(traj_time_sld[1, 1], left_year_val)
    Label(traj_time_sld[1, 4], right_year_val)

    # Placeholder store to control which trajectories are visible
    X = rs.inputs
    min_color_step = (1.0 / 0.05)
    init_weight = (1.0 / (size(X, 1) / min_color_step))

    # Group scenarios by type
    scen_groups::Dict{Symbol,BitVector} = ADRIA.analysis.scenario_types(rs.inputs)
    color_map = colors(scen_groups, init_weight)
    obs_color = Observable(color_map)

    seed_log = rs.seed_log[:, 1, :, :]

    # Trajectories
    series!(traj_display, years, tac_data; color=obs_color)

    # Color transparency for density plots
    # Note: Density plots currently cannot handle empty datasets
    #       as what might happen if user selects a region with no results.
    #       so instead we set alpha to 0.0 to hide it.
    cf_hist_alpha = Observable((:red, 0.5))
    ug_hist_alpha = Observable((:green, 0.5))
    g_hist_alpha = Observable((:blue, 0.5))

    has_cf = count(scen_groups[:counterfactual]) > 0
    has_ug = count(scen_groups[:unguided]) > 0
    has_g = count(scen_groups[:guided]) > 0

    # Density (TODO: Separate into own function)
    tac_scen_dist = dropdims(mean(tac_scens; dims=:timesteps); dims=:timesteps)
    obs_cf_scen_dist = Observable(tac_scen_dist[scen_groups[:counterfactual]])

    scen_hist = layout.scen_hist
    if has_cf
        density!(scen_hist, obs_cf_scen_dist; direction=:y, color=cf_hist_alpha)
    end

    obs_ug_scen_dist = Observable(tac_scen_dist[scen_groups[:unguided]])
    if has_ug
        density!(scen_hist, obs_ug_scen_dist; direction=:y, color=ug_hist_alpha)
    end

    obs_g_scen_dist = Observable(tac_scen_dist[scen_groups[:guided]])
    if has_g
        density!(scen_hist, obs_g_scen_dist; direction=:y, color=g_hist_alpha)
    end

    hidedecorations!(scen_hist)
    hidespines!(scen_hist)
    ylims!(scen_hist, 0.0, maximum(tac_scen_dist))

    ms = rs.model_spec
    intervention_components = ms[
        (ms.component .== "Intervention") .& (ms.fieldname .!= "guided"),
        [:name, :fieldname, :lower_bound, :upper_bound]
    ]
    interv_names = intervention_components.fieldname
    interv_idx = findall(x -> x in interv_names, names(X))

    # Add control grid
    # Controls for RCPs
    has_45 = count(X.RCP .== 45) > 0
    has_60 = count(X.RCP .== 60) > 0
    has_85 = count(X.RCP .== 85) > 0
    t_toggles = [
        Toggle(f; active=active) for
        active in [has_45, has_60, has_85, has_cf, has_ug, has_g]
    ]
    t_toggle_map = zip(
        t_toggles,
        ["RCP 4.5", "RCP 6.0", "RCP 8.5", "Counterfactual", "Unguided", "Guided"],
        [:black, :black, :black, :red, :green, :blue]
    )
    labels = [
        Label(f, "$l"; color=lift(x -> x ? c : :gray, t.active)) for
        (t, l, c) in t_toggle_map
    ]
    layout.controls[1:2, 1:2] = grid!(hcat(t_toggles, labels); tellheight=false)

    # Controls for guided type
    guide_toggle_map = zip(
        t_toggles[3:end],
        ["Counterfactual", "Unguided", "Guided"],
        [:red, :green, :blue]
    )

    # Controls for interventions
    interv_sliders = IntervalSlider[]
    interv_labels = []
    lc = layout.controls[3:6, 2] = GridLayout()
    for (i, v) in enumerate(eachrow(intervention_components))
        fn = v.name

        x = (v.lower_bound, v.upper_bound)
        l1 = Observable("$(round(x[1], digits=2))")
        l2 = Observable("$(round(x[2], digits=2))")
        push!(interv_sliders,
            IntervalSlider(
                lc[i, 2];
                range=LinRange(x[1], x[2], 10),
                startvalues=(x[1], x[2])
            )
        )

        push!(interv_labels, [l1, l2])

        Label(lc[i, 2], fn)
        Label(lc[i, 1], l1)
        Label(lc[i, 3], l2)
    end

    @info "Determining initial sensitivities"
    mean_tac_med = relative_sensitivities(X, Array(tac_scen_dist))
    mean_tac_med = mean_tac_med[interv_idx]

    mean_asv_med = relative_sensitivities(X, Array(asv_scen_dist))
    mean_asv_med = mean_asv_med[interv_idx]

    mean_juves_med = relative_sensitivities(X, Array(juves_scen_dist))
    mean_juves_med = mean_juves_med[interv_idx]

    # sample_cv = std(tac_scen_dist) ./ mean(tac_scen_dist)
    # cf_cv = std(tac_scen_dist[scen_groups[:counterfactual]]) ./ mean(tac_scen_dist[scen_groups[:counterfactual]])
    # ug_cv = std(tac_scen_dist[scen_groups[:unguided]]) ./ mean(tac_scen_dist[scen_groups[:unguided]])
    # g_cv = std(tac_scen_dist[scen_groups[:guided]]) ./ mean(tac_scen_dist[scen_groups[:guided]])

    ft_import = Axis(
        layout.importance[1, 1];
        xticks=([1, 2, 3], ["Mean TAC", "Mean ASV", "Mean Juveniles"]),
        yticks=(1:length(interv_names), intervention_components.name),
        title="Relative Importance"
    )
    ft_import.yreversed = true

    S_data = hcat(mean_tac_med, mean_asv_med, mean_juves_med)'
    sensitivities = Observable(S_data)
    heatmap!(ft_import, sensitivities)
    Colorbar(layout.importance[1, 2]; colorrange=(0.0, 1.0))

    scen_dist = tac_scen_dist
    hide_idx = falses(size(X, 1))
    show_idx = trues(size(X, 1))

    outcomes_ax = layout.outcomes
    probas = Observable(outcome_probability(scen_dist))
    barplot!(
        outcomes_ax,
        @lift($(probas).values);
        bar_labels=:y,
        direction=:y,
        flip_labels_at=@lift(maximum($(probas).values) * 0.9),
        color_over_bar=:white,
        grid=false,
        xticklabelsize=12
    )
    hideydecorations!(outcomes_ax)

    @info "Generating map"
    map_display = layout.map
    geodata = _get_geoms(rs.site_data)
    map_disp = create_map!(map_display, geodata, obs_mean_rc_sites, (:black, 0.05))

    curr_highlighted_sites = _get_seeded_sites(seed_log, (:), (:))
    obs_site_sel = Observable(FC(; features=geodata[curr_highlighted_sites]))
    obs_site_highlight = Observable((:lightgreen, 1.0))
    poly!(
        map_disp[1, 1],
        obs_site_sel;
        color=(:white, 0.0),
        strokecolor=obs_site_highlight,
        strokewidth=0.75,
        overdraw=true
    )

    # Image file for loading animation
    # loader_anim = load(LOADER)

    @info "Setting up interactions"
    function update_disp(
        time_val, tac_val, rcp45, rcp60, rcp85, c_tog, u_tog, g_tog, disp_vals...
    )
        # Display loading animation
        # load_anim_display = @async display_loader(traj_display[2, 1], loader_anim)

        # Convert time ranges to index values
        timespan = time_val[1]:time_val[2]

        show_idx .= trues(size(X, 1))

        # Update according to intervention slider values
        # Hide scenarios that do not meet selections based on selected intervention values
        # disp_vals = [i1_val, i2_val, i3_val, i4_val, i5_val, i6_val, i7_val, i8_val, i9_val, i10_val, i11_val]
        for (intv, bnds) in enumerate(interv_labels)
            bnds[1][] = "$(round(disp_vals[intv][1], digits=2))"
            bnds[2][] = "$(round(disp_vals[intv][2], digits=2))"

            show_idx .=
                show_idx .& (
                    (X[:, interv_names[intv]] .>= disp_vals[intv][1]) .&
                    (X[:, interv_names[intv]] .<= disp_vals[intv][2])
                )
        end

        # Hide/display scenario types
        if c_tog
            show_idx .= show_idx .| (X.guided .== -1.0)
        else
            show_idx .= show_idx .& (X.guided .!= -1.0)
        end

        if !u_tog
            show_idx .= show_idx .& (X.guided .!= 0.0)
        end
        if !g_tog
            show_idx .= show_idx .& (X.guided .<= 0.0)
        end

        if !rcp45
            show_idx .= show_idx .& (X.RCP .!= 45)
        end

        if !rcp60
            show_idx .= show_idx .& (X.RCP .!= 60)
        end

        if !rcp85
            show_idx .= show_idx .& (X.RCP .!= 85)
        end

        # Update hidden scenarios with inverse of show
        hide_idx .= Bool.(ones(Int64, length(hide_idx)) .⊻ show_idx)

        # Update map
        obs_mean_rc_sites[] = vec(
            mean(
                mean_rc_sites(; timesteps=timespan)[scenarios=show_idx];
                dims=(:scenarios, :timesteps)
            )
        )

        seeded_sites = _get_seeded_sites(seed_log, (:), show_idx)
        site_alpha = 1.0
        if seeded_sites != curr_highlighted_sites
            # Highlight seeded sites
            if any(seeded_sites .> 0.0) && any(show_idx)
                obs_site_sel[] = FC(; features=geodata[seeded_sites])
                site_alpha = 1.0
                curr_highlighted_sites .= seeded_sites
            else
                site_alpha = 0.0
            end
        elseif all(seeded_sites .== 0.0) || all(show_idx .== 0)
            site_alpha = 0.0
        end

        obs_site_highlight[] = (:lightgreen, site_alpha)

        # Update scenario density
        scen_dist = dropdims(
            mean(tac_scens(; timesteps=timespan); dims=:timesteps); dims=:timesteps
        )
        # Hide scenarios that were filtered out
        cf_dist = scen_dist[show_idx .& scen_groups[:counterfactual]]
        ug_dist = scen_dist[show_idx .& scen_groups[:unguided]]
        g_dist = scen_dist[show_idx .& scen_groups[:guided]]

        if c_tog && !isempty(cf_dist)
            obs_cf_scen_dist[] = cf_dist
            cf_hist_alpha[] = (:red, 0.5)
        else
            cf_hist_alpha[] = (:red, 0.0)
        end

        if u_tog && !isempty(ug_dist)
            obs_ug_scen_dist[] = ug_dist
            ug_hist_alpha[] = (:green, 0.5)
        else
            ug_hist_alpha[] = (:green, 0.0)
        end

        if g_tog && !isempty(g_dist)
            obs_g_scen_dist[] = g_dist
            g_hist_alpha[] = (:blue, 0.5)
        else
            g_hist_alpha[] = (:blue, 0.0)
        end

        # Update limits of density plot
        autolimits!(scen_hist)

        # Update visible trajectories
        # Determine level of transparency for each line (maximum of 0.5)
        color_weight = min((1.0 / (count(show_idx .> 0) / min_color_step)), 0.5)
        scenario_colors!(
            obs_color, color_map, scen_groups, color_weight, hide_idx, guide_toggle_map
        )

        # Update sensitivities (if there's enough samples...)
        if count(show_idx) > 16
            mean_tac_med = relative_sensitivities(X[show_idx, :], scen_dist[show_idx])[interv_idx]

            sel_asv_scens = dropdims(
                mean(asv_scens(; timesteps=timespan)[scenarios=show_idx]; dims=:timesteps);
                dims=:timesteps
            )
            mean_asv_med = relative_sensitivities(X[show_idx, :], sel_asv_scens)[interv_idx]

            sel_juves_scens = dropdims(
                mean(
                    juves_scens(; timesteps=timespan)[scenarios=show_idx]; dims=:timesteps
                );
                dims=:timesteps
            )
            mean_juves_med = relative_sensitivities(X[show_idx, :], sel_juves_scens)[interv_idx]
        else
            # Display nothing if no data is available
            mean_tac_med = fill(NaN, length(interv_idx))
            mean_asv_med = fill(NaN, length(interv_idx))
            mean_juves_med = fill(NaN, length(interv_idx))
        end

        S_data[1, :] .= mean_tac_med
        S_data[2, :] .= mean_asv_med
        S_data[3, :] .= mean_juves_med
        sensitivities[] = S_data

        # Update bar plot of outcome probability
        probas[] = outcome_probability(scen_dist[show_idx])
        return ylims!(layout.outcomes, minimum(probas[].values), maximum(probas[].values))

        # Clear loading animation
        # remove_loader(traj_display[2, 1], loader_anim_display)
    end

    # Trigger update only after some time since last interaction
    # TODO: Add update notification (spinner animation or something...)
    # up_timer = Timer(x -> x, 0.25)
    onany(time_slider.interval, tac_slider.interval,
        [t.active for t in t_toggles]...,
        [sld.interval for sld in interv_sliders]...
    ) do time_val, tac_val, rcp45, rcp60, rcp85, c_tog, u_tog, g_tog, intervs... # i1_val, i2_val, i3_val, i4_val, i5_val, i6_val, i7_val, i8_val, i9_val, i10_val, i11_val, i12_val

        # Update slider labels
        left_year_val[] = "$(Int(floor(time_val[1])))"
        right_year_val[] = "$(Int(ceil(time_val[2])))"
        tac_bot_val[] = tac_val[1]
        tac_top_val[] = tac_val[2]

        if @isdefined up_timer
            close(up_timer)
        end
        up_timer = Timer(
            x -> update_disp(
                time_val, tac_val, rcp45, rcp60, rcp85, c_tog, u_tog, g_tog, intervs...
            ),
            2
        )
    end

    @info "Displaying UI"
    gl_screen = display(f)
    # DataInspector()

    return wait(gl_screen)
    # close(up_timer)
end
function ADRIA.viz.explore(rs_path::String)
    return ADRIA.viz.explore(load_results(rs_path))
end

# function explore(rs::ADRIA.ResultSet)
#     layout = modeler_layout(size=(1920, 1080))

#     f = layout.figure
#     # controls = layout.controls
#     traj_display = layout.trajectory.temporal
#     traj_outcome_sld = layout.trajectory.outcome_slider
#     traj_time_sld = layout.trajectory.time_slider

#     scen_hist = layout.scen_hist
#     map_display = layout.map

#     interv_pcp_display = layout.interv_pcp
#     pair_display = layout.pairplot
#     outcome_pcp_display = layout.outcome_pcp

#     colsize!(f.layout, 1, Fixed(400))
#     colsize!(f.layout, 2, Fixed(1100))

#     color_map = scenario_colors(rs)
#     obs_color = Observable(color_map)

#     # n_visible_scenarios = Observable(size(rs.inputs, 1))
#     # seed_log = rs.seed_log[:, 1, :, :]

#     # Temporal controls
#     tac_scens = ADRIA.metrics.scenario_total_cover(rs)
#     # rc_scens = ADRIA.metrics.scenario_relative_cover(rs)
#     mean_tac_outcomes = vec(mean(tac_scens, dims=1))
#     # mean_rc_outcomes = vec(mean(rc_scens, dims=1))
#     tac_min_max = (minimum(tac_scens), maximum(tac_scens))

#     # tac_label = Label(traj_outcome_sld[1,1], "Mean TAC (m²)")  # , rotation = pi/2
#     num_steps = Int(ceil((tac_min_max[2] - tac_min_max[1]) + 1))
#     tac_slider = IntervalSlider(traj_outcome_sld[2, 1],
#         range=LinRange(floor(Int64, tac_min_max[1]) - 1, ceil(Int64, tac_min_max[2]) + 1, num_steps),
#         startvalues=tac_min_max,
#         horizontal=false
#         # width=350
#     )

#     # Dynamic label text for TAC slider
#     tac_bot_val = Observable(floor(tac_min_max[1]) - 1)
#     tac_top_val = Observable(ceil(tac_min_max[2]) + 1)

#     tac_bot = @lift("$(round($tac_bot_val / 1e6, digits=2))")
#     tac_top = @lift("$(round($tac_top_val / 1e6, digits=2))")
#     Label(traj_outcome_sld[1, 1], tac_top)
#     Label(traj_outcome_sld[3, 1], tac_bot)

#     # Time slider
#     years = timesteps(rs)
#     year_range = first(years), last(years)
#     time_slider = IntervalSlider(
#         traj_time_sld[1, 2],
#         range=LinRange(year_range[1], year_range[2], (year_range[2] - year_range[1]) + 1),
#         startvalues=year_range
#     )

#     # Dynamic label text for TAC slider
#     left_year_val = Observable("$(year_range[1])")
#     right_year_val = Observable("$(year_range[2])")
#     Label(traj_time_sld[1, 1], left_year_val)
#     Label(traj_time_sld[1, 3], right_year_val)

#     tac = ADRIA.metrics.scenario_total_cover(rs)
#     tac_data = Matrix(tac')

#     # asv = ADRIA.metrics.scenario_asv(rs)

#     # Histogram/Density plot
#     scen_types = scenario_type(rs)

#     scen_dist = dropdims(mean(tac, dims=:timesteps), dims=:timesteps)
#     # scen_dist = vec(mean(scen_dist, dims=:sites))

#     cf_scen_dist = scen_dist[scen_types.counterfactual]
#     ug_scen_dist = scen_dist[scen_types.unguided]
#     g_scen_dist = scen_dist[scen_types.guided]

#     obs_cf_scen_dist = Observable(cf_scen_dist)
#     obs_ug_scen_dist = Observable(ug_scen_dist)
#     obs_g_scen_dist = Observable(g_scen_dist)

#     # Color transparency for density plots
#     # Note: Density plots currently cannot handle empty datasets
#     #       as what might happen if user selects a region with no results.
#     #       so instead we set alpha to 0.0 to hide it.
#     cf_hist_alpha = Observable(0.5)
#     ug_hist_alpha = Observable(0.5)
#     g_hist_alpha = Observable(0.5)

#     # Get intervention/criteria inputs for each scenario
#     interv_criteria = ms[(ms.component.=="EnvironmentalLayer").|(ms.component.=="Intervention").|(ms.component.=="Criteria"), [:fieldname, :full_bounds]]
#     input_names = vcat(["RCP", interv_criteria.fieldname...])
#     in_pcp_data = normalize(Matrix(rs.inputs[:, input_names]))
#     # in_pcp_lines = Observable(in_pcp_data)

#     # Get mean outcomes for each scenario
#     outcome_pcp_data = hcat([
#         mean_tac_outcomes,
#         vec(mean(ADRIA.metrics.scenario_asv(rs), dims=1)),
#         vec(mean(ADRIA.metrics.scenario_rsv(rs), dims=1))
#     ]...)
#     disp_names = ["TAC", "ASV", "RSV"]

#     out_pcp_data = normalize(outcome_pcp_data)

#     # Specify interactive elements and behavior
#     # TODO: Lift on data to be plotted
#     #       Controls simply update transparency settings etc, and update the dataset to be plotted
#     #       All other elements simply update when the underlying dataset updates
#     # https://discourse.julialang.org/t/interactive-plot-with-makielayout/48843
#     onany(time_slider.interval, tac_slider.interval) do time_val, tac_val
#         # Update slider labels
#         tac_bot_val[] = tac_val[1]
#         tac_top_val[] = tac_val[2]

#         left_year_val[] = "$(Int(floor(time_val[1])))"
#         right_year_val[] = "$(Int(ceil(time_val[2])))"

#         # Trajectories
#         # tac_idx = (mean_tac_outcomes .>= tac_val[1]-0.5) .& (mean_tac_outcomes .<= tac_val[2]+0.5)

#         # Convert time ranges to index values
#         t_idx = Int(time_val[1] - (year_range[1]) + 1), Int(time_val[2] - (year_range[1]) + 1)

#         hide_idx = vec(all((tac_val[1] .<= tac_scens[t_idx[1]:t_idx[2], :] .<= tac_val[2]) .== 0, dims=1))
#         show_idx = Bool.(zeros(Int64, length(hide_idx)) .⊻ hide_idx)  # inverse of hide

#         scen_dist = dropdims(mean(tac[timesteps=t_idx[1]:t_idx[2]], dims=:timesteps), dims=:timesteps)

#         # Boolean index of scenarios to hide (inverse of tac_idx)
#         # hide_idx = Bool.(ones(Int64, length(tac_idx)) .⊻ tac_idx)
#         if !all(hide_idx .== 0)
#             # Hide scenarios that were filtered out
#             cf_dist = scen_dist[show_idx.&scen_types.counterfactual]
#             ug_dist = scen_dist[show_idx.&scen_types.unguided]
#             g_dist = scen_dist[show_idx.&scen_types.guided]
#         else
#             cf_dist = scen_dist[scen_types.counterfactual]
#             ug_dist = scen_dist[scen_types.unguided]
#             g_dist = scen_dist[scen_types.guided]
#         end

#         # Update scenario density plot
#         if !isempty(cf_dist)
#             obs_cf_scen_dist[] = cf_dist
#             cf_hist_alpha[] = 0.5
#         else
#             cf_hist_alpha[] = 0.0
#         end

#         if !isempty(ug_dist)
#             obs_ug_scen_dist[] = ug_dist
#             ug_hist_alpha[] = 0.5
#         else
#             ug_hist_alpha[] = 0.0
#         end

#         if !isempty(g_dist)
#             obs_g_scen_dist[] = g_dist
#             g_hist_alpha[] = 0.5
#         else
#             g_hist_alpha[] = 0.0
#         end

#         # Determine level of transparency for each line (maximum of 0.6)
#         min_step = (1 / 0.05)
#         color_weight = min((1.0 / (count(show_idx .> 0) / min_step)), 0.6)

#         obs_color[] = scenario_colors(rs, color_weight, hide_idx)
#     end

#     # Trajectories
#     # series!(traj_display, timesteps(rs), tac_data, color=@lift($obs_color[:]))
#     series!(traj_display, timesteps(rs), tac_data, color=obs_color)

#     # Legend(traj_display)  legend=["Counterfactual", "Unguided", "Guided"]
#     density!(scen_hist, @lift($obs_cf_scen_dist[:]), direction=:y, color=(:red, cf_hist_alpha))
#     density!(scen_hist, @lift($obs_ug_scen_dist[:]), direction=:y, color=(:green, ug_hist_alpha))
#     density!(scen_hist, @lift($obs_g_scen_dist[:]), direction=:y, color=(:blue, g_hist_alpha))

#     hidedecorations!(scen_hist)
#     hidespines!(scen_hist)

#     # TODO: Separate this out into own function
#     # Make temporary copy of GeoPackage as GeoJSON
#     tmpdir = mktempdir()
#     geo_fn = GDF.write(joinpath(tmpdir, "Aviz_$(rs.name).geojson"), rs.site_data; driver="geojson")
#     geodata = GeoMakie.GeoJSON.read(read(geo_fn))

#     # Get bounds to display
#     centroids = rs.site_centroids
#     lon = first.(centroids)
#     lat = last.(centroids)

#     # Display map
#     mean_rc_sites = mean(ADRIA.metrics.relative_cover(rs), dims=(:scenarios, :timesteps))
#     map_buffer = 0.005
#     spatial = GeoAxis(
#         map_display;  # any cell of the figure's layout
#         lonlims=(minimum(lon) - map_buffer, maximum(lon) + map_buffer),
#         latlims=(minimum(lat) - map_buffer, maximum(lat) + map_buffer),
#         xlabel="Long",
#         ylabel="Lat",
#         dest="+proj=latlong +datum=WGS84"
#     )

#     poly!(spatial, geodata, color=vec(mean_rc_sites), colormap=:plasma)
#     datalims!(spatial)

#     # Fill pairplot
#     # Get mean outcomes for each scenario
#     pairplot!(pair_display, outcome_pcp_data, disp_names)

#     # Parallel Coordinate Plot
#     pcp!(interv_pcp_display, in_pcp_data, input_names; color=@lift($obs_color[:]))
#     pcp!(outcome_pcp_display, out_pcp_data, disp_names; color=@lift($obs_color[:]))

#     gl_screen = display(f)
#     DataInspector()

#     wait(gl_screen)
# end
# function explore(rs_path::String)
#     explore(ADRIA.load_results(rs_path))
# end

end  # module

# Allow use from terminal if this file is run directly
# if abspath(PROGRAM_FILE) == @__FILE__
#     if "explore" in ARGS
#         rs_pkg = ARGS[2]
#         ADRIA.explore(rs_pkg)
#     end
# end
