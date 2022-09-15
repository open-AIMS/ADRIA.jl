module Aviz

using GLMakie
using GLMakie.GeometryBasics
using Statistics, Distributions
using GeoDataFrames, GeoInterface

using ADRIA


include("./plotting.jl")
include("./layout.jl")
include("./theme.jl")


"""Main entry point for app."""
function julia_main()::Cint
    main_menu()

    return 0
end

function main_menu()
    f = Figure()

    # img = load(assetpath("../assets/ADRIA_logo.png"))
    # logo = image(f[1,1], img)
    # hidedecorations!(logo)
    # hidespines!(logo)

    Label(f[1,1], "Enter ADRIA Result Set to analyze")
    rs_path_tb = Textbox(f[2, 1], placeholder="./Moore_RS")  # placeholder="Path to ADRIA Result Set"
    rs_path_tb.stored_string[] = "./Moore_RS"
    status_label = Label(f[3,1], "")

    launch_button = Button(f[4,1], label="Analyze")

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
                Aviz.gui_analysis(rs)
            end
        else
            rs_path_tb.bordercolor = :red
            status_label.text[] = "Invalid path"
        end
    end

    gl_screen = display(f)
    wait(gl_screen)
end


function gui_analysis(rs::ADRIA.ResultSet)
    layout = create_layout(resolution=(1920,1080))

    f = layout.figure
    controls = layout.controls
    traj_display = layout.trajectory
    scen_hist = layout.scen_hist
    map_display = layout.map

    interv_pcp_display = layout.interv_pcp
    pair_display = layout.pairplot
    outcome_pcp_display = layout.outcome_pcp

    colsize!(f.layout, 1, Fixed(400))
    colsize!(f.layout, 2, Fixed(1100))

    color_map = scenario_colors(rs)
    obs_color = Observable(color_map)

    # Controls
    ms = rs.model_spec
    # intervention_components = ms[ms.component .== "Intervention", [:fieldname, :full_bounds]]

    tac_scens = ADRIA.metrics.scenario_total_cover(rs)
    mean_tac_outcomes = vec(mean(tac_scens, dims=1))
    tac_min_max = (minimum(mean_tac_outcomes), maximum(mean_tac_outcomes))

    tac_label = Label(controls[1,1], "Mean TAC (m²)")
    tac_slider = IntervalSlider(controls[2,1], 
                                range=LinRange(floor(Int64, tac_min_max[1])-1, ceil(Int64, tac_min_max[2])+1, Int(5e5)),
                                startvalues=tac_min_max,
                                width=350)

    # Dynamic label text for TAC slider
    tac_lbl_txt = lift(tac_slider.interval) do intv
        string(round.(intv ./ 1e6, digits=4)) * "M"
    end
    Label(controls[2,1], tac_lbl_txt)

    scen_tac = ADRIA.metrics.scenario_total_cover(rs)
    tac_data = Matrix(scen_tac')
    tac_traj = Observable(tac_data)

    # tac = ADRIA.metrics.total_absolute_cover(rs);
    tac = ADRIA.metrics.scenario_total_cover(rs);

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
    interv_criteria = ms[(ms.component .== "EnvironmentalLayer") .| (ms.component .== "Intervention") .| (ms.component .== "Criteria"), [:fieldname, :full_bounds]]
    input_names = vcat(["RCP", interv_criteria.fieldname...])
    in_pcp_data = normalize(Matrix(rs.inputs[:, input_names]))
    # in_pcp_lines = Observable(in_pcp_data)


    # Get mean outcomes for each scenario
    outcome_pcp_data = hcat([
        mean_tac_outcomes, 
        vec(mean(ADRIA.metrics.scenario_relative_cover(rs), dims=1)),
        vec(mean(ADRIA.metrics.scenario_asv(rs), dims=1))
    ]...)
    disp_names = ["TAC", "RC", "ASV"]

    out_pcp_data = normalize(outcome_pcp_data)
    # out_pcp_lines = Observable(out_pcp_data)

    # Specify interactive elements and behavior
    lift(tac_slider.interval) do intv
        # Trajectories
        tac_idx = (mean_tac_outcomes .>= intv[1]-0.5) .& (mean_tac_outcomes .<= intv[2]+0.5)

        # Boolean index of scenarios to hide
        hide = Bool.(ones(Int64, length(tac_idx)) .⊻ tac_idx)

        t = copy(tac_data)
        out_pcp = copy(out_pcp_data)
        in_pcp = copy(in_pcp_data)
        s_dist = copy(scen_dist)

        if !all(hide .== 0)
            # Hide scenarios that were filtered out
            t[hide, :] .= NaN

            cf_dist = s_dist[tac_idx .& scen_types.counterfactual]
            ug_dist = s_dist[tac_idx .& scen_types.unguided]
            g_dist = s_dist[tac_idx .& scen_types.guided]

            in_pcp[hide, :] .= NaN
            out_pcp[hide, :] .= NaN
        else
            cf_dist = s_dist[scen_types.counterfactual]
            ug_dist = s_dist[scen_types.unguided]
            g_dist = s_dist[scen_types.guided]
        end

        tac_traj[] = t

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

        # Update PCPs
        # in_pcp_lines[] = in_pcp
        # out_pcp_lines[] = out_pcp
        # notify(out_pcp_lines)

        min_step = (1/0.05)
        color_weight = (1.0 / (count(tac_idx .> 0) / min_step))

        obs_color[] = scenario_colors(rs, color_weight, hide)
    end

    # Trajectories
    series!(traj_display, timesteps(rs), tac_data, color=@lift($obs_color[:]))  # , solid_color=(:blue, 0.1)

    # Legend(traj_display)  legend=["Counterfactual", "Unguided", "Guided"]
    density!(scen_hist, @lift($obs_cf_scen_dist[:]), direction=:y, color=(:red, @lift($cf_hist_alpha[])))
    density!(scen_hist, @lift($obs_ug_scen_dist[:]), direction=:y, color=(:green, @lift($ug_hist_alpha[])))
    density!(scen_hist, @lift($obs_g_scen_dist[:]), direction=:y, color=(:blue, @lift($g_hist_alpha[])))
    # density!(scen_hist, scen_dist, direction=:y)
    hidedecorations!(scen_hist)
    hidespines!(scen_hist)

    # Display map
    map_coords = GeoInterface.coordinates.(rs.site_data.geometry)
    plot_poly!.(map_display, map_coords)

    # Fill pairplot
    # Get mean outcomes for each scenario
    pairplot!(pair_display, outcome_pcp_data, disp_names)

    # Parallel Coordinate Plot
    pcp!(interv_pcp_display, in_pcp_data, input_names; color=@lift($obs_color[:]))
    pcp!(outcome_pcp_display, out_pcp_data, disp_names; color=@lift($obs_color[:]))

    gl_screen = display(f)
    DataInspector()

    wait(gl_screen);
end
function gui_analysis(rs_path::String)
    rs = ADRIA.load_results(rs_path)
    gui_analysis(rs)
end

end


# Allow use from terminal if this file is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    if "analyze" in ARGS
        rs_pkg = ARGS[2]
        Aviz.gui_analysis(rs_pkg)
    end
end
