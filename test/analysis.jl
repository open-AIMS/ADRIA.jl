using WGLMakie, GeoMakie, GraphMakie
using ADRIA
using ADRIA.metrics: total_absolute_cover
using Statistics

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

Makie.inline!(false)

"""Test larger scenario run with figure creation"""
function test_rs_w_fig(rs::ADRIA.ResultSet, scens::ADRIA.DataFrame)

    # Visualize results (in terms of absolute coral cover)
    s_tac = ADRIA.metrics.scenario_total_cover(rs)
    ADRIA.viz.scenarios(rs, s_tac)

    tac = ADRIA.metrics.total_absolute_cover(rs)
    rsv = ADRIA.metrics.relative_shelter_volume(rs)
    juves = ADRIA.metrics.relative_juveniles(rs)

    s_tac = ADRIA.metrics.scenario_total_cover(rs)
    s_rsv = ADRIA.metrics.scenario_rsv(rs)
    s_juves = ADRIA.metrics.scenario_relative_juveniles(rs)

    ## Visualization

    # Some shared options for the example plots below
    fig_opts = Dict(:size => (1600, 800))

    # Factors of Interest
    opts = Dict(
        :factors => [
            :RCP,
            :dhw_scenario,
            :wave_scenario,
            :guided,
            :N_seed_TA,
            :N_seed_CA,
            :fogging,
            :SRM,
            :a_adapt
        ]
    )

    ### Scenario outcomes

    fig_s_tac = ADRIA.viz.scenarios(
        rs, s_tac; fig_opts=fig_opts, axis_opts=Dict(:ylabel => "Scenario Total Cover")
    )
    # save("scenarios_tac.png", fig_s_tac)

    tf = Figure(; size=(1600, 600))  # resolution in pixels

    # Implicitly create a single figure with 2 columns
    ADRIA.viz.scenarios!(
        tf[1, 1],
        rs,
        s_tac;
        opts=Dict(:by_RCP => false, :legend => false),
        axis_opts=Dict(:title => "TAC [m²]")
    )
    ADRIA.viz.scenarios!(
        tf[1, 2],
        rs,
        s_juves;
        opts=Dict(:summarize => false),
        axis_opts=Dict(:title => "Juveniles [%]")
    )

    # tf  # display the figure
    # save("aviz_scenario.png", tf)  # save the figure to a file

    ### Intervention location selection - visualisation

    # TODO Fix this
    # Calculate frequencies with which each site was selected at each rank
    #rank_freq = ADRIA.decision.ranks_to_frequencies(
    #    rs.ranks[intervention=1];
    #    n_ranks=length(rs.ranks[intervention=1].sites),
    #    agg_func=x -> dropdims(sum(x; dims=:timesteps); dims=:timesteps),
    #)
    #
    ## Plot 1st rank frequencies as a colormap
    #rank_fig = ADRIA.viz.ranks_to_frequencies(
    #    rs, rank_freq, 1; fig_opts=Dict(:size => (1200, 800))
    #)
    ## save("single_rank_plot.png", rank_fig)
    #
    ## Plot 1st, 2nd and 3rd rank frequencies as an overlayed colormap
    #rank_fig = ADRIA.viz.ranks_to_frequencies(
    #    rs, rank_freq, [1, 2, 3]; fig_opts=Dict(:size => (1200, 800))
    #)
    # save("ranks_plot.png", rank_fig)

    ### PAWN sensitivity (heatmap overview)

    # Sensitivity (of mean scenario outcomes to factors)
    mean_s_tac = vec(mean(s_tac; dims=1))
    tac_Si = ADRIA.sensitivity.pawn(rs, mean_s_tac)
    pawn_fig = ADRIA.viz.pawn(
        tac_Si;
        opts,
        fig_opts
    )
    # save("pawn_si.png", pawn_fig)

    ### Temporal Sensitivity Analysis

    tsa_s = ADRIA.sensitivity.tsa(rs, s_tac)
    tsa_fig = ADRIA.viz.tsa(
        rs,
        tsa_s;
        opts,
        fig_opts
    )
    # save("tsa.png", tsa_fig)

    ### Convergence Analysis

    outcome = dropdims(
        mean(ADRIA.metrics.scenario_total_cover(rs); dims=:timesteps); dims=:timesteps
    )

    # Display convergence for specific factors of interest ("foi") within a single figure.
    # Bands represent the 95% confidence interval derived from the number of conditioning
    # points, the default for which is ten (i.e., 10 samples).
    # Due to the limited sample size, care should be taken when interpreting the figure.
    foi = [:dhw_scenario, :wave_scenario, :guided]
    Si_conv = ADRIA.sensitivity.convergence(scens, outcome, foi)
    ADRIA.viz.convergence(Si_conv, foi)

    # Convergence analysis of factors grouped by model component as a heat map
    components = [:EnvironmentalLayer, :Intervention, :Coral]
    Si_conv = ADRIA.sensitivity.convergence(rs, scens, outcome, components)
    ADRIA.viz.convergence(Si_conv, components; opts=Dict(:viz_type => :heatmap))

    ### Time Series Clustering

    # Extract metric from scenarios
    s_tac = ADRIA.metrics.scenario_total_cover(rs)

    # Cluster scenarios
    n_clusters = 4
    clusters = ADRIA.analysis.cluster_scenarios(s_tac, n_clusters)

    axis_opts = Dict{Symbol,Any}(
        :title => "Time Series Clustering with $n_clusters clusters",
        :ylabel => "TAC [m²]",
        :xlabel => "Timesteps [years]"
    )

    tsc_fig = ADRIA.viz.clustered_scenarios(
        s_tac, clusters; opts=Dict{Symbol,Any}(:summarize => true), fig_opts=fig_opts,
        axis_opts=axis_opts
    )

    # Save final figure
    # save("tsc.png", tsc_fig)

    ### Target clusters

    # Extract metric from scenarios
    asv = ADRIA.metrics.absolute_shelter_volume(rs)

    # Time series summarizing scenarios for each site
    asv_site_series = ADRIA.metrics.loc_trajectory(median, asv)

    # Cluster scenarios
    n_clusters = 6
    asv_clusters = ADRIA.analysis.cluster_scenarios(asv_site_series, n_clusters)

    # Target scenarios that belong to the two lowest value clusters
    lowest = x -> x .∈ [sort(x; rev=true)[1:2]]
    asv_target = ADRIA.analysis.find_scenarios(asv_site_series, asv_clusters, lowest)

    # Plot targeted scenarios
    axis_opts = Dict(:ylabel => "Absolute Shelter Volume", :xlabel => "Timesteps [years]")

    tsc_asc_fig = ADRIA.viz.clustered_scenarios(
        asv_site_series, asv_target; axis_opts=axis_opts, fig_opts=fig_opts
    )

    # Save final figure
    # save("tsc_asv.png", tsc_asc_fig)

    ### Multiple Time Series Clustering

    metrics::Vector{ADRIA.metrics.Metric} = [
        ADRIA.metrics.scenario_total_cover,
        ADRIA.metrics.scenario_asv,
        ADRIA.metrics.scenario_absolute_juveniles
    ]

    outcomes = ADRIA.metrics.scenario_outcomes(rs, metrics)
    n_clusters = 6

    # Clusters matrix
    outcomes_clusters::AbstractMatrix{Int64} = ADRIA.analysis.cluster_scenarios(
        outcomes, n_clusters
    )

    # Filter scenarios that belong to on of the 4 high value clusters for all outcomes
    highest_clusters(x) = x .∈ [sort(x; rev=true)[1:4]]
    robust_scens = ADRIA.analysis.find_scenarios(
        outcomes, outcomes_clusters, highest_clusters
    )

    ### Time Series Clustering Map

    # Extract metric from scenarios
    tac = ADRIA.metrics.total_absolute_cover(rs)

    # Get a timeseries summarizing the scenarios for each site
    tac_site_series = ADRIA.metrics.loc_trajectory(median, tac)

    # Cluster scenarios
    n_clusters = 6
    clusters = ADRIA.analysis.cluster_scenarios(tac_site_series, n_clusters)

    # Get a vector summarizing the scenarios and timesteps for each site
    tac_sites = ADRIA.metrics.per_loc(median, tac)

    # Plot figure
    tsc_map_fig = ADRIA.viz.map(rs, tac_sites, clusters)

    # Save final figure
    # save("tsc_map.png", tsc_map_fig)

    ### Rule Induction (using Series Clusters)

    # Find Time Series Clusters
    s_tac = ADRIA.metrics.scenario_total_cover(rs)
    num_clusters = 6
    clusters = ADRIA.analysis.cluster_scenarios(s_tac, num_clusters)

    # Target scenarios
    target_clusters = ADRIA.analysis.target_clusters(clusters, s_tac)

    # Select only desired features
    fields_iv =
        ADRIA.component_params(
            rs, [Intervention, FogCriteriaWeights, SeedCriteriaWeights]
        ).fieldname
    scenarios_iv = scens[:, fields_iv]

    # Use SIRUS algorithm to extract rules
    max_rules = 10
    rules_iv = ADRIA.analysis.cluster_rules(
        target_clusters, scenarios_iv, max_rules; remove_duplicates=true
    )
    rules_iv_duplicates = ADRIA.analysis.cluster_rules(
        target_clusters, scenarios_iv, max_rules; remove_duplicates=false
    )
    ADRIA.analysis.print_rules(rules_iv)
    ADRIA.analysis.print_rules(rules_iv_duplicates)

    # Plot scatters for each rule highlighting the area selected them
    rules_scatter_fig = ADRIA.viz.rules_scatter(
        rs,
        scenarios_iv,
        target_clusters,
        rules_iv;
        fig_opts=fig_opts,
        opts=opts
    )

    ADRIA.viz.rules_scatter(
        rs,
        scenarios_iv,
        target_clusters,
        rules_iv_duplicates;
        fig_opts=fig_opts,
        opts=opts
    )

    # Save final figure
    # save("rules_scatter.png", rules_scatter_fig)

    ### Regional Sensitivity Analysis

    foi = [
        :dhw_scenario,
        :wave_scenario,
        :N_seed_TA,
        :N_seed_CA,
        :fogging,
        :SRM
    ]

    tac_rs = ADRIA.sensitivity.rsa(rs, mean_s_tac; S=10)
    rsa_fig = ADRIA.viz.rsa(
        rs,
        tac_rs,
        foi;
        opts,
        fig_opts
    )
    # save("rsa.png", rsa_fig)

    ### Outcome mapping

    tf = Figure(; size=(1600, 1200))  # resolution in pixels

    # Indicate factor values that are in the top 50 percentile
    tac_om_50 = ADRIA.sensitivity.outcome_map(
        rs, mean_s_tac, x -> any(x .>= 0.5), foi; S=20
    )
    ADRIA.viz.outcome_map!(
        tf[1, 1],
        rs,
        tac_om_50,
        foi;
        axis_opts=Dict(
            :title => "Regions which lead to Top 50th Percentile Outcomes",
            :ylabel => "TAC [m²]"
        )
    )

    # Indicate factor values that are in the top 30 percentile
    tac_om_30 = ADRIA.sensitivity.outcome_map(
        rs, mean_s_tac, x -> any(x .>= 0.7), foi; S=20
    )
    ADRIA.viz.outcome_map!(
        tf[2, 1],
        rs,
        tac_om_30,
        foi;
        axis_opts=Dict(
            :title => "Regions which lead to Top 30th Percentile Outcomes",
            :ylabel => "TAC [m²]"
        ))
    # save("outcome_map.png", tf)

    return rs
end

test_rs_w_fig(TEST_RS, TEST_SCENS)
