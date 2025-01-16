using Test
using ADRIA.Distributions

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

@testset "site selection" begin
    # TODO: Complete tests with @tests

    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    p_tbl = ADRIA.param_table(dom)

    p_tbl[:, :depth_offset] .= 7.0

    # ranks = ADRIA.site_selection(dom, p_tbl, 1, 10, 1)
end

@testset "Unguided site selection" begin
    n_intervention_locs = 5
    pref_seed_sites = zeros(Int64, n_intervention_locs)
    pref_fog_sites = zeros(Int64, n_intervention_locs)
    seed_years = true
    fog_years = true
    max_cover = [0.0, 3000.0, 5000.0, 0.0, 0.0]
    depth_priority = collect(1:5)

    pref_seed_sites, pref_fog_sites = ADRIA.decision.unguided_site_selection(
        pref_seed_sites,
        pref_fog_sites,
        seed_years,
        fog_years,
        5,
        max_cover,
        depth_priority
    )

    # Check that only two sites are selected (the sites where k > 0.0)
    @test length(pref_seed_sites[pref_seed_sites .> 0]) == 2
    @test length(pref_fog_sites[pref_fog_sites .> 0]) == 2

    @test all([in(sid, [2, 3]) for sid in pref_seed_sites[pref_seed_sites .> 0]])
    @test all([in(sid, [2, 3]) for sid in pref_fog_sites[pref_fog_sites .> 0]])
end

@testset "Guided site selection without ADRIA ecological model" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    N = 2^3
    scens = ADRIA.sample_selection(dom, N)  # get scenario dataframe

    area_to_seed = 962.11  # Area of seeded corals in m^2.

    sum_cover = repeat(sum(dom.init_coral_cover; dims=1), size(scens, 1))
    ranks = ADRIA.decision.rank_locations(dom, scens, sum_cover, area_to_seed)

    @test length(ranks.scenarios) == sum(scens.guided .> 0) ||
        "Specified number of scenarios was not carried out."
    @test length(ranks.sites) == length(dom.loc_ids) ||
        "Ranks storage is not correct size for this domain."

    sel_sites = unique(ranks)
    sel_sites = sel_sites[sel_sites .!= 0.0]
    possible_ranks = collect(Float64, 1:(ADRIA.n_locations(dom) + 1.0))

    @test all([in(ss, possible_ranks) for ss in sel_sites]) || "Impossible rank assigned."
end

@testset "Test ranks line up with ordering" begin
    supported_methods = ADRIA.decision.mcda_methods()
    mcda_func = supported_methods[rand(1:length(supported_methods))]
    n_locs = 20

    S = ADRIA.decision.mcda_normalize(rand(Uniform(0, 1), n_locs, 6))

    weights = ADRIA.decision.mcda_normalize(rand(Uniform(0, 1), 6))
    n_site_int = 5
    loc_ids = collect(1:n_locs)
    S = hcat(loc_ids, S)

    rankings = Int64[loc_ids zeros(Int64, n_locs) zeros(Int64, n_locs)]

    prefsites, s_order = ADRIA.decision.rank_sites!(
        S, weights, rankings, n_site_int, mcda_func, 2
    )

    @test all([
        (rankings[rankings[:, 1] .== s_order[rank, 1], 2] .== rank)[1] for
        rank in 1:size(s_order, 1)
    ]) || "Ranking does not match mcda score ordering"
end

@testset "Test decision matrix spatial plotting" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    N = 2^3
    scens = ADRIA.sample_selection(dom, N)  # get scenario dataframe
    scen = scens[1, :]

    # Get seeding preferences
    seed_pref = ADRIA.decision.SeedPreferences(dom, scen)

    # Calculate criteria vectors
    # Cover
    sum_cover = vec(sum(dom.init_coral_cover; dims=1).data)
    # DHWS
    dhw_scens = dom.dhw_scens[:, :, Int64(scen["dhw_scenario"])]
    plan_horizon = Int64(scen["plan_horizon"])
    decay = 0.99 .^ (1:(plan_horizon + 1)) .^ 2
    dhw_projection = ADRIA.decision.weighted_projection(
        dhw_scens, 1, plan_horizon, decay, 75
    )
    # Connectivity
    area_weighted_conn = dom.conn.data .* ADRIA.site_k_area(dom)
    conn_cache = similar(area_weighted_conn)
    in_conn, out_conn, network = ADRIA.connectivity_strength(
        area_weighted_conn, sum_cover, conn_cache
    )

    # Create decision matrix
    seed_decision_mat = ADRIA.decision.decision_matrix(
        dom.loc_ids,
        seed_pref.names;
        seed_in_connectivity=in_conn,
        seed_out_connectivity=out_conn,
        seed_heat_stress=dhw_projection,
        seed_coral_cover=sum_cover
    )

    # Get results from applying MCDA algorithm
    crit_agg = ADRIA.decision.criteria_aggregated_scores(
        seed_pref, seed_decision_mat, mcda_funcs[1]
    )

    # Don't plot constant criteria
    is_const = Bool[length(x) == 1 for x in unique.(eachcol(seed_decision_mat.data))]

    # Plot normalized scores and criteria as map
    decision_mat_fig = ADRIA.viz.selection_criteria_map(
        dom, seed_decision_mat[criteria=.!is_const],
        crit_agg.scores ./ maximum(crit_agg.scores)
    )
end
