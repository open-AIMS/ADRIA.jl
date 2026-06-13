using Test
using ADRIA.Distributions

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

@testset "site selection" begin
    # TODO: Complete tests with @tests

    p_tbl = ADRIA.param_table(ADRIA_DOM_45)

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
    N = 2^3
    scens = ADRIA.sample_selection(ADRIA_DOM_45, N)  # get scenario dataframe

    area_to_seed = 962.11  # Area of seeded corals in m^2.

    sum_cover = repeat(sum(ADRIA_DOM_45.init_coral_cover; dims=1), size(scens, 1))
    ranks = ADRIA.decision.rank_locations(ADRIA_DOM_45, scens, sum_cover, area_to_seed)

    @test length(ranks.scenarios) == sum(scens.guided .> 0) ||
        "Specified number of scenarios was not carried out."
    @test length(ranks.sites) == length(ADRIA_DOM_45.loc_ids) ||
        "Ranks storage is not correct size for this domain."

    sel_sites = unique(ranks)
    sel_sites = sel_sites[sel_sites .!= 0.0]
    possible_ranks = collect(Float64, 1:(ADRIA.n_locations(ADRIA_DOM_45) + 1.0))

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
