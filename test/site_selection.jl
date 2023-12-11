using Distributions
using Test


if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Test_domain")
end

@testset "site selection" begin
    # TODO: Complete tests with @tests

    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    p_tbl = ADRIA.param_table(dom)

    p_tbl[:, :depth_offset] .= 7.0

    # ranks = ADRIA.site_selection(dom, p_tbl, 1, 10, 1)
end
@testset "MCDA variable constructor" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    criteria_df = ADRIA.sample_site_selection(dom, 1) # get scenario dataframe
    site_ids = collect(1:length(dom.site_ids))
    available_space = rand(Uniform(200, 30000), length(site_ids))

    area_to_seed = 962.11  # Area of seeded corals in m^2.

    dhw_scens = dom.dhw_scens[1, :, criteria_df.dhw_scenario[1]]
    wave_scens = dom.wave_scens[1, :, criteria_df.wave_scenario[1]]

    mcda_vars = ADRIA.decision.DMCDA_vars(
        dom,
        criteria_df[1, :],
        site_ids,
        available_space,
        area_to_seed,
        wave_scens,
        dhw_scens,
    )
    n_sites = length(mcda_vars.site_ids)
    @test (size(mcda_vars.conn, 1) == n_sites) && (size(mcda_vars.conn, 2) == n_sites) || "Connectivity input is incorrect size."
    @test length(mcda_vars.dam_prob) == n_sites || "Wave damage input is incorrect size."
    @test length(mcda_vars.heat_stress_prob) == n_sites || "Heat stress input is incorrect size."
    @test length(mcda_vars.leftover_space) == n_sites ||
        "Initial cover input is incorrect size."
end
@testset "Unguided site selection" begin
    n_intervention_sites = 5
    pref_seed_sites = zeros(Int64, n_intervention_sites)
    pref_fog_sites = zeros(Int64, n_intervention_sites)
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
        depth_priority,
    )

    # Check that only two sites are selected (the sites where k > 0.0)
    @test length(pref_seed_sites[pref_seed_sites .> 0]) == 2
    @test length(pref_fog_sites[pref_fog_sites .> 0]) == 2

    @test all([in(sid, [2, 3]) for sid in pref_seed_sites[pref_seed_sites .> 0]])
    @test all([in(sid, [2, 3]) for sid in pref_fog_sites[pref_fog_sites .> 0]])
end

@testset "Guided site selection without ADRIA ecological model" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    N = 2^3
    scens = ADRIA.sample_site_selection(dom, N)  # get scenario dataframe

    area_to_seed = 962.11  # Area of seeded corals in m^2.

    sum_cover = repeat(sum(dom.init_coral_cover, dims=1), size(scens, 1))
    ranks = ADRIA.decision.rank_locations(dom, scens, sum_cover, area_to_seed)

    @test length(ranks.scenarios) == sum(scens.guided .> 0) || "Specified number of scenarios was not carried out."
    @test length(ranks.sites) == length(dom.site_ids) || "Ranks storage is not correct size for this domain."

    sel_sites = unique(ranks)
    sel_sites = sel_sites[sel_sites.!=0.0]
    possible_ranks = collect(Float64, 1:ADRIA.n_locations(dom)+1.0)

    @test all([in(ss, possible_ranks) for ss in sel_sites]) || "Impossible rank assigned."
end

@testset "Test ranks line up with ordering" begin
    supported_methods = ADRIA.decision.mcda_methods()
    mcda_func = supported_methods[rand(1:length(supported_methods))]
    n_sites = 20

    S = ADRIA.decision.mcda_normalize(rand(Uniform(0, 1), n_sites, 6))

    weights = ADRIA.decision.mcda_normalize(rand(Uniform(0, 1), 6))
    n_site_int = 5
    site_ids = collect(1:n_sites)
    S = hcat(site_ids, S)

    rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    prefsites, s_order = ADRIA.decision.rank_sites!(
        S, weights, rankings, n_site_int, mcda_func, 2
    )

    @test all([(rankings[rankings[:, 1].==s_order[rank, 1], 2].==rank)[1] for rank in 1:size(s_order, 1)]) || "Ranking does not match mcda score ordering"

end
