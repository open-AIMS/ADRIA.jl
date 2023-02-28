using Distributions
using Test


if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")
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
    sum_cover = fill(0.1, nrow(criteria_df), nrow(dom.site_data))

    area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
    dhw_scens = dom.dhw_scens[1, :, criteria_df.dhw_scenario[1]]
    wave_scens = dom.wave_scens[1, :, criteria_df.wave_scenario[1]]

    mcda_vars = ADRIA.DMCDA_vars(dom, criteria_df[1, :], site_ids, sum_cover, area_to_seed, wave_scens, dhw_scens)
    n_sites = length(mcda_vars.site_ids)
    @test (length(mcda_vars.in_conn) == n_sites) && (length(mcda_vars.out_conn) == n_sites) || "Connectivity input is incorrect size."
    @test length(mcda_vars.dam_prob) == n_sites || "Wave damage input is incorrect size."
    @test length(mcda_vars.heat_stress_prob) == n_sites || "Heat stress input is incorrect size."
    @test length(mcda_vars.sum_cover) == n_sites || "Initial cover input is incorrect size."
end
@testset "Unguided site selection" begin
    n_intervention_sites = 5
    prefseedsites = zeros(Int, n_intervention_sites)
    prefshadesites = zeros(Int, n_intervention_sites)
    seed_years = true
    shade_years = true
    max_cover = [0.0, 3000.0, 5000.0, 0.0, 0.0]
    depth_priority = collect(1:5)

    prefseedsites, prefshadesites = ADRIA.unguided_site_selection(prefseedsites, prefshadesites, true, true, 5, max_cover, depth_priority)

    # Check that only two sites are selected (the sites where k > 0.0)
    @test length(prefseedsites[prefseedsites.>0]) == 2
    @test length(prefshadesites[prefshadesites.>0]) == 2

    @test all([in(sid, [2, 3]) for sid in prefseedsites[prefseedsites.>0]])
    @test all([in(sid, [2, 3]) for sid in prefshadesites[prefshadesites.>0]])
end

@testset "Guided site selection without ADRIA ecological model" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    N = 2^3
    criteria_df = ADRIA.sample_site_selection(dom, N)  # get scenario dataframe

    area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
    ts = 5  # time step to perform site selection at

    sum_cover = fill(0.1, N, ADRIA.n_locations(dom))
    ranks = ADRIA.run_site_selection(dom, criteria_df, sum_cover, area_to_seed, ts)

    @test size(ranks, 1) == sum(criteria_df.guided .> 0) || "Specified number of scenarios was not carried out."
    @test size(ranks, 2) == length(dom.site_ids) || "Ranks storage is not correct size for this domain."

    sel_sites = unique(ranks)
    sel_sites = sel_sites[sel_sites.!=0.0]
    possible_ranks = collect(Float64, 1:ADRIA.n_locations(dom)+1.0)

    @test all([in(ss, possible_ranks) for ss in sel_sites]) || "Impossible rank assigned."
end
