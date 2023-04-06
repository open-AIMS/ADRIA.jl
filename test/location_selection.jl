using Distributions
using Test


if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")
end

@testset "location selection" begin
    # TODO: Complete tests with @tests

    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    p_tbl = ADRIA.param_table(dom)

    p_tbl[:, :depth_offset] .= 7.0

    # ranks = ADRIA.location_selection(dom, p_tbl, 1, 10, 1)
end
@testset "MCDA variable constructor" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    criteria_df = ADRIA.sample_location_selection(dom, 1) # get scenario dataframe
    location_ids = collect(1:length(dom.location_ids))
    sum_cover = fill(0.1, nrow(criteria_df), nrow(dom.location_data))

    area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
    dhw_scens = dom.dhw_scens[1, :, criteria_df.dhw_scenario[1]]
    wave_scens = dom.wave_scens[1, :, criteria_df.wave_scenario[1]]

    mcda_vars = ADRIA.DMCDA_vars(dom, criteria_df[1, :], location_ids, sum_cover, area_to_seed, wave_scens, dhw_scens)
    n_locations = length(mcda_vars.location_ids)
    @test (length(mcda_vars.in_conn) == n_locations) && (length(mcda_vars.out_conn) == n_locations) || "Connectivity input is incorrect size."
    @test length(mcda_vars.dam_prob) == n_locations || "Wave damage input is incorrect size."
    @test length(mcda_vars.heat_stress_prob) == n_locations || "Heat stress input is incorrect size."
    @test length(mcda_vars.sum_cover) == n_locations || "Initial cover input is incorrect size."
end
@testset "Unguided location selection" begin
    n_intervention_locations = 5
    prefseedlocations = zeros(Int, n_intervention_locations)
    prefshadelocations = zeros(Int, n_intervention_locations)
    seed_years = true
    shade_years = true
    max_cover = [0.0, 3000.0, 5000.0, 0.0, 0.0]
    depth_priority = collect(1:5)

    prefseedlocations, prefshadelocations = ADRIA.unguided_location_selection(prefseedlocations, prefshadelocations, true, true, 5, max_cover, depth_priority)

    # Check that only two locations are selected (the locations where k > 0.0)
    @test length(prefseedlocations[prefseedlocations.>0]) == 2
    @test length(prefshadelocations[prefshadelocations.>0]) == 2

    @test all([in(sid, [2, 3]) for sid in prefseedlocations[prefseedlocations.>0]])
    @test all([in(sid, [2, 3]) for sid in prefshadelocations[prefshadelocations.>0]])
end

@testset "Guided location selection without ADRIA ecological model" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    N = 2^3
    criteria_df = ADRIA.sample_location_selection(dom, N)  # get scenario dataframe

    area_to_seed = 662.11  # area of seeded corals in m^2
    ts = 5  # time step to perform location selection at

    sum_cover = fill(0.1, N, ADRIA.n_locations(dom))
    ranks = ADRIA.run_location_selection(dom, criteria_df, sum_cover, area_to_seed, ts)

    @test size(ranks, 1) == sum(criteria_df.guided .> 0) || "Specified number of scenarios was not carried out."
    @test size(ranks, 2) == length(dom.location_ids) || "Ranks storage is not correct size for this domain."

    sel_locations = unique(ranks)
    sel_locations = sel_locations[sel_locations.!=0.0]
    possible_ranks = collect(Float64, 1:ADRIA.n_locations(dom)+1.0)

    @test all([in(ss, possible_ranks) for ss in sel_locations]) || "Impossible rank assigned."
end
