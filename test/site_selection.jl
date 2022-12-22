@testset "site selection" begin
    # TODO: Complete tests with @tests

    dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)
    p_tbl = ADRIA.param_table(dom)

    p_tbl[:, :depth_offset] .= 7.0

    # ranks = ADRIA.site_selection(dom, p_tbl, 1, 10, 1)
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
    criteria_df = ADRIA.sample(dom, 5) # get scenario dataframe
    criteria_df.dist_thresh .= 1.0
    area_to_seed = 1.5 * 10^-6 # area of seeded corals in km^2
    nreps = 30 # number of dhw and wave replicates you want to use
    ts = 5 # time step to perform site selection at
    scen = 5
    alg_ind = criteria_df.guided[scen]# MCDA algorithm to use (1-3)

    ranks = site_selection(dom, criteria_df, area_to_seed, ts, nreps, scen, alg_ind)

    # Check that only 4 sites make it through depth and heat/wave risk filter    
    @test size(ranks, 2) == 4 || "Sites which should have been filtered have still been ranked."
    @test size(ranks, 1) == nreps || "Specified number of replicates was not carried out."
    @test sum(ranks[:, :, [2, 3]]) .!= 0.0 || "No ranks assigned for any replicates."

end
