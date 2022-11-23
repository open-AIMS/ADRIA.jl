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
    depth_priority = [1 3 4 5]

    prefseedsites, prefshadesites = ADRIA.unguided_site_selection(prefseedsites, prefshadesites, true, true, 5, max_cover, depth_priority)

    # Check that only two sites are selected (the sites where k > 0.0)
    @test length(prefseedsites[prefseedsites.>0]) == 2
    @test length(prefshadesites[prefshadesites.>0]) == 2

    @test all([in(sid, [2, 3]) for sid in prefseedsites[prefseedsites.>0]])
    @test all([in(sid, [2, 3]) for sid in prefshadesites[prefshadesites.>0]])
end
