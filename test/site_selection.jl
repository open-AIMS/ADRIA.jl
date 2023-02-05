using Distributions


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
    dom = ADRIA.load_domain(joinpath(here, "Example_domain"))

    criteria_df = ADRIA.sample(ex_domain, 5)  # get scenario dataframe

    area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
    ts = 5  # time step to perform site selection at

    sumcover = 0.1 .* ones(5, size(dom.site_data, 1))
    ranks = run_site_selection(dom, criteria_df[criteria_df.guided.>0, :], sumcover, area_to_seed, ts)

    # Check that only 4 sites make it through depth and heat/wave risk filter  
    for nscen in 1:size(ranks, 1)
        @test all(ranks["scen_$nscen", 5, :] .== 0.0) || "Sites which should have been filtered in scenario $nscen have still been ranked."
    end

    @test size(ranks, 1) == sum(criteria_df.guided .> 0) || "Specified number of scenarios was not carried out."

    sel_sites = unique(ranks)
    sel_sites = sel_sites[sel_sites.!=0.0]
    site_ids = collect(Float64, 1:size(dom.site_data, 1))
    for ss in sel_sites
        @test in.(ss, site_ids) || "Sites outsite of site id set have been ranked."
    end

end
