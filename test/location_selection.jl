using Distributions
using Test
using NamedDims


if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const EXAMPLE_DOMAIN_PATH = joinpath(ADRIA_DIR, "examples", "Example_domain")
end

@testset "location selection" begin

    n_location_int = 5
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    scen = ADRIA.sample_location_selection(dom, 1)
    scen = NamedDimsArray(Vector(scen[1, :]), factors=names(scen))

    area_to_seed = 962.11

    location_ids = collect(1:length(dom.location_ids))
    criteria_store = ADRIA.create_criteria_store(location_ids, dom.mcda_criteria)
    tolerances = (iv__coral_space=(>, scen("iv__coral_space__tol") .* area_to_seed),
        iv__heat_stress=(>, 1 - scen("iv__heat_stress__tol")),
        iv__wave_stress=(>, 1 - scen("iv__wave_stress__tol")))

    dhw_scens = dom.dhw_scens
    wave_scens = dom.wave_scens
    criteria_store(:iv__wave_stress) .= ADRIA.env_stress_criteria(Array(dropdims(ADRIA.mean(wave_scens, dims=(:timesteps, :scenarios)) .+ ADRIA.var(wave_scens, dims=(:timesteps, :scenarios)), dims=:timesteps)))
    criteria_store(:iv__heat_stress) .= ADRIA.env_stress_criteria(Array(dropdims(ADRIA.mean(dhw_scens, dims=(:timesteps, :scenarios)) .+ ADRIA.var(dhw_scens, dims=(:timesteps, :scenarios)), dims=:timesteps)))

    int_logs = NamedDimsArray([true true false], scenarios=[1], log=[:seed, :fog, :shade])

    ranks = ADRIA.location_selection(criteria_store, dom.interventions, scen, tolerances, int_logs[scenarios=1],
        location_ids, dom.location_distances, dom.median_location_distance, n_location_int)

    @test all([in(r, collect(0:length(dom.location_ids)+1)) for r in ranks[:seed][:, 2]]) || "Some seed ranks outside of possible ranks (0-n_locations+1)"
    @test all([in(r, collect(0:length(dom.location_ids)+1)) for r in ranks[:fog][:, 2]]) || "Some shade ranks outside of possible ranks (0-n_locations+1)"
    @test (length(ranks[:seed][:, 2]) - length(unique(sort(ranks[:seed][:, 2])))) == sum(ranks[:seed][:, 2] .== 0) - 1 || " Some filtered locations not set to zero."

end
@testset "MCDA variable in Domain" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    n_locations = length(dom.location_ids)
    mcda_criteria = dom.mcda_criteria
    for crit in keys(mcda_criteria)
        @test (length(mcda_criteria[crit]) == n_locations) || "The criteria $crit has the wrong size."
    end
end
@testset "Unguided location selection" begin
    n_intervention_locations = 5
    int_logs = NamedDimsArray([true true false], scenarios=[1], log=[:seed, :fog, :shade])
    pref_locations = (seed=zeros(Int, n_intervention_locations), fog=shade = zeros(Int, n_intervention_locations), shade=zeros(Int, n_intervention_locations))
    max_cover = [0.0, 3000.0, 5000.0, 0.0, 0.0]
    depth_priority = collect(1:5)
    clusters = [1.0, 2.0, 3.0, 4.0, 5.0]
    pref_locations_new = ADRIA.unguided_location_selection(pref_locations, int_logs[scenarios=1], 5, max_cover, depth_priority, clusters)

    # Check that only two locations are selected (the locations where k > 0.0)
    @test length(pref_locations.seed[pref_locations.seed.>0]) == 2
    @test length(pref_locations.fog[pref_locations.fog.>0]) == 2

    @test all([in(sid, [2, 3]) for sid in pref_locations.seed[pref_locations.seed.>0]])
    @test all([in(sid, [2, 3]) for sid in pref_locations.shade[pref_locations.shade.>0]])
end

@testset "Guided location selection without ADRIA ecological model" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    N = 2^3
    criteria_df = ADRIA.sample_location_selection(dom, N)  # get scenario dataframe

    area_to_seed = 962.11  # area of seeded corals in m^2
    # define functions for tolerances
    f_coral_cover(param) = area_to_seed * param

    location_ids = collect(1:length(dom.location_ids))

    coral_cover = NamedDims.rename(repeat(sum(dom.init_coral_cover, dims=:species), size(criteria_df, 1)), (:scenarios, :locations))
    # initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of locations))
    tolerances = (iv__coral_space=(>, x -> f_coral_cover(x)),
        iv__heat_stress=(>, x -> 1 - x),
        iv__wave_stress=(>, x -> 1 - x))

    ranks = ADRIA.run_location_selection(dom, criteria_df, tolerances, coral_cover', target_seed_locations=location_ids, target_shade_locations=location_ids)

    @test size(ranks, 1) == sum(criteria_df.guided .> 0) || "Specified number of scenarios was not carried out."
    @test size(ranks, 2) == length(dom.location_ids) || "Ranks storage is not correct size for this domain."

    sel_locations = unique(ranks)
    sel_locations = sel_locations[sel_locations.!=0.0]
    possible_ranks = collect(Float64, 1:ADRIA.n_locations(dom)+1.0)

    @test all([in(ss, possible_ranks) for ss in sel_locations]) || "Impossible rank assigned."
end
