using YAXArrays
using ADRIA: DataCube, ZeroDataCube

@testset "Relative Shelter Volume" begin
    # Create scenario spec
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH)

    # Create scenario spec
    test_scens = ADRIA.sample(dom, 8)

    # Create dummy coral covers and make sum of values <= 1.0 (i.e., proportional to area)
    coral_cover = DataCube(rand(5, 35, 3, 1), (:timesteps, :species, :sites, :scenarios))
    coral_cover .= coral_cover ./ sum(coral_cover; dims=:species)

    k_area = Float64[70, 60, 50]  # in m²

    r_sv = ADRIA.metrics.relative_shelter_volume(
        coral_cover, k_area, DataFrame(test_scens[1, :])
    )
    @test all(0.0 .<= r_sv .<= 1.0)
    @test any(r_sv .>= 0.05)  # warn if all values ae very tiny values (catch Issue #91 : https://github.com/open-AIMS/ADRIA.jl/issues/91)

    # Test multi-scenario case
    coral_cover = DataCube(rand(5, 35, 3, 5), (:timesteps, :species, :sites, :scenarios))
    coral_cover .= coral_cover ./ sum(coral_cover; dims=:species)
    r_sv = ADRIA.metrics.relative_shelter_volume(
        coral_cover, k_area, DataFrame(test_scens[1:5, :])
    )

    @test all(0.0 .<= r_sv .<= 1.0) ||
        "Min CC: $(minimum(sum(coral_cover, dims=:species))); Max CC: $(maximum(sum(coral_cover, dims=:species))) | $((minimum(r_sv), maximum(r_sv)))"
    @test any(r_sv .>= 0.05)

    # Test zero value case
    coral_cover = ZeroDataCube((:timesteps, :species, :sites, :scenarios), (5, 35, 3, 5))
    r_sv = ADRIA.metrics.relative_shelter_volume(
        coral_cover, k_area, DataFrame(test_scens[1:5, :])
    )
    @test all(r_sv .== 0.0)

    # Maximum shelter volume case
    coral_cover = DataCube(rand(5, 35, 3, 5), (:timesteps, :species, :sites, :scenarios))
    coral_cover[species=24, sites=1:3] .= k_area'  # Coral type with maximum shelter density
    r_sv = ADRIA.metrics.relative_shelter_volume(
        coral_cover, k_area, DataFrame(test_scens[1:5, :])
    )
    @test all(r_sv .== 1.0) ||
        "Scenario with complete coral cover does not achieve max RSV | $(maximum(r_sv))"
end

# @testset "metric modifications" begin
#     dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)
#     sa = site_area(dom)
#     tac_metric = @extend_metric(TAC, total_absolute_cover, [sa])

#     @test tac_metric(rand(5, 26, 5, 5))
# end
