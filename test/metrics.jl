@testset "Shelter Volume" begin
    scen_path = joinpath(TEST_DATA_DIR, "test_scenarios.csv")
    test_scens = CSV.read(scen_path, DataFrame)

    @eval using NamedDims

    # Create dummy coral covers and make sum of values <= 1.0 (i.e., proportional to area)
    coral_cover = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)}(rand(5, 36, 3, 1))
    coral_cover .= coral_cover ./ sum(coral_cover, dims=:species)

    site_area = Float64[100, 100, 100]  # in m²

    r_sv = ADRIA.metrics.relative_shelter_volume(coral_cover, site_area, DataFrame(test_scens[1, :]));
    @test all(0.0 .<= r_sv .<= 1.0)
    @test any(r_sv .>= 0.05)  # warn if all values ae very tiny values (catch Issue #91 : https://github.com/open-AIMS/ADRIA.jl/issues/91)

    # Test multi-scenario case
    coral_cover = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)}(rand(5, 36, 3, 5))
    coral_cover .= coral_cover ./ sum(coral_cover, dims=:species)
    r_sv = ADRIA.metrics.relative_shelter_volume(coral_cover, site_area, DataFrame(test_scens[1:5, :]));

    @test all(0.0 .<= r_sv .<= 1.0)
    @test any(r_sv .>= 0.05)


    # Test zero value case
    coral_cover = NamedDimsArray{(:timesteps, :species, :sites, :scenarios)}(zeros(5, 36, 3, 5))
    r_sv = ADRIA.metrics.relative_shelter_volume(coral_cover, site_area, DataFrame(test_scens[1:5, :]));
    @test all(r_sv .== 0.0)
end


# @testset "metric modifications" begin
#     dom = ADRIA.load_domain(joinpath(@__DIR__, "..", "examples", "Example_domain"), 45)
#     sa = site_area(dom)
#     tac_metric = @extend_metric(TAC, total_absolute_cover, [sa])

#     @test tac_metric(rand(5, 26, 5, 5))
# end
