@testset "Shelter Volume" begin
    scen_path = joinpath(TEST_DATA_DIR, "test_scenarios.csv")
    test_scens = CSV.read(scen_path, DataFrame)

    r_sv = ADRIA.metrics.relative_shelter_volume(rand(5,5,5,5,5), [100, 100, 100, 100, 100], test_scens[1, :]);

    @test all(r_sv .<= 1.0)
    @test any(r_sv .>= 0.05)  # warn if all values are very tiny values (catch Issue #91 : https://github.com/open-AIMS/ADRIA.jl/issues/91)
end
