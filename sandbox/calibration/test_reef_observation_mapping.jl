using Test
include(joinpath(@__DIR__, "reef_observation_mapping.jl"))

@testset "Area-weighted reef observation contract" begin
    names = Dict("A" => "Reef A", "B" => "Reef B")
    mapping = reef_area_weights(
        ["s3", "s1", "s2"], ["B", "A", "A"], [5.0, 1.0, 3.0],
        names, ["Reef A", "Reef B"]
    )
    @test mapping["Reef A"].indices == [2, 3]
    @test mapping["Reef A"].weights ≈ [0.25, 0.75]
    @test mapping["Reef A"].area_m2 == 4.0
    densities = [100.0 2.0 6.0; 200.0 4.0 8.0]
    @test aggregate_reef_density(densities, mapping["Reef A"]) ≈ [5.0, 7.0]
    @test aggregate_reef_density(densities, mapping["Reef B"]) == [100.0, 200.0]
    @test_throws ErrorException reef_area_weights(
        ["s1", "s1"], ["A", "A"], [1.0, 2.0], names, ["Reef A"])
    @test_throws ErrorException reef_area_weights(
        ["s1"], ["A"], [0.0], names, ["Reef A"])
    @test_throws ErrorException reef_area_weights(
        ["s1"], ["missing"], [1.0], names, ["Reef A"])
    @test_throws ErrorException reef_area_weights(
        ["s1"], ["A"], [1.0], names, ["Reef C"])
end
