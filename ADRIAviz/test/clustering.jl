using ADRIAviz
using Random

@testset "clustering" begin
    @testset "scenarios(outcomes, clusters) returns Figure" begin
        # Create synthetic outcomes (timesteps × scenarios)
        Random.seed!(42)
        n_timesteps = 10
        n_scenarios = 12
        outcomes = rand(n_timesteps, n_scenarios)

        # Create cluster assignments (3 clusters)
        clusters = repeat([1, 2, 3, 1], 3)

        @test ADRIA.viz.scenarios(outcomes, clusters) isa Figure
    end

    @testset "clustered_scenarios(outcomes, clusters) returns Figure" begin
        Random.seed!(42)
        n_timesteps = 10
        n_scenarios = 12
        outcomes = rand(n_timesteps, n_scenarios)

        clusters = repeat([1, 2, 3, 1], 3)

        @test ADRIA.viz.clustered_scenarios(outcomes, clusters) isa Figure
    end
end
