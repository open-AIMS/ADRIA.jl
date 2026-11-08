using Test
using Random
using ADRIAanalysis

@testset "pawn" begin
    rng = Xoshiro(1)
    X = rand(rng, 200, 3)
    y = 5 .* X[:, 1] .+ 0.1 .* rand(rng, 200)
    Si = ADRIAanalysis.sensitivity.pawn(X, y, ["a", "b", "c"]; S=5)

    @test size(Si) == (3, 8)
    # Plain arrays reach callers (a static vector breaks Makie >= 0.24.15 plotting)
    @test Si.data isa Matrix{Float64}
    @test vec(Si.data[:, 4]) isa Vector{Float64}
    # The informative factor ranks first
    @test argmax(vec(Si.data[:, 4])) == 1
end
