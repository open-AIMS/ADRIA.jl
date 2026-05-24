using Test, ADRIAanalysis

@testset "Rule struct" begin
    @testset "sorts condition by feature name" begin
        r = ADRIAanalysis.Rule(
            Vector{Vector}([["z", :L, 0.5f0], ["a", :R, 0.3f0]]), [0.7, 0.3]
        )
        @test r.condition[1][1] == "a"
        @test r.condition[2][1] == "z"
    end

    @testset "direction symbol formats correctly in print_rules" begin
        r_L = ADRIAanalysis.Rule(Vector{Vector}([["feat", :L, 0.5f0]]), [0.7, 0.3])
        r_R = ADRIAanalysis.Rule(Vector{Vector}([["feat", :R, 0.5f0]]), [0.7, 0.3])
        @test_nowarn ADRIAanalysis.print_rules([r_L])
        @test_nowarn ADRIAanalysis.print_rules([r_R])
    end

    @testset "empty rule vector" begin
        empty_rules = ADRIAanalysis.Rule{Vector{Vector},Vector{Float64}}[]
        @test_nowarn ADRIAanalysis.print_rules(empty_rules)
        @test ADRIAanalysis.maximum_probability(empty_rules) == 0.0
    end
end

@testset "maximum_probability" begin
    r1 = ADRIAanalysis.Rule(Vector{Vector}([["a", :L, 0.1f0]]), [0.8, 0.2])
    r2 = ADRIAanalysis.Rule(Vector{Vector}([["b", :R, 0.3f0]]), [0.6, 0.4])
    @test ADRIAanalysis.maximum_probability([r1, r2]) ≈ 1.4
    @test ADRIAanalysis.maximum_probability([r1]) ≈ 0.8
end
