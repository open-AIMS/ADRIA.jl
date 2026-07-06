using ADRIA: DataCube
using ADRIAviz
using DataFrames
using Random

@testset "sensitivity" begin
    @testset "pawn(Si) returns Figure" begin
        # Fixture matches the structure returned by ADRIAanalysis.sensitivity.pawn():
        # dimensions are :factors (Symbol) × :Si (Symbol col names)
        factor_names = Symbol.(["Factor1", "Factor2", "Factor3"])
        Si_cols = [:min, :lb, :mean, :median, :ub, :max, :std, :cv]
        Si = DataCube(
            rand(length(factor_names), length(Si_cols));
            factors=factor_names, Si=Si_cols
        )

        @test ADRIA.viz.pawn(Si) isa Figure
    end

    @testset "rsa(X, y, foi) returns Figure" begin
        Random.seed!(1)
        X = DataFrame(; a=rand(30), b=rand(30), c=rand(30))
        y = rand(30)
        @test ADRIA.viz.rsa(X, y, (:a, :b)) isa Figure
    end

    @testset "outcome_map(X, y, factor) returns Figure" begin
        Random.seed!(2)
        X = DataFrame(; a=rand(30), b=rand(30))
        y = rand(30)
        @test ADRIA.viz.outcome_map(X, y, :a) isa Figure
    end

    @testset "outcome_map(X, y, factors) returns Figure" begin
        Random.seed!(3)
        X = DataFrame(; a=rand(30), b=rand(30), c=rand(30))
        y = rand(30)
        @test ADRIA.viz.outcome_map(X, y, [:a, :b, :c]) isa Figure
    end
end
