using Test
using DataFrames
using Statistics
using ADRIAanalysis

@testset "sensitivity.rsa" begin

    # ------------------------------------------------------------------
    # 1. Basic return schema
    # ------------------------------------------------------------------
    @testset "return schema" begin
        X = DataFrame(; a=1.0:10.0, b=rand(10), c=rand(10))
        y = collect(1.0:10.0)
        result = ADRIAanalysis.sensitivity.rsa(X, y)
        @test result isa DataFrame
        @test names(result) ==
            ["feature", "test", "statistic", "prob_superiority", "effect_size"]
        @test nrow(result) == 3
        @test eltype(result.feature) == Symbol
        @test eltype(result.test) == Symbol
        @test eltype(result.statistic) == Float64
        @test eltype(result.prob_superiority) == Float64
        @test eltype(result.effect_size) == Float64
        @test all(result.test .== :mann_whitney)
    end

    # ------------------------------------------------------------------
    # 2. Output is sorted descending by prob_superiority
    # ------------------------------------------------------------------
    @testset "sorted descending by prob_superiority" begin
        X = DataFrame(; a=1.0:20.0, b=rand(20), c=rand(20))
        y = collect(1.0:20.0)
        result = ADRIAanalysis.sensitivity.rsa(X, y)
        @test issorted(result.prob_superiority; rev=true)
    end

    # ------------------------------------------------------------------
    # 8. Input DataFrame not mutated
    # ------------------------------------------------------------------
    @testset "input DataFrame not mutated" begin
        X = DataFrame(; a=rand(30), b=rand(30))
        y = rand(30)
        X_copy = copy(X)
        ADRIAanalysis.sensitivity.rsa(X, y)
        @test X == X_copy
    end

    # ------------------------------------------------------------------
    # 9. Known-discriminating feature ranks first
    # ------------------------------------------------------------------
    @testset "discriminating factor ranks first" begin
        n = 100
        signal = vcat(zeros(50), ones(50))
        noise = rand(n)
        X = DataFrame(; signal=signal, noise=noise)
        y = signal .+ 0.01 .* randn(n)
        result = ADRIAanalysis.sensitivity.rsa(X, y)
        @test result.feature[1] == :signal
    end

    # ------------------------------------------------------------------
    # 4. prob_superiority and effect_size are in valid ranges
    # ------------------------------------------------------------------
    @testset "prob_superiority in [0,1], effect_size in [-1,1]" begin
        X = DataFrame(; a=rand(50), b=rand(50))
        y = rand(50)
        result = ADRIAanalysis.sensitivity.rsa(X, y)
        @test all(0.0 .<= result.prob_superiority .<= 1.0)
        @test all(-1.0 .<= result.effect_size .<= 1.0)
    end

    # ------------------------------------------------------------------
    # 5. Zero-variance column returns sentinel values with @warn
    # ------------------------------------------------------------------
    @testset "zero-variance column returns sentinel + warns" begin
        X = DataFrame(; flat=ones(20), varying=rand(20))
        y = rand(20)
        result = @test_logs (:warn, r"zero variance") ADRIAanalysis.sensitivity.rsa(
            X, y
        )
        flat_row = filter(:feature => f -> f == :flat, result)
        @test nrow(flat_row) == 1
        @test flat_row.statistic[1] == 0.0
        @test flat_row.prob_superiority[1] == 0.5
        @test flat_row.effect_size[1] == 0.0
    end

    # ------------------------------------------------------------------
    # 6. NaN input raises ArgumentError
    # ------------------------------------------------------------------
    @testset "NaN input raises ArgumentError" begin
        X = DataFrame(; a=[1.0, NaN, 3.0, 4.0])
        y = [1.0, 2.0, 3.0, 4.0]
        @test_throws ArgumentError ADRIAanalysis.sensitivity.rsa(X, y)
    end

    # ------------------------------------------------------------------
    # 7. All-true or all-false mask raises ArgumentError
    # ------------------------------------------------------------------
    @testset "degenerate mask raises ArgumentError" begin
        X = DataFrame(; a=rand(10))
        all_true = trues(10)
        all_false = falses(10)
        @test_throws ArgumentError ADRIAanalysis.sensitivity.rsa(X, all_true)
        @test_throws ArgumentError ADRIAanalysis.sensitivity.rsa(X, all_false)
    end

    # ------------------------------------------------------------------
    # 8. Pre-computed BitVector mask dispatch matches threshold dispatch
    # ------------------------------------------------------------------
    @testset "mask dispatch matches threshold dispatch" begin
        n = 80
        X = DataFrame(; a=rand(n), b=rand(n))
        y = rand(n)
        threshold = quantile(y, 0.9)
        mask = y .> threshold
        result_y = ADRIAanalysis.sensitivity.rsa(X, y; top_proportion=0.9)
        result_mask = ADRIAanalysis.sensitivity.rsa(X, mask)
        # Same rows, same order (both sorted by prob_superiority)
        @test result_y.feature == result_mask.feature
        @test result_y.statistic == result_mask.statistic
    end

    # ------------------------------------------------------------------
    # 9. top_proportion parameter controls threshold correctly
    # ------------------------------------------------------------------
    @testset "top_proportion alters selection group" begin
        n = 100
        X = DataFrame(; a=collect(1.0:n))
        y_sorted = collect(1.0:n)
        # With top_proportion=0.5 half the scenarios are selected;
        # with top_proportion=0.9 only 10% are. The statistic should differ.
        r50 = ADRIAanalysis.sensitivity.rsa(X, y_sorted; top_proportion=0.5)
        r90 = ADRIAanalysis.sensitivity.rsa(X, y_sorted; top_proportion=0.9)
        # Both should return a row for :a, but statistic values may differ
        @test r50.feature[1] == :a
        @test r90.feature[1] == :a
        @test r50.statistic[1] != r90.statistic[1]
    end
end
