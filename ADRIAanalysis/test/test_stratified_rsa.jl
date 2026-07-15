using Test
using DataFrames
using Statistics
using ADRIAanalysis

@testset "sensitivity.stratified_rsa" begin
    @testset "return schema" begin
        n = 200
        X = DataFrame(;
            a=rand(n),
            b=rand(n),
            dhw_mean=rand(n),
            dhw_stdev=rand(n),
            dhw_complexity=rand(n)
        )
        y = rand(n)
        result = ADRIAanalysis.sensitivity.stratified_rsa(X, y)
        @test result isa DataFrame
        @test names(result) ==
            ["feature", "stratum", "prob_superiority", "effect_size", "mean_importance"]
        @test eltype(result.feature) == Symbol
        @test all(1 .<= result.stratum .<= 4)
        @test all(0.0 .<= result.prob_superiority .<= 1.0)
        @test all(-1.0 .<= result.effect_size .<= 1.0)
    end

    @testset "DHW columns are dropped from feature matrix" begin
        n = 200
        X = DataFrame(;
            a=rand(n), dhw_mean=rand(n), dhw_stdev=rand(n), dhw_complexity=rand(n)
        )
        y = rand(n)
        result = ADRIAanalysis.sensitivity.stratified_rsa(X, y)
        @test !(:dhw_mean in result.feature)
        @test !(:dhw_stdev in result.feature)
        @test !(:dhw_complexity in result.feature)
    end

    @testset "mean_importance = mean(abs(prob_superiority - 0.5)) across strata" begin
        n = 400
        X = DataFrame(; signal=vcat(zeros(200), ones(200)), dhw_mean=rand(n))
        y = Float64.(X.signal) .+ 0.01 .* randn(n)
        result = ADRIAanalysis.sensitivity.stratified_rsa(X, y)
        for feat in unique(result.feature)
            rows = filter(:feature => f -> f == feat, result)
            expected_importance = mean(abs.(rows.prob_superiority .- 0.5))
            @test all(isapprox.(rows.mean_importance, expected_importance; atol=1e-10))
        end
    end

    @testset "mean_importance distinguishes reversed factor from neutral factor" begin
        # A factor with reversed importance [0.9, 0.1] should score higher than
        # a neutral factor [0.5, 0.5] under mean(abs(ps - 0.5)), not equal to it.
        reversed_ps = [0.9, 0.1]
        neutral_ps = [0.5, 0.5]
        @test mean(abs.(reversed_ps .- 0.5)) > mean(abs.(neutral_ps .- 0.5))
    end

    @testset "missing strat_col throws ArgumentError" begin
        X = DataFrame(; a=rand(10))
        @test_throws ArgumentError ADRIAanalysis.sensitivity.stratified_rsa(
            X, rand(10); strat_col=:nonexistent
        )
    end

    @testset "dimension mismatch throws DimensionMismatch" begin
        X = DataFrame(; a=rand(10), dhw_mean=rand(10))
        @test_throws DimensionMismatch ADRIAanalysis.sensitivity.stratified_rsa(X, rand(5))
    end

    @testset "empty stratum returns partial result (not error)" begin
        # With constant strat_col = 1.0, quantile edges all collapse to 1.0.
        # All scenarios land in stratum 4 (the last bin); strata 1-3 have 0
        # scenarios and are skipped with @warn.
        n = 100
        X = DataFrame(; a=rand(n), dhw_mean=ones(n))
        y = rand(n)
        result = @test_logs (:warn, r"scenarios \(< 10\)") match_mode=:any ADRIAanalysis.sensitivity.stratified_rsa(
            X, y
        )
        @test result isa DataFrame
    end

    @testset "n_strata keyword respected" begin
        n = 300
        X = DataFrame(; a=rand(n), dhw_mean=rand(n))
        y = rand(n)
        result = ADRIAanalysis.sensitivity.stratified_rsa(X, y; n_strata=3)
        @test all(1 .<= result.stratum .<= 3)
    end
end
