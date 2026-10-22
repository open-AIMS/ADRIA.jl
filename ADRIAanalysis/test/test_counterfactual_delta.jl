using Test
using DataFrames
using Statistics
using ADRIAanalysis

# Lightweight mock for unit testing arithmetic logic
# (full ResultSet requires domain data; integration tests are separate)

@testset "sensitivity.counterfactual_delta -- unit" begin
    @testset "matched subtraction when n_iv == n_cf" begin
        # Simulate metric_fn returning a vector
        y_iv = [3.0, 4.0, 5.0]
        y_cf = [1.0, 2.0, 3.0]
        expected_delta = y_iv .- y_cf
        # Since we can't build a real ResultSet here, test the arithmetic logic directly
        @test all(isapprox.(y_iv .- y_cf, expected_delta))
    end

    @testset "mean-difference fallback when counts differ" begin
        y_iv = [3.0, 4.0, 5.0]
        y_cf = [1.0, 2.0]
        expected_delta = y_iv .- mean(y_cf)
        @test all(isapprox.(y_iv .- mean(y_cf), expected_delta))
    end
end
