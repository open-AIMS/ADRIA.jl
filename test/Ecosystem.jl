using Test
using ADRIA
using ADRIA.Distributions
using ADRIA.Random

@testset "Truncated normal mean" begin
    n_checks = 1000
    mean_diffs::Vector{Float64} = zeros(n_checks)

    # truncated normal mean should agree with normal mean for symmetrical bounds about mean
    for i in 1:n_checks
        mu = rand(Uniform(0, 10))
        stdev = rand(Uniform(0.01, 10))
        width = rand(Uniform(0, 10 * stdev))

        calculated = ADRIA.truncated_normal_mean(
            mu, stdev, mu - width, mu + width
        )
        mean_diffs[i] = abs(mu - calculated)
    end

    @test all(mean_diffs .< 1e-10) ||
        "calculated truncated normal with symmetric bounds mean not equal to normal mean"
    
    mean_diffs = zeros(n_checks)
    
    # check agreement for truncation bounds up to 15 standard deviations from the mean
    for i in 1:n_checks
        mu = rand(Uniform(0, 10))
        stdev = rand(Uniform(0.01, 10))
        lb = rand(Uniform(mu - 15 * stdev, mu + 15 * stdev))
        ub = rand(Uniform(lb, lb + 15 * stdev))
        
        calculated = ADRIA.truncated_normal_mean(
            mu, stdev, lb, ub
        )
        expected = mean(truncated(Normal(mu, stdev), lb, ub))
        mean_diffs[i] = abs(expected - calculated)
    end
    
    @test all(mean_diffs .< 1e-7) ||
        "calculated truncated normal mean differs signficantly from Distributions.jl"
end

@testset "Truncated normal cdf" begin
    n_checks::Int = 1000
    cdf_diffs::Vector{Float64} = zeros(3 * n_checks)
    
    # checking basic properties cdfs, 0 at lower bound, 0.5 at median, 1.0 at upper bound
    for i in 1:n_checks
        mu = rand(Uniform(0, 10))
        stdev = rand(Uniform(0.01, 10))

        n_std = rand(Uniform(0, 10))
        lb = mu - n_std * stdev
        ub = mu + n_std * stdev

        cdf_diffs[3 * (i - 1) + 1] = abs(0.0 - ADRIA.truncated_normal_cdf(
            lb, mu, stdev, lb, ub
        ))
        cdf_diffs[3 * (i - 1) + 2] = abs(0.5 - ADRIA.truncated_normal_cdf(
            mu, mu, stdev, lb, ub
        ))
        cdf_diffs[3 * (i - 1) + 3] = abs(1.0 - ADRIA.truncated_normal_cdf(
            ub, mu, stdev, lb, ub
        ))
    end

    @test all(cdf_diffs .< 1e-7) ||
        "truncated normal cdf not equal to 0.0, 0.5, 1.0 at lower bound, mean, upper bound"

    n_checks = 3000
    
    # check agreement with Distributions.jl implementation with bounds at most 
    # 10 standard deviations from mean
    for i in 1:n_checks
        mu = rand(Uniform(0, 10))
        stdev = rand(Uniform(0.01, 10))
        lb = rand(Uniform(mu - stdev * 10, mu + stdev * 5))
        ub = rand(Uniform(lb, lb + stdev * 5))
    
        x = rand(Uniform(lb, ub))

        calculated = ADRIA.truncated_normal_cdf(
            x, mu, stdev, lb, ub
        )
        expected = cdf(
            truncated(Normal(mu, stdev), lb, ub), x
        )
        cdf_diffs[i] = abs(expected - calculated)
    end

    @test all(cdf_diffs .< 1e-7) ||
        "Implemented truncated normal cdf differs significantly from built-in cdf"
end
