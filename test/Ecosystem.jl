using Test
using ADRIA
using ADRIA.Distributions
using ADRIA.Random

@testset "Truncated normal mean" begin
    n_checks = 1000
    mean_diffs::Vector{Float64} = zeros(n_checks)

    # The truncted normal mean should agree with the normal mean for symmetrical bound about mean
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

    # The calculated truncated normal mean should agree with the Distributions.jl implementation
    for i in 1:n_checks
        mu = rand(Uniform(0, 10))
        stdev = rand(Uniform(0.01, 10))
        lb = rand(Uniform(mu - 6.0 * stdev, mu + 3.0 * stdev))
        ub = rand(Uniform(lb, lb + 3.0 * stdev))

        calculated = ADRIA.truncated_normal_mean(
            mu, stdev, lb, ub
        )
        expected = mean(truncated(Normal(mu, stdev), lb, ub))
        mean_diffs[i] = abs(expected - calculated)
    end

    @test all(mean_diffs .< 1e-4) ||
        "calculated truncated normal mean differs signficantly from Distributions.jl"
end

@testset "Truncated normal cdf" begin
    n_checks::Int = 1000
    cdf_diffs::Vector{Float64} = zeros(3 * n_checks)

    # The truncated normal cdf should return 0.0, 0.5, and 1.0 when evaluated at the lower
    # bound, median and upper bound respectively.
    for i in 1:n_checks
        mu = rand(Uniform(0, 10))
        stdev = rand(Uniform(0.01, 10))

        n_std = rand(Uniform(0, 10))
        lb = mu - n_std * stdev
        ub = mu + n_std * stdev

        cdf_diffs[3 * (i - 1) + 1] = abs(
            0.0 - ADRIA.truncated_normal_cdf(
                lb, mu, stdev, lb, ub
            )
        )
        cdf_diffs[3 * (i - 1) + 2] = abs(
            0.5 - ADRIA.truncated_normal_cdf(
                mu, mu, stdev, lb, ub
            )
        )
        cdf_diffs[3 * (i - 1) + 3] = abs(
            1.0 - ADRIA.truncated_normal_cdf(
                ub, mu, stdev, lb, ub
            )
        )
    end

    @test all(cdf_diffs .< 1e-7) ||
        "truncated normal cdf not equal to 0.0, 0.5, 1.0 at lower bound, mean, upper bound"

    n_checks = 3000

    # The truncated normal cdf should agree with the Distributions.jl implementation
    for i in 1:n_checks
        mu = rand(Uniform(0, 10))
        stdev = rand(Uniform(0.01, 10))
        lb = rand(Uniform(mu - stdev * 6.0, mu + stdev * 3.0))
        ub = rand(Uniform(lb, lb + stdev * 3.0))

        x = rand(Uniform(lb, ub))

        calculated = ADRIA.truncated_normal_cdf(
            x, mu, stdev, lb, ub
        )
        expected = cdf(
            truncated(Normal(mu, stdev), lb, ub), x
        )
        cdf_diffs[i] = abs(expected - calculated)
        if (cdf_diffs[i] > 1e-4)
            println("diff: $(cdf_diffs[i])")
            println("mu: $(mu), stdev: $(stdev), lb: $(lb), ub: $(ub)")
        end
    end

    @test all(cdf_diffs .< 1e-4) ||
        "Implemented truncated normal cdf differs significantly from built-in cdf"
end
