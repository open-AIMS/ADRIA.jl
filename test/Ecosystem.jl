using Test
using ADRIA
using ADRIA.Distributions

@testset "Truncated normal mean" begin

    means::Vector{Float64} = Vector{Float64}([0.0, 0.5, 1.5, 2.5, 4.0, 5.0])
    stdev::Vector{Float64} = Vector{Float64}([0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0])
    lower_bounds::Vector{Float64} = Vector{Float64}([0.0, 0.5, 1.5, 3.0, 7.0, 10.0])
    bound_widths::Vector{Float64} = Vector{Float64}([0.25, 0.5, 1.0, 2.0, 5.0, 10.0])

    n_checks = length(means) * length(stdev) * length(bound_widths)
    mean_diffs::Vector{Float64} = Vector{Float64}(undef, n_checks)

    # truncated normal mean should agree with normal mean for symmetrical bounds about mean
    ind::Integer = 1;
    for mu in means
        for std in stdev
            for width in bound_widths
                calculated = ADRIA.truncated_normal_mean(
                    mu, std, mu - width, mu + width
                )
                mean_diffs[ind] = abs(mu - calculated)
                ind += 1
            end
        end
    end

    @test all(mean_diffs .== 0.0) ||
        "calculated truncated normal mean not equal to normal mean for symmetric bounds"
    
    n_checks = length(means) * length(stdev) * length(bound_widths) * length(lower_bounds)
    mean_diffs = Vector{Float64}(undef, n_checks)
    
    ind = 1
    for mu in means
        for std in stdev
            for lower_b in lower_bounds
                for width in bound_widths
                    calculated = ADRIA.truncated_normal_mean(
                        mu, std, lower_b, lower_b + width
                    )
                    expected = mean(truncated(Normal(mu, std), lower_b, lower_b + width))
                    mean_diffs[ind] = abs(expected - calculated)
                    ind += 1
                end
            end
        end
    end

    @test all(mean_diffs .< 1e-7) ||
        "calculated truncated normal mean differs signficantly from Distributions.jl"
end

@testset "Truncated normal cdf" begin
    
    means::Vector{Float64} = Vector{Float64}([0.0, 0.5, 1.5, 2.5, 4.0, 5.0])
    stdevs::Vector{Float64} = Vector{Float64}([0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0])
    n_stdevs::Vector{Float64} = Vector{Float64}([1.0, 2.0, 3.0, 4.0, 5.0])

    n_checks::Int = 3 * length(means) * length(stdevs) * length(n_stdevs)
    cdf_diffs::Vector{Float64} = Vector{Float64}(undef, n_checks)
    
    ind::Int = 1
    for mu in means
        for std in stdevs
            for n_std in n_stdevs
                lower_bound = mu - n_std * std
                upper_bound = mu + n_std * std

                cdf_diffs[ind] = abs(0.0 - ADRIA.truncated_normal_cdf(
                    lower_bound, mu, std, lower_bound, upper_bound
                ))
                ind += 1
                cdf_diffs[ind] = abs(0.5 - ADRIA.truncated_normal_cdf(
                    mu, mu, std, lower_bound, upper_bound
                ))
                ind += 1
                cdf_diffs[ind] = abs(1.0 - ADRIA.truncated_normal_cdf(
                    upper_bound, mu, std, lower_bound, upper_bound
                ))
                ind += 1
            end
        end
    end

    @test all(cdf_diffs .< 1e-7) ||
        "truncated normal cdf not equal to 0.0, 0.5, 1.0 at lower bound, mean, upper bound"

    # evaluate cdf at proportion along domain of truncated normal
    eval_cdf_at::Vector{Float64} = Vector{Float64}([
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
    ])
    
    n_checks = length(means) * length(stdevs) * length(n_stdevs) * 
        length(n_stdevs) * 2 * length(eval_cdf_at)

    cdf_diffs = Vector{Float64}(undef, n_checks)
    
    ind = 1
    for mu in means
        for std in stdevs
            for upper_n_std in n_stdevs
                for lower_n_std in n_stdevs
                    lower_bound = mu - std * lower_n_std
                    upper_bound = lower_bound + std * upper_n_std
                    for p in eval_cdf_at
                        x = lower_bound + p * (upper_bound - lower_bound)
                        inbuilt = cdf(
                            truncated(Normal(mu, std), lower_bound, upper_bound), x
                        )
                        custom_cdf = ADRIA.truncated_normal_cdf(
                            x, mu, std, lower_bound, upper_bound
                        )
                        cdf_diffs[ind] = abs(custom_cdf - inbuilt)
                        ind += 1
                    end

                    lower_bound = mu + std * lower_n_std
                    upper_bound = lower_bound + std * upper_n_std
                    for p in eval_cdf_at
                        x = lower_bound + p * (upper_bound - lower_bound)
                        inbuilt = cdf(
                            truncated(Normal(mu, std), lower_bound, upper_bound), x
                        )
                        custom_cdf = ADRIA.truncated_normal_cdf(
                            x, mu, std, lower_bound, upper_bound
                        )
                        cdf_diffs[ind] = abs(custom_cdf - inbuilt)
                        ind += 1
                    end
                end
            end
        end
    end

    @test all(cdf_diffs .< 1e-7) ||
        "Implemented truncated normal cdf differs significantly from built in cdf"
end
