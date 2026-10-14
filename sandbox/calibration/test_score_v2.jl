using Test
using Statistics
include(joinpath(@__DIR__, "cots_cycle_metrics.jl"))
include(joinpath(@__DIR__, "score_v2.jl"))

# v1 total rebuilt from components with a substituted peak-count penalty (v1 default weights).
rebuild(s, pc) = 0.25 * s.matched_rmse + 0.001 * s.percent_bias_abs + 1.25 * pc +
    2.0 * s.peak_timing_penalty + 0.75 * s.period_penalty + 2.0 * s.amplitude_penalty +
    2.0 * s.flatline_penalty + 0.3 * s.lag_correlation_penalty

@testset "Frozen score v2 wrapper" begin
    @test score_v2_primary("North Direction Reef", "raw")
    @test !score_v2_primary("North Direction Reef", "smoothed_3y")
    for reef in ("Lizard Island Reef", "MacGillivray Reef", "Eyrie Reef")
        @test score_v2_primary(reef, "smoothed_3y") && !score_v2_primary(reef, "raw")
    end

    # Peak-count penalty: early peaks (< 2013 - 6) are dropped; an early peak cannot replace 2013.
    @test score_v2_peak_count_penalty([2000, 2014], [2013]) == (0.0, 1)
    @test score_v2_peak_count_penalty([2000], [2013]) == (1.0, 0)
    @test score_v2_peak_count_penalty([2007, 2014], [2013]) == (1.0, 2)   # 2007 is observable
    @test score_v2_peak_count_penalty(Int[], Int[]) == (0.0, 0)

    # Other reefs are identical to v1 whatever their peaks.
    @test score_v2_loss("North Direction Reef", 3.882, 0.0, [1999, 2013], [1995, 2013]) === 3.882
    @test score_v2_loss("Lizard Island Reef", 1.771, 1.0, [1990, 1997, 2012], [1997, 2013]) === 1.771

    # End to end on a synthetic Eyrie: triennial observations with only the 2013 peak,
    # model waves in 2000 and 2014.
    years = collect(1992:2024)
    sim = [exp(-((y - 2000) / 2.0)^2) + 0.6 * exp(-((y - 2014) / 2.0)^2) + 0.01 for y in years]
    obs_years = [1995, 1998, 2001, 2004, 2007, 2013, 2015, 2021, 2024]
    obs = [0.02, 0.13, 0.09, 0.00, 0.00, 0.85, 0.08, 0.00, 0.09]
    s = cots_cycle_score(years, sim, obs_years, obs; fixed_observation_scale=1.0)
    @test s.sim_peak_years == [2000, 2014] && s.obs_peak_years == [2013]
    @test s.peak_count_penalty == 1.0
    @test rebuild(s, s.peak_count_penalty) ≈ s.total_loss atol = 1e-12
    v2 = score_v2_loss("Eyrie Reef", s.total_loss, s.peak_count_penalty, s.sim_peak_years, s.obs_peak_years)
    @test v2 ≈ rebuild(s, 0.0) atol = 1e-12
    @test v2 ≈ s.total_loss - 1.25 atol = 1e-12

    # Control: no early model peak -> v2 equals v1 exactly.
    late = [0.6 * exp(-((y - 2014) / 2.0)^2) + 0.01 for y in years]
    s2 = cots_cycle_score(years, late, obs_years, obs; fixed_observation_scale=1.0)
    @test score_v2_loss("Eyrie Reef", s2.total_loss, s2.peak_count_penalty,
        s2.sim_peak_years, s2.obs_peak_years) === s2.total_loss

    # Early peak only: v2 penalises the missing 2013 wave (worse than v1).
    early = [exp(-((y - 2000) / 2.0)^2) + 0.01 for y in years]
    s3 = cots_cycle_score(years, early, obs_years, obs; fixed_observation_scale=1.0)
    @test score_v2_loss("Eyrie Reef", s3.total_loss, s3.peak_count_penalty,
        s3.sim_peak_years, s3.obs_peak_years) > s3.total_loss
end
