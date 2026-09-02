using Test
using Statistics
include("cots_cycle_metrics.jl")

years = collect(1985:2024)
obs = fill(0.02, length(years))
obs[years .== 1995] .= 1.0
obs[years .== 2010] .= 0.85
obs[years .== 1994] .= 0.45
obs[years .== 1996] .= 0.45
obs[years .== 2009] .= 0.35
obs[years .== 2011] .= 0.35

good = fill(0.02, length(years))
good[years .== 1995] .= 0.95
good[years .== 2010] .= 0.8
good[years .== 1994] .= 0.4
good[years .== 1996] .= 0.4
good[years .== 2009] .= 0.3
good[years .== 2011] .= 0.3

shifted_late = fill(0.02, length(years))
shifted_late[years .== 2002] .= 0.95
shifted_late[years .== 2017] .= 0.8

shifted_early = fill(0.02, length(years))
shifted_early[years .== 1992] .= 0.95
shifted_early[years .== 2007] .= 0.8
shifted_early[years .== 1991] .= 0.4
shifted_early[years .== 1993] .= 0.4
shifted_early[years .== 2006] .= 0.3
shifted_early[years .== 2008] .= 0.3

flat = fill(0.2, length(years))
one_peak = fill(0.02, length(years))
one_peak[years .== 1995] .= 1.0

@testset "COTS cycle metrics" begin
    obs_peaks = detect_cots_peaks(years, obs)
    @test [p.year for p in obs_peaks] == [1995, 2010]

    good_score = cots_cycle_score(years, good, years, obs)
    shifted_score = cots_cycle_score(years, shifted_late, years, obs)
    flat_score = cots_cycle_score(years, flat, years, obs)
    one_peak_score = cots_cycle_score(years, one_peak, years, obs)
    early_score = cots_cycle_score(years, shifted_early, years, obs)
    early_lag = lagged_correlation(years, shifted_early, years, obs; max_lag=6)

    @test good_score.n_obs_peaks == 2
    @test good_score.n_sim_peaks == 2
    @test good_score.peak_timing_penalty == 0.0
    @test good_score.period_penalty == 0.0
    @test good_score.best_lag_years == 0
    @test good_score.lag_correlation_penalty < 0.25
    @test early_lag.lag == 3
    @test early_lag.pearson > 0.8
    @test early_score.best_lag_years == 3
    @test flat_score.flatline_penalty == 1.0
    @test one_peak_score.peak_count_penalty > 0.0
    @test good_score.total_loss < shifted_score.total_loss
    @test good_score.total_loss < flat_score.total_loss
    @test one_peak_score.total_loss < flat_score.total_loss
end