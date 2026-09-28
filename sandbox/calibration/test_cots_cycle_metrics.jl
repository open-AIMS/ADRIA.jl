using Test
using Statistics
using CSV
using DataFrames
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

wrong_second_amplitude = copy(good)
wrong_second_amplitude[years .== 2010] .= 0.18
wrong_second_amplitude[years .== 2009] .= 0.08
wrong_second_amplitude[years .== 2011] .= 0.08
scaled_good = 12.0 .* good

@testset "COTS cycle metrics" begin
    obs_peaks = detect_cots_peaks(years, obs)
    @test [p.year for p in obs_peaks] == [1995, 2010]

    good_score = cots_cycle_score(years, good, years, obs)
    shifted_score = cots_cycle_score(years, shifted_late, years, obs)
    flat_score = cots_cycle_score(years, flat, years, obs)
    one_peak_score = cots_cycle_score(years, one_peak, years, obs)
    early_score = cots_cycle_score(years, shifted_early, years, obs)
    wrong_amplitude_score = cots_cycle_score(years, wrong_second_amplitude, years, obs)
    scaled_good_score = cots_cycle_score(years, scaled_good, years, obs)
    fixed_scale_score = cots_cycle_score(
        years, good, years, obs; fixed_observation_scale=0.5
    )
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
    @test good_score.amplitude_penalty < wrong_amplitude_score.amplitude_penalty
    @test good_score.total_loss < wrong_amplitude_score.total_loss
    @test isapprox(scaled_good_score.total_loss, good_score.total_loss; atol=1e-10)
    @test isapprox(scaled_good_score.observation_scale, good_score.observation_scale / 12.0)
    @test fixed_scale_score.observation_scale == 0.5
end

@testset "Observed peak regression" begin
    observations = CSV.read(joinpath(@__DIR__, "..", "data", "reef_cots.csv"), DataFrame)
    expected = Dict(
        "Lizard Isles" => [1998, 2013],
        "Macgillivray Reef" => [1996, 2015],
        "North Direction Island" => [1995, 2013],
        "Eyrie Reef" => [2013],
    )
    for (reef_name, expected_years) in expected
        reef = observations[observations.reef_name .== reef_name, :]
        yearly = combine(groupby(reef, :year), :cotsptow => mean => :value)
        peaks = detect_cots_peaks(Int.(yearly.year), Float64.(yearly.value))
        @test [peak.year for peak in peaks] == expected_years
    end
end

@testset "Calendar-aware surveys and scoring window" begin
    irregular_years = [1985, 1988, 1994, 1995, 1997, 2003, 2009, 2010, 2012, 2018]
    irregular_values = [0.02, 0.03, 0.45, 1.0, 0.2, 0.02, 0.35, 0.85, 0.15, 0.02]
    irregular_peaks = detect_cots_peaks(irregular_years, irregular_values)
    @test [peak.year for peak in irregular_peaks] == [1995, 2010]

    boundary_values = copy(irregular_values)
    boundary_values[1] = 2.0
    boundary_peaks = detect_cots_peaks(irregular_years, boundary_values)
    @test !(1985 in [peak.year for peak in boundary_peaks])

    extended_years = collect(1985:2026)
    extended_obs = vcat(obs, [2.0, 0.0])
    common_score = cots_cycle_score(years, good, extended_years, extended_obs)
    @test common_score.scoring_end_year == 2024
    @test !(2025 in common_score.obs_peak_years)

    duplicate_years, duplicate_values = yearly_series([1995, 1995, 1996], [0.8, 1.2, 0.4])
    @test duplicate_years == [1995, 1996]
    @test duplicate_values == [1.0, 0.4]
end
