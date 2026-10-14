using Test
include(joinpath(@__DIR__, "score_v3.jl"))

@testset "Frozen score v3 diver term" begin
    years = collect(2018:2024)
    adults = [0.1, 0.2, 1.0, 2.0, 4.0, 9.5, 3.0]
    # Exact match -> zero error; years absent from the model are skipped.
    @test diver_log_error(years, adults, [2020, 2023, 2026], [1.0, 9.5, 30.0]) == (0.0, 2)
    # A 2x error with the 0.5 offset: |log((A + 0.5) / (T + 0.5))|.
    e, n = diver_log_error(years, adults, [2022], [1.75])
    @test n == 1 && e ≈ abs(log(4.5 / 2.25))
    # Symmetric in over- and under-prediction on the log scale.
    @test diver_log_error([2020], [9.5], [2020], [4.5])[1] ≈ diver_log_error([2020], [4.5], [2020], [9.5])[1]
    # No overlap -> NaN; invalid targets throw.
    e0, n0 = diver_log_error(years, adults, [2030], [1.0])
    @test isnan(e0) && n0 == 0
    @test_throws ArgumentError diver_log_error(years, adults, [2020], [-1.0])
    @test_throws ArgumentError diver_log_error(years, adults, [2020], [NaN])
    # Combination and constants as preregistered.
    @test score_v3_loss(6.0, 2.0) == 6.0 + 1.5 * 2.0
    @test SCORE_V3_LAST_MANTA_YEAR == 2021 && SCORE_V3_FIRST_DIVER_YEAR == 2020 && SCORE_V3_OFFSET_HA == 0.5
end
