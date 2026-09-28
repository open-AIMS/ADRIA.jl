using Test
include(joinpath(@__DIR__, "calibration_parameters.jl"))

@testset "Calibration parameter modes" begin
    focused_names, focused_bounds = calibration_search_space("focused")
    expanded_names, expanded_bounds = calibration_search_space("expanded")
    pulse_names, pulse_bounds = calibration_search_space("expanded_pulse")
    @test length(focused_names) == length(focused_bounds) == 4
    @test length(expanded_names) == length(expanded_bounds) == 15
    @test length(pulse_names) == length(pulse_bounds) == 18
    @test pulse_names[(end - 2):end] == ["pulse_start", "pulse_duration", "pulse_relative_magnitude"]
    @test all(lower < upper for (lower, upper) in pulse_bounds)
    @test_throws ArgumentError calibration_search_space("invalid")
end
