"""Unit test for the labelled external-supply counterfactual (COTS_EXTERNAL_SUPPLY_WINDOW;
crash_mechanism_protocol.md, 2026-10-14). Run: julia --project=sandbox sandbox/calibration/test_supply_window.jl"""

using Test
using ADRIA

@testset "External-supply window counterfactual" begin
    @test ADRIA.parse_cots_supply_window("") === nothing
    window = ADRIA.parse_cots_supply_window("2017:2024:10")
    @test window == (2017, 2024, 10.0)
    @test ADRIA.cots_supply_window_factor(nothing, 2020) == 1.0          # unset: established path
    @test ADRIA.cots_supply_window_factor(window, 2016) == 1.0
    @test ADRIA.cots_supply_window_factor(window, 2017) == 10.0
    @test ADRIA.cots_supply_window_factor(window, 2024) == 10.0
    @test ADRIA.cots_supply_window_factor(window, 2025) == 1.0
    @test ADRIA.cots_supply_window_factor(ADRIA.parse_cots_supply_window("2019:2019:0"), 2019) == 0.0
    @test_throws ErrorException ADRIA.parse_cots_supply_window("2017:2024")
    @test_throws ErrorException ADRIA.parse_cots_supply_window("2024:2017:3")
    @test_throws ErrorException ADRIA.parse_cots_supply_window("2017:2024:-1")
    @test_throws ErrorException ADRIA.parse_cots_supply_window("2017:2024:Inf")
end
