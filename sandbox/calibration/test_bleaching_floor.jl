"""Unit test for the opt-in bleaching susceptibility floor (ADRIA_BLEACHING_DHW_FLOOR;
crash_mechanism_protocol.md, 2026-10-14). Run: julia --project=sandbox sandbox/calibration/test_bleaching_floor.jl"""

using Test
using ADRIA

function run_bleaching(eff_dhw::Vector{Float64}; kwargs...)
    n_groups, n_sizes, n_locs = 2, 3, length(eff_dhw)
    cover = fill(0.05, n_groups, n_sizes, n_locs)
    stdev = fill(1.5, n_groups, n_sizes)
    dist_t = fill(4.0, n_groups, n_sizes, n_locs)          # mean tolerance 4 DHW
    bleach = zeros(2, n_groups, n_sizes, n_locs)
    tol_ceil = fill(4.0 + ADRIA.HEAT_UB, n_groups, n_sizes, n_locs)
    ADRIA.bleaching_mortality!(cover, eff_dhw, stdev, dist_t, @view(bleach[1:2, :, :, :]),
        tol_ceil; kwargs...)
    return cover, dist_t
end

@testset "Bleaching susceptibility floor" begin
    dhw = [0.5, 2.0, 3.6, 4.0, 6.0, 8.0]
    default_cover, default_dist = run_bleaching(dhw)
    explicit_cover, explicit_dist = run_bleaching(dhw; heat_lb=ADRIA.HEAT_LB)
    @test default_cover == explicit_cover                    # default reproduces established path
    @test default_dist == explicit_dist
    @test all(default_cover[:, :, 2:3] .< 0.05)              # default: mortality at 2.0 and 3.6 DHW
    floor_cover, _ = run_bleaching(dhw; heat_lb=4.0)
    @test all(floor_cover[:, :, 1:4] .== 0.05)               # floor 4: none at or below 4 DHW
    @test all(floor_cover[:, :, 5:6] .< 0.05)                # mortality above the floor
    # Truncating the tolerance distribution from below renormalises its mass upwards, so
    # mortality above the floor is slightly LOWER than the default (0.992 vs 0.996 at 8 DHW
    # for mean 4, sd 1.5); a stronger high-DHW response comes from the bleaching scalar.
    @test all(floor_cover[:, :, 6] .>= default_cover[:, :, 6])
    @test all((0.05 .- floor_cover[:, :, 6]) .>= 0.95 .* (0.05 .- default_cover[:, :, 6]))
    @test_throws ArgumentError run_bleaching(dhw; heat_lb=0.5)
    # Low-DHW scaling (ADRIA_BLEACHING_LOW_DHW_SCALE): default identical; scaled below the
    # threshold; unchanged at or above it.
    unscaled_cover, unscaled_dist = run_bleaching(dhw; low_dhw_scale=1.0, low_dhw_threshold=4.0)
    @test unscaled_cover == default_cover && unscaled_dist == default_dist
    half_cover, _ = run_bleaching(dhw; low_dhw_scale=0.25, low_dhw_threshold=4.0)
    @test all(isapprox.(0.05 .- half_cover[:, :, 2:3], 0.25 .* (0.05 .- default_cover[:, :, 2:3]); rtol=1e-12))
    @test half_cover[:, :, 4:6] == default_cover[:, :, 4:6]
    zero_cover, _ = run_bleaching(dhw; low_dhw_scale=0.0, low_dhw_threshold=4.0)
    @test all(zero_cover[:, :, 1:3] .== 0.05) && zero_cover[:, :, 4:6] == default_cover[:, :, 4:6]
    @test_throws ArgumentError run_bleaching(dhw; low_dhw_scale=1.5)
end
