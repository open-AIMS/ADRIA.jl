using Test
using ADRIA
using ADRIA.Distributions

@testset "proportional adjustment" begin
    Y = rand(5, 36, 20)
    for i in axes(Y, 1)
        Y[i, :, :] .= Y[i, :, :] / sum(Y[i, :, :]; dims=1)
    end

    tmp = zeros(20)
    for i in axes(Y, 1)
        # No warning should be emitted when values are between 0 and 1
        Test.@test_nowarn ADRIA.proportional_adjustment!(Y[i, :, :], tmp)

        @test all(0.0 .<= Y[i, :, :] .<= 1.0)
    end
end

@testset "Coral Spec" begin
    coral_params = ADRIA.coral_spec().params

    bin_edge_diameters_cm2 = ADRIA.colony_mean_area(ADRIA.bin_edges())
    stored_colony_mean_areas = ADRIA.colony_mean_area(
        coral_params.mean_colony_diameter_m .* 100.0
    )

    # check colony areas in cm^2 are within bounds designated by bin edges
    for k in 1:6
        @test all(
            stored_colony_mean_areas[coral_params.class_id .== k] .>=
            bin_edge_diameters_cm2[k]
        ) ||
            "Some colony areas for size class $k are larger than the size class upper bound."
        @test all(
            stored_colony_mean_areas[coral_params.class_id .== k] .>=
            bin_edge_diameters_cm2[k]
        ) ||
            "Some colony areas for size class $k are smaller than the size class lower bound."
    end
end

@testset "Fecundity" begin
    n_locs = 216
    n_groups = 5
    n_sizes = 7
    fec_groups = zeros(n_groups, n_locs)

    fec_params = [
        0.0  0.0  5.17731e5  6.91902e5  1.05234e6  1.24466e6  1.46211e6;
        0.0  0.0  2.08702e5  2.14901e5  2.5563e5   2.58547e5  2.70719e5;
        0.0  0.0  1.79594e6  5.47688e6  1.2319e7   2.34875e7  3.79331e7;
        0.0  0.0  2.00452e5  2.44523e5  3.15532e5  4.16459e5  5.09551e5;
        0.0  0.0  2.12764e5  2.40148e5  3.38464e5  4.37654e5  5.08423e5
    ]

    C_cover_t = rand(Uniform(0.0, 0.1), n_groups, n_sizes, n_locs)

    habitable_areas = rand(Uniform(5, 1e7), 1, n_locs)

    ADRIA.fecundity_scope!(
        fec_groups, fec_params, C_cover_t, habitable_areas
    )

    @test any(fec_groups .> 1e8) ||
        "Fecundity is measured in m² and so should be a very large number"
    @test !any(fec_groups .< 0.0) || "Negative fecundity is not allowed"
end

@testset "Recruitment" begin
    n_locs = 334
    n_groups = 5

    total_loc_area = rand(n_locs) .* 1e6
    max_cover = rand(n_locs)
    avail_area = rand(n_locs)
    larval_pool = rand(n_groups, n_locs) .* 1e8

    recruits_per_m² = ADRIA.recruitment_rate(larval_pool, avail_area)
    abs_recruits = recruits_per_m² .* (avail_area .* max_cover .* total_loc_area)'

    @test any(abs_recruits .> 10^4) || "At least some recruitment values should be > 10,000"

    theoretical_max = ((avail_area .* max_cover .* total_loc_area)' * 51.8)
    for (i, rec) in enumerate(eachrow(abs_recruits))
        @test all(rec' .<= theoretical_max) ||
            "Species group $i exceeded maximum theoretical number of settlers"
    end
end
