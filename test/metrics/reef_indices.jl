
using ADRIA: metrics

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "reef_indices.jl" begin
    @testset "reduced_reef_tourism_index" begin
        test_metric(metrics.reef_tourism_index, (TEST_RS,))
    end

    @testset "scenario_reduced_rci" begin
        test_metric(metrics.scenario_reduced_rci, (TEST_RS,))
    end

    @testset "reduced_reef_tourism_index" begin
        test_metric(metrics.reef_tourism_index, (TEST_RS,))
    end

    @testset "scenario_reduced_rti" begin
        test_metric(metrics.scenario_reduced_rti, (TEST_RS,))
    end

    @testset "reef_fish_index" begin
        test_metric(metrics.reef_fish_index, (TEST_RS,))
    end

    @testset "scenario_rfi" begin
        test_metric(metrics.scenario_rfi, (TEST_RS,))
    end
end
