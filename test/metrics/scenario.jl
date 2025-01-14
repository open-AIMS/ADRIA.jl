using ADRIA: metrics

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "scenario.jl" begin
    _test_covers::Vector{YAXArray{Float64,4}} = mock_covers()
    _k_area::Vector{Float64} = k_area()

    # Metrics
    @testset "scenario_total_cover" begin
        tac = metrics.total_absolute_cover(TEST_RS)
        test_metric(metrics.scenario_total_cover, (tac,))
        test_metric(metrics.scenario_total_cover, (TEST_RS,))

        for cover in _test_covers
            tac = metrics.total_absolute_cover(metrics.relative_cover(cover), _k_area)
            test_metric(
                metrics.scenario_total_cover, (tac,)
            )
        end
    end

    @testset "scenario_relative_cover" begin
        test_metric(metrics.scenario_relative_cover, (TEST_RS,))
    end

    @testset "scenario_relative_juveniles" begin
        test_metric(metrics.scenario_relative_juveniles, (TEST_RS,))
        _coral_spec = ADRIA.to_coral_spec(ADRA.default_coral_params(), TEST_SCENS[1, :])

        for cover in _test_covers
            test_metric(
                metrics.scenario_relative_juveniles,
                (cover[:, :, :, 1], _coral_spec, _k_area)
            )
        end
    end

    @testset "scenario_absolute_juveniles" begin
        test_metric(metrics.scenario_absolute_juveniles, (TEST_RS,))
        _coral_spec = ADRIA.to_coral_spec(TEST_SCENS[1, :])

        for cover in _test_covers
            test_metric(
                metrics.scenario_absolute_juveniles,
                (cover[:, :, :, 1], _coral_spec, _k_area)
            )
        end
    end

    @testset "scenario_juvenile_indicator" begin
        test_metric(metrics.scenario_juvenile_indicator, (TEST_RS,))
        _coral_spec = ADRIA.to_coral_spec(TEST_SCENS[1, :])

        for cover in _test_covers
            test_metric(
                metrics.scenario_juvenile_indicator,
                (cover[:, :, :, 1], _coral_spec, _k_area)
            )
        end
    end

    @testset "scenario_asv" begin
        test_metric(metrics.scenario_asv, (TEST_RS,))

        for cover in _test_covers
            asv = metrics.absolute_shelter_volume(cover, _k_area, TEST_SCENS)
            test_metric(
                metrics.scenario_asv,
                (asv,)
            )
        end
    end

    @testset "scenario_rsv" begin
        test_metric(metrics.scenario_rsv, (TEST_RS,))

        # TODO
        # for cover in _test_covers
        #     rsv = metrics.relative_shelter_volume(
        #         cover[:, :, :, 1], _k_area, TEST_SCENS[1, :]
        #     )
        #     test_metric(metrics.scenario_rsv, (rsv,))
        # end
    end

    @testset "scenario_evenness" begin
        test_metric(metrics.scenario_evenness, (TEST_RS,))
        n_groups::Int64 = length(ADRIA.default_coral_spec().taxa_names)

        # TODO
        # for cover in _test_covers
        #     rltc = metrics.relative_loc_taxa_cover(cover[:, :, :, 1], _k_area, n_groups)
        #     ce = metrics.coral_evenness(rltc)
        #     test_metric(metrics.scenario_evenness, (ce,))
        # end
    end

    # TODO Helper functions
    # @testset "scenario_outcomes" begin end
    # @testset "scenario_trajectory" begin end
end
