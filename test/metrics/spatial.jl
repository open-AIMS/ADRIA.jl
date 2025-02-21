using ADRIA: metrics
using ADRIA: YAXArray

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "Spatially aggregated bootstrapped ensembled differences" begin
    test_metrics = [
        metrics.relative_cover,
        metrics.juvenile_indicator,
        metrics.absolute_shelter_volume,
        metrics.relative_juveniles,
        metrics.coral_evenness,
        metrics.relative_shelter_volume
    ]

    @testset "Values within bounds for guided/unguided diffs" begin
        for metric in test_metrics
            @testset "Default agg_metric (median)" begin
                guided_diff = metrics.ensemble_loc_difference(metric(TEST_RS), TEST_SCENS)
                unguided_diff = metrics.ensemble_loc_difference(
                    metric(TEST_RS), TEST_SCENS; diff_target=:unguided
                )

                for diff in [guided_diff, unguided_diff]
                    @test all(diff[1, :] .<= diff[2, :] .<= diff[3, :])
                end
            end

            @testset "Quantile as agg_metric" begin
                guided_diff = metrics.ensemble_loc_difference(
                    metric(TEST_RS), TEST_SCENS; agg_metric=0.3
                )
                unguided_diff = metrics.ensemble_loc_difference(
                    metric(TEST_RS), TEST_SCENS; diff_target=:unguided, agg_metric=0.7
                )

                for diff in [guided_diff, unguided_diff]
                    # 1 and 3 are the indexes for the lower and upper bounds of the
                    # confidence interval and 2 is the index for the aggregate metric
                    @test all(diff[1, :] .<= diff[2, :] .<= diff[3, :])
                    @test all(
                        diff[summary=At(:lower_bound)] .<= diff[summary=At(:agg_value)] .<=
                        diff[summary=At(:upper_bound)]
                    )
                end
            end
        end
    end
end
