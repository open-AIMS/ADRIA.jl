using ADRIA: metrics, metrics.cf_difference_loc
using ADRIA: YAXArray

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "cf_difference_loc" begin
    test_metrics = [
        metrics.relative_cover,
        metrics.juvenile_indicator,
        metrics.absolute_shelter_volume,
        metrics.relative_juveniles,
        metrics.coral_evenness,
        metrics.relative_shelter_volume
    ]

    for metric in test_metrics
        gd_diff::YAXArray{<:Real,2}, ug_diff::YAXArray{<:Real,2} = cf_difference_loc(
            metric(TEST_RS), TEST_SCENS
        )

        for diff in [gd_diff, ug_diff]
            @test all(diff[1, :] .<= diff[2, :] .<= diff[3, :])
        end
    end
end
