using WGLMakie, GeoMakie, GraphMakie
using ADRIA: viz
using ADRIA: YAXArray
using ADRIA: metrics

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

Makie.inline!(false)
@testset "Visualization of spatial.jl metrics" begin
    fig_opts = Dict(:size => (1600, 800))

    @testset "Bootstrapped ensemble difference to the counterfactual" begin
        guided_diff = metrics.ensemble_loc_difference(
            metrics.relative_cover(TEST_RS), TEST_SCENS
        )
        unguided_diff = metrics.ensemble_loc_difference(
            metrics.relative_cover(TEST_RS), TEST_SCENS; diff_target=:unguided
        )

        viz.map(
            TEST_RS, guided_diff[summary=At(:agg_value)]; diverging=true, fig_opts=fig_opts
        )
        viz.map(
            TEST_RS, unguided_diff[summary=At(:agg_value)]; diverging=true,
            fig_opts=fig_opts
        )
    end
end
