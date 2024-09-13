using WGLMakie, GeoMakie, GraphMakie
using ADRIA: viz
using ADRIA: YAXArray
using ADRIA: metrics, metrics.cf_difference_loc

if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

Makie.inline!(false)
@testset "spatial" begin
    fig_opts = Dict(:size => (1600, 800))

    @testset "diff_map" begin
        gd_diff::YAXArray{<:Real,2}, ug_diff::YAXArray{<:Real,2} = cf_difference_loc(
            metrics.relative_cover(TEST_RS), TEST_SCENS
        )

        viz.diff_map(TEST_RS, gd_diff[2, :]; fig_opts=fig_opts)
        viz.diff_map(TEST_RS, ug_diff[2, :]; fig_opts=fig_opts)
    end
end
