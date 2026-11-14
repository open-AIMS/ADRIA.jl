using Test
using ADRIA

@testset "Low-DHW bleaching reduction" begin
    @testset "reduction applies only below the threshold" begin
        n_groups, n_sizes, n_locs = 5, 7, 2
        eff_dhw = [3.0, 8.0]  # location 1 below a 6 DHW threshold, location 2 above
        function run_bleaching(; kwargs...)
            cover = fill(0.1, n_groups, n_sizes, n_locs)
            dist_t = fill(4.0, n_groups, n_sizes, n_locs)
            bleach = zeros(2, n_groups, n_sizes, n_locs)
            ADRIA.bleaching_mortality!(
                cover, eff_dhw, fill(1.5, n_groups, n_sizes), dist_t,
                @view(bleach[1:2, :, :, :]),
                fill(4.0 + ADRIA.HEAT_UB, n_groups, n_sizes, n_locs);
                kwargs...
            )
            return 0.1 .- cover  # cover lost to bleaching
        end
        base = run_bleaching()
        @test all(base .> 0.0)
        @test run_bleaching(; low_dhw_reduction=0.0, low_dhw_threshold=6.0) == base
        reduced = run_bleaching(; low_dhw_reduction=0.75, low_dhw_threshold=6.0)
        @test reduced[:, :, 1] ≈ 0.25 .* base[:, :, 1]
        @test reduced[:, :, 2] == base[:, :, 2]
        @test_throws ArgumentError run_bleaching(;
            low_dhw_reduction=1.5, low_dhw_threshold=6.0
        )
    end
end
