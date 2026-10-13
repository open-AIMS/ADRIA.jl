using Test
include(joinpath(@__DIR__, "connectivity_retention.jl"))

@testset "Source-conserving local retention" begin
    # Rows are sources, columns sinks. Source 1 exports mostly to 3; 1 and 2 form a group.
    M = [0.001 0.002 0.050 0.010;
         0.000 0.003 0.020 0.040;
         0.004 0.000 0.002 0.030;
         0.000 0.000 0.000 0.000]
    groups = [[1, 2]]
    out = retain_within_groups(M, groups, 0.5)
    @test sum(out; dims=2) ≈ sum(M; dims=2) atol = 1e-15        # source survival preserved
    @test all(>=(0.0), out)
    @test retained_share(out, [1, 2]) ≈ [0.5, 0.5]
    @test out[1, 1] / out[1, 2] ≈ M[1, 1] / M[1, 2]              # within-group pattern kept
    @test out[1, 3] / out[1, 4] ≈ M[1, 3] / M[1, 4]              # export pattern kept
    @test out[3:4, :] == M[3:4, :]                                # other sources untouched
    @test out[4, :] == zeros(4)                                   # zero row stays zero
    # Orientation: retention changes rows (sources), so column 1 inflow from source 3 is unchanged.
    @test out[3, 1] == M[3, 1]
    # Identity at the existing share.
    rho1 = retained_share(M, [1])[1]
    @test retain_within_groups(M, [[1]], rho1)[1, :] ≈ M[1, :]
    # Singleton groups (all-reef control) and a source without a within-group path.
    single = retain_within_groups(M, [[1], [2], [3]], 0.9)
    @test sum(single; dims=2) ≈ sum(M; dims=2) atol = 1e-15
    @test [single[i, i] / sum(single[i, :]) for i in 1:3] ≈ fill(0.9, 3)
    @test_throws ArgumentError retain_within_groups(M, [[1, 2], [2]], 0.5)
    @test_throws ArgumentError retain_within_groups(M, groups, 1.5)
    @test_throws ArgumentError retain_within_groups(-M, groups, 0.5)
end
