
if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "Log element types are Float64 on load" begin
    @test eltype(TEST_RS.ranks) == Float64
    @test eltype(TEST_RS.mc_log) == Float64
    @test eltype(TEST_RS.seed_log) == Float64
    @test eltype(TEST_RS.shading_log) == Float64
    @test eltype(TEST_RS.coral_dhw_tol_log) == Float64

    # coral_cover_log: enabled via log_cover=true in test/config.toml
    @test !isnothing(TEST_RS.coral_cover_log)
    @test eltype(TEST_RS.coral_cover_log) == Float64

    # Verify UInt16 unscaling: cover values must be in [0, 1]
    cover_data = Array(TEST_RS.coral_cover_log)
    @test all(0.0 .<= cover_data .<= 1.0)
end

@testset "Combine result sets" begin
    combined_rs = ADRIA.combine_results(TEST_RS, TEST_RS)
    @test combined_rs isa ADRIA.ResultSet
end
