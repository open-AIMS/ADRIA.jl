
if !@isdefined(TEST_RS)
    const TEST_DOM, TEST_N_SAMPLES, TEST_SCENS, TEST_RS = test_rs()
end

@testset "Combine result sets" begin
    combined_rs = ADRIA.combine_results(TEST_RS, TEST_RS)
    @test combined_rs isa ADRIA.ResultSet
end
