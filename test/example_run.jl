using Test
using ADRIA

@testset "Full example run" begin
    rs = try
        # Should default to full example with figure creation
        # when running full test suite
        TEST_RS
    catch
        # Otherwise use the smaller example run
        test_small_spec_rs()
    end

    @test typeof(rs) <: ADRIA.ResultSet
end
