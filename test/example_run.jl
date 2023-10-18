using Test
using ADRIA

@testset "Full example run" begin
    rs = try
        TEST_RS
    catch
        test_rs()
    end

    @test typeof(rs) <: ADRIA.ResultSet
end
