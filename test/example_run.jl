using Test
using ADRIA

@testset "Full example run" begin
    rs = get_example_resultset()

    @test typeof(rs) <: ADRIA.ResultSet

    # Test RMEDomain loading
    dom = test_reefmod_engine_domain()
    @test typeof(dom) <: ADRIA.RMEDomain
end
