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

    # # Test RMEDomain loading
    reefmod_dom = ADRIA.load_domain(
        ReefModDomain,
        joinpath(TEST_DATA_DIR, "Reefmod_test_domain"),
        "45"
    )
    @test typeof(reefmod_dom) <: ADRIA.ReefModDomain

    rme_dom = ADRIA.load_domain(RMEDomain, joinpath(TEST_DATA_DIR, "RME_test_domain"), "45")
    @test typeof(rme_dom) <: ADRIA.RMEDomain
end
