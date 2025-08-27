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

    n_samples = 2^5

    @testset "ReefModDomain run" begin
        # Test RMEDomain loading
        reefmod_dom = ADRIA.load_domain(
            ReefModDomain,
            joinpath(TEST_DATA_DIR, "Reefmod_test_domain"),
            "45"
        )
        @test typeof(reefmod_dom) <: ADRIA.ReefModDomain
        reefmod_samples = ADRIA.sample(reefmod_dom, n_samples)
    end
    @testset "RMEDomain run" begin
        rme_dom = ADRIA.load_domain(
            RMEDomain, joinpath(TEST_DATA_DIR, "RME_test_domain"), "45"
        )
        @test typeof(rme_dom) <: ADRIA.RMEDomain
        rme_samples = ADRIA.sample(rme_dom, n_samples)

        rme_dom = ADRIA.load_domain(
            RMEDomain, joinpath(TEST_DATA_DIR, "RME_test_domain"), "45";
            timeframe=(2008, 2022)
        )
    end
end
