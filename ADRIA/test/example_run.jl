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

        half = length(reefmod_dom.loc_ids) ÷ 2
        target_locs = reefmod_dom.loc_ids[1:half]
        ADRIA.set_seed_target_locations!(
            reefmod_dom,
            [(weight=TEST_TARGET_WEIGHT_1, target_locs=target_locs),
                (
                    weight=TEST_TARGET_WEIGHT_2,
                    target_locs=reefmod_dom.loc_ids[(half + 1):end]
                )]
        )
        ADRIA.set_mc_target_locations!(
            reefmod_dom,
            [(weight=TEST_TARGET_WEIGHT_1, target_locs=target_locs),
                (
                    weight=TEST_TARGET_WEIGHT_2,
                    target_locs=reefmod_dom.loc_ids[(half + 1):end]
                )]
        )
        @test reefmod_dom.seed_target_locations[1].weight == TEST_TARGET_WEIGHT_1
        @test reefmod_dom.seed_target_locations[2].weight == TEST_TARGET_WEIGHT_2
        @test reefmod_dom.mc_target_locations[1].weight == TEST_TARGET_WEIGHT_1
        @test reefmod_dom.mc_target_locations[2].weight == TEST_TARGET_WEIGHT_2
    end
    @testset "RMEDomain run" begin
        rme_dom = ADRIA.load_domain(
            RMEDomain, joinpath(TEST_DATA_DIR, "RME_test_domain"), "45"
        )
        @test typeof(rme_dom) <: ADRIA.RMEDomain
        rme_samples = ADRIA.sample(rme_dom, n_samples)

        half = length(rme_dom.loc_ids) ÷ 2
        target_locs = rme_dom.loc_ids[1:half]
        ADRIA.set_seed_target_locations!(
            rme_dom,
            [(weight=TEST_TARGET_WEIGHT_1, target_locs=target_locs),
                (weight=TEST_TARGET_WEIGHT_2, target_locs=rme_dom.loc_ids[(half + 1):end])]
        )
        ADRIA.set_mc_target_locations!(
            rme_dom,
            [(weight=TEST_TARGET_WEIGHT_1, target_locs=target_locs),
                (weight=TEST_TARGET_WEIGHT_2, target_locs=rme_dom.loc_ids[(half + 1):end])]
        )
        @test rme_dom.seed_target_locations[1].weight == TEST_TARGET_WEIGHT_1
        @test rme_dom.seed_target_locations[2].weight == TEST_TARGET_WEIGHT_2
        @test rme_dom.mc_target_locations[1].weight == TEST_TARGET_WEIGHT_1
        @test rme_dom.mc_target_locations[2].weight == TEST_TARGET_WEIGHT_2

        rme_dom = ADRIA.load_domain(
            RMEDomain, joinpath(TEST_DATA_DIR, "RME_test_domain"), "45";
            timeframe=(2008, 2022)
        )
    end
end
