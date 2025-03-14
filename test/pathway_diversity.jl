using ADRIA

@testset "option_similarity returns correct values" begin
    @testset "for standard domain" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
        seleceted_locations_1 = dom.loc_data.reef_siteid[4:8]
        seleceted_locations_2 = dom.loc_data.reef_siteid[6:end]

        result = ADRIA.analysis.option_similarity(
            seleceted_locations_1, seleceted_locations_2, dom.loc_data
        )
        @test isapprox(result, 0.85, atol=0.01)
    end

    @testset "for ReefMod Engine domain" begin
        dom = ADRIA.load_domain(ADRIA.RMEDomain, TEST_REEFMOD_ENGINE_DOMAIN_PATH, "45")
        seleceted_locations_1 = dom.loc_data.UNIQUE_ID[5:15]
        seleceted_locations_2 = dom.loc_data.UNIQUE_ID[10:end]

        result = ADRIA.analysis.option_similarity(
            seleceted_locations_1, seleceted_locations_2, dom.loc_data
        )
        @test isapprox(result, 0.62, atol=0.01)
    end
end
