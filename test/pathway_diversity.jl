using ADRIA

@testset "distance_selected_locations returns float" begin
    @testset "for standard domain" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH)

        seleceted_locations_1 = ADRIA.sample(dom.loc_data.reef_siteid,3)
        seleceted_locations_2 = ADRIA.sample(dom.loc_data.reef_siteid,3)

        result = ADRIA.analysis.distance_selected_locations(seleceted_locations_1, seleceted_locations_2, dom.loc_data)
        @test result isa Float64
    end

    @testset "for ReefMod Engine domain" begin
        dom  = ADRIA.load_domain(ADRIA.RMEDomain, TEST_REEFMOD_ENGINE_DOMAIN_PATH, "45")
        seleceted_locations_1 = ADRIA.sample(dom.loc_data.UNIQUE_ID,3)
        seleceted_locations_2 = ADRIA.sample(dom.loc_data.UNIQUE_ID,3)

        result = ADRIA.analysis.distance_selected_locations(seleceted_locations_1, seleceted_locations_2, dom.loc_data)
        @test result isa Float64
    end
end
