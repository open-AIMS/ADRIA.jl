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

@testset "Absolute cost" begin
    @testset "_distance_port" begin
        # List of main ports in the GBR
        # ref: https://elibrary.gbrmpa.gov.au/jspui/retrieve/5391c720-b846-4fae-93c7-0ff53f829ca2/Ports%20and%20Shipping%20Information%20sheet-29May2013.pdf
        ports = Dict(
           :quintell_beach => (143.5444588,-12.8437284),
           :cooktown => (145.2475149,-15.4610395),
           :cairns => (145.7808046,-16.9213696),
           :townsville => (146.8332585,-19.2529659),
           :abbot_point => (148.0948982,-19.9220318),
           :hay_point => (149.2737503,-21.2934747),
           :gladstone => (151.2993586,-23.8740713),
        )
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)

        distance_port = ADRIA.analysis._distance_port(dom.loc_data, ports)
        @test isapprox(distance_port, 137573.76, atol=0.01)
    end

    @testset "_dispersion" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)

        dispersion = ADRIA.analysis._dispersion(dom.loc_data)
        @test isapprox(dispersion, 5906.83, atol=0.01)
    end

    @testset "absolute cost" begin
        dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
        selected_locations = dom.loc_data.reef_siteid

        # List of main ports in the GBR
        # ref: https://elibrary.gbrmpa.gov.au/jspui/retrieve/5391c720-b846-4fae-93c7-0ff53f829ca2/Ports%20and%20Shipping%20Information%20sheet-29May2013.pdf
        ports = Dict(
           :quintell_beach => (143.5444588,-12.8437284),
           :cooktown => (145.2475149,-15.4610395),
           :cairns => (145.7808046,-16.9213696),
           :townsville => (146.8332585,-19.2529659),
           :abbot_point => (148.0948982,-19.9220318),
           :hay_point => (149.2737503,-21.2934747),
           :gladstone => (151.2993586,-23.8740713),
        )

        absolute_cost = ADRIA.analysis.absolute_cost(selected_locations, dom.loc_data, ports)

        @test isapprox(absolute_cost, 0.65, atol=0.01)
    end
end
