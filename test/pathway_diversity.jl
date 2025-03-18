using ADRIA

@testset "switching_probability returns correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
    past_locations = dom.loc_data.reef_siteid[1:3]
    min_locs = 5  # select at least 10 locations
    valid_locations = dom.loc_data.reef_siteid

    supported_methods = ADRIA.decision.mcda_methods()
    mcda_method = supported_methods[rand(1:length(supported_methods))]

    seed_pref_names = collect(fieldnames(ADRIA.SeedCriteriaWeights))
    decision_matrix = ADRIA.decision_matrix(
        valid_locations,
        seed_pref_names,
        rand(length(valid_locations), length(seed_pref_names))
    )

    result = ADRIA.analysis.switching_probability(
        past_locations, decision_matrix, dom.loc_data, mcda_method, min_locs
    )
    @test sum(values(result)) == 1.0
    @test keys(result) == keys(ADRIA.analysis.OPTION_WEIGHS)

    option = first(keys(ADRIA.analysis.OPTION_WEIGHS))
    new_result = ADRIA.analysis.switching_probability(past_locations, decision_matrix,
        dom.loc_data, mcda_method,
        min_locs, option)
    @test new_result isa Float64
end

@testset "distance_selected_locations returns correct values" begin
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

@testset "cost_index return correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    selected_locations = dom.loc_data.reef_siteid

    @testset "when list only dataframe and no optional argument is passed" begin
        cost_index = ADRIA.analysis.cost_index(dom.loc_data)
        @test isapprox(cost_index, 0.76, atol=0.01)
    end

    @testset "when list with locations and no optional argument is passed" begin
        cost_index = ADRIA.analysis.cost_index(selected_locations, dom.loc_data)
        @test isapprox(cost_index, 0.76, atol=0.01)
    end

    @testset "when new weigh is passed" begin
        cost_index = ADRIA.analysis.cost_index(selected_locations, dom.loc_data; weigh=0.4)
        @test isapprox(cost_index, 0.79, atol=0.01)
    end

    @testset "when new normalization is passed" begin
        cost_index = ADRIA.analysis.cost_index(selected_locations, dom.loc_data;
            max_distance_port=500000.0, max_dispersion=40000.0)
        @test isapprox(cost_index, 0.25, atol=0.01)
    end
end

@testset "_distance_port return correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)

    distance_port = ADRIA.analysis._distance_port(dom.loc_data)
    @test isapprox(distance_port, 137573.76, atol=0.01)

    ports = ADRIA.analysis.PORTS
    filtred_ports_key = [:cairns, :townsville, :gladstone]
    filtred_ports_dict = Dict(k => ports[k] for k in filtred_ports_key)
    distance_port = ADRIA.analysis._distance_port(dom.loc_data, filtred_ports_dict)
    @test isapprox(distance_port, 227809.14, atol=0.01)
end

@testset "_dispersion return correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)

    dispersion = ADRIA.analysis._dispersion(dom.loc_data)
    @test isapprox(dispersion, 8593.79, atol=0.01)
end
