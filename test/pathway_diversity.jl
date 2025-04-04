using ADRIA

@testset "switching_probability returns correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH)
    past_option = :heat_stress
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
    option_names = ADRIA.analysis.option_seed_preference().option_name

    @testset "when only necessary inputs are passed" begin
        result = ADRIA.analysis.switching_probability(
            past_option, decision_matrix, dom.loc_data, mcda_method, min_locs
        )

        @test isapprox(sum(result.probability), 1.0, atol=0.01)
        @test result.option_name == option_names
    end

    @testset "when option is passed" begin
        option = first(option_names)
        result = ADRIA.analysis.switching_probability(past_option, decision_matrix,
            dom.loc_data, mcda_method, min_locs, option)
        @test result isa Float64
    end

    @testset "when ports dataframe is passed" begin
        filtered_ports = ADRIA.analysis._ports()[1:3, :]
        result = ADRIA.analysis.switching_probability(past_option, decision_matrix,
            dom.loc_data, mcda_method, min_locs; ports=filtered_ports)
        @test isapprox(sum(result.probability), 1.0, atol=0.01)
        @test result.option_name == option_names
    end

    @testset "option_seed_preference return correctly" begin
        options = ADRIA.analysis.option_seed_preference(; include_weights=true)

        @test all(typeof.(options[:, :preference]) .== ADRIA.decision.SeedPreferences)

        criteria_names = collect(fieldnames(ADRIA.SeedCriteriaWeights))
        assumed_names = [:seed_heat_stress, :seed_wave_stress, :seed_in_connectivity,
            :seed_out_connectivity, :seed_depth, :seed_coral_cover, :seed_cluster_diversity,
            :seed_geographic_separation]

        msg = """
        The number/order of MCDA criteria has changed. Change the definition of the option_seed_preference to follow
        the new criteria."""

        @test criteria_names == assumed_names || msg
    end
end

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

@testset "cost_index return correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    selected_locations = dom.loc_data.reef_siteid
    ports = ADRIA.analysis._ports()

    @testset "when list only dataframe and no optional argument is passed" begin
        cost_index = ADRIA.analysis.cost_index(dom.loc_data, ports)
        @test isapprox(cost_index, 0.76, atol=0.01)
    end

    @testset "when list with locations and no optional argument is passed" begin
        cost_index = ADRIA.analysis.cost_index(selected_locations, dom.loc_data, ports)
        @test isapprox(cost_index, 0.76, atol=0.01)
    end

    @testset "when new weight is passed" begin
        cost_index = ADRIA.analysis.cost_index(
            selected_locations, dom.loc_data, ports; weight=0.4
        )
        @test isapprox(cost_index, 0.79, atol=0.01)
    end

    @testset "when new normalization is passed" begin
        cost_index = ADRIA.analysis.cost_index(selected_locations, dom.loc_data, ports;
            max_distance_port=500000.0, max_dispersion=40000.0)
        @test isapprox(cost_index, 0.25, atol=0.01)
    end
end

@testset "_distance_port return correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    ports = ADRIA.analysis._ports()

    distance_port = ADRIA.analysis._distance_port(dom.loc_data, ports)
    @test isapprox(distance_port, 137573.76, atol=0.01)

    selected_ports = [:cairns, :townsville, :gladstone]
    filtred_ports = ports[ports.name .∈ [selected_ports], :]
    distance_port = ADRIA.analysis._distance_port(dom.loc_data, filtred_ports)
    @test isapprox(distance_port, 227809.14, atol=0.01)
end

@testset "_dispersion return correct values" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)

    dispersion = ADRIA.analysis._dispersion(dom.loc_data)
    @test isapprox(dispersion, 8593.79, atol=0.01)
end
