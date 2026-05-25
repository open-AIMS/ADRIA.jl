using ADRIA

const DOM = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
const RME_DOM = ADRIA.load_domain(ADRIA.RMEDomain, joinpath(TEST_DATA_DIR, "RME_test_domain"), "45")

@testset "switching_probability returns correct values" begin
    past_option = :heat_stress
    min_locs = 5  # select at least 10 locations
    valid_locations = DOM.loc_data.reef_siteid

    supported_methods = ADRIA.decision.mcda_methods()
    mcda_method = supported_methods[rand(1:length(supported_methods))]

    seed_pref_names = collect(fieldnames(ADRIA.SeedCriteriaWeights))
    decision_matrix = ADRIA.decision.decision_matrix(
        valid_locations,
        seed_pref_names,
        rand(length(valid_locations), length(seed_pref_names))
    )
    option_names = ADRIA.analysis.option_seed_preference().option_name

    @testset "when only necessary inputs are passed" begin
        result = ADRIA.analysis.switching_probability(
            past_option, decision_matrix, DOM.loc_data, mcda_method, min_locs
        )

        @test isapprox(sum(result.probability), 1.0, atol=0.01)
        @test result.option_name == option_names
    end

    @testset "when option is passed" begin
        option = first(option_names)
        result = ADRIA.analysis.switching_probability(past_option, decision_matrix,
            DOM.loc_data, mcda_method, min_locs, option)
        @test result isa Float64
    end

    @testset "when ports dataframe is passed" begin
        ports = ADRIA.analysis._ports()
        selected_ports = [:cairns, :townsville, :gladstone]
        filtered_ports = ports[ports.name .∈ [selected_ports], :]
        result = ADRIA.analysis.switching_probability(past_option, decision_matrix,
            DOM.loc_data, mcda_method, min_locs; ports=filtered_ports)
        @test isapprox(sum(result.probability), 1.0, atol=0.01)
        @test result.option_name == option_names
    end

    @testset "when custom weights are passed" begin
        result = ADRIA.analysis.switching_probability(
            past_option, decision_matrix, DOM.loc_data, mcda_method, min_locs;
            weights=(0.5, 0.3, 0.2)
        )
        @test isapprox(sum(result.probability), 1.0, atol=0.01)
        @test result.option_name == option_names
    end

    @testset "option_seed_preference return correctly" begin
        options = ADRIA.analysis.option_seed_preference(; include_weights=true)

        @test all(typeof.(options[:, :preference]) .== ADRIA.decision.SeedPreferences)

        criteria_names = collect(fieldnames(ADRIA.SeedCriteriaWeights))
        assumed_names = [:seed_heat_stress, :seed_in_connectivity,
            :seed_out_connectivity, :seed_depth, :seed_coral_cover, :seed_cluster_diversity,
            :seed_geographic_separation, :seed_coral_diversity]

        msg = """
        The number/order of MCDA criteria has changed. Change the definition of the option_seed_preference to follow
        the new criteria."""

        @test criteria_names == assumed_names || msg
    end
end

@testset "option_similarity returns correct values" begin
    @testset "for standard domain" begin
        seleceted_locations_1 = DOM.loc_data[4:8, :]
        seleceted_locations_2 = DOM.loc_data[6:end, :]

        result = ADRIA.analysis.option_similarity(
            seleceted_locations_1, seleceted_locations_2
        )
        @test isapprox(result, 0.99, atol=0.01)
    end

    @testset "for ReefMod Engine domain" begin
        seleceted_locations_1 = RME_DOM.loc_data[1:6, :]
        seleceted_locations_2 = RME_DOM.loc_data[4:end, :]

        result = ADRIA.analysis.option_similarity(
            seleceted_locations_1, seleceted_locations_2
        )
        @test isapprox(result, 0.98, atol=0.01)
    end
end


@testset "_distance_port return correct values" begin
    ports = ADRIA.analysis._ports()

    distance_port = ADRIA.analysis._distance_port(DOM.loc_data, ports)
    @test isapprox(distance_port, 142084.75, atol=0.01)

    selected_ports = [:cairns, :townsville, :gladstone]
    filtred_ports = ports[ports.name .∈ [selected_ports], :]
    distance_port = ADRIA.analysis._distance_port(DOM.loc_data, filtred_ports)
    @test isapprox(distance_port, 230851.47, atol=0.01)
end

@testset "_dispersion return correct values" begin
    dispersion = ADRIA.analysis._dispersion(DOM.loc_data)
    @test isapprox(dispersion, 8593.79, atol=0.01)
end

@testset "distance_port_score return correct values" begin
    ports = ADRIA.analysis._ports()
    locs1 = DOM.loc_data[1:5, :]
    locs2 = DOM.loc_data[6:end, :]

    score = ADRIA.analysis.distance_port_score(locs1, locs2, ports)
    @test 0.0 <= score <= 1.0

    @testset "equal locations return 0.5" begin
        score_eq = ADRIA.analysis.distance_port_score(locs1, locs1, ports)
        @test isapprox(score_eq, 0.5, atol=1e-10)
    end

    @testset "both empty return 0.5" begin
        empty_locs = DOM.loc_data[[], :]
        score_empty = ADRIA.analysis.distance_port_score(empty_locs, empty_locs, ports)
        @test isapprox(score_empty, 0.5, atol=1e-10)
    end
end

@testset "dispersion_score return correct values" begin
    locs1 = DOM.loc_data[1:5, :]
    locs2 = DOM.loc_data[6:end, :]

    score = ADRIA.analysis.dispersion_score(locs1, locs2)
    @test 0.0 <= score <= 1.0

    @testset "equal locations return 0.5" begin
        score_eq = ADRIA.analysis.dispersion_score(locs1, locs1)
        @test isapprox(score_eq, 0.5, atol=1e-10)
    end
end
