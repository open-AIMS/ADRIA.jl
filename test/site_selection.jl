using Distributions
using Combinatorics, IterTools
using Test
using ADRIA.Distributions

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

function test_site_ranks(
	weights_set::Vector{Vector{Float64}},
	A::Matrix{Float64},
	rankings::Matrix{Int64},
	n_site_int::Int64,
	mcda_func::Function,
	inv::Float64,
)
    S = ADRIA.decision.mcda_normalize(A)
    S[:, 1] .= A[:, 1]
	criteria_names = [
		"heat stress",
		"wave stress",
		"median depth",
		"coral cover space",
		"in connectivity",
		"out connectivity",
	]
    for weights in weights_set
        crit_inds = findall(weights .> 0.0)
        criteria = vec(
            sum(
                weights[crit_inds]' .*
                ADRIA.decision.mcda_normalize(A[:, 2:end][:, crit_inds]);
                dims=2,
            ),
        )
        site_max = get_site_order(criteria, site_ids)
        prefsites, s_order = ADRIA.decision.rank_sites!(
            S, weights, rankings, n_site_int, mcda_func, inv
        )

        # check that 5 worst sites aren't in those selected
        names_temp = criteria_names[crit_inds]
        names_string = string(["$(name), " for name in names_temp[1:(end - 1)]]...)

        @test !any(in.(Int.(site_max[(n_site_int + 1):end, 1]), prefsites)) || string(
            "For the ",
            names_string,
            "and $(names_temp[end]) criteria, some of the worst sites were still selected during site selection.",
        )
        @test any(prefsites .== 5) & any(prefsites .== 6) || string(
            "For the ",
            names_string,
            "and $(names_temp[end]) criteria, the best sites (5 and 6) were not selected.",
        )
    end
end

function test_mcda_funcs(rankings, S, weights, mcda_funcs, n_site_int)
    for mf in mcda_funcs
        prefsites, s_order = ADRIA.decision.rank_sites!(
            S, weights, rankings, n_site_int, mf, 2
        )
        @test in(5, prefsites) .& in(6, prefsites) ||
            "The best overall sites were not chosen by method $mf"
    end
end

function get_test_decision_matrix(dom)
    cover = sum(dom.init_coral_cover; dims=:species)[species=1]
    leftover_space = 1 .- cover
    k_area = dom.site_data.area .* dom.site_data.k
    dhw_av = ADRIA.decision.summary_stat_env(dom.dhw_scens, (:timesteps, :scenarios))
    wave_av = ADRIA.decision.summary_stat_env(dom.wave_scens, (:timesteps, :scenarios))
    depth_med = dom.site_data.depth_med

    TP_data = ADRIA.connectivity_strength(
        dom.TP_data .* ADRIA.site_k_area(dom), collect(cover), dom.TP_data
    )

    site_ids = dom.site_data.site_id

    heat_stress =
        1 .- vec((dhw_av .- minimum(dhw_av)) ./ (maximum(dhw_av) - minimum(dhw_av)))
    wave_stress =
        1 .- vec((wave_av .- minimum(wave_av)) ./ (maximum(wave_av) - minimum(wave_av)))
    space_area = leftover_space .* k_area
    in_conn = TP_data.in_conn
    out_conn = TP_data.out_conn

    A = hcat(site_ids, heat_stress, wave_stress, depth_med, space_area, in_conn, out_conn)

	return A
    return A, criteria_names
end

@testset "site selection" begin
    # TODO: Complete tests with @tests

    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    p_tbl = ADRIA.param_table(dom)

    p_tbl[:, :depth_offset] .= 7.0

    # ranks = ADRIA.site_selection(dom, p_tbl, 1, 10, 1)
end
@testset "MCDA variable constructor" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
    criteria_df = ADRIA.param_table(dom)  # Get single scenario dataframe

    site_ids = collect(1:length(dom.site_ids))
    available_space = rand(Uniform(200, 30000), length(site_ids))

    area_to_seed = 962.11  # Area of seeded corals in m^2.

    dhw_scens = dom.dhw_scens[1, :, criteria_df.dhw_scenario[1]]
    wave_scens = dom.wave_scens[1, :, criteria_df.wave_scenario[1]]

    mcda_vars = ADRIA.decision.DMCDA_vars(
        dom,
        criteria_df[1, :],
        site_ids,
        available_space,
        area_to_seed,
        wave_scens,
        dhw_scens,
    )
    n_sites = length(mcda_vars.site_ids)
    @test (size(mcda_vars.conn, 1) == n_sites) && (size(mcda_vars.conn, 2) == n_sites) || "Connectivity input is incorrect size."
    @test length(mcda_vars.dam_prob) == n_sites || "Wave damage input is incorrect size."
    @test length(mcda_vars.heat_stress_prob) == n_sites || "Heat stress input is incorrect size."
    @test length(mcda_vars.leftover_space) == n_sites ||
        "Initial cover input is incorrect size."
end
@testset "Unguided site selection" begin
    n_intervention_sites = 5
    pref_seed_sites = zeros(Int64, n_intervention_sites)
    pref_fog_sites = zeros(Int64, n_intervention_sites)
    seed_years = true
    fog_years = true
    max_cover = [0.0, 3000.0, 5000.0, 0.0, 0.0]
    depth_priority = collect(1:5)

    pref_seed_sites, pref_fog_sites = ADRIA.decision.unguided_site_selection(
        pref_seed_sites,
        pref_fog_sites,
        seed_years,
        fog_years,
        5,
        max_cover,
        depth_priority,
    )

    # Check that only two sites are selected (the sites where k > 0.0)
    @test length(pref_seed_sites[pref_seed_sites .> 0]) == 2
    @test length(pref_fog_sites[pref_fog_sites .> 0]) == 2

    @test all([in(sid, [2, 3]) for sid in pref_seed_sites[pref_seed_sites .> 0]])
    @test all([in(sid, [2, 3]) for sid in pref_fog_sites[pref_fog_sites .> 0]])
end

@testset "Guided site selection without ADRIA ecological model" begin
    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
    site_ids = collect(1:length(dom.site_data.site_id))
	A = get_test_decision_matrix(dom)
    n_sites = length(site_ids)
    n_site_int = 5

    rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    # Test 0.0 or 1.0 weights in all combinations
    weights = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    weights_choice = collect(0.0:0.1:1.0)

    for num_crit in 1:6
        weights[num_crit] = 1.0
        weights_set = collect(distinct(permutations(weights, 6)))
        test_site_ranks(
            weights_set, A, rankings, criteria_names, n_site_int, mcda_funcs[1], 1
        )
    end
    weights = rand(100, 6)
    weights = weights ./ sum(weights; dims=2)
    for ww in eachrow(weights)
        weights_set = collect(distinct(permutations(ww, 6)))
        test_site_ranks(
            weights_set, A, rankings, criteria_names, n_site_int, mcda_funcs[1], 1
        )
    end

end

@testset "Test ranks line up with ordering" begin
    supported_methods = ADRIA.decision.mcda_methods()
    mcda_func = supported_methods[rand(1:length(supported_methods))]
    n_sites = 20

    S = ADRIA.decision.mcda_normalize(rand(Uniform(0, 1), n_sites, 6))

    weights = ADRIA.decision.mcda_normalize(rand(Uniform(0, 1), 6))
    n_site_int = 5
    site_ids = collect(1:n_sites)
    S = hcat(site_ids, S)

    rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    prefsites, s_order = ADRIA.decision.rank_sites!(
        S, weights, rankings, n_site_int, mcda_func, 2
    )

    @test all([(rankings[rankings[:, 1].==s_order[rank, 1], 2].==rank)[1] for rank in 1:size(s_order, 1)]) || "Ranking does not match mcda score ordering"

end

@testset "Test each mcda method is working" begin
    mcda_funcs = ADRIA.decision.mcda_methods()

    dom = ADRIA.load_domain(EXAMPLE_DOMAIN_PATH, 45)
	A = get_test_decision_matrix(dom)

    site_ids = dom.site_data.site_id
    n_sites = length(site_ids)
    n_site_int = 5

    rankings = Int64[site_ids zeros(Int64, n_sites) zeros(Int64, n_sites)]

    S = ADRIA.decision.mcda_normalize(A)
    S[:, 1] .= A[:, 1]

    weights = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0]

    test_mcda_funcs(rankings, S, weights, mcda_funcs, n_site_int)
end