using Test
using Distributions
using ADRIA.decision:
    mcda_normalize,
    create_decision_matrix,
    create_seed_matrix,
    create_fog_matrix,
    priority_location_criteria,
    priority_zones_criteria

function create_test_decision_matrix(n_sites::Int64, risk_tol::Float64)
    # Combine decision criteria into decision matrix A
    site_ids = collect(1:n_sites)
    area = rand(Uniform(200.0, 1000.0), n_sites)
    centr_out = rand(Uniform(0.0, 1.0), n_sites)
    centr_in = rand(Uniform(0.0, 1.0), n_sites)
    dam_prob = rand(Uniform(0.0, 1.0), n_sites)
    heat_stress_prob = rand(Uniform(0.0, 1.0), n_sites)

    zones = [fill("Good", n_sites - 2)..., fill("Bad", 2)...]
    site_depth = rand(Uniform(5.0, 10.0), n_sites)
    strong_pred = rand(1:n_sites, n_sites)

    # Dummy priority predecessors
    predec = priority_location_criteria(strong_pred, [1, 2], site_ids)
    zones_criteria = priority_zones_criteria(strong_pred, zones, ["Good"], site_ids)
    Main.@infiltrate
    prop_cover = [0.3, 0.5, 0.9, 0.6, 0.0]
    max_cover = [0.8, 0.75, 0.95, 0.7, 0.0]

    k_area = max_cover .* area
    leftover_space = (max_cover .- prop_cover) .* area
    A, filtered = create_decision_matrix(
        collect(1:n_sites),
        centr_in,
        centr_out,
        leftover_space,
        dam_prob,
        heat_stress_prob,
        site_depth,
        predec,
        zones_criteria,
        risk_tol,
    )
    return A, filtered, k_area
end
@testset "Create decision matrix" begin
    n_sites = 5
    risk_tol = 0.8

    A, filtered, k_area = create_test_decision_matrix(n_sites, risk_tol)

    @test !any(isnan.(A)) || "NaNs found in decision matrix"
    @test !any(isinf.(A)) || "Infs found in decision matrix"

    @test A[end, 6] == 0.0 || "Site with 0 max cover should be ignored but was not"
end

@testset "MCDA seed matrix creation" begin
    wtconseedout, wtconseedin, wt_waves, wt_heat, wt_predec_seed, wt_zones_seed, wt_lo_cover, wt_depth_seed = rand(
        Uniform(0.0, 1.0),
        8,
    )

    A, filtered, k_area = create_test_decision_matrix(5, 1.0)
    min_area = 20.0
    SE, wse = create_seed_matrix(
        A,
        min_area,
        wtconseedin,
        wtconseedout,
        wt_waves,
        wt_heat,
        wt_predec_seed,
        wt_zones_seed,
        wt_lo_cover,
        wt_depth_seed,
    )

    @test (sum(filtered)) == size(A, 1) || "Site where heat stress > risk_tol not filtered out"
    @test size(SE, 1) == sum(A[:, 8] .> min_area) ||
        "Sites where space available < min_area not filtered out"
end

@testset "MCDA fog matrix creation" begin
    wt_conn_fog, wt_waves, wt_heat, wt_predec_fog, wt_zones_fog, wt_hi_cover = rand(
        Uniform(0.0, 1.0), 6
    )

    n_sites = 5
    A, filtered, k_area = create_test_decision_matrix(n_sites, 0.8)
    SH, wsh = create_fog_matrix(
        A,
        k_area[filtered],
        wt_conn_fog,
        wt_waves,
        wt_heat,
        wt_predec_fog,
        wt_zones_fog,
        wt_hi_cover,
    )

    @test maximum(SH[:, 8]) ==
          (maximum(k_area[convert(Vector{Int64}, A[:, 1])] .- A[:, 8])) ||
        "Largest site with most coral area should have highest score"
end

@testset "MCDA normalisation" begin
    # randomised weights
    w = vec(rand(Uniform(0, 1), 7))

    # randomised decision matrix
    A = zeros(5, 7)
    A[:, 1] = [1.0, 2.0, 3.0, 4.0, 5.0]
    A[:, 2:6] = rand(Uniform(0, 1), (5, 5))
    A[:, 7] = rand(Uniform(100, 1000), (5, 1))

    norm_A = mcda_normalize(A[:, 2:end])
    norm_w = mcda_normalize(w)

    @test all((sqrt.(sum(norm_A .^ 2, dims=1)) .- 1.0) .< 0.0001) || "Decision matrix normalization not giving column sums = 1."
    @test (sum(norm_w) - 1.0) <= 0.001 || "MCDA weights not summing to one."
end

@testset "Priority predecessor and zones criteria" begin
    mcda_funcs = ADRIA.decision.mcda_methods()
    n_site_int = 5

    # Create simplistic set of 8 sites 
    # Site 1 is the strongest predecessor for 2, 3, 4
    # Site 5 is the strongest predecessor for sites 4, 5, 6
    site_ids = collect(1:8)
    n_sites = length(site_ids)
    strong_pred = [0; 1; 1; 1; 0; 5; 5; 5]

    # If 2,3,4 are the priority sites, site 1 should be selected first
    priority_sites = [2, 3, 4]

    priority_sites_crit = ADRIA.decision.priority_location_criteria(
        strong_pred, priority_sites, site_ids
    )
    @test findall(priority_sites_crit .> 0)... == 1 ||
        "The priority zones criteria is valuing sites which are not the strongest predecessor for priority sites"
    zones = [
        "DarkBlue"
        "DarkBlue"
        "DarkBlue"
        "Pink"
        "DarkBlue"
        "DarkBlue"
        "DarkBlue"
        "Pink"
    ]

    priority_zones_crit = ADRIA.decision.priority_zones_criteria(
        strong_pred, zones, ["Pink"], site_ids
    )
    @test all(findall(priority_zones_crit .> 0) .== [1, 4, 5, 8]) ||
        "Some sites which are priority zones or strongest predecessor to priority zones have zero valued zone criteria."

    zones = [
        "Black"
        "Black"
        "Black"
        "Pink"
        "DarkBlue"
        "Black"
        "Black"
        "Pink"
    ]

    priority_zones_crit = ADRIA.decision.priority_zones_criteria(
        strong_pred, zones, ["Pink", "DarkBlue"], site_ids
    )
    @test findall(
        priority_zones_crit .== minimum(priority_zones_crit[priority_zones_crit .!= 0.0])
    )[1] == 1 ||
        "Strongest predecessor site to priority zone in non-priority zone should have smallest non-zero value."

    @test all(findall(priority_zones_crit .== maximum(priority_zones_crit)) == [4, 8]) ||
        "Highest priority zone sites should have highest value in the zone criteria."
end
