using Test
using ADRIA.Distributions
using ADRIA: mcda_normalize, create_decision_matrix, create_seed_matrix, create_fog_matrix


@testset "Create decision matrix" begin

    # Dummy data to create decision matrix from
    n_sites = 5
    area = [1000.0, 800.0, 600.0, 200.0, 200.0]
    centr_out = [1.0, 0.5, 0.5, 0.5, 0.5]
    centr_in = [0.1, 0.0, 0.3, 0.1, 0.1]
    dam_prob = [0.05, 0.1, 0.1, 0.5, 0.0]
    heat_stress_prob = [0.05, 0.1, 0.1, 0.5, 0.0]
    zones_criteria = [1.0, 1.33, 0.333, 1.0, 0.333]
    site_depth = [5.0, 5.2, 7.1, 6.3, 8.6]

    # Dummy priority predecessors
    priority_sites = zeros(n_sites)
    predec = zeros(n_sites, 3)
    predec[:, 1:2] .= rand(n_sites, 2)
    predprior = predec[in.(predec[:, 1], [priority_sites']), 2]
    predprior = [x for x in predprior if !isnan(x)]
    predec[predprior, 3] .= 1.0
    risk_tol = 0.8

    prop_cover = [0.3, 0.5, 0.9, 0.6, 0.0]
    max_cover = [0.8, 0.75, 0.95, 0.7, 0.0]
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

    @test !any(isnan.(A)) || "NaNs found in decision matrix"
    @test !any(isinf.(A)) || "Infs found in decision matrix"

    @test A[end, 6] == 0.0 || "Site with 0 max cover should be ignored but was not"
end


@testset "MCDA seed matrix creation" begin
    wtconseedout, wtconseedin, wt_waves, wt_heat, wt_predec_seed, wt_zones_seed, wt_lo_cover, wt_depth_seed = [
        1.0, 1.0, 0.7, 1.0, 0.6, 0.6, 0.6, 0.7
    ]

    # Combine decision criteria into decision matrix A
    n_sites = 5
    area = [1000.0, 800.0, 600.0, 200.0, 200.0]
    centr_out = [1.0, 0.5, 0.5, 0.5, 0.5]
    centr_in = [0.1, 0.0, 0.3, 0.1, 0.1]
    dam_prob = [0.05, 0.1, 0.1, 0.5, 0.0]
    heat_stress_prob = [0.05, 0.1, 0.1, 0.5, 0.0]
    zones_criteria = [1.0, 1.33, 0.333, 1.0, 0.333]
    site_depth = [5.0, 5.2, 7.1, 6.3, 8.6]

    # Dummy priority predecessors
    priority_sites = zeros(n_sites)
    predec = zeros(n_sites, 3)
    predec[:, 1:2] .= rand(n_sites, 2)
    predprior = predec[in.(predec[:, 1], [priority_sites']), 2]
    predprior = [x for x in predprior if !isnan(x)]
    predec[predprior, 3] .= 1.0
    risk_tol = 0.8

    sum_cover = [0.3, 0.5, 0.9, 0.6, 0.0]
    max_cover = [0.8, 0.75, 0.95, 0.7, 0.0]

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
        0.8,
    )

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
    @test size(SE, 1) == n_sites - 2 ||
        "Sites where space available < min_area not filtered out"
end

@testset "MCDA fog matrix creation" begin
    wt_conn_fog, wt_waves, wt_heat, wt_predec_fog, wt_zones_fog, wt_hi_cover = [
        1.0, 0.7, 1.0, 0.6, 0.6, 0.6
    ]

    # Combine decision criteria into decision matrix A
    n_sites = 5
    area = [1000.0, 800.0, 600.0, 200.0, 200.0]
    centr_out = [1.0, 0.5, 0.5, 0.5, 0.5]
    centr_in = [0.1, 0.0, 0.3, 0.1, 0.1]
    dam_prob = [0.05, 0.1, 0.1, 0.5, 0.0]
    heat_stress_prob = [0.05, 0.1, 0.1, 0.5, 0.0]
    zones_criteria = [1.0, 1.33, 0.333, 1.0, 0.333]
    site_depth = [5.0, 5.2, 7.1, 6.3, 8.6]

    # Dummy priority predecssors
    priority_sites = zeros(n_sites)
    predec = zeros(n_sites, 3)
    predec[:, 1:2] .= rand(n_sites, 2)
    predprior = predec[in.(predec[:, 1], [priority_sites']), 2]
    predprior = [x for x in predprior if !isnan(x)]
    predec[predprior, 3] .= 1.0

    prop_cover = [0.75, 0.5, 0.3, 0.7, 0.0]
    max_cover = [0.8, 0.75, 0.6, 0.77, 0.0]
    area_max_cover = max_cover .* area
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
        0.8,
    )

    SH, wsh = create_fog_matrix(
        A,
        k_area[filtered],
        area_max_cover[filtered],
        wt_conn_fog,
        wt_waves,
        wt_heat,
        wt_predec_fog,
        wt_zones_fog,
        wt_hi_cover,
    )

    @test maximum(SH[:, 8]) ==
          (maximum(area_max_cover[convert(Vector{Int64}, A[:, 1])] .- A[:, 8])) ||
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
