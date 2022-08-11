using Test
import ADRIA: mcda_normalize, create_decision_matrix, create_seed_matrix, create_shade_matrix

using Infiltrator
@testset "Create decision matrix" begin

    # Dummy data to create decision matrix from
    n_sites = 5
    area = [1000.0, 800.0, 600.0, 200.0, 200.0]
    centr = [1.0, 0.5, 0.5, 0.5, 0.5]
    damprob = [0.05, 0.1, 0.1, 0.5, 0.0]
    heatstressprob = [0.05, 0.1, 0.1, 0.5, 0.0]

    # Dummy priority predecessors
    prioritysites = zeros(n_sites)
    predec = zeros(n_sites, 3)
    predec[:, 1:2] .= rand(n_sites, 2)
    predprior = predec[in.(predec[:, 1], [prioritysites']), 2]
    predprior = [x for x in predprior if !isnan(x)]
    predec[predprior, 3] .= 1.0
    risktol = 0.8

    sumcover = [0.3, 0.5, 0.9, 0.6, 0.0]
    maxcover = [0.8, 0.75, 0.95, 0.7, 0.0]

    A = create_decision_matrix(1:n_sites, centr, sumcover, maxcover, area, damprob, heatstressprob, predec, risktol)

    @test all(0.0 .<= A[:, 2:end] .<= 1.0) || "`A` decision matrix out of bounds"

    @test !any(isnan.(A)) || "NaNs found in decision matrix"
    @test !any(isinf.(A)) || "Infs found in decision matrix"

    @test A[end, 6] == 0.0 || "Site with 0 max cover should be ignored but was not"
end


@testset "MCDA seed matrix creation" begin
    wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover = [1.0, 0.7, 1.0, 0.6, 0.6]

    # Combine decision criteria into decision matrix A
    n_sites = 5
    area = [1000.0, 800.0, 600.0, 200.0, 200.0]
    centr = [1.0, 0.5, 0.5, 0.5, 0.5]
    damprob = [0.05, 0.1, 0.1, 0.5, 0.0]
    heatstressprob = [0.05, 0.1, 0.1, 0.5, 0.0]

    # Dummy priority predecssors
    prioritysites = zeros(n_sites)
    predec = zeros(n_sites, 3)
    predec[:, 1:2] .= rand(n_sites, 2)
    predprior = predec[in.(predec[:, 1], [prioritysites']), 2]
    predprior = [x for x in predprior if !isnan(x)]
    predec[predprior, 3] .= 1.0

    sumcover = [0.3, 0.5, 0.9, 0.6, 0.0]
    maxcover = [0.8, 0.75, 0.6, 0.77, 0.0]

    A = create_decision_matrix(1:n_sites, centr, sumcover, maxcover, area, damprob, heatstressprob, predec, 0.8)

    SE = zeros(size(A, 1), 6)
    SE = create_seed_matrix(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)

    @test size(SE, 1) == (size(A, 1) - 2) || "Site where cover > carrying capacity not filtered out"
    @test maximum(SE[:, 1]) != maximum(A[:, 1]) || "Last site should be filtered out due to no space"
    @test SE[1, 6] == 1.0 || "Largest site with lots of space should have highest score"

    # After normalization, all entries for seeding decision matrix should be ∈ [0,1]
    @test all(0.0 .<= SE[:, 2:end] .<= 1.0) || "Seeding Decision matrix values out of bounds"
end

@testset "MCDA shade matrix creation" begin
    wtconshade, wtwaves, wtheat, wtpredecshade, wthicover = [1.0, 0.7, 1.0, 0.6, 0.6]

    # Combine decision criteria into decision matrix A
    n_sites = 5
    area = [1000.0, 800.0, 600.0, 200.0, 200.0]
    centr = [1.0, 0.5, 0.5, 0.5, 0.5]
    damprob = [0.05, 0.1, 0.1, 0.5, 0.0]
    heatstressprob = [0.05, 0.1, 0.1, 0.5, 0.0]
    
    @infiltrate
    # Dummy priority predecssors
    prioritysites = zeros(n_sites)
    predec = zeros(n_sites, 3)
    predec[:, 1:2] .= rand(n_sites, 2)
    predprior = predec[in.(predec[:, 1], [prioritysites']), 2]
    predprior = [x for x in predprior if !isnan(x)]
    predec[predprior, 3] .= 1.0

    sumcover = [0.75, 0.5, 0.3, 0.7, 0.0]
    maxcover = [0.8, 0.75, 0.6, 0.77, 0.0]

    A = create_decision_matrix(1:n_sites, centr, sumcover, maxcover, area, damprob, heatstressprob, predec, 0.8)

    SH = zeros(size(A, 1), 6)
    SH = create_shade_matrix(SH, A, wtconshade, wtwaves, wtheat, wtpredecshade, wthicover)
    
    @test SH[4, 6] == 1.0 || "Largest site with most coral should have highest score"

    # After normalization, all entries for seeding decision matrix should be ∈ [0,1] 
    @test all(0.0 .<= SH[:, 2:end] .<= 1.0) || "Shading Decision matrix values out of bounds"
end
