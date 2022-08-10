using Test
import ADRIA: mcda_normalize, create_decision_matrix, filter_seed_sites


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

    # There should be nothing in `A` that is Inf or NaN
    @test !any(isnan.(A))
    @test !any(isinf.(A))

    # Site with known 0 max cover should be 0
    @test A[end, 6] == 0.0

    # After normalization, all entries should be in [0,1]
    @test !any(A .> 1.0)
    @test !any(A .< 0.0)
end


@testset "MCDA seed site filter" begin
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
    SE = filter_seed_sites(SE, A, wtconseed, wtwaves, wtheat, wtpredecseed, wtlocover)

    # Last site should be filtered out due to no space
    # Third site should be filtered out due to cover >carrying capacity
    @test size(SE, 1) == (size(A, 1) - 2)
    @test maximum(SE[:, 1]) != maximum(A[:, 1])

    # Largest site with plenty of space should have highest score
    @test SE[1, 6] == 1.0
end
