using Test
using DataFrames

@testset "sensitivity._cramers_v" begin

    # ------------------------------------------------------------------
    # 1. Perfect association: category is fully determined by group
    #    membership -- Cramer's V should be exactly 1.0.
    #
    #    Contingency table (levels x groups), levels = [A, B],
    #    groups = [selected=true, selected=false]:
    #
    #        true  false
    #    A |  20     0
    #    B |   0    20
    #
    #    n = 40, row_sums = [20, 20], col_sums = [20, 20]
    #    expected count in every cell = row*col/n = 20*20/40 = 10
    #    chi2 = 4 * (20-10)^2/10 = 4 * 10 = 40
    #    V = sqrt(chi2 / (n * (min(2,2) - 1))) = sqrt(40 / 40) = 1.0
    # ------------------------------------------------------------------
    @testset "perfect association gives V ≈ 1" begin
        col_vals = vcat(fill("A", 20), fill("B", 20))
        selection_mask = vcat(trues(20), falses(20))

        chi2, V = ADRIAanalysis.sensitivity._cramers_v(col_vals, selection_mask)

        @test isapprox(chi2, 40.0; atol=1e-8)
        @test isapprox(V, 1.0; atol=1e-8)
    end

    # ------------------------------------------------------------------
    # 2. No association: category counts split identically between the
    #    two groups regardless of mask -- Cramer's V should be exactly 0.
    #
    #    Contingency table (levels x groups), levels = [A, B, C]:
    #
    #        true  false
    #    A |  10     10
    #    B |  10     10
    #    C |  10     10
    #
    #    n = 60, row_sums = [20, 20, 20], col_sums = [30, 30]
    #    expected count in every cell = row*col/n = 20*30/60 = 10 = observed
    #    chi2 = 0
    #    V = sqrt(0 / (n * (min(3,2) - 1))) = 0.0
    # ------------------------------------------------------------------
    @testset "no association gives V ≈ 0" begin
        col_vals = vcat(
            repeat(["A", "B", "C"], 10),
            repeat(["A", "B", "C"], 10)
        )
        selection_mask = vcat(trues(30), falses(30))

        chi2, V = ADRIAanalysis.sensitivity._cramers_v(col_vals, selection_mask)

        @test isapprox(chi2, 0.0; atol=1e-8)
        @test isapprox(V, 0.0; atol=1e-8)
    end

    # ------------------------------------------------------------------
    # 3. Partial association: an independently-derived reference value
    #    from a known 2x2 table, distinct from the two extremes above.
    #
    #    Contingency table (levels x groups), levels = [A, B]:
    #
    #        true  false
    #    A |  15     5
    #    B |   5    15
    #
    #    n = 40, row_sums = [20, 20], col_sums = [20, 20]
    #    expected count in every cell = 20*20/40 = 10
    #    chi2 = 4 * (15-10)^2/10 = 4 * 2.5 = 10
    #    V = sqrt(10 / (40 * 1)) = sqrt(0.25) = 0.5
    # ------------------------------------------------------------------
    @testset "partial association matches hand-computed reference" begin
        col_vals = vcat(fill("A", 15), fill("B", 5), fill("A", 5), fill("B", 15))
        selection_mask = vcat(trues(20), falses(20))

        chi2, V = ADRIAanalysis.sensitivity._cramers_v(col_vals, selection_mask)

        @test isapprox(chi2, 10.0; atol=1e-8)
        @test isapprox(V, 0.5; atol=1e-8)
    end

    # ------------------------------------------------------------------
    # 4. Degenerate case: a single unique category value returns the
    #    documented sentinel (0.0, 0.0) without erroring, with a @warn.
    # ------------------------------------------------------------------
    @testset "degenerate single-level column returns sentinel + warns" begin
        col_vals = fill("A", 20)
        selection_mask = vcat(trues(10), falses(10))

        chi2, V = @test_logs (:warn, r"fewer than 2 unique") ADRIAanalysis.sensitivity._cramers_v(
            col_vals, selection_mask
        )

        @test chi2 == 0.0
        @test V == 0.0
    end
end
