using Test
using DataFrames
using Random

@testset "sensitivity.rsa -- method dispatch" begin

    # Categorical column is stored as small integer codes (Float64-valued), the
    # same representation ADRIA scenario tables use for nominal factors (e.g.
    # `mcda_method`) -- this matters below: forcing `method=:mann_whitney`
    # must still be able to compute (Float64.(col_vals) on a genuine string
    # column would throw MethodError, which is not what a "uniform override"
    # test should be exercising).

    # ------------------------------------------------------------------
    # 1. :auto dispatch routes per-column based on colmetadata "ptype"
    # ------------------------------------------------------------------
    @testset "auto dispatch: tagged categorical -> cramers_v, untagged -> mann_whitney" begin
        Random.seed!(1)
        n = 60
        X = DataFrame(;
            cont_col=rand(n),
            cat_col=Float64.(rand(1:3, n))
        )
        selection_mask = vcat(trues(30), falses(30))

        colmetadata!(X, :cat_col, "ptype", "unordered categorical"; style=:note)

        result = ADRIAanalysis.sensitivity.rsa(X, selection_mask)

        cat_row = only(filter(:feature => f -> f == :cat_col, result))
        cont_row = only(filter(:feature => f -> f == :cont_col, result))

        @test cat_row.test == :cramers_v
        @test cont_row.test == :mann_whitney
    end

    @testset "auto dispatch: explicitly-tagged continuous column still uses mann_whitney" begin
        Random.seed!(2)
        n = 60
        X = DataFrame(;
            cont_col=rand(n),
            cat_col=Float64.(rand(1:3, n))
        )
        selection_mask = vcat(trues(30), falses(30))

        colmetadata!(X, :cont_col, "ptype", "continuous"; style=:note)
        colmetadata!(X, :cat_col, "ptype", "unordered categorical"; style=:note)

        result = ADRIAanalysis.sensitivity.rsa(X, selection_mask)

        cat_row = only(filter(:feature => f -> f == :cat_col, result))
        cont_row = only(filter(:feature => f -> f == :cont_col, result))

        @test cat_row.test == :cramers_v
        @test cont_row.test == :mann_whitney
    end

    # ------------------------------------------------------------------
    # 2. method= override forces uniform application, ignoring ptype tags
    # ------------------------------------------------------------------
    @testset "method=:mann_whitney forces every row to mann_whitney" begin
        Random.seed!(3)
        n = 60
        X = DataFrame(;
            cont_col=rand(n),
            cat_col=Float64.(rand(1:3, n))
        )
        selection_mask = vcat(trues(30), falses(30))
        colmetadata!(X, :cat_col, "ptype", "unordered categorical"; style=:note)

        result = ADRIAanalysis.sensitivity.rsa(X, selection_mask; method=:mann_whitney)
        @test all(result.test .== :mann_whitney)
    end

    @testset "method=:cramers_v forces every row to cramers_v" begin
        Random.seed!(4)
        n = 60
        X = DataFrame(;
            cont_col=rand(n),
            cat_col=Float64.(rand(1:3, n))
        )
        selection_mask = vcat(trues(30), falses(30))
        colmetadata!(X, :cat_col, "ptype", "unordered categorical"; style=:note)

        result = ADRIAanalysis.sensitivity.rsa(X, selection_mask; method=:cramers_v)
        @test all(result.test .== :cramers_v)
    end

    # ------------------------------------------------------------------
    # 3. Invalid method throws ArgumentError
    # ------------------------------------------------------------------
    @testset "unknown method throws ArgumentError" begin
        X = DataFrame(; a=rand(10))
        selection_mask = vcat(trues(5), falses(5))
        @test_throws ArgumentError ADRIAanalysis.sensitivity.rsa(
            X, selection_mask; method=:bogus
        )
    end

    # ------------------------------------------------------------------
    # 4. Backwards-compatible fallback for plain, untagged DataFrames
    #    containing a column whose *eltype* looks non-numeric (the
    #    `looks_nominal` heuristic checks `eltype(col) <: Real`, which a
    #    genuine `Vector{String}` column would also trip -- but a truly
    #    nominal/string-valued column is NOT actually Mann-Whitney
    #    computable (`Float64.(col_vals)` throws `MethodError` on
    #    strings), so it can't be used to test the "does not throw"
    #    fallback path. `Union{Missing,Int}` reproduces the same
    #    eltype-heuristic trigger (`eltype <: Real` is false for a Union
    #    type) while still being numeric and Float64-convertible when no
    #    values are actually missing, so it is the realistic case where
    #    the documented "informational only, falls back to mann_whitney"
    #    behaviour actually holds without erroring.
    # ------------------------------------------------------------------
    @testset "untagged DataFrame falls back to mann_whitney everywhere, with @warn" begin
        Random.seed!(5)
        n = 40
        nom_col = Union{Int,Missing}[rand(1:3) for _ = 1:n]
        X = DataFrame(;
            cont_col=rand(n),
            nom_col=nom_col
        )
        selection_mask = vcat(trues(20), falses(20))

        # No colmetadata attached at all (plain DataFrame)
        @test isempty(DataFrames.colmetadatakeys(X))
        # Sanity check: the eltype heuristic is actually triggered
        @test !(eltype(X.nom_col) <: Real)

        # match_mode=:any: nom_col's random draws may or may not also trip the
        # ">20% tied values" warning depending on RNG stream, so only assert
        # the ptype warning is present rather than requiring an exact log set.
        result = @test_logs (:warn, r"No column \"ptype\" metadata found") match_mode = :any ADRIAanalysis.sensitivity.rsa(
            X, selection_mask
        )

        @test result isa DataFrame
        @test all(result.test .== :mann_whitney)
    end
end
