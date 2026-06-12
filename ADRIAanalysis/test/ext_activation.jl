if get(ENV, "ADRIA_TEST_RULES_EXT", "0") != "1"
    @info "Skipping extension activation tests (set ADRIA_TEST_RULES_EXT=1)"
else
    using Test, ADRIAanalysis, SIRUS, MLJ
    import DataFrames: DataFrame
    import StableRNGs: StableRNG

    @testset "ADRIAanalysisRulesExt activation" begin
        ext = Base.get_extension(ADRIAanalysis, :ADRIAanalysisRulesExt)
        @test !isnothing(ext)

        @test hasmethod(ADRIAanalysis.cluster_rules,
            (ADRIA.ResultSet, Vector{Int64}, DataFrame, Vector{Symbol}, Int64))
        @test hasmethod(ADRIAanalysis.rules, (SIRUS.StableRules{Int64},))

        @testset "rules() returns ADRIAanalysis.Rule objects" begin
            n = 60
            rng = StableRNG(42)
            scens = DataFrame(; x1=rand(rng, n), x2=rand(rng, n))
            clusters = vcat(ones(Int64, 30), zeros(Int64, 30))
            model = StableRulesClassifier(; max_rules=2, rng=StableRNG(1))
            mach = machine(model, scens, clusters)
            MLJ.fit!(mach)
            rv = ADRIAanalysis.rules(mach.fitresult)

            @test rv isa Vector
            @test all(r -> r isa ADRIAanalysis.Rule, rv)
            @test all(r -> parentmodule(typeof(r)) === ADRIAanalysis, rv)
        end

        @testset "_remove_duplicates deduplication" begin
            r1 = ADRIAanalysis.Rule(Vector{Vector}([["x", :L, 0.3f0]]), [0.9, 0.1])
            r2 = ADRIAanalysis.Rule(Vector{Vector}([["x", :L, 0.9f0]]), [0.7, 0.3])
            deduped = ext._remove_duplicates([r1, r2])
            @test length(deduped) == 1
            @test deduped[1].consequent[1] >= 0.9
        end

        @testset "_condition internal fields smoke test" begin
            n = 60
            rng = StableRNG(7)
            scens = DataFrame(; x1=rand(rng, n), x2=rand(rng, n))
            clusters = vcat(ones(Int64, 30), zeros(Int64, 30))
            model = StableRulesClassifier(; max_rules=3, rng=StableRNG(7))
            mach = machine(model, scens, clusters)
            MLJ.fit!(mach)
            @test_nowarn ADRIAanalysis.rules(mach.fitresult)
        end
    end
end
