if get(ENV, "ADRIA_TEST_RULES_EXT", "0") != "1"
    @info "Skipping integration tests (set ADRIA_TEST_RULES_EXT=1)"
else
    using Test, ADRIA, ADRIAanalysis, SIRUS, MLJ, DataFrames

    domain_path = get(ENV, "ADRIA_TEST_DOMAIN_PATH",
        joinpath(@__DIR__, "..", "..", "test", "data", "Test_domain"))

    if !isdir(domain_path)
        @warn "No domain path (ADRIA_TEST_DOMAIN_PATH) -- integration tests skipped"
    else
        @testset "Integration: cluster_rules end-to-end" begin
            dom = ADRIA.load_domain(domain_path, "45")
            scens = ADRIA.sample(dom, 64)
            rs = ADRIA.run_scenarios(dom, scens, "45")
            s_tac = ADRIA.metrics.scenario_total_cover(rs)
            clusters_4 = ADRIA.analysis.cluster_scenarios(s_tac, 4)
            target = ADRIAanalysis.target_clusters(clusters_4, s_tac)
            factors = ADRIA.component_params(rs, [ADRIA.Intervention]).fieldname

            rv = ADRIAanalysis.cluster_rules(rs, target, scens, factors, 6)
            @test rv isa Vector
            @test all(r -> r isa ADRIAanalysis.Rule, rv)
            @test all(r -> parentmodule(typeof(r)) === ADRIAanalysis, rv)
            @test_nowarn ADRIAanalysis.print_rules(rv)
            @test ADRIAanalysis.maximum_probability(rv) > 0.0

            # BitVector overload
            bv = BitVector(target .== 1)
            rv_bv = ADRIAanalysis.cluster_rules(rs, bv, scens, factors, 6)
            @test rv_bv isa Vector{<:ADRIAanalysis.Rule}

            # Constant-only factors must raise ArgumentError
            ms = ADRIA.model_spec(rs)
            const_factors = ms[ms.is_constant .== true, :fieldname]
            if !isempty(const_factors)
                @test_throws ArgumentError ADRIAanalysis.cluster_rules(
                    rs, target, scens, const_factors, 3)
            end
        end
    end
end
