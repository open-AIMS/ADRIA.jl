using Test
using ADRIA

if !@isdefined(TEST_DOMAIN_PATH)
    const TEST_DOMAIN_PATH = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")
end

@testset "COTS submodel in scenario runs" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
    ms = ADRIA.model_spec(dom)
    @test all(ms[ms.component .== "CotsParams", :is_constant])  # never sampled by default
    scens = ADRIA.sample(dom, 4)
    @test all(scens.cots_on .== 0)

    @testset "no COTS logs when every scenario has COTS off" begin
        rs_off = ADRIA.run_scenarios(dom, scens[1:2, :], "45")
        @test isnothing(rs_off.cots_pop_log)
        @test isnothing(rs_off.cots_condition_log)
        @test rs_off.logs == NamedTuple()
    end

    # Scenarios 1 and 2 run COTS (with juvenile predation on in scenario 2); 3 and 4 do not
    scens.cots_on[1:2] .= 1
    scens.cots_juvenile_predation[2] = 1
    rs = ADRIA.run_scenarios(dom, scens, "45")

    @testset "COTS logs are written for COTS scenarios only" begin
        @test !isnothing(rs.cots_pop_log)
        pop = rs.cots_pop_log
        @test size(pop, :scenarios) == 4
        @test size(pop, :age_class) == 3
        @test any(pop[scenarios = 1].data .> 0.0)
        @test any(pop[scenarios = 2].data .> 0.0)
        @test all(pop[scenarios = 3].data .== 0.0)
        @test all(pop[scenarios = 4].data .== 0.0)
        @test all(isfinite, pop.data)
        @test all(rs.cots_condition_log[scenarios = 3].data .== 0.0)
        # The logs are submodel logs of the result set; the old properties are deprecated
        @test keys(rs.logs) == (:cots,)
        @test keys(rs.logs.cots) == (:condition, :pop)
        @test rs.logs.cots.pop === rs.cots_pop_log
        @test rs.logs.cots.condition === rs.cots_condition_log
    end

    @testset "result sets with COTS logs combine" begin
        combined = ADRIA.combine_results(rs, rs)
        @test size(combined.cots_pop_log, :scenarios) == 8
        @test size(combined.cots_condition_log, :scenarios) == 8
        @test combined.cots_pop_log[scenarios = 5].data ==
            rs.cots_pop_log[scenarios = 1].data
    end

    @testset "result sets with and without COTS logs combine, filling the missing log with NaN" begin
        rs_off = ADRIA.run_scenarios(dom, scens[3:4, :], "45")
        for (a, b, on_range, off_range) in
            ((rs, rs_off, 1:4, 5:6), (rs_off, rs, 3:6, 1:2))
            combined = ADRIA.combine_results(a, b)
            @test size(combined.cots_pop_log, :scenarios) == 6
            @test size(combined.cots_condition_log, :scenarios) == 6
            @test all(isnan, combined.cots_pop_log[scenarios = off_range].data)
            @test all(isnan, combined.cots_condition_log[scenarios = off_range].data)
            @test combined.cots_pop_log[scenarios = on_range].data == rs.cots_pop_log.data
            # Same factors in both sources: the scenario table is the plain stack
            @test combined.inputs == vcat(a.inputs, b.inputs)
            # One model specification per source, in order, and the source of each scenario
            @test length(ADRIA.model_specs(combined)) == 2
            @test isequal(ADRIA.model_specs(combined)[1], a.model_spec)
            @test isequal(ADRIA.model_specs(combined)[2], b.model_spec)
            @test ADRIA.scenario_sources(combined) ==
                [fill(1, size(a.inputs, 1)); fill(2, size(b.inputs, 1))]
            @test isequal(ADRIA.model_spec(combined), a.model_spec)
            # Merged results can be merged again, keeping every source
            again = ADRIA.combine_results(combined, rs)
            @test length(ADRIA.model_specs(again)) == 3
            @test ADRIA.scenario_sources(again) ==
                [ADRIA.scenario_sources(combined); fill(3, 4)]
        end
        # Sets that are not merged have a single source
        @test length(ADRIA.model_specs(rs)) == 1
        @test ADRIA.scenario_sources(rs) == fill(1, 4)
    end

    @testset "COTS runs from the data package's inputs" begin
        @isdefined(make_cots_test_package) ||
            include(joinpath(@__DIR__, "..", "io", "cots_test_package.jl"))
        pkg_dom = ADRIA.load_domain(
            make_cots_test_package(joinpath(mktempdir(), "Test_domain")), "45"
        )
        # COTS off; COTS on (package connectivity and starting state); COTS on with
        # separated recruitment, drawn annual matrices and outside outbreak supply
        pkg_scens = scens[[3, 1, 1, 1, 1, 1], :]
        pkg_scens.cots_separated_recruitment[3:6] .= 1
        pkg_scens.cots_larval_fecundity[3:6] .= 20.0
        pkg_scens.cots_larval_conn_mode[3] = 2
        pkg_scens.cots_larval_conn_seed[3] = 20260930
        pkg_scens.cots_external_outbreak_production[3] = 1000.0
        # Remaining dispersal paths, all with separated recruitment: the mean reef
        # connectivity (4), annual matrices in turn (5), and drawn matrices without
        # outside outbreak supply (6)
        pkg_scens.cots_larval_conn_mode[5] = 1
        pkg_scens.cots_larval_conn_mode[6] = 2
        pkg_scens.cots_larval_conn_seed[6] = 20260930
        pkg_rs = ADRIA.run_scenarios(pkg_dom, pkg_scens, "45")

        pop = pkg_rs.cots_pop_log
        for i in 4:6
            @test any(pop[scenarios = i].data .> 0.0)
        end
        # Each dispersal path gives its own trajectory: non-separated (2), mean (4),
        # annual in turn (5), drawn (6), and drawn with outside outbreak supply (3)
        paths = [2, 4, 5, 6, 3]
        for (k, i) in enumerate(paths), j in paths[(k + 1):end]
            @test pop[scenarios = i].data != pop[scenarios = j].data
        end
        @test all(isfinite, pop.data)
        @test all(pop[scenarios = 1].data .== 0.0)
        @test any(pop[scenarios = 3].data .> 0.0)
        @test pop[scenarios = 2].data != rs.cots_pop_log[scenarios = 1].data
        @test pop[scenarios = 3].data != pop[scenarios = 2].data
        # COTS inputs leave coral-only scenarios unchanged
        for (name, cube) in pairs(pkg_rs.outcomes)
            @test cube[scenarios = 1].data == rs.outcomes[name][scenarios = 3].data
        end
    end

    @testset "unsupported COTS settings fail before any scenario runs" begin
        bad = scens[1:2, :]
        bad.cots_on .= [1, 1]
        # Separated recruitment needs annual reef connectivity, which Test_domain lacks
        bad.cots_separated_recruitment[2] = 1
        err = try
            ADRIA.validate_cots_settings(dom, bad)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("Scenario 2", err.msg)
        @test occursin("annual reef connectivity", err.msg)
        @test_throws ErrorException ADRIA.run_scenarios(dom, bad, "45")

        # Outside outbreak supply needs separated recruitment
        bad.cots_separated_recruitment[2] = 0
        bad.cots_external_outbreak_production[2] = 10.0
        @test_throws "needs separated recruitment" ADRIA.validate_cots_settings(dom, bad)

        # Scenarios with COTS off are not checked
        bad.cots_on[2] = 0
        @test isnothing(ADRIA.validate_cots_settings(dom, bad))
    end

    # Last, as it changes the domain: a COTS factor can be opened up for sampling
    @testset "COTS factors can be opened up for sampling" begin
        ADRIA.set_factor_bounds!(dom; cots_juvenile_predation_capacity=(10.0, 40.0, 20.0))
        ms2 = ADRIA.model_spec(dom)
        @test !only(ms2[ms2.fieldname .== :cots_juvenile_predation_capacity, :is_constant])
    end
end
