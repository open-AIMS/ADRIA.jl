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
    end

    @testset "result sets with COTS logs combine" begin
        combined = ADRIA.combine_results(rs, rs)
        @test size(combined.cots_pop_log, :scenarios) == 8
        @test size(combined.cots_condition_log, :scenarios) == 8
        @test combined.cots_pop_log[scenarios = 5].data ==
            rs.cots_pop_log[scenarios = 1].data
    end

    @testset "COTS runs from the data package's inputs" begin
        @isdefined(make_cots_test_package) ||
            include(joinpath(@__DIR__, "..", "io", "cots_test_package.jl"))
        pkg_dom = ADRIA.load_domain(
            make_cots_test_package(joinpath(mktempdir(), "Test_domain")), "45"
        )
        # COTS off; COTS on (package connectivity and starting state); COTS on with
        # separated recruitment, drawn annual matrices and outside outbreak supply
        pkg_scens = scens[[3, 1, 1], :]
        pkg_scens.cots_separated_recruitment[3] = 1
        pkg_scens.cots_larval_fecundity[3] = 20.0
        pkg_scens.cots_larval_conn_mode[3] = 2
        pkg_scens.cots_larval_conn_seed[3] = 20260930
        pkg_scens.cots_external_outbreak_production[3] = 1000.0
        pkg_rs = ADRIA.run_scenarios(pkg_dom, pkg_scens, "45")

        pop = pkg_rs.cots_pop_log
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

    # Last, as it changes the domain: a COTS factor can be opened up for sampling
    @testset "COTS factors can be opened up for sampling" begin
        ADRIA.set_factor_bounds!(dom; cots_juvenile_predation_capacity=(10.0, 40.0, 20.0))
        ms2 = ADRIA.model_spec(dom)
        @test !only(ms2[ms2.fieldname .== :cots_juvenile_predation_capacity, :is_constant])
    end
end
