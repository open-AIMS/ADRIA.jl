using Test
using ADRIA
using ADRIA.CSV, ADRIA.DataFrames

@isdefined(make_cots_test_package) || include("cots_test_package.jl")

@testset "Optional COTS inputs in the data package" begin
    @testset "a package without COTS inputs has none" begin
        inputs = ADRIA.load_cots_inputs(TEST_DOMAIN_PATH, ["reef_01"], collect(2025:2099))
        @test all(isnothing(getfield(inputs, f)) for f in fieldnames(ADRIA.CotsInputs))
    end

    pkg = make_cots_test_package(joinpath(mktempdir(), "Test_domain"))
    ids = names(CSV.read(joinpath(pkg, "cots_connectivity", "cots_conn.csv"), DataFrame))[2:end]
    timeframe = collect(2025:2099)
    inputs = ADRIA.load_cots_inputs(pkg, ids, timeframe)

    @testset "every input loads" begin
        @test size(inputs.conn) == (10, 10)
        @test inputs.conn[1, 1] == 0.3
        @test inputs.conn[1, 2] == 0.05

        f = inputs.reef_forcing
        @test f.years == [2001, 2002]
        @test f.site_to_reef == [fill(1, 5); fill(2, 5)]
        @test f.reef_habitable_area_ha == [50.0, 80.0]
        @test f.reef_connectivity[2] == [0.1 0.2; 0.2 0.1]
        @test f.outbreak_years == [2026, 2027]
        @test size(f.outbreak_supply) == (2, 2, 2)
        @test f.outbreak_supply[2, 1, 2] ≈ 4e-4  # reef B, 2026, hydrodynamic year 2002

        @test inputs.initial_state[3, :] == [3.0, 0.5, 2.0]
        @test size(inputs.larval_survival) == (75, 10)
        @test inputs.larval_survival[2, 1] == 0.5  # 2026 is the second timestep
        @test inputs.larval_survival[6, 2] == 0.0  # 2030
        @test count(!=(1.0), inputs.larval_survival) == 2  # 1990 is outside the timeframe
        @test inputs.settlement_modifier == [fill(0.8, 9); 0.5]
    end

    @testset "invalid inputs are rejected" begin
        dir = mktempdir()
        bad(fn, df) = (CSV.write(joinpath(dir, fn), df); joinpath(dir, fn))
        @test_throws ErrorException ADRIA.load_cots_initial_state(
            bad(
                "init.csv",
                DataFrame(site_id=ids[1:9], juvenile_ha=1.0, subadult_ha=1.0,
                    adult_ha=1.0)
            ), ids
        )
        @test_throws ErrorException ADRIA.load_cots_initial_state(
            bad(
                "init.csv",
                DataFrame(site_id=ids, juvenile_ha=-1.0, subadult_ha=1.0,
                    adult_ha=1.0)
            ), ids
        )
        @test_throws ErrorException ADRIA.load_cots_settlement_modifier(
            bad("settle.csv", DataFrame(site_id=[ids; "elsewhere"], multiplier=0.5)), ids
        )
        @test_throws ErrorException ADRIA.load_cots_settlement_modifier(
            bad("settle.csv", DataFrame(site_id=ids, multiplier=1.5)), ids
        )
        @test_throws ErrorException ADRIA.load_cots_larval_survival(
            bad("surv.csv", DataFrame(site_id=["elsewhere"], year=2026, multiplier=0.5)),
            ids, timeframe
        )
    end

    @testset "annual matrices are used in turn or drawn reproducibly" begin
        f = inputs.reef_forcing
        @test [ADRIA.larval_conn_index(f, t, 1, 0) for t in 1:3] == [1, 2, 1]
        drawn = [ADRIA.larval_conn_index(f, t, 2, 20260930) for t in 1:20]
        @test all(in((1, 2)), drawn) && length(unique(drawn)) == 2
        @test drawn == [ADRIA.larval_conn_index(f, t, 2, 20260930) for t in 1:20]
        @test_throws ArgumentError ADRIA.larval_conn_index(f, 1, 0, 0)
    end

    @testset "outside outbreak supply" begin
        f = inputs.reef_forcing
        supply = ADRIA.outbreak_supply_by_location(f, 2, 2026, 1000.0)
        @test supply ≈ [fill(0.2, 5); fill(0.4, 5)]
        @test ADRIA.outbreak_supply_by_location(f, 0, 2026, 1000.0) ≈
            [fill(0.15, 5); fill(0.3, 5)]  # mean over the annual matrices
        @test all(iszero, ADRIA.outbreak_supply_by_location(f, 1, 2030, 1000.0))
    end

    @testset "inputs act on the COTS state" begin
        params = ADRIA.COTSMod.COTSParams(; separated_recruitment=true)
        models = ADRIA.initialize_cots(10, params; rng=ADRIA.Random.Xoshiro(1))
        ADRIA.set_initial_state!(models, inputs.initial_state)
        @test models[3].N == [3.0, 0.5, 2.0]

        for m in models
            m.last_local_fecundity = 2.0
            m.last_settlement_gate = 1.0
        end
        ADRIA.apply_larval_survival!(models, inputs.larval_survival[2, :])
        @test models[1].last_local_fecundity == 1.0
        @test models[2].last_local_fecundity == 2.0
        ADRIA.apply_settlement_modifier!(models, inputs.settlement_modifier)
        @test models[10].last_settlement_gate == 0.5
        @test models[3].N[1] == 3.0 * 0.8
    end

    @testset "the domain carries the inputs" begin
        dom = ADRIA.load_domain(pkg, "45")
        @test !isnothing(dom.cots_inputs.reef_forcing)
        @test ADRIA.cots_connectivity(dom) === dom.cots_inputs.conn
        @test dom.cots_inputs.initial_state == inputs.initial_state
    end
end
