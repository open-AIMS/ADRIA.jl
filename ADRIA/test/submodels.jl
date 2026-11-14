using Test
using ADRIA
using ADRIA: ADRIAinterfaces, DataFrames

# Mock submodels are shared with the ADRIAinterfaces package's own tests
@isdefined(MockSubmodels) || include(
    joinpath(pkgdir(ADRIA), "..", "ADRIAinterfaces", "test", "MockSubmodels.jl")
)
using .MockSubmodels

if !@isdefined(TEST_DOMAIN_PATH)
    const TEST_DOMAIN_PATH = joinpath(pkgdir(ADRIA), "test", "data", "Test_domain")
end

# Relative cover (timesteps x locations) of a run_model/run_scenario result, as stored
function _rel_cover(res)
    x = ADRIA.relative_cover(res.raw)
    return Float32.(hasproperty(x, :data) ? x.data : x)
end

@testset "Submodel hooks in the run loop" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
    scens = ADRIA.sample(dom, 4)
    row = scens[2, :]
    tf = length(dom.env_layer_md.timeframe)

    @testset "no submodels attached by default" begin
        @test isempty(ADRIA.attached_submodels(dom))
        res = ADRIA.run_scenario(dom, row)
        @test res.aux_logs == (;)
    end

    base = ADRIA.run_scenario(dom, row)

    @testset "attached submodels are copied with the domain's spec store" begin
        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d2, NoOpMock())
        @test length(ADRIA.attached_submodels(d2)) == 1
        @test length(ADRIA.attached_submodels(deepcopy(d2))) == 1
        @test length(copy(d2.spec_store).submodels) == 1
        @test isempty(ADRIA.attached_submodels(dom))
    end

    @testset "a submodel that does nothing leaves the run bit-identical" begin
        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d2, NoOpMock())
        res = ADRIA.run_scenario(d2, row)
        @test isequal(res.raw, base.raw)
        @test res.aux_logs == (;)
    end

    @testset "hooks are called in order, once per timestep" begin
        calls = Tuple{Symbol,Symbol,Int}[]
        a = RecordingMock(:a, calls)
        b = RecordingMock(:b, calls)
        off = RecordingMock(:off, calls; on=false)
        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        for m in (a, off, b)
            ADRIA._attach_submodel!(d2, m)
        end
        res = ADRIA.run_scenario(d2, row)

        # Attach, then (direct call) validate and prepare in-line, then per-scenario init
        # for the active submodels only, then steps in attach order followed by records
        @test calls[1:3] == [(:a, :load_inputs, 0), (:off, :load_inputs, 0),
            (:b, :load_inputs, 0)]
        @test calls[4:9] == [(:a, :validate, 0), (:off, :validate, 0), (:b, :validate, 0),
            (:a, :prepare, 0), (:off, :prepare, 0), (:b, :prepare, 0)]
        @test calls[10:11] == [(:a, :init_state, 0), (:b, :init_state, 0)]
        per_step = calls[12:end]
        @test length(per_step) == 4 * (tf - 1)
        for (k, t) in enumerate(2:tf)
            @test per_step[(4 * (k - 1) + 1):(4 * k)] ==
                [(:a, :step!, t), (:b, :step!, t), (:a, :record!, t), (:b, :record!, t)]
        end
        @test !any(c -> c[1] == :off && c[2] in (:init_state, :step!, :record!), calls)

        # Observers do not change the run, and see the cover of the run itself
        @test isequal(res.raw, base.raw)
        @test length(a.seen) == tf - 1 && a.seen == b.seen
        @test all(>=(0.0), a.seen)
    end

    @testset "a submodel reducing cover: decrease-only contract" begin
        frac = 0.2
        seen_calls = Tuple{Symbol,Symbol,Int}[]
        obs = RecordingMock(:obs, seen_calls)
        obs0 = RecordingMock(:obs0, seen_calls)

        # frac = 0: multiplies by 1.0, identical to no submodel
        d0 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d0, ShrinkMock(0.0))
        r0 = ADRIA.run_scenario(d0, row)
        @test isequal(r0.raw, base.raw)

        d1 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d1, ShrinkMock(frac))
        ADRIA._attach_submodel!(d1, obs)
        r1 = ADRIA.run_scenario(d1, row)
        @test all(isfinite, r1.raw)
        @test all(>=(0.0), r1.raw)
        # Cover never exceeds the unreduced run's at the first step (same starting cover)
        @test sum(r1.raw[2, :, :, :]) < sum(base.raw[2, :, :, :])
        @test sum(r1.raw[2, :, :, :]) ≈ (1 - frac) * sum(base.raw[2, :, :, :]) rtol = 1e-6
        # The reduction reaches the next step: its cover reflects it
        @test sum(r1.raw[3, :, :, :]) < sum(base.raw[3, :, :, :])
        # The observer after the reducer saw the reduced cover (and the reducer ran first)
        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d2, obs0)
        ADRIA.run_scenario(d2, row)
        @test obs.seen[1] ≈ (1 - frac) * obs0.seen[1] rtol = 1e-6

        # Logs of the reducer are returned among the results
        @test keys(r1.aux_logs) == (:mock_cover, :mock_calls)
        @test r1.aux_logs.mock_calls[2:tf] == collect(1.0:(tf - 1))
        @test r1.aux_logs.mock_cover[2:tf] == collect(2.0:tf)
    end

    @testset "a submodel that adds cover trips the survival assertion" begin
        d1 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d1, GrowMock(1.5))
        @test_throws AssertionError ADRIA.run_scenario(d1, row)
    end

    @testset "run_scenario equals the matching slice of run_scenarios" begin
        calls = Tuple{Symbol,Symbol,Int}[]
        d1 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d1, ShrinkMock(0.1))
        rs = ADRIA.run_scenarios(d1, scens, "45")
        for i in (1, 3)
            single = ADRIA.run_scenario(d1, scens[i, :])
            @test _rel_cover(single) == rs.outcomes[:relative_cover][scenarios = i].data
        end

        # and the reduction is visible against a domain without the submodel
        rs0 = ADRIA.run_scenarios(dom, scens, "45")
        @test sum(rs.outcomes[:relative_cover][scenarios = 1].data) <
            sum(rs0.outcomes[:relative_cover][scenarios = 1].data)
    end
end
