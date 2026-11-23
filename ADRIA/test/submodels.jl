using Test
using ADRIA
using ADRIA: ADRIAinterfaces, DataFrames
using Random

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

# The submodel phase of the run loop (snapshot, step, decrease-only check, record) over
# timesteps 2:tf, as `run_model` calls it. `slots_ref[]` is read inside, so with a
# `Ref{Tuple}` the slots are only known as an abstract `Tuple` (as in `run_model`), with a
# concretely typed `Ref` they are known exactly.
function _hook_phase_loop!(slots_ref, before, view, tf)
    for tstep in 2:tf
        slots = slots_ref[]
        copyto!(before, ADRIAinterfaces.cover(view))
        ADRIA._step_submodels!(
            slots, ADRIAinterfaces.AfterDisturbances(), view, tstep, 2000 + tstep
        )
        ADRIA._check_decrease_only(before, ADRIAinterfaces.cover(view), tstep)
        ADRIA._record_submodels!(slots, tstep)
    end
    return nothing
end

@testset "Submodel hooks in the run loop" begin
    dom = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
    scens = ADRIA.sample(dom, 4)
    row = scens[2, :]
    tf = length(dom.env_layer_md.timeframe)

    @testset "no submodels attached by default" begin
        @test isempty(ADRIA.attached_submodels(dom))
    end

    base = ADRIA.run_scenario(dom, row)
    @test base.aux_logs == (;)

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

    @testset "an increase that the survival ratio cannot see is caught per cell" begin
        # Cells holding no cover: survival stays 1 whatever is added, so only the
        # per-cell snapshot check can catch it
        d1 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d1, FillMock(1e-6))
        err = try
            ADRIA.run_scenario(d1, row)
            nothing
        catch e
            e
        end
        @test err isa AssertionError
        @test occursin("increased (or invalidated) coral cover", err.msg)

        # An increase that only offsets a disturbance loss (survival 0.9 <= 1 would pass)
        before = fill(0.5, 1, 1, 1)
        @test ADRIA._check_decrease_only(before, fill(0.5, 1, 1, 1), 2) === nothing
        @test ADRIA._check_decrease_only(before, fill(0.4, 1, 1, 1), 2) === nothing
        @test_throws AssertionError ADRIA._check_decrease_only(
            before, fill(0.5 + 1e-12, 1, 1, 1), 2
        )
        @test_throws AssertionError ADRIA._check_decrease_only(
            before, fill(NaN, 1, 1, 1), 2
        )
        @test_throws AssertionError ADRIA._check_decrease_only(
            before, fill(-0.1, 1, 1, 1), 2
        )
    end

    @testset "run_scenario equals the matching slice of run_scenarios" begin
        d1 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d1, ShrinkMock(0.1))
        rs = ADRIA.run_scenarios(d1, scens, "45")
        single = ADRIA.run_scenario(d1, scens[3, :])
        @test _rel_cover(single) == rs.outcomes[:relative_cover][scenarios = 3].data

        # and the reduction is visible against the run without the submodel (scenario 2)
        @test sum(rs.outcomes[:relative_cover][scenarios = 2].data) <
            sum(_rel_cover(base))
    end

    @testset "the context carries the domain's connectivity matrix" begin
        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        m = ConnMock()
        ADRIA._attach_submodel!(d2, m)
        a = only(ADRIA.attached_submodels(d2))
        ctx = ADRIA._submodel_context(d2, a)
        @test ctx isa ADRIAinterfaces.SubmodelContext
        @test ctx.connectivity === d2.conn.data
        @test size(ctx.connectivity) == (ctx.n_locs, ctx.n_locs)
        @test eltype(ctx.connectivity) === Float64
        @test isconcretetype(typeof(ctx))

        # Reachable from the hooks (inference/allocation: see the ADRIAinterfaces tests)
        st = ADRIAinterfaces.init_state(m, ctx, (;), Random.Xoshiro(1))
        @test st.conn === d2.conn.data
        @inferred ADRIAinterfaces.init_state(m, ctx, (;), Random.Xoshiro(1))
        cov = zeros(ctx.n_groups, ctx.n_sizes, ctx.n_locs)
        v = ADRIAinterfaces.CoralStepView(cov, zeros(ctx.n_locs), ctx.location_k)
        step_conn(m, st, ctx, v) = ADRIAinterfaces.step!(
            m, ADRIAinterfaces.AfterDisturbances(), st, ctx, v, 2, 2001
        )
        step_conn(m, st, ctx, v)
        @test (@allocated step_conn(m, st, ctx, v)) == 0
        @test st.total == 2 * sum(d2.conn.data)

        # Attached and run: the run is unchanged
        res = ADRIA.run_scenario(d2, row)
        @test isequal(res.raw, base.raw)
    end

    @testset "a guided reactive scenario runs with a submodel attached" begin
        # Reactive strategies use run_model locals (cover losses etc.) next to the hook
        # phase; a name clash there broke attached runs once
        sc = ADRIA.sample_guided(dom, 4)
        for f in ("iv_CAq_strategy", "iv_Fog_strategy", "iv_LvM_strategy")
            sc[!, f] .= 2.0
        end
        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d2, ShrinkMock(0.1))
        res = ADRIA.run_scenario(d2, sc[2, :])
        @test all(isfinite, res.raw)
        @test res.aux_logs.mock_calls[tf] == tf - 1
    end

    @testset "the submodel phase allocates nothing per step after warm-up" begin
        # The hook calls of the time loop on their own, with slots built as `run_model`
        # builds them for a seeded scenario row. Concrete slot types: no allocation at
        # all. Abstract `Tuple` (as in `run_model`): the dynamic dispatch boxes a few
        # arguments, a small bounded amount per step (48 B on Julia 1.12)
        old_rng = copy(Random.default_rng())
        Random.seed!(7)
        srow = ADRIA.sample(dom, 4)[2, :]
        copy!(Random.default_rng(), old_rng)

        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d2, ShrinkMock(0.1))
        ps = ADRIA.DataCube(Vector(srow); factors=names(srow))
        caches = ADRIA._prepare_submodels(d2)
        slots = ADRIA._init_submodels(
            d2, ps, collect(ps.factors.val), caches, Random.Xoshiro(1)
        )
        @test length(slots) == 1 && isconcretetype(typeof(slots))

        cg = d2.coral_growth
        cov = rand(Random.Xoshiro(2), cg.n_groups, cg.n_sizes, cg.n_locs)
        before = similar(cov)
        view = ADRIAinterfaces.CoralStepView(cov, zeros(cg.n_locs), ADRIA.loc_k_area(d2))

        concrete = Ref(slots)
        abstract = Ref{Tuple}(slots)
        _hook_phase_loop!(concrete, before, view, tf)  # warm-up (compile)
        _hook_phase_loop!(abstract, before, view, tf)
        @test (@allocated _hook_phase_loop!(concrete, before, view, tf)) == 0
        @test (@allocated _hook_phase_loop!(abstract, before, view, tf)) <= 256 * (tf - 1)

        # The hooks really ran: every step reduced the cover and was recorded
        @test slots[1].state.nsteps == 4 * (tf - 1)
        @test slots[1].bufs.mock_calls[tf] == 4 * (tf - 1)
    end

    @testset "attaching a submodel adds little allocation to a run" begin
        # Loose whole-run guard (the per-step check is the test above). Difference, not an
        # absolute value: a run allocates a lot regardless. ShrinkMock so that `step!` and
        # `record!` really run; the extra (slots, logs, per-cell snapshot, in-line
        # validate/prepare) must stay a small fraction of the run's allocation
        d2 = ADRIA.load_domain(TEST_DOMAIN_PATH, "45")
        ADRIA._attach_submodel!(d2, ShrinkMock(0.0))  # runs, leaves cover unchanged
        ADRIA.run_scenario(d2, row)  # compile
        detached = @allocated ADRIA.run_scenario(dom, row)
        attached = @allocated ADRIA.run_scenario(d2, row)
        @test attached - detached < 0.02 * detached
    end

    @testset "fixed-seed checksum of a short run (detached)" begin
        # Total of the raw cover array of scenario 2 of a seeded 4-scenario sample on
        # Test_domain RCP 45 (75 timesteps). Generated from the current code (7490170f
        # plus Stage 0): a change here means simulation results changed, which needs a
        # deliberate decision; update the value only for an intended result change.
        # Not exact on purpose: the last digits differ across Julia versions (1.11.9 and
        # 1.12.7 give 163.98244844000132, 1.13.4 gives 163.9824484400013); rtol 1e-10 is
        # far above that and far below any real result change. The regression harness
        # (test/regression) is the exact, same-version gate.
        old_rng = copy(Random.default_rng())
        Random.seed!(42)
        s42 = ADRIA.sample(dom, 4)
        copy!(Random.default_rng(), old_rng)
        res = ADRIA.run_scenario(dom, s42[2, :])
        @test size(res.raw) == (75, 5, 7, 10)
        @test isapprox(sum(res.raw), 163.9824484400013; rtol=1e-10)
    end
end
