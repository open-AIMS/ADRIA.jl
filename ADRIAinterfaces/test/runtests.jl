using Test
using Random
using ADRIAinterfaces
import ADRIAinterfaces: step!
include("MockSubmodels.jl")
using .MockSubmodels

function make_ctx(;
    n_locs=4, tf=5, inputs=nothing, cache=nothing,
    conn=[Float64(i == j) for i in 1:n_locs, j in 1:n_locs]
)
    return SubmodelContext(
        n_locs, 2, 3, tf, collect(2000:(2000 + tf - 1)), [:a, :b],
        fill(10.0, n_locs), conn, inputs, cache
    )
end

@testset "ADRIAinterfaces" begin
    @testset "Aqua" begin
        using Aqua
        Aqua.test_all(ADRIAinterfaces)
    end

    @testset "no-op defaults" begin
        sm = NoOpMock()
        ctx = make_ctx()
        @test submodel_name(sm) === :NoOpMock
        @test phases(sm) == ()
        @test !has_phase(sm, AfterDisturbances())
        @test parameters(sm) == ParamSpec[]
        @test load_inputs(sm, "somewhere") === nothing
        @test validate(sm, nothing, nothing) === nothing
        @test prepare(sm, ctx) === nothing
        @test active(sm, (;))
        @test seed_excluded(sm) == Symbol[]
        @test log_specs(sm, ctx) == LogSpec[]
        @test init_state(sm, ctx, (;), Xoshiro(1)) === nothing

        cov = fill(1.0, 2, 3, 4)
        v = CoralStepView(cov, zeros(4), fill(10.0, 4))
        @test step!(sm, AfterDisturbances(), nothing, ctx, v, 2, 2001) === nothing
        @test cov == fill(1.0, 2, 3, 4)   # default step! leaves cover untouched
        @test record!(sm, (;), nothing, 2) === nothing
    end

    @testset "types" begin
        p = ParamSpec(:m_rate, 0.5, (0, 1))
        @test p.bounds === (0.0, 1.0) && p.dist === :uniform && p.dist_params == [0.0, 1.0]
        @test_throws ArgumentError ParamSpec(:m_rate, 2.0, (0.0, 1.0))
        @test_throws ArgumentError ParamSpec(:m_rate, 0.5, (1.0, 0.0))

        # Additive fields: defaults keep the old behaviour, kind and ptype derive from each
        # other
        @test p.label == "m_rate" && p.kind === :real && p.ptype == "real" && p.constant
        @test ParamSpec(:m, 1, (0, 2); ptype="integer").kind === :integer
        b = ParamSpec(:m_on, 0, (0, 1); kind=:bool)
        @test b.ptype == "categorical" && b.levels == [0.0, 1.0] && b.dist === :categorical
        c = ParamSpec(:m_mode, 2, (0, 2); kind=:categorical, levels=[0, 1, 2], label="Mode")
        @test c.label == "Mode" && c.levels == [0.0, 1.0, 2.0]
        @test_throws ArgumentError ParamSpec(:m_mode, 0, (0, 2); kind=:categorical)
        @test_throws ArgumentError ParamSpec(
            :m_mode, 5, (0, 9); kind=:categorical, levels=[0, 1]
        )
        @test_throws ArgumentError ParamSpec(:m_x, 0, (0, 1); kind=:other)
        @test_throws ArgumentError ParamSpec(:m_x, 0, (0, 1); kind=:real, ptype="integer")
        @test_throws ArgumentError ParamSpec(:m_x, 0.5, (0, 1); kind=:bool)

        s = LogSpec(:l, [:timesteps, :locations], [5, 4]; units="u")
        @test s.size == [5, 4] && s.units == "u" && isempty(s.labels)
        @test_throws ArgumentError LogSpec(:l, [:timesteps], [5, 4])

        bufs = alloc_logs([s, LogSpec(:m, [:timesteps], [5])])
        @test keys(bufs) == (:l, :m)
        @test size(bufs.l) == (5, 4) && all(iszero, bufs.m)
        @test alloc_logs(LogSpec[]) == (;)

        ctx = make_ctx(; inputs=(; x=1))
        ctx2 = with_cache(ctx, [1.0])
        @test ctx2.cache == [1.0] && ctx2.inputs === ctx.inputs
        @test ctx2 isa
            SubmodelContext{@NamedTuple{x::Int},Vector{Float64},Matrix{Float64}}
        @test ctx2.connectivity === ctx.connectivity   # shared, not copied
        @test size(ctx.connectivity) == (ctx.n_locs, ctx.n_locs)

        cov = zeros(2, 3, 4)
        v = CoralStepView(cov, zeros(4), fill(10.0, 4))
        @test cover(v) === cov && length(dhw(v)) == 4 && location_k(v)[1] == 10.0
        @test v isa CoralStepView{Array{Float64,3},Vector{Float64}}
    end

    @testset "mock contract" begin
        sm = ShrinkMock(0.25)
        @test has_phase(sm, AfterDisturbances())
        @test active(sm, (; mock_on=1.0)) && !active(sm, (; mock_on=0.0))
        ctx = make_ctx(; inputs=load_inputs(sm, "p"))
        ctx = with_cache(ctx, prepare(sm, ctx))
        @test ctx.cache.n_locs == 4
        bufs = alloc_logs(log_specs(sm, ctx))
        state = init_state(sm, ctx, (; mock_on=1.0), Xoshiro(1))

        cov = fill(2.0, 2, 3, 4)
        v = CoralStepView(cov, zeros(4), fill(10.0, 4))
        step!(sm, AfterDisturbances(), state, ctx, v, 2, 2001)
        record!(sm, bufs, state, 2)
        @test all(==(1.5), cov)
        @test bufs.mock_calls[2] == 1 && bufs.mock_cover[2] == 2
    end

    @testset "type stability and allocations" begin
        sm = ShrinkMock(0.01)
        ctx = make_ctx()
        bufs = alloc_logs(log_specs(sm, ctx))
        state = init_state(sm, ctx, (; mock_on=1.0), Xoshiro(1))
        cov = fill(1.0, 2, 3, 4)
        v = CoralStepView(cov, zeros(4), fill(10.0, 4))

        @inferred step!(sm, AfterDisturbances(), state, ctx, v, 2, 2001)
        @inferred record!(sm, bufs, state, 2)
        @inferred init_state(sm, ctx, (; mock_on=1.0), Xoshiro(1))
        @inferred has_phase(sm, AfterDisturbances())
        @inferred cover(v)

        # The connectivity matrix is reachable through the context, type-stable and free
        cm = ConnMock()
        cctx = make_ctx(; conn=reshape(collect(1.0:16.0), 4, 4))
        cstate = init_state(cm, cctx, (;), Xoshiro(1))
        @inferred step!(cm, AfterDisturbances(), cstate, cctx, v, 2, 2001)
        get_conn(c) = c.connectivity
        @inferred get_conn(cctx)
        run_conn(cm, cstate, cctx, v) =
            step!(cm, AfterDisturbances(), cstate, cctx, v, 2, 2001)
        run_conn(cm, cstate, cctx, v)
        @test (@allocated run_conn(cm, cstate, cctx, v)) == 0
        @test cstate.total == 3 * sum(1.0:16.0)   # three step! calls above

        # Wrap in functions so the measurement does not include global-scope overhead.
        run_step(sm, state, ctx, v) = step!(sm, AfterDisturbances(), state, ctx, v, 2, 2001)
        run_rec(sm, bufs, state) = record!(sm, bufs, state, 2)
        run_step(sm, state, ctx, v);
        run_rec(sm, bufs, state)  # compile
        @test (@allocated run_step(sm, state, ctx, v)) == 0
        @test (@allocated run_rec(sm, bufs, state)) == 0

        # Defaults are free too
        n = NoOpMock()
        run_noop(n, ctx, v) = step!(n, AfterDisturbances(), nothing, ctx, v, 2, 2001)
        run_noop(n, ctx, v)
        @test (@allocated run_noop(n, ctx, v)) == 0
    end
end
