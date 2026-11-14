using Test
using Random
using ADRIAinterfaces
import ADRIAinterfaces: step!
include("MockSubmodels.jl")
using .MockSubmodels

function make_ctx(; n_locs=4, tf=5, inputs=nothing, cache=nothing)
    return SubmodelContext(
        n_locs, 2, 3, tf, collect(2000:(2000 + tf - 1)), [:a, :b],
        fill(10.0, n_locs), inputs, cache
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
        @test ctx2 isa SubmodelContext{@NamedTuple{x::Int},Vector{Float64}}

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

        # Wrap in functions so the measurement does not include global-scope overhead.
        run_step(sm, state, ctx, v) = step!(sm, AfterDisturbances(), state, ctx, v, 2, 2001)
        run_rec(sm, bufs, state) = record!(sm, bufs, state, 2)
        run_step(sm, state, ctx, v); run_rec(sm, bufs, state)  # compile
        @test (@allocated run_step(sm, state, ctx, v)) == 0
        @test (@allocated run_rec(sm, bufs, state)) == 0

        # Defaults are free too
        n = NoOpMock()
        run_noop(n, ctx, v) = step!(n, AfterDisturbances(), nothing, ctx, v, 2, 2001)
        run_noop(n, ctx, v)
        @test (@allocated run_noop(n, ctx, v)) == 0
    end
end
