# Mock submodels used to test the interface contract, in this package and in ADRIA.
# Include this file from a scope where `ADRIAinterfaces` is bound (e.g. Main).
# TEST SUPPORT ONLY: not part of the package API.
module MockSubmodels

using ..ADRIAinterfaces
using Random
import ..ADRIAinterfaces:
    phases, step!, record!, init_state, log_specs, parameters, active, prepare,
    load_inputs, validate, seed_excluded, submodel_name

export NoOpMock, ShrinkMock, RecordingMock, GrowMock, ConnMock, MockState

"""Implements nothing: everything comes from the no-op defaults."""
struct NoOpMock <: Submodel end

mutable struct MockState
    nsteps::Int
    last_tstep::Int
end

"""
Reduces all cover by the fraction `frac` at `AfterDisturbances`, and logs the total cover
it saw (`:mock_cover`, `tf` x 1) and the number of calls (`:mock_calls`, `tf`).
Parameter `mock_on` (default 1.0) switches it off per scenario when 0.
"""
struct ShrinkMock <: Submodel
    frac::Float64
end
phases(::ShrinkMock) = (AfterDisturbances(),)
submodel_name(::ShrinkMock) = :shrink_mock
parameters(::ShrinkMock) = [ParamSpec(:mock_on, 1.0, (0.0, 1.0); ptype="integer")]
active(::ShrinkMock, params::NamedTuple) = params.mock_on == 1.0
seed_excluded(::ShrinkMock) = [:mock_on]
load_inputs(::ShrinkMock, path::AbstractString) = (; path=String(path))
prepare(::ShrinkMock, ctx::SubmodelContext) = (; n_locs=ctx.n_locs)
function log_specs(::ShrinkMock, ctx::SubmodelContext)
    return [
        LogSpec(:mock_cover, [:timesteps], [ctx.tf]; units="m2"),
        LogSpec(:mock_calls, [:timesteps], [ctx.tf])
    ]
end
function init_state(::ShrinkMock, ctx::SubmodelContext, ::NamedTuple, ::Random.AbstractRNG)
    return MockState(0, 0)
end
function step!(
    sm::ShrinkMock, ::AfterDisturbances, state::MockState, ::SubmodelContext,
    view::CoralStepView, tstep::Int, year::Int
)::Nothing
    c = cover(view)
    @inbounds for i in eachindex(c)
        c[i] *= (1.0 - sm.frac)
    end
    state.nsteps += 1
    state.last_tstep = tstep
    return nothing
end
function record!(::ShrinkMock, bufs::NamedTuple, state::MockState, tstep::Int)::Nothing
    bufs.mock_calls[tstep] = state.nsteps
    bufs.mock_cover[tstep] = state.last_tstep
    return nothing
end

"""
Records every hook call (in order) into the shared `calls` vector, and the cover total seen
at `AfterDisturbances` into `seen`. Does not touch cover. `on=false` makes it inactive.
"""
struct RecordingMock <: Submodel
    tag::Symbol
    calls::Vector{Tuple{Symbol,Symbol,Int}}  # (tag, hook, tstep)
    seen::Vector{Float64}
    on::Bool
end
RecordingMock(tag::Symbol, calls; on::Bool=true) = RecordingMock(tag, calls, Float64[], on)
active(m::RecordingMock, ::NamedTuple) = m.on
phases(::RecordingMock) = (AfterDisturbances(),)
function load_inputs(m::RecordingMock, ::AbstractString)
    push!(m.calls, (m.tag, :load_inputs, 0))
    return nothing
end
function validate(m::RecordingMock, inputs, scenarios)::Nothing
    push!(m.calls, (m.tag, :validate, 0))
    return nothing
end
function prepare(m::RecordingMock, ctx::SubmodelContext)
    push!(m.calls, (m.tag, :prepare, 0))
    return nothing
end
function init_state(
    m::RecordingMock, ctx::SubmodelContext, ::NamedTuple, ::Random.AbstractRNG
)
    push!(m.calls, (m.tag, :init_state, 0))
    return nothing
end
function step!(
    m::RecordingMock, ::AfterDisturbances, state, ctx::SubmodelContext,
    view::CoralStepView, tstep::Int, year::Int
)::Nothing
    push!(m.calls, (m.tag, :step!, tstep))
    push!(m.seen, sum(cover(view)))
    return nothing
end
function record!(m::RecordingMock, bufs::NamedTuple, state, tstep::Int)::Nothing
    push!(m.calls, (m.tag, :record!, tstep))
    return nothing
end

"""Violates the decrease-only contract: multiplies cover by `factor` (> 1)."""
struct GrowMock <: Submodel
    factor::Float64
end
phases(::GrowMock) = (AfterDisturbances(),)
function step!(
    m::GrowMock, ::AfterDisturbances, state, ::SubmodelContext, view::CoralStepView,
    tstep::Int, year::Int
)::Nothing
    cover(view) .*= m.factor
    return nothing
end

"""
Reads the connectivity matrix from the context: `init_state` returns a state holding a
reference to it and `step!` accumulates its sum into `state.total` (no allocation). Records
nothing. Used to check that the host hands the domain's connectivity to submodels.
"""
struct ConnMock <: Submodel end
mutable struct ConnState{M<:AbstractMatrix{Float64}}
    conn::M
    total::Float64
end
phases(::ConnMock) = (AfterDisturbances(),)
init_state(::ConnMock, ctx::SubmodelContext, ::NamedTuple, ::Random.AbstractRNG) =
    ConnState(ctx.connectivity, 0.0)
function step!(
    ::ConnMock, ::AfterDisturbances, state::ConnState, ctx::SubmodelContext,
    ::CoralStepView, tstep::Int, year::Int
)::Nothing
    state.total += sum(ctx.connectivity)
    return nothing
end

end # module
