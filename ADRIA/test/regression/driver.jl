# Pinned baseline driver: runs a small fixed-seed ensemble on an ADRIA checkout.
#
#   julia --project=<checkout>/ADRIA -t4 driver.jl <checkout>/ADRIA <outdir>
#
# Writes one result store per case under <outdir>/<case>/ (ADRIA_OUTPUT_DIR is set per
# case) and <outdir>/META (Julia version and thread count of this process). Nothing is
# written inside the checkout. Cases:
#   mixed    guided + unguided scenarios with interventions, COTS off and on (plus juvenile
#            predation on for some rows): exercises outcomes, ranks, intervention logs and
#            the COTS logs.
#   cots_off the same guided scenarios, all with COTS off (no COTS logs in the store).
#   attached 4 of the guided scenarios with COTS off and 4 with COTS on, on a domain with
#            `ShrinkMock(0.1)` attached (defined below, so the reference and the candidate
#            run identical mock code whatever the checkout's test support holds): exercises
#            the submodel hook path (step, decrease-only check, record) and the propagation
#            of the reduced cover. Needs `_attach_submodel!` and the ADRIAinterfaces hook
#            API in the checkout (present since the hook plumbing landed).
#
# Per-stage variants: if a stage renames accessors, copy this file to e.g. driver_stage1.jl
# and point ADRIA_REGRESSION_DRIVER at it (see README.md). Keep the ensemble identical.
using Random
using ADRIA
using ADRIA: ADRIAinterfaces

"""
Pinned copy of `ShrinkMock` from ADRIAinterfaces/test/MockSubmodels.jl (as of the hook
plumbing, 331795c88): reduces all cover by `frac` at `AfterDisturbances`. Kept here, not
included from the checkout, so that a change to the checkout's mocks cannot change the
`attached` case. Uses only interface functions present since the hook plumbing landed;
all of them are qualified, so the checkout's exports do not matter.
"""
module PinnedShrinkMock
using Random
using ADRIA: ADRIAinterfaces
const AI = ADRIAinterfaces

struct ShrinkMock <: AI.Submodel
    frac::Float64
end
mutable struct MockState
    nsteps::Int
    last_tstep::Int
end
AI.phases(::ShrinkMock) = (AI.AfterDisturbances(),)
AI.submodel_name(::ShrinkMock) = :shrink_mock
AI.parameters(::ShrinkMock) = [AI.ParamSpec(:mock_on, 1.0, (0.0, 1.0); ptype="integer")]
AI.active(::ShrinkMock, params::NamedTuple) = params.mock_on == 1.0
AI.seed_excluded(::ShrinkMock) = [:mock_on]
AI.load_inputs(::ShrinkMock, path::AbstractString) = (; path=String(path))
AI.prepare(::ShrinkMock, ctx::AI.SubmodelContext) = (; n_locs=ctx.n_locs)
function AI.log_specs(::ShrinkMock, ctx::AI.SubmodelContext)
    return [
        AI.LogSpec(:mock_cover, [:timesteps], [ctx.tf]; units="m2"),
        AI.LogSpec(:mock_calls, [:timesteps], [ctx.tf])
    ]
end
function AI.init_state(
    ::ShrinkMock, ::AI.SubmodelContext, ::NamedTuple, ::Random.AbstractRNG
)
    return MockState(0, 0)
end
function AI.step!(
    sm::ShrinkMock, ::AI.AfterDisturbances, state::MockState, ::AI.SubmodelContext,
    view::AI.CoralStepView, tstep::Int, year::Int
)::Nothing
    c = AI.cover(view)
    @inbounds for i in eachindex(c)
        c[i] *= (1.0 - sm.frac)
    end
    state.nsteps += 1
    state.last_tstep = tstep
    return nothing
end
function AI.record!(::ShrinkMock, bufs::NamedTuple, state::MockState, tstep::Int)::Nothing
    bufs.mock_calls[tstep] = state.nsteps
    bufs.mock_cover[tstep] = state.last_tstep
    return nothing
end
end # module
using .PinnedShrinkMock: ShrinkMock

const N_GUIDED = 8
const N_UNGUIDED = 4
const SEED = 20261116
const SHRINK_FRAC = 0.1

function build_cases(dom)
    Random.seed!(SEED)
    guided = ADRIA.sample_guided(dom, N_GUIDED)
    unguided = ADRIA.sample_unguided(dom, N_UNGUIDED)

    off = copy(guided)
    off.cots_on .= 0

    on = copy(guided)
    on.cots_on .= 1
    on.cots_juvenile_predation[1:2:end] .= 1

    unguided_on = copy(unguided)
    unguided_on.cots_on .= 1

    mixed = vcat(off, on, unguided, unguided_on)
    attached = vcat(off[1:4, :], on[1:4, :])
    return Dict("mixed" => mixed, "cots_off" => off, "attached" => attached)
end

function main(checkout::String, outdir::String)
    mkpath(outdir)
    # Apply the default configuration now: otherwise the first run_scenarios call does it
    # and overrides ADRIA_OUTPUT_DIR (the cwd must not hold a config.toml; the runner uses
    # an empty scratch dir)
    ADRIA.setup()
    dom = ADRIA.load_domain(joinpath(checkout, "test", "data", "Test_domain"), "45")
    cases = build_cases(dom)
    dom_attached = ADRIA.load_domain(
        joinpath(checkout, "test", "data", "Test_domain"), "45"
    )
    ADRIA._attach_submodel!(dom_attached, ShrinkMock(SHRINK_FRAC))
    domains = Dict("attached" => dom_attached)
    for name in sort!(collect(keys(cases)))
        case_dir = joinpath(outdir, name)
        mkpath(case_dir)
        withenv("ADRIA_OUTPUT_DIR" => case_dir) do
            t = @elapsed ADRIA.run_scenarios(
                get(domains, name, dom), cases[name], "45"; show_progress=false
            )
            println(
                "case $(name): $(size(cases[name], 1)) scenarios, $(round(t; digits=1)) s"
            )
        end
    end
    # Run metadata of this process (run_regression.sh appends the driver and the source):
    # recorded here, so it is the Julia and thread count that produced the results
    write(
        joinpath(outdir, "META"),
        "julia_version=$(VERSION)\nthreads=$(Threads.nthreads())\n"
    )
    return nothing
end

main(abspath(ARGS[1]), abspath(ARGS[2]))
