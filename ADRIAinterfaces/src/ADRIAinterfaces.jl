"""
    ADRIAinterfaces

Contracts between the ADRIA run loop and optional submodels (models that run alongside the
coral model, for example a predator or a disease model).

This package holds contracts only: no ADRIA types, no model code and no simulation loop.
ADRIA depends on it and calls the hooks defined here at fixed points of a scenario run; a
submodel package (or an ADRIA-side adapter for it) depends on it and implements them.

# Call order (per scenario run)

1. Once per run: [`validate`](@ref), then [`prepare`](@ref) (read-only; the result is shared
   by all scenarios and threads). [`load_inputs`](@ref) is called once when the submodel is
   attached to a domain.
2. Once per scenario: [`active`](@ref) decides whether the submodel takes part. If so,
   [`log_specs`](@ref) is queried, log buffers are allocated, and [`init_state`](@ref) is
   called.
3. Every timestep, in each [`Phase`](@ref) the submodel lists in [`phases`](@ref):
   [`step!`](@ref), then [`record!`](@ref).

Within a phase, the order is the one documented on the phase type: core model processes
first, then submodels in attach order, then the core model carries on.

Every hook has a no-op default, so a submodel implements only what it needs. State returned
by [`init_state`](@ref) and the context/view types are concretely typed so that hooks can be
type-stable and allocation free inside the timestep loop; implementations should keep their
state concrete too.
"""
module ADRIAinterfaces

using Random: AbstractRNG

export Submodel, Phase, AfterDisturbances
export ParamSpec, LogSpec, SubmodelContext, CoralStepView
export phases, has_phase, step!, record!, init_state, log_specs, validate, prepare
export active, seed_excluded, load_inputs, parameters, submodel_name
export cover, dhw, location_k, with_cache, alloc_logs

"""Version of this interface. Breaking changes to a hook signature bump the major version;
adding phases or optional hooks does not."""
const SUBMODEL_INTERFACE_VERSION = v"0.1.0"

# --- Types -------------------------------------------------------------------------------

"""
    Submodel

Supertype of all submodels. Subtype it and implement the hooks you need; the generic
functions below have no-op defaults. Instances are immutable configuration: anything that
changes during a scenario belongs in the state returned by [`init_state`](@ref).
"""
abstract type Submodel end

"""
    Phase

A point in the yearly cycle of the core model at which submodels are stepped. A submodel
lists the phases it takes part in with [`phases`](@ref) and implements
[`step!`](@ref) for each of them. New phases can be added without breaking this interface.
"""
abstract type Phase end

"""
    AfterDisturbances()

The only phase currently called. It happens in each timestep after the core model's
disturbances (bleaching and cyclone mortality) and before survival is propagated to the
coral groups. Order within the phase: core disturbances, then submodels in attach order,
then survival propagation.

At this phase a submodel may only reduce coral cover (see [`CoralStepView`](@ref)).
"""
struct AfterDisturbances <: Phase end

"""
    ParamSpec(name, default, bounds; ptype="real", dist=:uniform, dist_params=[bounds...],
              description="")

Declaration of one scenario parameter (factor) of a submodel. Plain data: the host decides
how to turn it into a sampled/constant factor.

- `name::Symbol`: unique and prefixed with the submodel name, e.g. `:mymodel_rate`.
- `default::Float64` and `bounds::Tuple{Float64,Float64}` (lower, upper). Equal bounds
  declare a constant.
- `ptype::String`: `"real"`, `"integer"` or `"categorical"`.
- `dist::Symbol`: name of the sampling distribution (`:uniform`, `:triang`, ...), with
  parameters `dist_params`.
"""
struct ParamSpec
    name::Symbol
    default::Float64
    bounds::Tuple{Float64,Float64}
    ptype::String
    dist::Symbol
    dist_params::Vector{Float64}
    description::String
end
function ParamSpec(
    name::Symbol, default::Real, bounds::Tuple{<:Real,<:Real};
    ptype::AbstractString="real", dist::Symbol=:uniform,
    dist_params::AbstractVector{<:Real}=Float64[bounds[1], bounds[2]],
    description::AbstractString=""
)
    lo, hi = Float64(bounds[1]), Float64(bounds[2])
    lo <= hi || throw(ArgumentError("$name: lower bound $lo exceeds upper bound $hi"))
    lo <= default <= hi ||
        throw(ArgumentError("$name: default $default is outside bounds ($lo, $hi)"))
    return ParamSpec(
        name, Float64(default), (lo, hi), String(ptype), dist,
        Float64.(dist_params), String(description)
    )
end

"""
    LogSpec(name, dims, size; labels=Dict(), units="", description="")

Declaration of one per-scenario output log of a submodel. The host allocates a
`Float64` array of `size` for it, passes it to [`record!`](@ref) and returns it among the
run results under `name`. A scenario axis is appended by the host when results are stored.

`dims` names each axis (e.g. `[:timesteps, :age_class, :locations]`); `labels` optionally
maps a dim to its labels.
"""
struct LogSpec
    name::Symbol
    dims::Vector{Symbol}
    size::Vector{Int}
    labels::Dict{Symbol,Vector{String}}
    units::String
    description::String
end
function LogSpec(
    name::Symbol, dims::AbstractVector{Symbol}, size::AbstractVector{<:Integer};
    labels::AbstractDict=Dict{Symbol,Vector{String}}(), units::AbstractString="",
    description::AbstractString=""
)
    length(dims) == length(size) ||
        throw(ArgumentError("LogSpec $name: dims and size must have equal length"))
    return LogSpec(
        name, collect(dims), collect(Int, size),
        Dict{Symbol,Vector{String}}(labels), String(units), String(description)
    )
end

"""
    SubmodelContext{I,C}

What a submodel may see about the domain and run; built by the host (once per scenario run)
and passed to [`init_state`](@ref), [`log_specs`](@ref) and [`step!`](@ref). Read-only by
contract. `inputs::I` is the result of [`load_inputs`](@ref) and `cache::C` the result of
[`prepare`](@ref); both concrete type parameters so field access is type-stable.

Fields: `n_locs`, `n_groups`, `n_sizes`, `tf` (number of timesteps), `timeframe` (year of
each timestep), `group_names` (functional group names, so submodels map groups by name and
not by position), `location_k` (habitable area per location), `inputs`, `cache`.
"""
struct SubmodelContext{I,C}
    n_locs::Int
    n_groups::Int
    n_sizes::Int
    tf::Int
    timeframe::Vector{Int}
    group_names::Vector{Symbol}
    location_k::Vector{Float64}
    inputs::I
    cache::C
end

"""
    with_cache(ctx::SubmodelContext, cache)

Copy of `ctx` with `cache` set (type parameter `C` follows).
"""
function with_cache(ctx::SubmodelContext, cache)
    return SubmodelContext(
        ctx.n_locs, ctx.n_groups, ctx.n_sizes, ctx.tf, ctx.timeframe, ctx.group_names,
        ctx.location_k, ctx.inputs, cache
    )
end

"""
    CoralStepView(cover, dhw, location_k)

The only window a [`step!`](@ref) hook has on the coral model at a phase.

- `cover(v)`: coral cover (`n_groups x n_sizes x n_locs`, as absolute area, the same units as
  `location_k`). **Mutable in place, decrease only**: a hook may reduce any entry (down to
  zero) and must never increase one. The host propagates the reduction to the coral
  groups through the survival ratio, and asserts that survival is `<= 1`, so an increase
  fails there.
- `dhw(v)`: degree heating weeks of this timestep per location. Read-only by contract.
- `location_k(v)`: habitable area per location. Read-only by contract.

Immutable; the fields are references to host arrays, so building one does not copy.
"""
struct CoralStepView{A<:AbstractArray{Float64,3},V<:AbstractVector{Float64}}
    cover::A
    dhw::V
    location_k::V
end

cover(v::CoralStepView) = v.cover
dhw(v::CoralStepView) = v.dhw
location_k(v::CoralStepView) = v.location_k

# --- Hooks (all with no-op defaults) -------------------------------------------------------

"""
    submodel_name(sm::Submodel)::Symbol

Short identifier, used in messages. Default: the type name.
"""
submodel_name(sm::Submodel)::Symbol = nameof(typeof(sm))

"""
    phases(sm::Submodel)::Tuple

Phases (instances of [`Phase`](@ref) subtypes) in which `sm` is stepped. Default `()`: never
stepped.
"""
phases(::Submodel)::Tuple = ()

"""
    has_phase(sm::Submodel, phase::Phase)::Bool

Whether `phase` is listed in `phases(sm)`.
"""
has_phase(sm::Submodel, phase::Phase)::Bool = any(p -> p isa typeof(phase), phases(sm))

"""
    parameters(sm::Submodel)::Vector{ParamSpec}

Scenario parameters the submodel declares. Default: none.
"""
parameters(::Submodel)::Vector{ParamSpec} = ParamSpec[]

"""
    load_inputs(sm::Submodel, data_path::AbstractString)

Called once when the submodel is attached to a domain, with the path of the domain data
package. Returns the submodel's static inputs (stored as `ctx.inputs`; must be treated as
read-only afterwards). Default: `nothing`.
"""
load_inputs(::Submodel, ::AbstractString) = nothing

"""
    validate(sm::Submodel, inputs, scenarios)::Nothing

Called once per run, before any scenario runs, with the scenario table (a `DataFrame`-like
object with one column per factor). Throw an error to reject a run, e.g. when parameters
need inputs the domain does not provide. Must not mutate its arguments. Default: no-op.
"""
validate(::Submodel, inputs, scenarios)::Nothing = nothing

"""
    prepare(sm::Submodel, ctx::SubmodelContext)

Called once per run (after [`validate`](@ref)) with a context whose `cache` is `nothing`.
Returns a cache (stored as `ctx.cache`) that is **shared read-only** by all scenarios and
threads, e.g. a sparse connectivity matrix. Default: `nothing`.
"""
prepare(::Submodel, ::SubmodelContext) = nothing

"""
    active(sm::Submodel, params::NamedTuple)::Bool

Whether the submodel takes part in a scenario. `params` holds this submodel's declared
parameters ([`parameters`](@ref)) by name, with the scenario's values. When `false`, no hook
is called and no log is allocated for that scenario. Default: `true`.
"""
active(::Submodel, params::NamedTuple)::Bool = true

"""
    seed_excluded(sm::Submodel)::Vector{Symbol}

Names of factors that must not contribute to the scenario's random seed (so that switching
or tuning the submodel does not reshuffle the coral model's random draws). Default: none.
"""
seed_excluded(::Submodel)::Vector{Symbol} = Symbol[]

"""
    log_specs(sm::Submodel, ctx::SubmodelContext)::Vector{LogSpec}

Output logs the submodel wants. Default: none.
"""
log_specs(::Submodel, ::SubmodelContext)::Vector{LogSpec} = LogSpec[]

"""
    alloc_logs(specs::AbstractVector{LogSpec})::NamedTuple

Zero-filled `Float64` buffers, one per spec, keyed by `LogSpec.name`.
"""
function alloc_logs(specs::AbstractVector{LogSpec})::NamedTuple
    return NamedTuple{Tuple(s.name for s in specs)}(Tuple(zeros(s.size...) for s in specs))
end

"""
    init_state(sm::Submodel, ctx::SubmodelContext, params::NamedTuple, rng::AbstractRNG)

Called once per active scenario. Returns the mutable per-scenario state handed to
[`step!`](@ref) and [`record!`](@ref). Return a concretely typed value. `params` holds the
submodel's declared parameters by name. `rng` is the scenario's random stream (interim: the
coral model's; a dedicated stream per submodel is planned, so do not rely on draw
positions). Default: `nothing`.
"""
init_state(::Submodel, ::SubmodelContext, ::NamedTuple, ::AbstractRNG) = nothing

"""
    step!(sm::Submodel, phase::Phase, state, ctx::SubmodelContext, view::CoralStepView,
          tstep::Int, year::Int)::Nothing

Advance the submodel for one timestep at `phase` (only called for phases in
[`phases`](@ref)). May mutate `state`, and `cover(view)` as the phase's contract allows
(decrease only at [`AfterDisturbances`](@ref)); must not mutate anything else, including
`ctx`. Called in attach order. Default: no-op.
"""
function step!(
    ::Submodel, ::Phase, state, ::SubmodelContext, ::CoralStepView, tstep::Int, year::Int
)::Nothing
    return nothing
end

"""
    record!(sm::Submodel, bufs::NamedTuple, state, tstep::Int)::Nothing

Called once per timestep for every active submodel, after all submodels' [`step!`](@ref)
calls at [`AfterDisturbances`](@ref) of that timestep (so: all steps in attach order, then all
records in attach order), whether or not the submodel lists any phase. Writes the timestep's
values into `bufs`, the buffers allocated for this submodel from [`log_specs`](@ref).
Must not mutate coral cover. Default: no-op.
"""
record!(::Submodel, bufs::NamedTuple, state, tstep::Int)::Nothing = nothing

end
