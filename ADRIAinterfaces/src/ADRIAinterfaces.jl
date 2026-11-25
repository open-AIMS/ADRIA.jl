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
    ParamSpec(name, default, bounds; label=string(name), kind=nothing, levels=Float64[],
              ptype=nothing, dist=nothing, dist_params=[bounds...], constant=true,
              description="")

Declaration of one scenario parameter (factor) of a submodel. Plain data (strings, symbols
and numbers only): the host decides how to turn it into a sampled/constant factor (ADRIA
does so in `ADRIA.spec_rows`).

- `name::Symbol`: unique and prefixed with the submodel name, e.g. `:mymodel_rate`.
- `label::String`: human-readable name shown in the parameter table (defaults to `name`).
- `default::Float64` and `bounds::Tuple{Float64,Float64}` (lower, upper). Equal bounds
  declare a constant.
- `kind::Symbol`: `:real`, `:integer`, `:bool` or `:categorical`. Derived from `ptype` when
  not given (`"real"`, `"integer"`, `"categorical"`), `:real` if neither is given.
- `levels::Vector{Float64}`: the values of a `:bool` (`[0, 1]`, filled in) or
  `:categorical` factor (required).
- `ptype::String`: `"real"`, `"integer"` or `"categorical"`. Derived from `kind` when not
  given (`:bool` is `"categorical"`).
- `dist::Symbol`: sampling distribution, with parameters `dist_params` (the full range):
  `:uniform` (lo, hi), `:triang` (lo, hi, mode), `:discrete_ordered` (lo, hi, step),
  `:discrete_uniform` (lo, hi) or `:categorical` (the host uses `levels`, `dist_params` is
  ignored). Defaults
  to `:categorical` for `:bool`/`:categorical` kinds and `:uniform` otherwise.
- `constant::Bool`: if `true` (the default) the factor is held at `default` when sampling
  and `dist_params` is only the range it may be widened to.

Note: the struct grew from 7 to 11 fields, so the positional constructor changed. Use the
keyword constructor. Version 0.1.0 of this package has not been tagged or released (it is
covered by the 0.2.0-dev window), so this is not a breaking release change.
"""
struct ParamSpec
    name::Symbol
    default::Float64
    bounds::Tuple{Float64,Float64}
    ptype::String
    dist::Symbol
    dist_params::Vector{Float64}
    description::String
    label::String
    kind::Symbol
    levels::Vector{Float64}
    constant::Bool
end

const _PARAM_KINDS = (:real, :integer, :bool, :categorical)
const _PTYPE_KIND = Dict(
    "real" => :real, "integer" => :integer, "categorical" => :categorical
)

function ParamSpec(
    name::Symbol, default::Real, bounds::Tuple{<:Real,<:Real};
    ptype::Union{Nothing,AbstractString}=nothing, dist::Union{Nothing,Symbol}=nothing,
    dist_params::AbstractVector{<:Real}=Float64[bounds[1], bounds[2]],
    description::AbstractString="", label::AbstractString=String(name),
    kind::Union{Nothing,Symbol}=nothing, levels::AbstractVector{<:Real}=Float64[],
    constant::Bool=true
)
    lo, hi = Float64(bounds[1]), Float64(bounds[2])
    lo <= hi || throw(ArgumentError("$name: lower bound $lo exceeds upper bound $hi"))
    lo <= default <= hi ||
        throw(ArgumentError("$name: default $default is outside bounds ($lo, $hi)"))

    if kind === nothing
        kind = ptype === nothing ? :real : get(_PTYPE_KIND, String(ptype), nothing)
        kind === nothing && throw(ArgumentError("$name: unknown ptype \"$ptype\""))
    end
    kind in _PARAM_KINDS ||
        throw(ArgumentError("$name: kind must be one of $_PARAM_KINDS, got :$kind"))
    if ptype !== nothing && get(_PTYPE_KIND, String(ptype), kind) !== kind &&
        !(kind === :bool && ptype == "categorical")
        throw(ArgumentError("$name: ptype \"$ptype\" does not match kind :$kind"))
    end
    kind === :bool && !(default in (0, 1)) &&
        throw(ArgumentError("$name: a :bool factor has default 0 or 1"))
    ptype === nothing && (ptype = if kind === :real
        "real"
    elseif kind === :integer
        "integer"
    else
        "categorical"
    end)

    levels = Float64.(levels)
    if kind === :bool
        isempty(levels) && (levels = [0.0, 1.0])
        levels == [0.0, 1.0] ||
            throw(ArgumentError("$name: a :bool factor has levels [0, 1]"))
    elseif kind === :categorical
        length(levels) >= 2 ||
            throw(ArgumentError("$name: a :categorical factor needs at least 2 levels"))
        default in levels ||
            throw(ArgumentError("$name: default $default is not one of the levels"))
    end
    if dist === nothing
        dist = kind in (:bool, :categorical) ? :categorical : :uniform
    end

    return ParamSpec(
        name, Float64(default), (lo, hi), String(ptype), dist,
        Float64.(dist_params), String(description), String(label), kind, levels,
        constant
    )
end

"""
    LogSpec(name, dims, size; labels=Dict(), units="", description="")

Declaration of one per-scenario output log of a submodel. The host allocates a
`Float64` array of `size` for it, passes it to [`record!`](@ref) and returns it among the
run results, grouped by submodel: `submodel_logs.<submodel_name>.<name>`, so `name` only
needs to be unique within the submodel. When the host stores logs in a result set (today
only the built-in COTS logs; the logs of attached submodels are returned by a run but not
stored yet), it appends a scenario axis.

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
    SubmodelContext{I,C,M}

What a submodel may see about the domain and run; built by the host (once per scenario run)
and passed to [`init_state`](@ref), [`log_specs`](@ref) and [`step!`](@ref). Read-only by
contract. `inputs::I` is the result of [`load_inputs`](@ref) and `cache::C` the result of
[`prepare`](@ref); both concrete type parameters so field access is type-stable.
`connectivity::M` is the domain's connectivity matrix (`M <: AbstractMatrix{Float64}`;
source locations in rows, sink locations in columns, `n_locs x n_locs`, as used by the
core model, i.e. not normalised to transition probabilities). It is shared with the host,
not copied: read-only by contract.

Fields: `n_locs`, `n_groups`, `n_sizes`, `tf` (number of timesteps), `timeframe` (year of
each timestep), `group_names` (functional group names, so submodels map groups by name and
not by position), `location_k` (habitable area per location), `connectivity`, `inputs`, `cache`.
"""
struct SubmodelContext{I,C,M<:AbstractMatrix{Float64}}
    n_locs::Int
    n_groups::Int
    n_sizes::Int
    tf::Int
    timeframe::Vector{Int}
    group_names::Vector{Symbol}
    location_k::Vector{Float64}
    connectivity::M
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
        ctx.location_k, ctx.connectivity, ctx.inputs, cache
    )
end

"""
    CoralStepView(cover, dhw, location_k)

The only window a [`step!`](@ref) hook has on the coral model at a phase.

- `cover(v)`: coral cover (`n_groups x n_sizes x n_locs`) as a **fraction of each
  location's habitable area** (not absolute area; multiply by `location_k` for area).
  **Mutable in place, decrease only**: a hook may reduce any entry (down to zero) and must
  never increase one. The host snapshots cover before the submodel phase and, when any
  submodel is attached, checks per cell that the value after the phase is `<=` the value
  before (exact comparison) and raises an error otherwise, so an increase of any size is
  caught, including one that only offsets a disturbance loss or fills an empty cell. The
  reduction is then propagated to the coral groups through the survival ratio.
- `dhw(v)`: degree heating weeks of this timestep per location. Read-only by contract.
- `location_k(v)`: habitable area per location (absolute area). Read-only by contract.

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

Name of the submodel. It identifies the submodel in messages and provenance records and
names the group of its logs: `submodel_logs.<name>` among the results of a run and
`rs.logs.<name>` in a result set. It must be a `Symbol` that is a valid identifier, unique
among the submodels attached to a domain (define it per instance to attach two instances
of one type), and not `core`, `cots` or the name of a core log array; the host checks this
when the submodel is attached. Default: the type name.
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

Output logs the submodel wants. Their names must be unique within the submodel (logs are
grouped by [`submodel_name`](@ref)). Default: none.
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
