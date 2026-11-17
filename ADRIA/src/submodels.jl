# Plumbing that lets the run loop call the hooks of ADRIAinterfaces on the submodels attached
# to a domain. Interim: there is no public attach API yet (see `_attach_submodel!`), and
# COTS does not go through this path. With no submodel attached, none of this runs and the
# run is unchanged (no extra RNG draws, no extra allocations).

using TOML: TOML
using ADRIAinterfaces: ADRIAinterfaces
using ADRIAinterfaces: Submodel, SubmodelContext, CoralStepView, AfterDisturbances, LogSpec

"""
    AttachedSubmodel(model::Submodel, inputs)

A submodel attached to a domain together with the result of `load_inputs` for that domain.
Held in `domain.spec_store.submodels`, in attach order.
"""
struct AttachedSubmodel
    model::Submodel
    inputs::Any
end

"""
    attached_submodels(domain::Domain)::Vector{AttachedSubmodel}

Submodels attached to `domain`, in attach order (empty by default).
"""
attached_submodels(domain::Domain)::Vector{AttachedSubmodel} = domain.spec_store.submodels

"""
Names a submodel may not have: its name is the group its logs are stored and loaded under
(`rs.logs.<name>`). `core` is kept for the core logs; `cots` is the group of the built-in
COTS logs, until COTS is attached like any other submodel.
"""
const RESERVED_SUBMODEL_NAMES = (:core, :cots)

"""
    _reserved_submodel_names()::Tuple

`RESERVED_SUBMODEL_NAMES` and the names of the core log arrays (`CORE_LOG_NAMES`), which
share the log group with the submodel subgroups.
"""
_reserved_submodel_names()::Tuple = (RESERVED_SUBMODEL_NAMES..., Symbol.(CORE_LOG_NAMES)...)

"""
    _check_submodel_name(domain::Domain, model::Submodel)::Symbol

Name of `model` (`ADRIAinterfaces.submodel_name`), which names the group of its logs in the
result set. Throws an `ArgumentError` if the name is not a `Symbol`, is not a valid Julia
identifier, is reserved (`core`, `cots` or the name of a core log, see
[`_reserved_submodel_names`](@ref)) or is already used by a submodel attached to `domain`.
Two instances of one model type share a name unless the type defines `submodel_name` per
instance, so such instances cannot be attached together otherwise.
"""
function _check_submodel_name(domain::Domain, model::Submodel)::Symbol
    name = ADRIAinterfaces.submodel_name(model)
    name isa Symbol || throw(
        ArgumentError(
            "`submodel_name` of a `$(typeof(model))` returned the $(typeof(name)) " *
            "`$(repr(name))`; it must return a Symbol"
        )
    )

    Base.isidentifier(name) || throw(
        ArgumentError(
            "Submodel name `$name` is not a valid identifier; it names the group of the " *
            "submodel's logs in the result set"
        )
    )
    reserved = _reserved_submodel_names()
    name in reserved && throw(
        ArgumentError(
            "Submodel name `$name` is reserved (reserved names: $(join(reserved, ", ")))"
        )
    )

    any(a -> ADRIAinterfaces.submodel_name(a.model) == name, attached_submodels(domain)) &&
        throw(
            ArgumentError(
                "A submodel named `$name` is already attached; attached submodels need " *
                "distinct names (define `submodel_name` for the instance)"
            )
        )

    return name
end

"""
    _attach_submodel!(domain::Domain, model::Submodel)::Domain

Interim mechanism (not the final attach API): attach `model` to `domain` after the
already attached submodels. Calls `ADRIAinterfaces.load_inputs(model, dpkg_path)` once.
Attached submodels are copied with the domain's `SpecStore` (`copy`/`deepcopy`). The model's
name must be valid and unique, and not reserved: see [`_check_submodel_name`](@ref).
"""
function _attach_submodel!(domain::Domain, model::Submodel)::Domain
    _check_submodel_name(domain, model)
    inputs = ADRIAinterfaces.load_inputs(model, String(domain.env_layer_md.dpkg_path))
    push!(domain.spec_store.submodels, AttachedSubmodel(model, inputs))
    return domain
end

"""Remove all submodels attached with `_attach_submodel!`."""
function _detach_submodels!(domain::Domain)::Domain
    empty!(domain.spec_store.submodels)
    return domain
end

# Context for submodel `a`, before `prepare` (cache is `nothing`).
function _submodel_context(domain::Domain, a::AttachedSubmodel)
    return SubmodelContext(
        domain.coral_growth.n_locs,
        domain.coral_growth.n_groups,
        domain.coral_growth.n_sizes,
        length(domain.env_layer_md.timeframe),
        Vector{Int}(domain.env_layer_md.timeframe),
        functional_group_names(),
        loc_k_area(domain),
        domain.conn.data,
        a.inputs,
        nothing
    )
end

"""
    _validate_submodels(domain::Domain, scens::DataFrame)::Nothing

Call `validate` of each attached submodel (once per run, before any scenario runs).
"""
function _validate_submodels(domain::Domain, scens::DataFrame)::Nothing
    for a in attached_submodels(domain)
        ADRIAinterfaces.validate(a.model, a.inputs, scens)
    end
    return nothing
end

"""
    _prepare_submodels(domain::Domain)::Vector{Any}

Call `prepare` of each attached submodel (once per run). The caches are shared, read-only,
by all scenarios and threads.
"""
function _prepare_submodels(domain::Domain)::Vector{Any}
    return Any[
        ADRIAinterfaces.prepare(a.model, _submodel_context(domain, a))
        for a in attached_submodels(domain)
    ]
end

# Factors excluded from the random seed by the attached submodels.
function _submodel_seed_excluded(domain::Domain)::Vector{String}
    subs = attached_submodels(domain)
    isempty(subs) && return String[]
    return String[
        string(f) for a in subs for f in ADRIAinterfaces.seed_excluded(a.model)
    ]
end

# Scenario's values of the parameters a submodel declares (declared default if the scenario
# has no such factor).
function _submodel_params(model::Submodel, param_set::YAXArray, factor_names)::NamedTuple
    specs = ADRIAinterfaces.parameters(model)
    names = Tuple(s.name for s in specs)
    vals = Tuple(
        string(s.name) in factor_names ? Float64(param_set[At(string(s.name))]) :
        s.default for s in specs
    )
    return NamedTuple{names}(vals)
end

# Per-scenario state of one active submodel. Concretely typed, so that a tuple of slots
# is stepped without dynamic dispatch.
struct _SubmodelSlot{M<:Submodel,C,S,B<:NamedTuple}
    model::M
    ctx::C
    state::S
    bufs::B
end

"""
    _init_submodels(domain, param_set, factor_names, caches, rng)::Tuple

Per scenario: for each attached submodel that is `active` for this scenario, build its
context, allocate its logs and call `init_state`. Returns a tuple of slots (empty when
nothing is attached, which costs nothing).
"""
function _init_submodels(
    domain::Domain, param_set::YAXArray, factor_names, caches::Union{Nothing,Vector{Any}},
    rng
)::Tuple
    subs = attached_submodels(domain)
    isempty(subs) && return ()
    slots = []
    for (i, a) in enumerate(subs)
        params = _submodel_params(a.model, param_set, factor_names)
        ADRIAinterfaces.active(a.model, params) || continue
        ctx = ADRIAinterfaces.with_cache(_submodel_context(domain, a), caches[i])
        specs = ADRIAinterfaces.log_specs(a.model, ctx)
        _check_unique_log_names(ADRIAinterfaces.submodel_name(a.model), specs)
        bufs = ADRIAinterfaces.alloc_logs(specs)
        state = ADRIAinterfaces.init_state(a.model, ctx, params, rng)
        push!(slots, _SubmodelSlot(a.model, ctx, state, bufs))
    end

    return Tuple(slots)
end

"""
    _check_unique_log_names(submodel::Symbol, specs::AbstractVector{LogSpec})::Nothing

Throw an `ArgumentError` if `specs` declares a log name twice. Log names need only be
unique within a submodel: logs are grouped by submodel name.
"""
function _check_unique_log_names(submodel::Symbol, specs::AbstractVector{LogSpec})::Nothing
    names = [s.name for s in specs]
    length(unique(names)) == length(names) || throw(
        ArgumentError("Submodel `$submodel` declares logs with the same name: $names")
    )

    return nothing
end

"""
    _submodel_logs(slots::Tuple)::NamedTuple

Logs of the active submodels of a run that have any, grouped by submodel name:
`(name = (log = buffer, ...), ...)`. Returned by `run_model` as `submodel_logs`; not to be
confused with the `submodel_logs` of the collected scenario results (the COTS logs written
to the result store, see `_cots_log_results`).
"""
function _submodel_logs(slots::Tuple)::NamedTuple
    logging = filter(s -> !isempty(s.bufs), collect(Any, slots))
    isempty(logging) && return (;)

    names = Tuple(ADRIAinterfaces.submodel_name(s.model) for s in logging)
    return NamedTuple{names}(Tuple(s.bufs for s in logging))
end

"""
    _step_submodels!(slots, phase, view, tstep, year)::Nothing

Phase hook position: after core disturbances (bleaching, cyclones), before survival is
propagated. Submodels step in attach order, then record.

Function barrier: `slots` is only known as an abstract `Tuple` in `run_model`, so each call
to these hooks is a dynamic dispatch. The callee is compiled for the concrete slot types
(once per submodel combination) and is type-stable inside, and the time loop is not
recompiled per combination. The dispatch boxes some arguments (tens of bytes per step).
`@noinline` is not what makes the barrier; it keeps the calls out-of-line, as does the
`@noinline` on `_check_decrease_only`, which keeps its error path out of `run_model`.
"""
@noinline function _step_submodels!(
    slots::Tuple, phase::ADRIAinterfaces.Phase, view::CoralStepView, tstep::Int, year::Int
)::Nothing
    foreach(slots) do s
        if ADRIAinterfaces.has_phase(s.model, phase)
            ADRIAinterfaces.step!(s.model, phase, s.state, s.ctx, view, tstep, year)
        end
    end
    return nothing
end

@noinline function _record_submodels!(slots::Tuple, tstep::Int)::Nothing
    foreach(s -> ADRIAinterfaces.record!(s.model, s.bufs, s.state, tstep), slots)
    return nothing
end

"""
    _check_decrease_only(before, after, tstep)

Enforce the submodel contract that cover may only decrease: per cell, `after <= before`
(exact), where `before` is the cover snapshot taken just before the submodel phase.
Only called when submodels are attached. The comparison is exact on purpose, so a
submodel that rescales cover must not round a value up by even one ulp. NaN and negative
values also fail.

Caveat (Stage 3a, COTS onto the hooks): if the host ever handed over a negative value (it
currently does not), this check fails whatever the submodel does, since clamping it up to
`0.0` is an increase and leaving it is negative. Today COTSMod's
`clamp!(C_cover_t, 0.0, 1.0)` silently repairs such values on the unchecked COTS path; on
the hook path the same input would raise here instead.
"""
@noinline function _check_decrease_only(
    before::AbstractArray{Float64,3}, after::AbstractArray{Float64,3}, tstep::Int
)::Nothing
    # Fast path: branch-free reduction (vectorises); same condition as the loop below
    ok = true
    @inbounds @simd for i in eachindex(before, after)
        a = after[i]
        ok &= (a <= before[i]) & !(a < 0.0)
    end
    ok && return nothing

    # Slow path, only on a violation: find the first offending cell and report it
    @inbounds for i in eachindex(before, after)
        if !(after[i] <= before[i]) || after[i] < 0.0
            throw(
                AssertionError(
                    "Submodel increased (or invalidated) coral cover at timestep $(tstep) (cell " *
                    "$(Tuple(CartesianIndices(before)[i]))): $(before[i]) -> $(after[i]). " *
                    "Submodels may only reduce cover."
                )
            )
        end
    end
    return nothing
end

"""
    _manifest_field(name::AbstractString, field::AbstractString)::Union{String,Nothing}

Value of `field` in the entry of package `name` in the active environment's manifest, or
`nothing` when the package or the field is not recorded there or the manifest cannot be
read. Uses only `TOML`, already a dependency.
"""
function _manifest_field(name::AbstractString, field::AbstractString)::Union{String,Nothing}
    try
        manifest = Base.project_file_manifest_path(Base.active_project())
        isnothing(manifest) && return nothing
        entries = get(TOML.parsefile(manifest), "deps", Dict())
        haskey(entries, name) || return nothing
        val = get(first(entries[name]), field, nothing)
        return isnothing(val) ? nothing : String(val)
    catch
        return nothing
    end
end

"""
    _pkg_revision(name::AbstractString)::Union{String,Nothing}

Revision of package `name` that was asked for, as recorded in the active environment's
manifest (`repo-rev`): a commit SHA, or a branch or tag name if that is what was added (so
it does not pin the code by itself; see `_pkg_tree_hash`). `nothing` when it is
not recorded (registered or path packages) or the manifest cannot be read.
"""
_pkg_revision(name::AbstractString)::Union{String,Nothing} =
    _manifest_field(name, "repo-rev")

"""
    _pkg_tree_hash(name::AbstractString)::Union{String,Nothing}

Git tree hash of the installed source of package `name` (`git-tree-sha1` in the active
environment's manifest), which identifies the code exactly; `nothing` when it is not
recorded (path or `dev`ed packages) or the manifest cannot be read.
"""
_pkg_tree_hash(name::AbstractString)::Union{String,Nothing} =
    _manifest_field(name, "git-tree-sha1")

"""
    _model_provenance(model::Submodel)::Dict{String,Any}

Provenance of one attached submodel: its name, its type as `<module>.<TypeName>` (without
type parameters, so it does not depend on what the session imported), the module that
defines the type and that module's package version (`nothing` when the module is not a
package, e.g. `Main`; note that `pkgversion` of a package *extension* module may also be
`nothing`, so a model living in an extension reports no version). `revision` and
`tree_hash` are the manifest `repo-rev` and `git-tree-sha1` of the package the module
belongs to (that of its root module, so a submodule of a package resolves to the
package), if any.
"""
function _model_provenance(model::Submodel)::Dict{String,Any}
    T = typeof(model)
    mod = parentmodule(T)
    ver = pkgversion(mod)
    pkg = string(nameof(Base.moduleroot(mod)))
    return Dict{String,Any}(
        "name" => string(ADRIAinterfaces.submodel_name(model)),
        "type" => string(mod) * "." * string(nameof(T)),
        "module" => string(mod),
        "version" => isnothing(ver) ? nothing : string(ver),
        "revision" => isnothing(ver) ? nothing : _pkg_revision(pkg),
        "tree_hash" => isnothing(ver) ? nothing : _pkg_tree_hash(pkg)
    )
end

"""
    _provenance_record(domain::Domain, scens::DataFrame)::Dict{String,Any}

What produced a run, written to the result store: ADRIA (with a leading `v`, like the
`ADRIA_VERSION` attr; other versions are plain), Julia and submodel interface versions;
`models`, the attached submodels in attach order followed by the built-in COTS model
(`COTSMod`, with its package version, manifest revision and tree hash) when any scenario
has `cots_on = 1`; and `overrides` (always empty: there is no mechanism to override
submodel parameters except the scenario table, which is stored in `inputs`).
Deterministic for a given run apart from the versions.
"""
function _provenance_record(domain::Domain, scens::DataFrame)::Dict{String,Any}
    models = Dict{String,Any}[
        _model_provenance(a.model) for a in attached_submodels(domain)
    ]

    if "cots_on" in names(scens) && any(==(1), scens.cots_on)
        cots_ver = pkgversion(COTSMod)
        push!(
            models,
            Dict{String,Any}(
                "name" => "COTS", "type" => "builtin", "module" => "COTSMod",
                "version" => isnothing(cots_ver) ? nothing : string(cots_ver),
                "revision" => _pkg_revision("COTSMod"),
                "tree_hash" => _pkg_tree_hash("COTSMod")
            )
        )
    end

    return Dict{String,Any}(
        "adria_version" => "v" * string(pkgversion(@__MODULE__)),
        "julia_version" => string(VERSION),
        "interface_version" => string(ADRIAinterfaces.SUBMODEL_INTERFACE_VERSION),
        "models" => models,
        "overrides" => Any[]
    )
end
