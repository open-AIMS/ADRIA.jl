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
    _attach_submodel!(domain::Domain, model::Submodel)::Domain

Interim mechanism (not the final attach API): attach `model` to `domain` after the
already attached submodels. Calls `ADRIAinterfaces.load_inputs(model, dpkg_path)` once.
Attached submodels are copied with the domain's `SpecStore` (`copy`/`deepcopy`).
"""
function _attach_submodel!(domain::Domain, model::Submodel)::Domain
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
        bufs = ADRIAinterfaces.alloc_logs(ADRIAinterfaces.log_specs(a.model, ctx))
        state = ADRIAinterfaces.init_state(a.model, ctx, params, rng)
        push!(slots, _SubmodelSlot(a.model, ctx, state, bufs))
    end
    slots_t = Tuple(slots)
    names = Tuple(k for s in slots_t for k in keys(s.bufs))
    length(unique(names)) == length(names) ||
        error("Submodels declare logs with the same name: $names")
    return slots_t
end

_aux_logs(slots::Tuple)::NamedTuple = merge((s.bufs for s in slots)..., (;))

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
    _pkg_revision(name::AbstractString)::Union{String,Nothing}

Git revision of package `name` as recorded in the active environment's manifest
(`repo-rev`), or `nothing` when it is not recorded there (registered or path packages) or
the manifest cannot be read. Uses only `TOML`, already a dependency.
"""
function _pkg_revision(name::AbstractString)::Union{String,Nothing}
    try
        manifest = Base.project_file_manifest_path(Base.active_project())
        isnothing(manifest) && return nothing
        entries = get(TOML.parsefile(manifest), "deps", Dict())
        haskey(entries, name) || return nothing
        rev = get(first(entries[name]), "repo-rev", nothing)
        return isnothing(rev) ? nothing : String(rev)
    catch
        return nothing
    end
end

"""
    _model_provenance(model::Submodel)::Dict{String,Any}

Provenance of one attached submodel: its name, its type as `<module>.<TypeName>` (without
type parameters, so it does not depend on what the session imported), the module that
defines the type and that module's package version (`nothing` when the module is not a
package, e.g. `Main`; note that `pkgversion` of a package *extension* module may also be
`nothing`, so a model living in an extension reports no version). `revision` is the
manifest `repo-rev` of the package of that name, if any.
"""
function _model_provenance(model::Submodel)::Dict{String,Any}
    T = typeof(model)
    mod = parentmodule(T)
    ver = pkgversion(mod)
    return Dict{String,Any}(
        "name" => string(ADRIAinterfaces.submodel_name(model)),
        "type" => string(mod) * "." * string(nameof(T)),
        "module" => string(mod),
        "version" => isnothing(ver) ? nothing : string(ver),
        "revision" => isnothing(ver) ? nothing : _pkg_revision(string(nameof(mod)))
    )
end

"""
    _provenance_record(domain::Domain, scens::DataFrame)::Dict{String,Any}

What produced a run, written to the result store: ADRIA (with a leading `v`, like the
`ADRIA_VERSION` attr; other versions are plain), Julia and submodel interface versions;
`models`, the attached submodels in attach order followed by the built-in COTS model
(`COTSMod`, with its package version and manifest revision) when any scenario has
`cots_on = 1`; and `overrides` (always empty: there is no mechanism to override submodel
parameters except the scenario table, which is stored in `inputs`). Deterministic for a
given run apart from the versions.
"""
function _provenance_record(domain::Domain, scens::DataFrame)::Dict{String,Any}
    models = Dict{String,Any}[
        _model_provenance(a.model) for a in attached_submodels(domain)
    ]
    if "cots_on" in names(scens) && any(==(1), scens.cots_on)
        push!(
            models,
            Dict{String,Any}(
                "name" => "COTS", "type" => "builtin", "module" => "COTSMod",
                "version" => string(pkgversion(COTSMod)),
                "revision" => _pkg_revision("COTSMod")
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
