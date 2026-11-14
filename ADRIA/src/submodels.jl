# Plumbing that lets the run loop call the hooks of ADRIAinterfaces on the submodels attached
# to a domain. Interim: there is no public attach API yet (see `_attach_submodel!`), and
# COTS does not go through this path. With no submodel attached, none of this runs and the
# run is unchanged (no extra RNG draws, no extra allocations).

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

# Phase hook position: after core disturbances (bleaching, cyclones), before survival is
# propagated. Submodels step in attach order, then record.
@inline function _step_submodels!(
    slots::Tuple, phase::ADRIAinterfaces.Phase, view::CoralStepView, tstep::Int, year::Int
)::Nothing
    foreach(slots) do s
        if ADRIAinterfaces.has_phase(s.model, phase)
            ADRIAinterfaces.step!(s.model, phase, s.state, s.ctx, view, tstep, year)
        end
    end
    return nothing
end

@inline function _record_submodels!(slots::Tuple, tstep::Int)::Nothing
    foreach(s -> ADRIAinterfaces.record!(s.model, s.bufs, s.state, tstep), slots)
    return nothing
end
