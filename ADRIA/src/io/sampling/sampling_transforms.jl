# sampling_transforms.jl
#
# Post-sampling zeroing and reparameterization transforms, applied to already-
# drawn samples rather than to the pre-sampling spec (contrast with
# sampling_dependencies.jl, which fixes spec columns before sampling). This is
# needed for dependencies that can only be evaluated per-sample — e.g.
# gating on a sampled parameter's value (`fogging > 0`) rather than a
# top-level setting known in advance — and for reparameterizations
# (gamma-to-Dirichlet, mcda_normalize) that only make sense on drawn values.
#
# Byte-for-byte fidelity to what `adjust_samples` previously did is maintained,
# with one exception: seed_strategy/fog_strategy/mc_strategy are EXCLUDED from
# the :seed_group/:fog_group/:mc_group substring matches (commit d5871840,
# 2025-12-15, fixed a regression where including them here conflicted with
# :strategy_group in sampling_dependencies.jl, which already owns them with
# fix_to=-1.0 for CF).
#
# Reuses _CONDITIONS from sampling_dependencies.jl via broadcast; no re-import
# needed since sampling_dependencies.jl is included first.
#
# ---------------------------------------------------------------------------
# How to add a rule
#
# There are four mechanisms in this file. Pick the first one that fits,
# ordered from simplest to most involved.
#
# 1. Single parent, single child/group -> _TRANSFORM_DEPENDENCIES
#    Use when a child is zeroed based on ONE other sampled column, evaluated
#    row-wise (e.g. "zero fog_group when fogging == 0"). Add a row:
#      (parent=:col, child=:target, op=:eq, value=0.0, negate=true, fix_to=0.0)
#    `child` can be a bare spec column (fixed directly) or a group name
#    resolved via `_transform_group_columns` (prefix match / extras, see
#    _TRANSFORM_GROUP_* below). See the negate-semantics note further down
#    for how `op`+`negate` combine into the "active" check.
#
# 2. Multiple parents combined -> _TRANSFORM_GATE_DEFS + _TRANSFORM_GATE_RULES
#    Use when "active" depends on MORE than one column at once (e.g. "any of
#    the five N_seed_* columns is > 0"). Two steps:
#      a. Add a combiner function to _TRANSFORM_GATE_COMBINERS if none of the
#         existing ones (`:any_gt0`, `:any_reactive`) fit. It takes the
#         parents' columns as a sub-DataFrame and returns a BitVector/Vector{Bool}
#         of per-row "active".
#      b. Add a row to _TRANSFORM_GATE_DEFS naming the parents and combiner,
#         then a row to _TRANSFORM_GATE_RULES pointing a child at that gate
#         name, with a fix_to. Gate masks are computed once per call in
#         `_apply_transform_dependencies!` and reused by `_apply_transform_calls!`
#         (see `masks` dict); don't recompute a gate you already defined here.
#
# 3. A child needs new group membership -> _TRANSFORM_GROUP_PREFIXES /
#    _TRANSFORM_GROUP_EXTRA_MEMBERS / _TRANSFORM_GROUP_EXCLUDED_MEMBERS
#    Only needed if your rule's `child` (in 1 or 2 above) is NOT a single spec
#    column but should expand to several. Membership resolves as:
#    (columns matching the group's prefix, if any) + (explicit extras, if
#    any) - (explicit exclusions, if any). You rarely need all three; most
#    groups use only one. See `_transform_group_columns` docstring.
#
# 4. A reparameterization/normalization call, gated by activity ->
#    _TRANSFORM_CALLS
#    Use when the rule isn't zeroing but instead re-deriving values in place
#    (currently only `mcda_normalize`). Add a row with either a direct
#    `parent`/`op`/`value`/`negate` condition, OR `gate=:some_gate_name` to
#    reuse a mask already computed in step 2 above (mutually exclusive: set
#    exactly one of {parent-based condition, gate}). `columns` must be a key
#    into _CRITERIA_WEIGHT_GROUPS (add one if the target isn't already there).
#
# Before choosing 1-4: if the rule depends ONLY on `guided` and is knowable
# BEFORE sampling (not on a value that only exists after a row is drawn), it
# belongs in sampling_dependencies.jl instead, see that file's own "How to
# add a rule" note. Everything in this file fires strictly AFTER sampling,
# because the condition depends on the sampled value of another column
# (`fogging > 0`) rather than a top-level setting fixed for the whole call.
#
# Whichever mechanism you use, double check `_apply_transforms!`'s docstring
# below: the ordering between its numbered steps is load-bearing (e.g.
# gamma_to_dirichlet must run before the zeroing rules that depend on it).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Broadcast gate combiners and gate defs (post-sampling, vector masks)
# ---------------------------------------------------------------------------

const _TRANSFORM_GATE_COMBINERS = Dict{Symbol,Function}(
    :any_gt0 => cols -> vec(any(Matrix(cols) .> 0; dims=2)),
    :any_reactive =>
        cols -> reduce(
            .|, (is_reactive(cols[:, c]) for c in propertynames(cols))
        )
)

const _TRANSFORM_GATE_DEFS = [
    (name=:any_seeding,
        parents=[:N_seed_TA, :N_seed_CA, :N_seed_CNA, :N_seed_SM, :N_seed_LM],
        combine=:any_gt0),
    (name=:any_reactive, parents=[:seed_strategy, :fog_strategy, :mc_strategy],
        combine=:any_reactive)
]

# ---------------------------------------------------------------------------
# Single-parent deferred zeroing rules
#
# negate semantics follow the same convention as _PARAM_DEPENDENCIES: for each
# row, `active = op(parent, value)`, then, if `negate=true`, inverted to
# `active = !active`. `negate` exists because a rule's natural condition is
# sometimes stated as the INACTIVE case rather than the active one, and
# inverting it in-place is simpler than picking a different `op`. The child is
# then zeroed on rows where `!active` (fix_mask = .!active), regardless of
# whether that flip came from `op` or from `negate`.
#
# e.g. if fogging is disabled for a sample, every fog-related parameter for
# that sample is meaningless and should be zeroed too. :fog_group encodes
# this with op=:eq, value=0.0, negate=true:
#   active = .!(fogging .== 0) = (fogging .!= 0)
#   fix_mask = .!active = (fogging .== 0)
# i.e. fog_group is zeroed on rows where fogging == 0.
# ---------------------------------------------------------------------------

const _TRANSFORM_DEPENDENCIES = [
    (parent=:fogging, child=:fog_group, op=:eq, value=0.0, negate=true, fix_to=0.0),
    (parent=:N_mc_settlers, child=:mc_group, op=:eq, value=0.0, negate=true, fix_to=0.0),
    (parent=:SRM, child=:shade_group, op=:eq, value=0.0, negate=true, fix_to=0.0),
    (
        parent=:seed_strategy,
        child=:seed_deployment_freq,
        op=:is_reactive,
        value=nothing,
        negate=true,
        fix_to=0.0
    ),
    (
        parent=:mc_strategy,
        child=:mc_deployment_freq,
        op=:is_reactive,
        value=nothing,
        negate=true,
        fix_to=0.0
    ),
    (
        parent=:fog_strategy,
        child=:fog_deployment_freq,
        op=:is_periodic,
        value=nothing,
        negate=false,
        fix_to=0.0
    )
]

# ---------------------------------------------------------------------------
# Gate rules for compound-parent rows
# ---------------------------------------------------------------------------

const _TRANSFORM_GATE_RULES = [
    (child=:seed_group, gate=:any_seeding, fix_to=0.0),
    (child=:a_adapt, gate=:any_seeding, fix_to=0.0),
    (child=:reactive_group, gate=:any_reactive, fix_to=0.0)
]

# ---------------------------------------------------------------------------
# Live substring-based group membership for post-sampling groups
#
# seed_strategy/fog_strategy/mc_strategy are EXCLUDED from
# :seed_group/:fog_group/:mc_group. :strategy_group (pre-sampling,
# sampling_dependencies.jl) already owns these with fix_to=-1.0 for CF;
# re-touching them here with fix_to=0.0 would reproduce the d5871840
# regression (see file header).
#
# Possible improvement: these groups are currently maintained by hand
# (prefix/extras/exclusions) rather than derived from the model spec's own
# component definitions (e.g. SeedCriteriaWeights, FogCriteriaWeights,
# component_params groupings used elsewhere in this file/_apply_transforms!).
# If the component definitions already encode which fields belong together,
# deriving _TRANSFORM_GROUP_* from them instead of hand-maintaining prefixes
# would remove a source of drift between the two representations.
# ---------------------------------------------------------------------------

const _TRANSFORM_GROUP_PREFIXES = Dict{Symbol,String}(
    :seed_group => "seed_",
    :fog_group => "fog_",
    :mc_group => "mc_",
    :shade_group => "shade_"
)

const _TRANSFORM_GROUP_EXTRA_MEMBERS = Dict{Symbol,Vector{Symbol}}(
    :reactive_group => [
        :reactive_absolute_threshold,
        :reactive_loss_threshold,
        :reactive_min_cover_remaining,
        :reactive_response_delay,
        :reactive_cooldown_period
    ]
)

const _TRANSFORM_GROUP_EXCLUDED_MEMBERS = Dict{Symbol,Vector{Symbol}}(
    :seed_group => [:seed_strategy],
    :fog_group => [:fog_strategy],
    :mc_group => [:mc_strategy]
)

"""
    _transform_group_columns(samples::DataFrame, group::Symbol)::Vector{Symbol}

Return the column names in `samples` that belong to `group` for post-sampling
zeroing. Membership is derived by substring match plus explicit extras, minus
the exclusions in `_TRANSFORM_GROUP_EXCLUDED_MEMBERS`.
"""
function _transform_group_columns(samples::DataFrame, group::Symbol)::Vector{Symbol}
    cols = Symbol[]
    if haskey(_TRANSFORM_GROUP_PREFIXES, group)
        prefix = _TRANSFORM_GROUP_PREFIXES[group]
        append!(cols, filter(n -> contains(String(n), prefix), propertynames(samples)))
    end
    if haskey(_TRANSFORM_GROUP_EXTRA_MEMBERS, group)
        append!(cols, _TRANSFORM_GROUP_EXTRA_MEMBERS[group])
    end
    excluded = get(_TRANSFORM_GROUP_EXCLUDED_MEMBERS, group, Symbol[])
    cols = setdiff(cols, excluded)
    return unique(cols)
end

# ---------------------------------------------------------------------------
# Post-sampling application of the pre-sampling guided dependency table
#
# `_resolve_conditional_spec!` (sampling_dependencies.jl) can only fix
# `_PARAM_DEPENDENCIES` columns pre-sampling when `guided` is a single known
# constant for the whole call (as in sample_cf/sample_unguided/sample_guided/
# sample_selection). Plain `sample(dom, n)` samples `guided` as its own
# dimension, so each row can realize a different regime and the pre-sampling
# resolver structurally cannot apply. `_apply_guided_dependencies!` re-applies
# the same table row-wise, after the fact, using each row's realized `guided`
# value. For calls where `guided` was already fixed pre-sampling, this is a
# no-op (every row is already at its `fix_to` value).
#
# Rule order matters for :strategy_group, which has two rows targeting the
# same child (CF -> -1.0, then guided<=0 -> periodic). The second row's mask
# (guided<=0) is a superset of the first's (guided==-1), so once a row is
# fixed by an earlier rule for a given child, later rules for that same child
# must not re-fix it — mirrored here via `fixed_rows`, analogous to the
# `haskey(known, rule.child)` guard in `_resolve_conditional_spec!`.
# ---------------------------------------------------------------------------

"""
    _apply_guided_dependencies!(samples::DataFrame)::Nothing

Zero/sentinel `_PARAM_DEPENDENCIES` columns row-wise based on each row's
realized `guided` value. No-op if `samples` has no `:guided` column (e.g.
component-scoped sampling that excludes Intervention).
"""
function _apply_guided_dependencies!(samples::DataFrame)::Nothing
    :guided in propertynames(samples) || return nothing

    guided = samples.guided
    fixed_rows = Dict{Symbol,BitVector}()
    for rule in _PARAM_DEPENDENCIES
        already = get!(fixed_rows, rule.child, falses(length(guided)))
        active = _CONDITIONS[rule.op].(guided, rule.value)
        rule.negate && (active = .!active)
        fix_mask = .!active .& .!already
        fixed_rows[rule.child] = already .| fix_mask
        any(fix_mask) || continue

        cols = get(_GROUP_MEMBERS, rule.child, [rule.child])
        cols = filter(c -> c in propertynames(samples), cols)
        isempty(cols) && continue
        samples[fix_mask, cols] .= rule.fix_to
    end
    return nothing
end

# ---------------------------------------------------------------------------
# mcda_normalize gating rules
#
# fog_weights deliberately recomputes fogging > 0 directly rather than
# reusing the :fog_group fix_mask, preserving an asymmetry present in the
# original implementation.
# ---------------------------------------------------------------------------

const _CRITERIA_WEIGHT_GROUPS = Dict{Symbol,Vector{Symbol}}(
    :seed_weights => [
        :seed_heat_stress, :seed_wave_stress, :seed_in_connectivity,
        :seed_out_connectivity, :seed_depth, :seed_coral_cover,
        :seed_cluster_diversity, :seed_geographic_separation
    ],
    :fog_weights => [
        :fog_heat_stress, :fog_wave_stress, :fog_in_connectivity,
        :fog_out_connectivity, :fog_depth, :fog_coral_cover,
        :fog_cluster_diversity, :fog_geographic_separation
    ],
    :mc_weights => [
        :mc_heat_stress, :mc_wave_stress, :mc_in_connectivity,
        :mc_out_connectivity, :mc_depth, :mc_coral_cover,
        :mc_cluster_diversity, :mc_geographic_separation
    ]
)

const _TRANSFORM_CALLS = [
    # fog_weights: deliberately NOT gate-based (see section header above).
    (transform=:mcda_normalize, columns=:fog_weights,
        parent=:fogging, op=:gt, value=0.0, negate=false, gate=nothing),
    (transform=:mcda_normalize, columns=:seed_weights,
        parent=nothing, op=nothing, value=nothing, negate=false, gate=:any_seeding),
    (transform=:mcda_normalize, columns=:mc_weights,
        parent=:N_mc_settlers, op=:eq, value=0.0, negate=true, gate=nothing)
]

# ---------------------------------------------------------------------------
# _apply_transform_dependencies!
# ---------------------------------------------------------------------------

"""
    _apply_transform_dependencies!(samples::DataFrame)::Dict{Symbol,BitVector}

Apply single-parent and gate-based post-sampling zeroing rules.
Returns a dict of gate-name => fix_mask, reused by `_apply_transform_calls!`.
"""
function _apply_transform_dependencies!(samples::DataFrame)::Dict{Symbol,BitVector}
    masks = Dict{Symbol,BitVector}()
    sample_cols = propertynames(samples)

    for rule in _TRANSFORM_DEPENDENCIES
        rule.parent in sample_cols || continue
        active = _CONDITIONS[rule.op].(samples[:, rule.parent], rule.value)
        rule.negate && (active = .!active)
        fix_mask = .!active
        any(fix_mask) || continue
        cols =
            rule.child in sample_cols ?
            [rule.child] :
            _transform_group_columns(samples, rule.child)
        isempty(cols) && continue
        samples[fix_mask, cols] .= rule.fix_to
    end

    for rule in _TRANSFORM_GATE_RULES
        gdef = only(filter(g -> g.name == rule.gate, _TRANSFORM_GATE_DEFS))
        all(p -> p in sample_cols, gdef.parents) || continue
        active = _TRANSFORM_GATE_COMBINERS[gdef.combine](samples[:, gdef.parents])
        fix_mask = .!active
        masks[rule.gate] = fix_mask
        any(fix_mask) || continue
        cols =
            rule.child in sample_cols ?
            [rule.child] :
            _transform_group_columns(samples, rule.child)
        isempty(cols) && continue
        samples[fix_mask, cols] .= rule.fix_to
    end

    return masks
end

# ---------------------------------------------------------------------------
# _apply_transform_calls!
# ---------------------------------------------------------------------------

"""
    _apply_transform_calls!(samples::DataFrame, masks::Dict{Symbol,BitVector})::Nothing

Apply mcda_normalize to criteria weight columns, gated per-intervention-type.
`masks` is the gate-name => fix_mask dict returned by `_apply_transform_dependencies!`.
"""
function _apply_transform_calls!(
    samples::DataFrame, masks::Dict{Symbol,BitVector}
)::Nothing
    :guided in propertynames(samples) || return nothing
    guided_active = samples.guided .> 0
    sample_cols = propertynames(samples)
    for rule in _TRANSFORM_CALLS
        rule.gate === nothing && !(rule.parent in sample_cols) && continue
        cols = _CRITERIA_WEIGHT_GROUPS[rule.columns]
        all(c -> c in sample_cols, cols) || continue
        parent_active = if rule.gate !== nothing
            .!masks[rule.gate]
        else
            active = _CONDITIONS[rule.op].(samples[:, rule.parent], rule.value)
            rule.negate ? .!active : active
        end
        mask = parent_active .& guided_active
        any(mask) || continue
        samples[mask, cols] .= mcda_normalize(samples[mask, cols])
    end
    return nothing
end

# ---------------------------------------------------------------------------
# _apply_transforms! — assembled post-sampling transform function
# ---------------------------------------------------------------------------

"""
    _apply_transforms!(spec::DataFrame, samples::DataFrame)::DataFrame

Post-sampling transforms: reparameterization + conditional zeroing.
Replaces the body of `adjust_samples`.

Ordering is load-bearing:
1. `_apply_guided_dependencies!` — row-wise re-application of the pre-sampling
   guided dependency table, needed for callers (plain `sample(dom, n)`) where
   `guided` was sampled rather than fixed. Runs first so later steps see
   already-zeroed intervention/criteria/strategy columns for inactive rows.
2. `floor.(a_adapt_ref)` — no ordering dependency vs step 1 (zeroed rows floor
   to 0.0 either way).
3. `gamma_to_dirichlet` on weight columns — must precede the zeroing rules
   below, which depend on the reparameterized values.
4. `_apply_transform_dependencies!` — single-parent + gate zeroing rules.
5. `_apply_transform_calls!` — mcda_normalize gates.
6. `seed_wave_stress` zeroing — must run LAST: mcda_normalize normalizes seed
   weights (including seed_wave_stress) to sum=1 first; zeroing
   seed_wave_stress afterward deliberately leaves the remaining 7 weights
   summing to <1 for wave_scenario==0 rows, matching production behaviour.
"""
function _apply_transforms!(spec::DataFrame, samples::DataFrame)::DataFrame
    _apply_guided_dependencies!(samples)

    sample_cols = propertynames(samples)

    if :a_adapt_ref in sample_cols
        samples[:, :a_adapt_ref] .= floor.(samples[:, :a_adapt_ref])
    end

    # Must precede the zeroing rules below (see docstring above, point 2).
    if :guided in sample_cols
        guided_mask = samples.guided .> 0
        if any(guided_mask)
            seed_weights = component_params(spec, SeedCriteriaWeights)
            fog_weights = component_params(spec, FogCriteriaWeights)
            mc_weights = component_params(spec, MCCriteriaWeights)
            for wf in (seed_weights.fieldname, fog_weights.fieldname, mc_weights.fieldname)
                wf = filter(c -> c in sample_cols, wf)
                isempty(wf) && continue
                samples[guided_mask, wf] .= gamma_to_dirichlet(
                    Matrix(samples[guided_mask, wf])
                )
            end
        end
    end

    masks = _apply_transform_dependencies!(samples)
    _apply_transform_calls!(samples, masks)

    # Must run LAST (see docstring above, point 5). seed_wave_stress is a
    # SeedCriteriaWeights member already normalised by mcda_normalize above;
    # zeroing it here deliberately does NOT renormalise the remaining 7 seed
    # weights back to sum=1, matching production behaviour.
    if :wave_scenario in sample_cols && :seed_wave_stress in sample_cols
        no_wave = _CONDITIONS[:eq].(samples.wave_scenario, 0.0)
        samples[no_wave, :seed_wave_stress] .= 0.0
    end

    if size(unique(Matrix(samples); dims=1), 1) < nrow(samples)
        perc = "$(@sprintf("%.3f", (1.0 - (nrow(unique(samples)) / nrow(samples))) * 100.0))%"
        @warn "Non-unique samples created: $perc of the samples are duplicates."
    end

    return samples
end
