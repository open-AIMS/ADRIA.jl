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
# with one exception: iv_CAq_strategy/iv_Fog_strategy/iv_LvM_strategy are EXCLUDED from
# the :CAq_group/:Fog_group/:LvM_group substring matches (commit d5871840,
# 2025-12-15, fixed a regression where including them here conflicted with
# :strategy_group in sampling_dependencies.jl, which already owns them with
# fix_to=-1.0 for CF). Fog_strategy/LvM_strategy are instead gated on their own
# intervention's activity via the dedicated :Fog_strategy_gate/:LvM_strategy_gate
# rules below (fix_to=periodic when fogging/iv_LvM_N_settlers == 0), which check
# the CF sentinel before touching the column so they don't reintroduce that
# regression. iv_CAq_strategy is deliberately left out of this: it isn't gated
# on :any_CAq here, so a coral-aquaculture-inactive scenario can still draw
# iv_CAq_strategy=reactive independently.
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
#    row-wise (e.g. "zero Fog_group when fogging == 0"). Add a row:
#      (parent=:col, child=:target, op=:eq, value=0.0, negate=true, fix_to=0.0)
#    `child` can be a bare spec column (fixed directly) or a group name
#    resolved via `_transform_group_columns` (prefix match / extras, see
#    _TRANSFORM_GROUP_* below). See the negate-semantics note further down
#    for how `op`+`negate` combine into the "active" check.
#
# 2. Multiple parents combined -> _TRANSFORM_GATE_DEFS + _TRANSFORM_GATE_RULES
#    Use when "active" depends on MORE than one column at once (e.g. "any of
#    the five iv_CAq_N_* columns is > 0"). Two steps:
#      a. Add a combiner function to _TRANSFORM_GATE_COMBINERS if none of the
#         existing ones (`:any_gt0`, `:any_reactive`, `:strategy_gate`) fit. It takes the
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
        ),
    # cols[:,1] is the intervention-activity column (e.g. fogging,
    # iv_LvM_N_settlers), cols[:,2] is the *_strategy column itself. "active"
    # (i.e. leave alone) when the intervention is on, OR when the strategy
    # column already holds the CF sentinel (-1.0, set pre-sampling by
    # :strategy_group / re-applied by _apply_guided_dependencies!); this
    # gate must never overwrite that sentinel.
    :strategy_gate => cols -> (cols[:, 1] .!= 0) .| (cols[:, 2] .== -1.0)
)

const _TRANSFORM_GATE_DEFS = [
    (name=:any_CAq,
        parents=[:iv_CAq_N_TA, :iv_CAq_N_CA, :iv_CAq_N_CNA, :iv_CAq_N_SM, :iv_CAq_N_LM],
        combine=:any_gt0),
    (name=:any_reactive, parents=[:iv_CAq_strategy, :iv_Fog_strategy, :iv_LvM_strategy],
        combine=:any_reactive),
    (name=:Fog_strategy_gate, parents=[:iv_Fog, :iv_Fog_strategy], combine=:strategy_gate),
    (name=:LvM_strategy_gate, parents=[:iv_LvM_N_settlers, :iv_LvM_strategy],
        combine=:strategy_gate)
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
# that sample is meaningless and should be zeroed too. :Fog_group encodes
# this with op=:eq, value=0.0, negate=true:
#   active = .!(fogging .== 0) = (fogging .!= 0)
#   fix_mask = .!active = (fogging .== 0)
# i.e. Fog_group is zeroed on rows where fogging == 0.
# ---------------------------------------------------------------------------

const _TRANSFORM_DEPENDENCIES = [
    (parent=:iv_Fog, child=:Fog_group, op=:eq, value=0.0, negate=true, fix_to=0.0),
    (
        parent=:iv_LvM_N_settlers,
        child=:LvM_group,
        op=:eq,
        value=0.0,
        negate=true,
        fix_to=0.0
    ),
    (parent=:iv_Shd, child=:Shd_group, op=:eq, value=0.0, negate=true, fix_to=0.0),
    (
        parent=:iv_CAq_strategy,
        child=:iv_CAq_deployment_freq,
        op=:is_reactive,
        value=nothing,
        negate=true,
        fix_to=0.0
    ),
    (
        parent=:iv_LvM_strategy,
        child=:iv_LvM_deployment_freq,
        op=:is_reactive,
        value=nothing,
        negate=true,
        fix_to=0.0
    ),
    (
        parent=:iv_Fog_strategy,
        child=:iv_Fog_deployment_freq,
        op=:is_periodic,
        value=nothing,
        negate=false,
        fix_to=0.0
    ),
    # iv_CAq_a_adapt_ref (adaptation reference period) only matters when assisted adaptation is
    # actually applied -- zero it whenever iv_CAq_a_adapt itself is 0, whether that's a direct
    # sampled draw or an explicit override (e.g. a sensitivity analysis forcing it to 0).
    (
        parent=:iv_CAq_a_adapt,
        child=:iv_CAq_a_adapt_ref,
        op=:eq,
        value=0.0,
        negate=true,
        fix_to=0.0
    )
]

# ---------------------------------------------------------------------------
# Gate rules for compound-parent rows
# ---------------------------------------------------------------------------

const _TRANSFORM_GATE_RULES = [
    (child=:CAq_group, gate=:any_CAq, fix_to=0.0),
    # :a_adapt_group (not just :a_adapt) so a_adapt_ref is also zeroed here -- the
    # single-parent rule above runs BEFORE this gate loop, so it can't see a_adapt
    # get zeroed by the no-coral-aquaculture gate itself; this covers that ordering gap.
    (child=:a_adapt_group, gate=:any_CAq, fix_to=0.0),
    # Must run BEFORE :reactive_group's :any_reactive rule below: that rule
    # reads :Fog_strategy/:LvM_strategy, and needs to see them already fixed
    # to periodic here for fogging/iv_LvM_N_settlers-inactive rows, otherwise an
    # orphaned reactive draw on an inactive intervention keeps reactive_group
    # alive as pure sampling noise (see file header).
    (child=:iv_Fog_strategy, gate=:Fog_strategy_gate, fix_to=DECISION_STRATEGY[:periodic]),
    (child=:iv_LvM_strategy, gate=:LvM_strategy_gate, fix_to=DECISION_STRATEGY[:periodic]),
    (child=:reactive_group, gate=:any_reactive, fix_to=0.0)
]

# ---------------------------------------------------------------------------
# Live substring-based group membership for post-sampling groups
#
# iv_CAq_strategy/iv_Fog_strategy/iv_LvM_strategy are EXCLUDED from
# :CAq_group/:Fog_group/:LvM_group. :strategy_group (pre-sampling,
# sampling_dependencies.jl) already owns these with fix_to=-1.0 for CF;
# re-touching them here with fix_to=0.0 would reproduce the d5871840
# regression (see file header).
#
# Behaviour note (intervention-parameter rename, 2026): iv_CAq_devices_per_m2 (formerly
# CAq_devices_per_m2) now starts with the :CAq_group prefix "iv_CAq_", so it is newly
# swept into :CAq_group's post-sampling zeroing (fixed to 0 whenever :any_CAq is
# false). Previously "seed_" did not match "CAq_devices_per_m2" as a substring (the
# character after "seed" was "i", not "_"), so this column was NOT zeroed when coral aquaculture
# was inactive. This is a disclosed, intentional side effect of the rename, not a
# regression: the column is genuinely meaningless without coral aquaculture, so zeroing it when
# inactive is arguably a latent-bug fix rather than a behaviour change worth guarding
# against.
#
# Possible improvement: these groups are currently maintained by hand
# (prefix/extras/exclusions) rather than derived from the model spec's own
# component definitions (e.g. CAqCriteriaWeights, FogCriteriaWeights,
# component_params groupings used elsewhere in this file/_apply_transforms!).
# If the component definitions already encode which fields belong together,
# deriving _TRANSFORM_GROUP_* from them instead of hand-maintaining prefixes
# would remove a source of drift between the two representations.
# ---------------------------------------------------------------------------

const _TRANSFORM_GROUP_PREFIXES = Dict{Symbol,String}(
    :CAq_group => "iv_CAq_",
    :Fog_group => "iv_Fog_",
    :LvM_group => "iv_LvM_",
    :Shd_group => "iv_Shd_"
)

const _TRANSFORM_GROUP_EXTRA_MEMBERS = Dict{Symbol,Vector{Symbol}}(
    :reactive_group => [
        :reactive_absolute_threshold,
        :reactive_loss_threshold,
        :reactive_min_cover_remaining,
        :reactive_response_delay
    ],
    :a_adapt_group => [:iv_CAq_a_adapt, :iv_CAq_a_adapt_ref]
)

const _TRANSFORM_GROUP_EXCLUDED_MEMBERS = Dict{Symbol,Vector{Symbol}}(
    :CAq_group => [:iv_CAq_strategy],
    :Fog_group => [:iv_Fog_strategy],
    :LvM_group => [:iv_LvM_strategy]
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
# Fog_weights deliberately recomputes fogging > 0 directly rather than
# reusing the :Fog_group fix_mask, preserving an asymmetry present in the
# original implementation.
# ---------------------------------------------------------------------------

const _CRITERIA_WEIGHT_GROUPS = Dict{Symbol,Vector{Symbol}}(
    :CAq_weights => [
        :iv_CAq_heat_stress, :iv_CAq_wave_stress, :iv_CAq_in_connectivity,
        :iv_CAq_out_connectivity, :iv_CAq_depth, :iv_CAq_coral_cover,
        :iv_CAq_cluster_diversity, :iv_CAq_geographic_separation
    ],
    :Fog_weights => [
        :iv_Fog_heat_stress, :iv_Fog_wave_stress, :iv_Fog_in_connectivity,
        :iv_Fog_out_connectivity, :iv_Fog_depth, :iv_Fog_coral_cover,
        :iv_Fog_cluster_diversity, :iv_Fog_geographic_separation
    ],
    :LvM_weights => [
        :iv_LvM_heat_stress, :iv_LvM_wave_stress, :iv_LvM_in_connectivity,
        :iv_LvM_out_connectivity, :iv_LvM_depth, :iv_LvM_coral_cover,
        :iv_LvM_cluster_diversity, :iv_LvM_geographic_separation
    ]
)

const _TRANSFORM_CALLS = [
    # Fog_weights: deliberately NOT gate-based (see section header above).
    (transform=:mcda_normalize, columns=:Fog_weights,
        parent=:iv_Fog, op=:gt, value=0.0, negate=false, gate=nothing),
    (transform=:mcda_normalize, columns=:CAq_weights,
        parent=nothing, op=nothing, value=nothing, negate=false, gate=:any_CAq),
    (transform=:mcda_normalize, columns=:LvM_weights,
        parent=:iv_LvM_N_settlers, op=:eq, value=0.0, negate=true, gate=nothing)
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
5. Floor each `*_revisit_cadence` column at its paired `*_deployment_freq`
   column (`max(cadence, freq)`) — must run AFTER step 4, since step 4 is what
   settles each `*_deployment_freq` column to either its real sampled value or
   `0.0` (via three different mechanisms across coral aquaculture/fog/larval methods — the
   `any_CAq` gate rule, and single-parent `_TRANSFORM_DEPENDENCIES` rows
   for fog/larval methods — all applied inside step 4). Without this floor, sensitivity
   analyses over the raw cadence factor would have a large inert region
   whenever `deployment_freq` already exceeds it (cadence has zero effect on
   model behaviour in that region).
6. `_apply_transform_calls!` — mcda_normalize gates.
7. `iv_CAq_wave_stress` zeroing — must run LAST: mcda_normalize normalizes coral aquaculture
   weights (including iv_CAq_wave_stress) to sum=1 first; zeroing
   iv_CAq_wave_stress afterward deliberately leaves the remaining 7 weights
   summing to <1 for wave_scenario==0 rows, matching production behaviour.
"""
function _apply_transforms!(spec::DataFrame, samples::DataFrame)::DataFrame
    _apply_guided_dependencies!(samples)

    sample_cols = propertynames(samples)

    if :iv_CAq_a_adapt_ref in sample_cols
        samples[:, :iv_CAq_a_adapt_ref] .= floor.(samples[:, :iv_CAq_a_adapt_ref])
    end

    # Must precede the zeroing rules below (see docstring above, point 2).
    if :guided in sample_cols
        guided_mask = samples.guided .> 0
        if any(guided_mask)
            CAq_weights = component_params(spec, CAqCriteriaWeights)
            Fog_weights = component_params(spec, FogCriteriaWeights)
            LvM_weights = component_params(spec, LvMCriteriaWeights)
            for wf in (CAq_weights.fieldname, Fog_weights.fieldname, LvM_weights.fieldname)
                wf = filter(c -> c in sample_cols, wf)
                isempty(wf) && continue
                samples[guided_mask, wf] .= gamma_to_dirichlet(
                    Matrix(samples[guided_mask, wf])
                )
            end
        end
    end

    masks = _apply_transform_dependencies!(samples)

    # Must run after _apply_transform_dependencies! (see docstring above, point 5).
    for prefix in ("iv_CAq", "iv_Fog", "iv_LvM")
        freq_col = Symbol("$(prefix)_deployment_freq")
        cadence_col = Symbol("$(prefix)_revisit_cadence")
        (freq_col in sample_cols && cadence_col in sample_cols) || continue
        samples[:, cadence_col] .= max.(samples[:, cadence_col], samples[:, freq_col])
    end

    _apply_transform_calls!(samples, masks)

    # Must run LAST (see docstring above, point 5). iv_CAq_wave_stress is a
    # CAqCriteriaWeights member already normalised by mcda_normalize above;
    # zeroing it here deliberately does NOT renormalise the remaining 7 coral aquaculture
    # weights back to sum=1, matching production behaviour.
    if :wave_scenario in sample_cols && :iv_CAq_wave_stress in sample_cols
        no_wave = _CONDITIONS[:eq].(samples.wave_scenario, 0.0)
        samples[no_wave, :iv_CAq_wave_stress] .= 0.0
    end

    if size(unique(Matrix(samples); dims=1), 1) < nrow(samples)
        perc = "$(@sprintf("%.3f", (1.0 - (nrow(unique(samples)) / nrow(samples))) * 100.0))%"
        @warn "Non-unique samples created: $perc of the samples are duplicates."
    end

    return samples
end
