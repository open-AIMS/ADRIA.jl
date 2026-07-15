# sampling_dependencies.jl
#
# Many sampled parameters are only meaningful under certain scenario settings
# (e.g. seeding/fogging deployment parameters only matter when an intervention
# is actually enabled, and decision-strategy weights only matter under
# :guided). Sampling every parameter unconditionally wastes sample budget on
# combinations that don't affect model output and dilutes sensitivity/PAWN
# analyses with dead dimensions.
#
# This file defines a declarative parent→child dependency table
# (_PARAM_DEPENDENCIES) plus resolution logic that, given a known top-level
# setting (e.g. `guided`), fixes every dependent parameter that is rendered
# inactive to a constant value before sampling — rather than sampling it and
# discarding the result. Parameter groups (_GROUP_MEMBERS) let a single rule
# fix many physical spec columns at once.
#
# Only rows whose parent is :guided can fire here, since that is the only
# value known at pre-sampling resolve time. Dependencies on other settings are
# resolved later in sampling_transforms.jl.

using ADRIA.decision: is_reactive, DECISION_STRATEGY

# ---------------------------------------------------------------------------
# Op registry (scalar conditions only; broadcast use is in
# sampling_transforms.jl)
# ---------------------------------------------------------------------------

const _CONDITIONS = Dict{Symbol,Function}(
    :gt => (x, t) -> x > t,
    :lt => (x, t) -> x < t,
    :ge => (x, t) -> x >= t,
    :le => (x, t) -> x <= t,
    :eq => (x, t) -> x == t,
    :neq => (x, t) -> x != t,
    :is_reactive => (x, _) -> is_reactive(x),
    :is_periodic => (x, _) -> x == DECISION_STRATEGY[:periodic]
)

# ---------------------------------------------------------------------------
# Dependency table. Rows with a parent other than :guided are structurally
# unreachable at pre-sampling resolve time; they are handled post-sampling in
# sampling_transforms.jl instead.
#
# Both :strategy_group rows use negate=false. The CF row's fix_to is -1.0,
# not 0.0, to distinguish it from the :intervention_group columns (see
# _GROUP_MEMBERS below).
#
# Each row is a NamedTuple with fields:
#   parent  — spec column whose (known) value the condition is evaluated against
#   child   — spec column, or a key into _GROUP_MEMBERS, to fix when inactive
#   op      — key into _CONDITIONS; op(parent_value, value) determines whether
#             the child is considered "active"
#   value   — threshold/comparison value passed to `op`
#   negate  — if true, invert the boolean result of `op(parent_value, value)`
#             before it is used as the activity check. Needed when the
#             natural predicate is stated as the INACTIVE condition rather
#             than the active one — e.g. :strategy_group's CF row wants
#             "active unless guided == -1", which is `op=:neq` on -1.0
#             directly (negate=false); but if a rule were instead phrased as
#             "inactive when guided == -1", you'd write `op=:eq, value=-1.0,
#             negate=true` to get the same activity result.
#   fix_to  — constant the child is fixed to when NOT active
#
# Rows are grouped by parent (all rows here have parent=:guided — see note
# above on why other parents are structurally unreachable at this stage).
#
# To add a new rule:
#   1. If it depends on `guided` and is knowable pre-sampling, add a row here.
#      Otherwise it belongs in sampling_transforms.jl instead.
#   2. If `child` covers more than one physical spec column, add an entry to
#      _GROUP_MEMBERS mapping the child name to those column names.
#   3. Run `_validate_dependencies` (called from the sampling entry point) to
#      catch unknown parents/ops/children, cycles, and conflicting fix_to
#      values on shared columns.
#
# e.g. if :my_param only matters under a guided strategy, it should be
# fixed to 0.0 for every other regime (CF, unguided). That is:
# Example — fix :my_param to 0.0 whenever guided <= 0 (i.e. not guided):
#   (parent=:guided, child=:my_param, op=:gt, value=0.0, negate=false, fix_to=0.0)
# ---------------------------------------------------------------------------

const _PARAM_DEPENDENCIES = [
    # parent=:guided
    (
        parent=:guided,
        child=:intervention_group,
        op=:neq,
        value=-1.0,
        negate=false,
        fix_to=0.0
    ),
    (parent=:guided, child=:criteria_weights, op=:gt, value=0.0, negate=false, fix_to=0.0),
    # mcda_method is only meaningful in the guided regime (guided > 0); fix_to=0.0
    # is a sentinel meaning "no MCDA method selected", never a valid method index,
    # since decision.mcda_methods() is indexed from 1.
    (parent=:guided, child=:mcda_method, op=:gt, value=0.0, negate=false, fix_to=0.0),
    (parent=:guided, child=:plan_horizon, op=:neq, value=0.0, negate=false, fix_to=0.0),
    # fix_to must match :intervention_group's fix_to=0.0 above, since
    # `projection_confidence` is also a member of that group (see
    # _GROUP_MEMBERS) and _validate_dependencies rejects two child symbols
    # disagreeing on a fix_to for the same physical column. The value is
    # cosmetic either way, not correctness-affecting: `plan_horizon`'s own row
    # above already collapses the decay window to a single element whenever
    # this row fires, and `1^exponent == 1` for any exponent -- so this
    # factor's fixed value has no effect on `build_decay`'s output once
    # `plan_horizon` is inactive.
    (
        parent=:guided,
        child=:projection_confidence,
        op=:neq,
        value=0.0,
        negate=false,
        fix_to=0.0
    ),
    (
        parent=:guided,
        child=:depth_thresholds,
        op=:neq,
        value=-1.0,
        negate=false,
        fix_to=0.0
    ),
    (parent=:guided, child=:strategy_group, op=:neq, value=-1.0, negate=false, fix_to=-1.0),
    (
        parent=:guided,
        child=:strategy_group,
        op=:gt,
        value=0.0,
        negate=false,
        fix_to=DECISION_STRATEGY[:periodic]
    )
]

# ---------------------------------------------------------------------------
# Group membership (explicit Dict; pre-sampling groups only)
#
# :seed_strategy/:fog_strategy/:mc_strategy appear ONLY in :strategy_group,
# not in :intervention_group. :strategy_group fixes them to -1.0 for CF and
# to DECISION_STRATEGY[:periodic] for unguided; including them in
# :intervention_group (fix_to=0.0) would produce a conflicting fix_to for the
# same column in _fix_child!.
#
# Possible improvement: _GROUP_MEMBERS is a hand-maintained list of column
# names per group. If the model spec's own component definitions (the same
# ones used with component_params in sampling_transforms.jl) already encode
# equivalent groupings, deriving these entries from those component
# definitions instead of hand-maintaining them here would remove a source of
# drift between the two representations.
# ---------------------------------------------------------------------------

const _GROUP_MEMBERS = Dict{Symbol,Vector{Symbol}}(
    :intervention_group => [
        :N_seed_TA, :N_seed_CA, :N_seed_CNA, :N_seed_SM, :N_seed_LM,
        :N_mc_settlers, :seeding_devices_per_m2, :min_iv_locations,
        :mc_min_iv_locations, :fogging, :SRM, :a_adapt, :a_adapt_ref,
        :seed_years, :shade_years, :fog_years, :plan_horizon, :projection_confidence,
        :seed_deployment_freq, :fog_deployment_freq, :shade_deployment_freq,
        :mc_deployment_freq, :seed_year_start, :shade_year_start,
        :fog_year_start, :mc_year_start, :mc_years,
        :mcb_albedo, :mcb_duration, :mcb_deployment_freq,
        # NOTE: :seed_strategy/:fog_strategy/:mc_strategy are intentionally
        # excluded here — they belong exclusively to :strategy_group (see above).
        :reactive_absolute_threshold, :reactive_loss_threshold,
        :reactive_min_cover_remaining, :reactive_response_delay,
        :reactive_cooldown_period
    ],
    :criteria_weights => [
        :seed_heat_stress, :seed_wave_stress, :seed_in_connectivity,
        :seed_out_connectivity, :seed_depth, :seed_coral_cover,
        :seed_cluster_diversity, :seed_geographic_separation,
        :fog_heat_stress, :fog_wave_stress, :fog_in_connectivity,
        :fog_out_connectivity, :fog_depth, :fog_coral_cover,
        :fog_cluster_diversity, :fog_geographic_separation,
        :mc_heat_stress, :mc_wave_stress, :mc_in_connectivity,
        :mc_out_connectivity, :mc_depth, :mc_coral_cover,
        :mc_cluster_diversity, :mc_geographic_separation
    ], :depth_thresholds => [:depth_min, :depth_offset],
    :strategy_group => [:seed_strategy, :fog_strategy, :mc_strategy]
)

# ---------------------------------------------------------------------------
# _fix_child!  (no-op if already fixed, to the same value or a caller-chosen one)
# ---------------------------------------------------------------------------

"""
    _fix_child!(spec::DataFrame, child::Symbol, fix_to::Real)::Nothing

Fix a spec column (or all members of a named group) to `fix_to` as a constant.

- No-op if the column is already `is_constant=true` with a value `isequal` to `fix_to`.
- Also a no-op (with a `@debug` note) if the column is already `is_constant=true` with a
  *different* value — `fix_to` here is a dead-dimension sentinel (the column doesn't
  matter under the current regime), and a column the caller already fixed via
  `fix_factor!`/`set_factor_bounds!` has already achieved that "no longer varies" goal,
  just to a real value instead of the sentinel. Genuine rule-vs-rule disagreements
  (two `_PARAM_DEPENDENCIES` rows wanting different fix_to values for the same column)
  are caught statically by `_validate_dependencies`, not here.
- Otherwise sets `is_constant=true`, `val`, `lower_bound`, `upper_bound` to `fix_to`,
  and `dist_params` to a same-length tuple filled with `fix_to`.

Columns absent from `spec` (e.g. domain-filtered fields) are silently skipped.
"""
function _fix_child!(spec::DataFrame, child::Symbol, fix_to::Real)::Nothing
    physical_cols = get(_GROUP_MEMBERS, child, [child])
    for col in physical_cols
        row_mask = spec.fieldname .== col
        any(row_mask) || continue
        is_const = spec[row_mask, :is_constant][1]
        cur_val = spec[row_mask, :val][1]
        if is_const
            if !isequal(cur_val, fix_to)
                @debug "Column :$col is already fixed to $cur_val; leaving as-is " *
                    "instead of overwriting with dependency sentinel $fix_to"
            end
            continue
        end
        existing_dp = spec[row_mask, :dist_params][1]
        new_dp = Tuple(fill(Float64(fix_to), length(existing_dp)))
        spec[row_mask, :val] .= Float64(fix_to)
        spec[row_mask, :lower_bound] .= Float64(fix_to)
        spec[row_mask, :upper_bound] .= Float64(fix_to)
        spec[row_mask, :dist_params] .= [new_dp]
        spec[row_mask, :is_constant] .= true
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Fixed-point resolution loop (gate-free; only parent→child edges)
# ---------------------------------------------------------------------------

"""
    _resolve_conditional_spec!(spec::DataFrame, fixed::NamedTuple)::DataFrame

Pre-sampling conditional spec resolution.

Fix spec columns/groups whose activation condition can be determined from `fixed`
(typically `(guided=value,)` for one of CF/unguided/guided). Runs a fixed-point
loop until no new rules fire. Rules whose parent is not in `fixed` are silently
skipped — they belong to the post-sampling `_apply_transforms!` step.

# Example
Called once per scenario regime, before sampling, with the known top-level
`guided` value; every parameter rendered inactive under that regime is fixed
to a constant so it is never sampled.
```julia
_resolve_conditional_spec!(spec, (guided=-1.0,))  # counterfactual
_resolve_conditional_spec!(spec, (guided=0.0,))   # unguided
_resolve_conditional_spec!(spec, (guided=1.0,))   # guided (top-split only)
```
"""
function _resolve_conditional_spec!(spec::DataFrame, fixed::NamedTuple)::DataFrame
    known = Dict{Symbol,Any}(pairs(fixed))
    changed = true
    while changed
        changed = false
        for rule in _PARAM_DEPENDENCIES
            haskey(known, rule.child) && continue
            haskey(known, rule.parent) || continue
            active = _CONDITIONS[rule.op](known[rule.parent], rule.value)
            rule.negate && (active = !active)
            if !active
                _fix_child!(spec, rule.child, rule.fix_to)
                known[rule.child] = rule.fix_to
                changed = true
            end
        end
    end
    return spec
end

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

"""
    _check_dependency_acyclic(rules)::Nothing

Check that the parent→child edges in `rules` form a DAG. Errors if a cycle is found.
"""
function _check_dependency_acyclic(rules)::Nothing
    adj = Dict{Symbol,Vector{Symbol}}()
    for rule in rules
        push!(get!(adj, rule.parent, Symbol[]), rule.child)
    end

    colors = Dict{Symbol,Int}()  # 0=unvisited, 1=in-stack, 2=done

    function dfs(node)
        get(colors, node, 0) == 2 && return nothing
        if get(colors, node, 0) == 1
            throw(ArgumentError("Cycle detected in _PARAM_DEPENDENCIES involving :$node"))
        end
        colors[node] = 1
        for child in get(adj, node, Symbol[])
            dfs(child)
        end
        colors[node] = 2
    end

    all_nodes = unique(
        vcat(Symbol[r.parent for r in rules], Symbol[r.child for r in rules])
    )
    for node in all_nodes
        get(colors, node, 0) == 0 && dfs(node)
    end
    return nothing
end

"""
    _validate_dependencies(spec::DataFrame)::Nothing

Validate `_PARAM_DEPENDENCIES` and `_GROUP_MEMBERS` against `spec`.

Checks:
- Every rule's parent is a known spec column.
- Every rule uses a registered op from `_CONDITIONS`.
- Every rule's child is either a spec column or a defined group name.
- The parent→child graph is acyclic.
- No two rules fix the same physical column to *different* values; same-value
  overlaps are permitted (no-op branch).
"""
function _validate_dependencies(spec::DataFrame)::Nothing
    # spec_cols = set of parameter field names (rows of spec), not the DataFrame's column names
    spec_cols = Set(spec.fieldname)

    for rule in _PARAM_DEPENDENCIES
        rule.parent in spec_cols ||
            throw(
                ArgumentError(
                    "Dependency rule for :$(rule.child) references unknown parent " *
                    ":$(rule.parent)"
                )
            )
        haskey(_CONDITIONS, rule.op) ||
            throw(
                ArgumentError(
                    "Dependency rule for :$(rule.child) uses unknown op :$(rule.op)"
                )
            )
        rule.child in spec_cols || rule.child in keys(_GROUP_MEMBERS) ||
            throw(
                ArgumentError(
                    "Dependency rule references unknown child :$(rule.child) " *
                    "(not a spec column or a defined group name)"
                )
            )
    end

    _check_dependency_acyclic(_PARAM_DEPENDENCIES)

    # Duplicate-child detection on expanded physical columns.
    # Rows with the SAME child symbol are guarded by the `haskey(known, rule.child)`
    # check in _resolve_conditional_spec! — only one can fire per call, so different
    # fix_to values are permitted (they represent mutually-exclusive regimes, e.g. the
    # two :strategy_group rows for CF and unguided).
    # Only flag when DIFFERENT child symbols expand to the same physical column with
    # disagreeing fix_to values — those lack the haskey guard and are genuine conflicts.
    column_owners = Dict{Symbol,Vector{Any}}()
    for rule in _PARAM_DEPENDENCIES
        cols = get(_GROUP_MEMBERS, rule.child, [rule.child])
        for c in cols
            push!(get!(column_owners, c, []), (rule.child, rule.fix_to))
        end
    end
    for (col, owners) in column_owners
        length(owners) <= 1 && continue
        # Group by child symbol; only check cross-child conflicts
        by_child = Dict{Symbol,Set}()
        for (child_sym, ft) in owners
            push!(get!(by_child, child_sym, Set()), ft)
        end
        child_syms = collect(keys(by_child))
        length(child_syms) <= 1 && continue   # all from the same child — guarded, skip
        # Multiple different child symbols touch this column — check they agree
        fix_tos = [only(fts) for fts in values(by_child) if length(fts) == 1]
        length(fix_tos) == length(by_child) ||
            throw(
                ArgumentError(
                    "Column :$col is fixed by multiple child symbols, " *
                    "some of which have ambiguous fix_to values"
                )
            )
        length(unique(isequal, fix_tos)) == 1 ||
            throw(
                ArgumentError(
                    "Column :$col is targeted by rules with different child symbols " *
                    "($(child_syms)) and disagreeing fix_to values: $(fix_tos)"
                )
            )
    end

    return nothing
end
