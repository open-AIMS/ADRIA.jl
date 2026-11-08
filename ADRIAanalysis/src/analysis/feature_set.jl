"""
    _filter_constants(scens::DataFrame)::DataFrame

Filter out features/factors that do not vary.
"""
function _filter_constants(scens::DataFrame)::DataFrame
    varying_cols = [
        i for (i, col) in enumerate(eachcol(scens))
              if !all(==(first(col)), col)
    ]

    return scens[:, varying_cols]
end

"""
    _iv_CAq_stats(rs::ResultSet)::Tuple{Vector}

Extract total and average deployment, currently in terms of proportional increase to cover
relative to the locations' carrying capacity.

# Arguments
- `rs` : ResultSet holding scenario outcomes

# Returns
DataFrames of mean, total and event-count deployment for each coral group
"""
function _iv_log_stats(
    logs::YAXArray; prefix::String=""
)::Tuple{DataFrame,DataFrame,DataFrame}
    deployed_corals = logs.coral_id

    n = length(deployed_corals)
    mean_deployment = Vector{Vector{Float64}}(undef, n)
    total_deployment = Vector{Vector{Float64}}(undef, n)
    n_events = Vector{Vector{Float64}}(undef, n)
    for (i, c_id) in enumerate(deployed_corals)
        group_log = logs[:, c_id, :, :]
        ts = sum(group_log; dims=:timesteps)
        mean_deployment[i] = mean(ts; dims=:locations).data[:]
        total_deployment[i] = sum(ts; dims=:locations).data[:]
        # A deployment event is a timestep where this group received a nonzero volume at
        # any location -- counted straight off the log rather than derived from the
        # CAq/LvM schedule levers, since a reactive strategy's event timing depends on
        # simulated cover, not just the levers.
        per_timestep_total = sum(group_log; dims=:locations)
        n_events[i] = Float64.(sum(per_timestep_total .> 0; dims=:timesteps).data[:])
    end

    col_names = string.(ADRIA.functional_group_names())
    μ = DataFrame(hcat(mean_deployment...), "$(prefix)deployed_volume_mean_" .* col_names)
    T = DataFrame(hcat(total_deployment...), "$(prefix)volume_total_" .* col_names)
    N = DataFrame(hcat(n_events...), "$(prefix)n_events_" .* col_names)

    return μ, T, N
end

"""
    _iv_event_count(logs::YAXArray)::Vector{Float64}

Count deployment events per scenario: timesteps where ANY group received a nonzero
deployment at any location. Used for larval methods, whose target (`iv_LvM_N_settlers`)
is a single scalar spanning all groups rather than a per-group lever, so events must be
counted across groups together rather than per-group as in `_iv_log_stats`.
"""
function _iv_event_count(logs::YAXArray)::Vector{Float64}
    any_group_total = sum(logs; dims=(:coral_id, :locations))
    return Float64.(sum(any_group_total .> 0; dims=:timesteps).data[:])
end

"""
    _iv_unique_locations(logs::YAXArray)::Vector{Float64}

Number of distinct locations that received any deployment (any group, any timestep) over
the whole run, per scenario: the cumulative spatial footprint. Unlike `n_loc_*_mean`
(locations per deployment year, set mostly by the minimum-locations lever), this also
reflects how much the selection rotates between reefs across years (revisit cadence,
state-responsive criteria such as `iv_CAq_juvenile_share`).
"""
function _iv_unique_locations(logs::YAXArray)::Vector{Float64}
    per_loc_total = sum(logs; dims=(:timesteps, :coral_id))
    return Float64.(sum(per_loc_total .> 0; dims=:locations).data[:])
end

# Maps each coral aquaculture-deployable functional group (as named by
# `ADRIA.functional_group_names()`) to its per-event deployment-target lever. Only 5 of
# the 6 defined functional groups are currently deployable -- `arborescent_Acropora` has
# no `iv_CAq_N_*` lever.
const _CAQ_TARGET_LEVER = Dict(
    "tabular_Acropora" => "iv_CAq_N_TA",
    "corymbose_Acropora" => "iv_CAq_N_CA",
    "corymbose_non_Acropora" => "iv_CAq_N_CNA",
    "small_massives" => "iv_CAq_N_SM",
    "large_massives" => "iv_CAq_N_LM"
)

"""
    dhw_spatial_features(rs::ResultSet)::DataFrame

Compute per-scenario spatial-heterogeneity features from time-mean DHW exposure across
locations: the coefficient of variation of per-site mean DHW, and the gap between the
median site and the coolest available refugia (10th percentile site). These characterize
how much spatial contrast in heat exposure a scenario's `dhw_scenario` draw offers --
i.e. whether there is any "cooler refugia" for guided site selection to target.

# Arguments
- `rs` : ResultSet holding scenario outcomes

# Returns
DataFrame with `dhw_site_cv` and `dhw_refugia_gap` columns, one row per scenario in
`rs.inputs`. Rows whose `dhw_scenario` is a no-heat-stress counterfactual (`< 1`, e.g. a
row sampled with `dhw_scenario = 0`) get `0.0` for both features.
"""
function dhw_spatial_features(rs::ResultSet)::DataFrame
    rcp_ids = collect(keys(rs.dhw_stats))
    rcp_id = first(rcp_ids)

    # dims: (scenarios, locations) -- time-mean DHW per site, per dhw_scenario draw
    site_means = rs.dhw_stats[rcp_id][stat = At("mean")]

    idx = Int64.(rs.inputs.dhw_scenario)
    n = length(idx)
    cv = Vector{Float64}(undef, n)
    refugia_gap = Vector{Float64}(undef, n)
    for (i, s_idx) in enumerate(idx)
        # `dhw_scenario < 1` marks a no-heat-stress counterfactual (e.g. a row sampled
        # with `dhw_scenario = 0`): it has no entry in the per-trajectory DHW stats, and
        # by definition there is no spatial heat-stress contrast, so both features are 0.
        if s_idx < 1
            cv[i] = 0.0
            refugia_gap[i] = 0.0
            continue
        end
        sites = site_means[s_idx, :].data[:]
        μ = mean(sites)
        cv[i] = μ == 0 ? 0.0 : std(sites; mean=μ) / μ
        refugia_gap[i] = median(sites) - quantile(sites, 0.1)
    end

    return DataFrame(; dhw_site_cv=cv, dhw_refugia_gap=refugia_gap)
end

"""
    feature_set(rs::ResultSet)::DataFrame

Extract a feature set from results for analysis purposes.

In addition to the raw realized-deployment columns (`n_loc_*_mean`,
`*_volume_mean_*`, `*_volume_total_M_*`), a `*_effort` counterpart is added for
each: a min-max normalization of that column computed over the *current
scenario set only* (e.g. `n_loc_CAq_mean_effort = (n_loc_CAq_mean .-
minimum(n_loc_CAq_mean)) ./ (maximum(n_loc_CAq_mean) - minimum(n_loc_CAq_mean))`).
This is analogous in spirit to `intervention_effort` (`performance.jl`), but
normalizes *actual simulated deployment* rather than pre-simulation sampled
targets. The resulting 0-1 scale is **relative to the most/least effort seen
among the scenarios explored in this particular analysis** -- it is NOT a
universal/absolute effort scale and is not comparable across different
scenario sets or studies. Columns that are constant within the current
scenario set are skipped (no `_effort` counterpart is emitted for them).

Scenarios whose `dhw_scenario` is a no-heat-stress counterfactual (`< 1`, e.g.
a row sampled with `dhw_scenario = 0`) have no entry in the per-trajectory DHW
statistics; their DHW-derived columns (`dhw_mean`, `dhw_stdev`,
`dhw_complexity`, `dhw_site_cv`, `dhw_refugia_gap`) are set to `0.0`.

Also adds a target-vs-realised comparison for each intervention: `CAq_total_target_coral_M`/
`LvM_total_target_coral_M` (each group's per-event target lever times how many times it
actually deployed, counted from the log) and `CAq_deployment_delta_M`/
`LvM_deployment_delta_M` (realised minus target; negative means the plan under-delivered,
whether from skipped revisits or a deployment throttled by site capacity).

`n_unique_loc_CAq`/`n_unique_loc_LvM` count the distinct locations deployed to over the
whole run (cumulative footprint), as opposed to `n_loc_*_mean` (per deployment year).
"""
function feature_set(rs::ResultSet)::DataFrame
    scens = copy(rs.inputs)

    rcp_ids = collect(keys(rs.dhw_stats))
    rcp_id = first(rcp_ids)
    if length(rcp_ids) > 1
        @warn "Multiple RCPs found, assigning stats for first id: $(rcp_id)"
    end

    # Add DHW statistics
    dhw_stat = mean(rs.dhw_stats[rcp_id]; dims=:locations)
    dhw_means = dhw_stat[stat = At("mean")].data[:]
    dhw_stdevs = dhw_stat[stat = At("std")].data[:]
    dhw_complexities = dhw_stat[stat = At("complexity")].data[:]

    # `dhw_scenario < 1` marks a no-heat-stress counterfactual (e.g. a row sampled with
    # `dhw_scenario = 0`): no entry exists for it in the per-trajectory DHW stats arrays,
    # and its DHW-derived features are zero by definition. `safe_idx` keeps the lookup
    # in bounds for those rows; `ifelse` then zeroes them.
    idx = Int64.(scens.dhw_scenario)
    valid = idx .>= 1
    safe_idx = ifelse.(valid, idx, 1)
    insertcols!(
        scens, 2,
        :dhw_mean => ifelse.(valid, dhw_means[safe_idx], zero(eltype(dhw_means))),
        :dhw_stdev => ifelse.(valid, dhw_stdevs[safe_idx], zero(eltype(dhw_stdevs))),
        :dhw_complexity =>
            ifelse.(valid, dhw_complexities[safe_idx], zero(eltype(dhw_complexities)))
    )
    colmetadata!(scens, :dhw_mean, "ptype", "continuous"; style=:note)
    colmetadata!(scens, :dhw_mean, "label", "Mean DHW"; style=:note)
    colmetadata!(scens, :dhw_stdev, "ptype", "continuous"; style=:note)
    colmetadata!(scens, :dhw_stdev, "label", "DHW standard deviation"; style=:note)
    colmetadata!(scens, :dhw_complexity, "ptype", "continuous"; style=:note)
    colmetadata!(scens, :dhw_complexity, "label", "DHW time series complexity"; style=:note)

    # Add DHW spatial-heterogeneity features (is there any cooler refugia to target?)
    dhw_spatial = dhw_spatial_features(rs)
    DataFrames.hcat!(scens, dhw_spatial)
    colmetadata!(scens, :dhw_site_cv, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :dhw_site_cv, "label",
        "Spatial CV of mean site-level DHW"; style=:note
    )
    colmetadata!(scens, :dhw_refugia_gap, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :dhw_refugia_gap, "label",
        "Median minus 10th-percentile site DHW"; style=:note
    )

    # Add indicators of deployments
    CAq_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :caq)
    Fog_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :fog)
    LvM_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :lvm)

    # Only attach mean of deployment effort
    insertcols!(
        scens,
        :n_loc_CAq_mean => CAq_stats[stats = At(:mean)].data[:],
        :n_loc_Fog_mean => Fog_stats[stats = At(:mean)].data[:],
        :n_loc_LvM_mean => LvM_stats[stats = At(:mean)].data[:]
    )
    colmetadata!(scens, :n_loc_CAq_mean, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :n_loc_CAq_mean, "label", "Mean coral aquaculture locations"; style=:note
    )
    colmetadata!(scens, :n_loc_Fog_mean, "ptype", "continuous"; style=:note)
    colmetadata!(scens, :n_loc_Fog_mean, "label", "Mean fogged locations"; style=:note)
    colmetadata!(scens, :n_loc_LvM_mean, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :n_loc_LvM_mean, "label", "Mean larval methods locations"; style=:note
    )

    # Cumulative footprint: distinct locations deployed to over the run
    scens.n_unique_loc_CAq = _iv_unique_locations(rs.iv_CAq_log)
    scens.n_unique_loc_LvM = _iv_unique_locations(rs.iv_LvM_log)
    colmetadata!(scens, :n_unique_loc_CAq, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :n_unique_loc_CAq, "label",
        "Distinct coral aquaculture locations over the run"; style=:note
    )
    colmetadata!(scens, :n_unique_loc_LvM, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :n_unique_loc_LvM, "label",
        "Distinct larval methods locations over the run"; style=:note
    )

    # Replace `depth_offset` with maximum depth
    scens.depth_max = scens.depth_min .+ scens.depth_offset
    colmetadata!(scens, :depth_max, "ptype", "ordered categorical"; style=:note)
    colmetadata!(scens, :depth_max, "label", "Maximum depth"; style=:note)
    scens = scens[:, Not(:depth_offset)]

    # Transform aggregate deployment totals into units of millions
    CAq_volume_mean, CAq_volume_total, CAq_n_events = _iv_log_stats(
        rs.iv_CAq_log; prefix="CAq_"
    )
    CAq_volume_total_M = DataFrame(
        Matrix(CAq_volume_total) ./ 1e6,
        replace.(names(CAq_volume_total), "volume_total_" => "volume_total_M_")
    )
    DataFrames.hcat!(scens, CAq_volume_mean)
    for col in names(CAq_volume_mean)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens,
            col,
            "label",
            "Mean coral aquaculture deployment volume ($col)";
            style=:note
        )
    end
    DataFrames.hcat!(scens, CAq_volume_total_M)
    for col in names(CAq_volume_total_M)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, col, "label", "Total coral aquaculture deployment volume ($col)";
            style=:note
        )
    end
    scens.CAq_total_deployed_coral_M = vec(sum(Matrix(CAq_volume_total_M); dims=2))
    colmetadata!(scens, :CAq_total_deployed_coral_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :CAq_total_deployed_coral_M, "label",
        "Total coral aquaculture deployment (millions)"; style=:note
    )

    # Target cumulative deployment: each group's per-event target (`iv_CAq_N_*`, still
    # present in `scens` at this point) times how many times that group actually
    # deployed (`CAq_n_events`, from the log -- not the schedule levers, so it reflects
    # a reactive strategy's actual event count, not a periodic-schedule assumption).
    # `CAq_deployment_delta_M` is realised minus target: negative means the plan
    # under-delivered against what was intended (revisits skipped, or a deployment event
    # throttled by site capacity), positive means it over-delivered.
    CAq_groups = filter(
        g -> haskey(_CAQ_TARGET_LEVER, g), string.(ADRIA.functional_group_names())
    )
    CAq_target_cols = [
        Float64.(scens[!, _CAQ_TARGET_LEVER[g]]) .* CAq_n_events[!, "CAq_n_events_$(g)"] ./
        1e6
        for g in CAq_groups
    ]
    CAq_target_total_M = DataFrame(
        hcat(CAq_target_cols...), "CAq_target_volume_total_M_" .* CAq_groups
    )
    DataFrames.hcat!(scens, CAq_target_total_M)
    for col in names(CAq_target_total_M)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, col, "label", "Target coral aquaculture deployment volume ($col)";
            style=:note
        )
    end
    scens.CAq_total_target_coral_M = vec(sum(Matrix(CAq_target_total_M); dims=2))
    colmetadata!(scens, :CAq_total_target_coral_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :CAq_total_target_coral_M, "label",
        "Total target coral aquaculture deployment (millions)"; style=:note
    )
    scens.CAq_deployment_delta_M =
        scens.CAq_total_deployed_coral_M .- scens.CAq_total_target_coral_M
    colmetadata!(scens, :CAq_deployment_delta_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :CAq_deployment_delta_M, "label",
        "Realised minus target coral aquaculture deployment (millions)"; style=:note
    )

    # Per-group event counts aren't used for LvM's target (`iv_LvM_N_settlers` is a single
    # scalar spanning all groups, not a per-group lever like CAq's `iv_CAq_N_*`) -- see
    # `_iv_event_count` below.
    LvM_volume_mean, LvM_volume_total, _ = _iv_log_stats(rs.iv_LvM_log; prefix="LvM_")
    LvM_volume_total_M = DataFrame(
        Matrix(LvM_volume_total) ./ 1e6,
        replace.(names(LvM_volume_total), "volume_total_" => "volume_total_M_")
    )
    DataFrames.hcat!(scens, LvM_volume_mean)
    for col in names(LvM_volume_mean)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, col, "label", "Mean larval methods deployment volume ($col)"; style=:note
        )
    end
    DataFrames.hcat!(scens, LvM_volume_total_M)
    for col in names(LvM_volume_total_M)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, col, "label", "Total larval methods deployment volume ($col)";
            style=:note
        )
    end
    scens.LvM_total_deployed_coral_M = vec(sum(Matrix(LvM_volume_total_M); dims=2))
    colmetadata!(scens, :LvM_total_deployed_coral_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :LvM_total_deployed_coral_M, "label",
        "Total larval methods deployment (millions)"; style=:note
    )

    # Target cumulative deployment for LvM, mirroring CAq's above -- `iv_LvM_N_settlers`
    # (a single target spanning all groups) times how many times LvM actually deployed.
    LvM_n_events = _iv_event_count(rs.iv_LvM_log)
    scens.LvM_total_target_coral_M =
        Float64.(scens.iv_LvM_N_settlers) .* LvM_n_events ./ 1e6
    colmetadata!(scens, :LvM_total_target_coral_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :LvM_total_target_coral_M, "label",
        "Total target larval methods deployment (millions)"; style=:note
    )
    scens.LvM_deployment_delta_M =
        scens.LvM_total_deployed_coral_M .- scens.LvM_total_target_coral_M
    colmetadata!(scens, :LvM_deployment_delta_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :LvM_deployment_delta_M, "label",
        "Realised minus target larval methods deployment (millions)"; style=:note
    )

    # Add normalized "actual effort" columns: min-max normalization of realized
    # deployment columns, computed over the current scenario set (relative, not
    # absolute/universal -- see docstring). Columns constant within this scenario
    # set are skipped (min == max would divide by zero).
    effort_source_cols = filter(
        c -> c in names(scens),
        vcat(
            ["n_loc_CAq_mean", "n_loc_Fog_mean", "n_loc_LvM_mean",
                "n_unique_loc_CAq", "n_unique_loc_LvM"],
            names(CAq_volume_mean), names(CAq_volume_total_M),
            names(LvM_volume_mean), names(LvM_volume_total_M)
        )
    )
    for col in effort_source_cols
        vals = Float64.(scens[!, col])
        lo, hi = extrema(vals)
        hi == lo && continue
        effort_col = Symbol(col, "_effort")
        scens[!, effort_col] = (vals .- lo) ./ (hi - lo)
        colmetadata!(scens, effort_col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, effort_col, "label",
            "Normalized realized effort ($col); relative to this scenario set only, " *
            "not a universal/absolute effort scale"; style=:note
        )
    end

    # Remove `dhw_scenario` as Scenario IDs are not very informative for analyses
    scens = scens[:, Not(:dhw_scenario)]

    # Remove correlated features
    # Remove coral aquaculture deployment target values as `iv_CAq_N_*` factors indicate
    # maximum (desired) deployment effort, not actual simulated deployment
    scens = scens[:, .!contains.(names(scens), "iv_CAq_N")]

    # Set missing values to 0
    for col in eachcol(scens)
        replace!(col, missing => 0.0)
    end

    return _filter_constants(scens)
end
