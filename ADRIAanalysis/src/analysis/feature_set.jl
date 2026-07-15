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
    _seeding_stats(rs::ResultSet)::Tuple{Vector}

Extract total and average deployment, currently in terms of proportional increase to cover
relative to the locations' carrying capacity.

# Arguments
- `rs` : ResultSet holding scenario outcomes

# Returns
DataFrames of mean and total deployment for each coral group
"""
function _iv_log_stats(logs::YAXArray; prefix::String="")::Tuple{DataFrame,DataFrame}
    deployed_corals = logs.coral_id

    n = length(deployed_corals)
    mean_deployment = Vector{Vector{Float64}}(undef, n)
    total_deployment = Vector{Vector{Float64}}(undef, n)
    for (i, c_id) in enumerate(deployed_corals)
        ts = sum(logs[:, c_id, :, :]; dims=:timesteps)
        mean_deployment[i] = mean(ts; dims=:locations).data[:]
        total_deployment[i] = sum(ts; dims=:locations).data[:]
    end

    col_names = string.(ADRIA.functional_group_names())
    μ = DataFrame(hcat(mean_deployment...), "$(prefix)deployed_volume_mean_" .* col_names)
    T = DataFrame(hcat(total_deployment...), "$(prefix)volume_total_" .* col_names)

    return μ, T
end

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
`rs.inputs`.
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
scenario set only* (e.g. `n_loc_seed_mean_effort = (n_loc_seed_mean .-
minimum(n_loc_seed_mean)) ./ (maximum(n_loc_seed_mean) - minimum(n_loc_seed_mean))`).
This is analogous in spirit to `intervention_effort` (`performance.jl`), but
normalizes *actual simulated deployment* rather than pre-simulation sampled
targets. The resulting 0-1 scale is **relative to the most/least effort seen
among the scenarios explored in this particular analysis** -- it is NOT a
universal/absolute effort scale and is not comparable across different
scenario sets or studies. Columns that are constant within the current
scenario set are skipped (no `_effort` counterpart is emitted for them).
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

    idx = Int64.(scens.dhw_scenario)
    insertcols!(
        scens, 2,
        :dhw_mean => dhw_means[idx],
        :dhw_stdev => dhw_stdevs[idx],
        :dhw_complexity => dhw_complexities[idx]
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
    seed_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :seed)
    fog_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :fog)
    mc_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :mc)

    # Only attach mean of deployment effort
    insertcols!(
        scens,
        :n_loc_seed_mean => seed_stats[stats = At(:mean)].data[:],
        :n_loc_fog_mean => fog_stats[stats = At(:mean)].data[:],
        :n_loc_mc_mean => mc_stats[stats = At(:mean)].data[:]
    )
    colmetadata!(scens, :n_loc_seed_mean, "ptype", "continuous"; style=:note)
    colmetadata!(scens, :n_loc_seed_mean, "label", "Mean seeded locations"; style=:note)
    colmetadata!(scens, :n_loc_fog_mean, "ptype", "continuous"; style=:note)
    colmetadata!(scens, :n_loc_fog_mean, "label", "Mean fogged locations"; style=:note)
    colmetadata!(scens, :n_loc_mc_mean, "ptype", "continuous"; style=:note)
    colmetadata!(scens, :n_loc_mc_mean, "label", "Mean moving-coral locations"; style=:note)

    # Replace `depth_offset` with maximum depth
    scens.depth_max = scens.depth_min .+ scens.depth_offset
    colmetadata!(scens, :depth_max, "ptype", "ordered categorical"; style=:note)
    colmetadata!(scens, :depth_max, "label", "Maximum depth"; style=:note)
    scens = scens[:, Not(:depth_offset)]

    # Transform aggregate deployment totals into units of millions
    seed_volume_mean, seed_volume_total = _iv_log_stats(rs.seed_log; prefix="seed_")
    seed_volume_total_M = DataFrame(
        Matrix(seed_volume_total) ./ 1e6,
        replace.(names(seed_volume_total), "volume_total_" => "volume_total_M_")
    )
    DataFrames.hcat!(scens, seed_volume_mean)
    for col in names(seed_volume_mean)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(scens, col, "label", "Mean seed deployment volume ($col)"; style=:note)
    end
    DataFrames.hcat!(scens, seed_volume_total_M)
    for col in names(seed_volume_total_M)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, col, "label", "Total seed deployment volume ($col)"; style=:note
        )
    end
    scens.seed_total_deployed_coral_M = vec(sum(Matrix(seed_volume_total_M); dims=2))
    colmetadata!(scens, :seed_total_deployed_coral_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :seed_total_deployed_coral_M, "label",
        "Total seeded coral deployment (millions)"; style=:note
    )

    mc_volume_mean, mc_volume_total = _iv_log_stats(rs.mc_log; prefix="mc_")
    mc_volume_total_M = DataFrame(
        Matrix(mc_volume_total) ./ 1e6,
        replace.(names(mc_volume_total), "volume_total_" => "volume_total_M_")
    )
    DataFrames.hcat!(scens, mc_volume_mean)
    for col in names(mc_volume_mean)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, col, "label", "Mean moving-coral deployment volume ($col)"; style=:note
        )
    end
    DataFrames.hcat!(scens, mc_volume_total_M)
    for col in names(mc_volume_total_M)
        colmetadata!(scens, col, "ptype", "continuous"; style=:note)
        colmetadata!(
            scens, col, "label", "Total moving-coral deployment volume ($col)"; style=:note
        )
    end
    scens.mc_total_deployed_coral_M = vec(sum(Matrix(mc_volume_total_M); dims=2))
    colmetadata!(scens, :mc_total_deployed_coral_M, "ptype", "continuous"; style=:note)
    colmetadata!(
        scens, :mc_total_deployed_coral_M, "label",
        "Total moving-coral deployment (millions)"; style=:note
    )

    # Add normalized "actual effort" columns: min-max normalization of realized
    # deployment columns, computed over the current scenario set (relative, not
    # absolute/universal -- see docstring). Columns constant within this scenario
    # set are skipped (min == max would divide by zero).
    effort_source_cols = filter(
        c -> c in names(scens),
        vcat(
            ["n_loc_seed_mean", "n_loc_fog_mean", "n_loc_mc_mean"],
            names(seed_volume_mean), names(seed_volume_total_M),
            names(mc_volume_mean), names(mc_volume_total_M)
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
    # Remove seed deployment target values as `N_seed_*` factors indicate
    # maximum (desired) deployment effort, not actual simulated deployment
    scens = scens[:, .!contains.(names(scens), "N_seed")]

    # Set missing values to 0
    for col in eachcol(scens)
        replace!(col, missing => 0.0)
    end

    return _filter_constants(scens)
end
