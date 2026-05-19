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

    mean_deployment = [
        mean(sum(logs[:, c_id, :, :]; dims=:timesteps); dims=:locations).data[:]
        for c_id in deployed_corals
    ]

    total_deployment = [
        sum(sum(logs[:, c_id, :, :]; dims=:timesteps); dims=:locations).data[:]
        for c_id in deployed_corals
    ]

    col_names = ["fg_$(i)" for i in string.(collect(deployed_corals))]
    μ = DataFrame(hcat(mean_deployment...), "$(prefix)deployed_volume_mean_" .* col_names)
    T = DataFrame(hcat(total_deployment...), "$(prefix)volume_total_" .* col_names)

    return μ, T
end

"""
    feature_set(rs::ResultSet)::DataFrame

Extract a feature set from results for analysis purposes.
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
    dhw_means = dhw_stat[stat=At("mean")].data[:]
    dhw_stdevs = dhw_stat[stat=At("std")].data[:]

    insertcols!(scens, 2, :dhw_mean => -99.0, :dhw_stdev => -99.0)
    for (i, r) in enumerate(eachrow(scens))
        scens[i, :dhw_mean] = dhw_means[Int64(r.dhw_scenario)]
        scens[i, :dhw_stdev] = dhw_stdevs[Int64(r.dhw_scenario)]
    end

    @assert all(scens.dhw_mean .> -99.0) "Unknown DHW scenario found, check rows: $(findall(scens.dhw_mean .== -99.0))"
    @assert all(scens.dhw_stdev .> -99.0) "Unknown DHW scenario found, check rows: $(findall(scens.dhw_mean .== -99.0))"

    # Add indicators of deployments
    seed_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :seed)
    fog_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :fog)
    mc_stats = ADRIA.decision.deployment_summary_stats(rs.ranks, :mc)

    # Only attach mean of deployment effort
    insertcols!(
        scens,
        :n_loc_seed_mean => seed_stats[stats=At(:mean)].data[:],
        :n_loc_fog_mean => fog_stats[stats=At(:mean)].data[:],
        :n_loc_mc_mean => mc_stats[stats=At(:mean)].data[:]
    )

    # Replace `depth_offset` with maximum depth
    scens.depth_max = scens.depth_min .+ scens.depth_offset
    scens = scens[:, Not(:depth_offset)]

    seed_volume_mean, seed_volume_total = _iv_log_stats(rs.seed_log; prefix="seed_")
    DataFrames.hcat!(scens, seed_volume_mean)
    DataFrames.hcat!(scens, seed_volume_total)
    scens.seed_total_deployed_coral = sum.(eachrow(seed_volume_total))

    mc_volume_mean, mc_volume_total = _iv_log_stats(rs.mc_log; prefix="mc_")
    DataFrames.hcat!(scens, mc_volume_mean)
    DataFrames.hcat!(scens, mc_volume_total)
    scens.mc_total_deployed_coral = sum.(eachrow(mc_volume_total))

    # Remove `dhw_scenario` as Scenario IDs are not very informative for analyses
    scens = scens[:, Not(:dhw_scenario)]

    # Remove correlated features
    # Remove desired seed deployment targets
    # N_seed factors indicate maximum deployment effort, not actual simulated deployment
    scens = scens[:, .!contains.(names(scens), "N_seed")]

    # Set missing values to 0
    for col in eachcol(scens)
        replace!(col, missing => 0.0)
    end

    return _filter_constants(scens)
end
