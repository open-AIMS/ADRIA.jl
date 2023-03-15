using ADRIA: mcda_normalize, DMCDA_vars


function zones_criteria(zones, priority_zones, strong_pred, site_ids)
    n_sites = size(zones, 1)
    zone_ids = intersect(priority_zones, unique(zones))
    zone_weights = mcda_normalize(collect(length(zone_ids):-1:1))
    zone_preds = zeros(n_sites, 1)
    zone_sites = zeros(n_sites, 1)

    for k in axes(zone_ids, 1)
        # find sites which are strongest predecessors of sites in the zone
        zone_preds_temp = strong_pred[zones.==zone_ids[k]]
        for s in unique(zone_preds_temp)
            # for each predecessor site, add zone_weights* (no. of zone sites the site is a strongest predecessor for)
            zone_preds[site_ids.==s] .= zone_preds[site_ids.==s] .+ (zone_weights[k]) .* sum(zone_preds_temp .== s)
        end
        #     # add zone_weights for sites in the zone (whether a strongest predecessor of a zone or not)
        zone_sites[zones.==zone_ids[k]] .= zone_weights[k]
    end

    # # add weights for strongest predecessors and zones to get zone criteria
    zones_criteria = zone_preds .+ zone_sites
    return zones_criteria
end

function priority_predecessor_criteria(strong_pred, priority_sites)
    # # work out which priority predecessors are connected to priority sites
    predec::Array{Float64} = zeros(length(strong_pred), 3)
    predec[:, 1:2] .= strong_pred
    predprior = predec[in.(predec[:, 1], [priority_sites']), 2]
    predprior = [x for x in predprior if !isnan(x)]

    predec[predprior, 3] .= 1.0
    return predec
end

function coral_cover_criteria(site_data, coral_cover)
    max_area = site_data.k .* site_data.area
    coral_cover_area = site_data.area .* coral_cover'
    return max.(coral_cover_area, 0.0), max.(max_area .- coral_cover_area, 0.0)
end

function env_stress_criteria(env_stress)
    return maximum(env_stress) != 0.0 ? (env_stress .- minimum(env_stress)) ./ (maximum(env_stress) - minimum(env_stress)) : zeros(Float64, size(env_stress))
end

function connectivity_criteria(conn, sum_cover, area)
    cov_area = conn .* sum_cover .* area
    return maximum(cov_area) != 0.0 ? cov_area / maximum(cov_area) : zeros(Float64, size(cov_area))
end

function initialize_mcda(domain, param_set, sim_params, site_data, depth_priority, init_sum_cover, area_to_seed)
    n_sites = length(site_data.site_id)
    rankings = [depth_priority zeros(Int, length(depth_priority)) zeros(Int, length(depth_priority))]

    # initialize weights
    weights = create_weights_df(param_set("coral_cover_high"), param_set("coral_cover_low"),
        param_set("in_seed_connectivity"), param_set("out_seed_connectivity"),
        param_set("shade_connectivity"), param_set("heat_stress"), param_set("wave_stress"),
        (param_set("shade_priority"), "predec_shade_wt"),
        (param_set("seed_priority"), "predec_seed_wt"),
        (param_set("zone_seed"), "zones_seed_wt"),
        (param_set("zone_shade"), "zones_shade_wt"))

    # initialize thresholds
    thresholds = create_thresholds_df([param_set("coral_cover_tol") .* area_to_seed, "gt"],
        [param_set("deployed_coral_risk_tol"), "lt"],
        [param_set("deployed_coral_risk_tol"), "lt"])

    # initialize criteria
    zones = zones_criteria(site_data.zone_type, sim_params.priority_zones, domain.strong_pred, collect(1:n_sites))
    predec = priority_predecessor_criteria(domain.strong_pred, sim_params.priority_sites)
    coral_cover, coral_space = coral_cover_criteria(site_data, init_sum_cover)
    heat_stress = zeros(1, n_sites)
    wave_stress = zeros(1, n_sites)

    # Prep site selection
    mcda_vars = DMCDA_vars(domain, sim_params.seed_criteria_names, sim_params.shade_criteria_names,
        param_set("use_dist"), domain.median_site_distance - domain.median_site_distance * param_set("dist_thresh"),
        Int(param_set("top_n")), weights, thresholds)

    criteria_df = create_criteria_df(depth_priority, coral_cover, coral_space,
        domain.in_conn, domain.out_conn, heat_stress, wave_stress,
        ("zones", zones), ("predec", predec))

    return rankings, mcda_vars, criteria_df

end

function update_criteria_df!(criteria_df, wave_stress, heat_stress, in_conn, out_conn, site_area, site_coral_cover, site_data, depth_priority)
    criteria_df[:, "wave_stress"] .= env_stress_criteria(wave_stress)[depth_priority]
    criteria_df[:, "heat_stress"] .= env_stress_criteria(heat_stress)[depth_priority]
    criteria_df[:, "connectivity_in"] .= connectivity_criteria(in_conn, site_coral_cover, site_area)[depth_priority]
    criteria_df[:, "connectivity_out"] .= connectivity_criteria(out_conn, site_coral_cover, site_area)[depth_priority]
    coral_cover, coral_space = coral_cover_criteria(site_data, site_coral_cover)
    criteria_df[:, "coral_cover"] .= coral_cover[depth_priority]
    criteria_df[:, "coral_space"] .= coral_space[depth_priority]
    return criteria_df
end